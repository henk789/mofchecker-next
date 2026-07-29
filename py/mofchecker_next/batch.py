"""Batch validation API.

Run MOFChecker diagnostics over many structures efficiently -- useful for
validating generated structures (e.g. from a diffusion model).

- Accepts pymatgen ``Structure``/``IStructure``, ASE ``Atoms``, or CIF path
  (``str``/``Path``) interchangeably -- or a mix.
- Builds each structure's graph once (via the ``MOFChecker`` class) and reuses it
  across all checks.
- Parallelizes across structures with ``multiprocessing`` (CPU-bound Python;
  structures are independent).
- Lets you pick a subset of descriptors.

Example
-------
>>> from mofchecker_next.batch import check_structures
>>> results = check_structures(list_of_atoms_or_structures_or_paths, n_workers=16)
>>> bad = [r for r in results if r["has_atomic_overlaps"]]
"""

from __future__ import annotations

import concurrent.futures as cf
import json
import os
import queue
import signal
import subprocess
import sys
import threading
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Sequence

from mofchecker_next.core import DEFAULT_DESCRIPTORS as ALL_DESCRIPTORS  # full set (metadata + hashes + diagnostics)
from mofchecker_next.core import LEGACY_DEFAULT_DESCRIPTORS as LEGACY_ALL_DESCRIPTORS
from mofchecker_next.core import NEXT_DEFAULT_DESCRIPTORS as NEXT_ALL_DESCRIPTORS
from mofchecker_next.checks.charge_oms import METALS as DEFAULT_METALS
from mofchecker_next.core import MOFChecker, normalize_structure
from mofchecker_next.profiles import (
    ADIT_PRESENCE_FLAGS,
    ADIT_PROBLEM_FLAGS,
    DEFAULT_MODE,
    MODES,
    profile_manifest,
    resolve_profile,
)

# Default batch descriptor set: the validation-relevant diagnostics (fast path).
# It omits the metadata/symmetry/hash descriptors (which add SpacegroupAnalyzer
# and graph-hash cost); request those explicitly or use ``ALL_DESCRIPTORS``.
DEFAULT_DESCRIPTORS = (
    "has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_stray_atom", "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
    "has_high_charges",
)
NEXT_DEFAULT_DESCRIPTORS = DEFAULT_DESCRIPTORS + (
    "max_abs_eqeq_charge", "eqeq_charge_sum", "eqeq_expected_total_charge",
    "eqeq_charge_residual", "eqeq_charge_threshold",
)
LEGACY_DEFAULT_DESCRIPTORS = (
    "has_carbon", "has_hydrogen", "has_metal", "has_atomic_overlaps",
    "has_overcoordinated_c", "has_overcoordinated_n", "has_overcoordinated_h",
    "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "has_high_charges",
)

# Composite validity = ADiT / Mofasa "Validity rate (all passed)" (Mofasa paper
# Appendix G, Table 4): a structure is valid iff all 3 presence checks are True
# and all 12 problem checks are False. NB this differs from the older native
# composite -- it adds has_hydrogen + has_geometrically_exposed_metal and
# DROPS has_3d_connected_graph.
PROBLEM_FLAGS = ADIT_PROBLEM_FLAGS
PRESENCE_FLAGS = ADIT_PRESENCE_FLAGS


def _input_id(obj, index: int) -> str:
    if isinstance(obj, (str, Path)):
        return Path(obj).name
    return str(index)


def _input_protocol(obj, corrected: bool) -> dict:
    if isinstance(obj, (str, Path)):
        input_type = "cif"
        normalization = "CifParser(frac_tolerance=0,site_tolerance=0,primitive=False)" if corrected else "compatibility-parser"
    elif hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_positions"):
        input_type, normalization = "ase.Atoms", "AseAtomsAdaptor"
    else:
        input_type, normalization = "pymatgen.Structure", "preserve-cell-and-sites"
    return {
        "input_type": input_type,
        "normalization": normalization,
        "primitive": False if corrected else None,
        "symmetrized": False if corrected else None,
        "coordinate_snapping": False if corrected else None,
    }


def _result_metadata(profile, mode: str, obj, names, metals, method, distance_scale, clash_scale) -> dict:
    identity = profile_manifest(mode)
    return {
        "mode": mode,
        "resolved_profile": profile.id,
        "profile_provisional": profile.provisional,
        "implementation": {key: identity[key] for key in (
            "package_version", "source_sha", "source_dirty", "source_tree_sha256",
            "rust_extension_sha256",
        )},
        "composite_name": profile.composite_name,
        "composite_problem_flags": list(profile.problem_flags),
        "composite_presence_flags": list(profile.presence_flags),
        "report_only_fields": list(profile.report_only_fields),
        "input_protocol": _input_protocol(obj, profile.corrected_input),
        "evaluation_settings": {
            "method": method,
            "metals": sorted(str(metal) for metal in (DEFAULT_METALS if metals is None else metals)),
            "distance_scale": float(distance_scale),
            "clash_scale": float(clash_scale),
            "requested_descriptors": list(names),
            "total_charge": profile.total_charge,
            "eqeq_threshold": profile.eqeq_threshold,
            "eqeq_charge_sum_tolerance": profile.eqeq_charge_sum_tolerance,
        },
    }


def _atom_evidence(checker, field: str) -> tuple[list[dict], dict]:
    properties = {
        "has_atomic_overlaps": checker.get_overlapping_indices,
        "has_overcoordinated_c": lambda: checker.overvalent_c_indices,
        "has_overcoordinated_n": lambda: checker.overcoordinated_n_indices,
        "has_overcoordinated_h": lambda: checker.overvalent_h_indices,
        "has_stray_atom": lambda: [i for component in checker.stray_atom_indices for i in component],
        "has_lone_molecule": lambda: [i for component in checker.lone_molecule_indices for i in component],
        "has_high_charges": lambda: [
            i for i, charge in enumerate(checker.eqeq_charges)
            if abs(charge) > checker.eqeq_charge_threshold
        ],
    }
    if field not in properties:
        return [], {}
    indices = sorted(set(int(i) for i in properties[field]()))
    atoms = [{"index": index, "image": [0, 0, 0]} for index in indices]
    values = {}
    if field == "has_atomic_overlaps":
        values["contacts"] = [
            {
                "atoms": [
                    {"index": atom.index, "image": list(atom.image)}
                    for atom in diagnostic.atoms
                ],
                **diagnostic.values,
            }
            for diagnostic in checker.overlap_diagnostics
        ]
    elif field in {"has_overcoordinated_c", "has_overcoordinated_n"}:
        values["neighbors"] = {
            str(index): [
                {
                    "index": int(site.index),
                    "image": list(getattr(site, "jimage", (0, 0, 0))),
                    "distance_angstrom": float(getattr(site, "dist", 0.0)),
                }
                for site in checker.graph.get_connected_sites(index)
            ]
            for index in indices
        }
    elif field in {"has_stray_atom", "has_lone_molecule"}:
        values["components"] = [list(map(int, component)) for component in (
            checker.stray_atom_indices if field == "has_stray_atom" else checker.lone_molecule_indices
        )]
    elif field == "has_high_charges":
        values = {
            "threshold": checker.eqeq_charge_threshold,
            "max_abs_charge": checker.max_abs_eqeq_charge,
            "charges": list(checker.eqeq_charges),
        }
    return atoms, values


def _corrected_check_results(result: dict, checker, profile) -> dict:
    checks = {}
    errors = result.get("errors", {})
    rules = {
        "input_sanity": "finite non-singular ordered structure with real elements",
        "has_carbon": "application scope contains carbon",
        "has_metal": "application scope contains a configured metal",
        "has_atomic_overlaps": "pair distance is below scaled covalent-radius cutoff",
        "has_overcoordinated_c": "carbon graph coordination exceeds four after exclusions",
        "has_overcoordinated_n": "nitrogen graph coordination exceeds four after exclusions",
        "has_overcoordinated_h": "hydrogen has more than one neighbor within 1.1 angstrom",
        "has_stray_atom": "finite detached component contains one atom",
        "has_lone_molecule": "finite detached component contains multiple atoms",
        "has_high_charges": "true-element EQeq max absolute partial charge exceeds threshold",
    }
    for field in (*profile.presence_flags, *profile.problem_flags):
        if field in errors:
            unsupported = "outside EQeq's parameter table" in errors[field]
            checks[field] = {
                "status": "unsupported" if unsupported else "indeterminate",
                "severity": "error",
                "error": errors[field],
                "rule": rules[field], "atoms": [], "values": {},
            }
            continue
        value = result.get(field)
        if not isinstance(value, bool):
            checks[field] = {
                "status": "indeterminate", "severity": "error",
                "error": "required check missing or non-Boolean", "rule": rules[field],
                "atoms": [], "values": {},
            }
            continue
        is_problem = field in profile.problem_flags
        failed = value if is_problem else not value
        atoms, values = _atom_evidence(checker, field)
        checks[field] = {
            "status": "fail" if failed else "pass",
            "severity": "error" if failed else "info",
            "rule": rules[field],
            "atoms": atoms,
            "values": values,
        }
    return checks


def _descriptor_names(profile, descriptors: Sequence[str] | None) -> list[str]:
    defaults = {
        "legacy": LEGACY_DEFAULT_DESCRIPTORS,
        "modern": DEFAULT_DESCRIPTORS,
        "next": NEXT_DEFAULT_DESCRIPTORS,
    }[profile.descriptor_set]
    return list(descriptors) if descriptors is not None else list(defaults)


def _profile_statuses(result: dict, profile) -> dict[str, str]:
    statuses = {}
    errors = result.get("errors", {})
    for name, kind, fields in profile.status_groups:
        if (
            result.get("error")
            or any(field in errors for field in fields)
            or any(field not in result or not isinstance(result[field], bool) for field in fields)
        ):
            statuses[name] = "indeterminate"
        elif kind == "problem":
            statuses[name] = "fail" if any(result[field] for field in fields) else "pass"
        else:
            statuses[name] = "pass" if all(result[field] for field in fields) else "fail"
    return statuses


def _failed_result(
    obj, *, mode: str, message: str, descriptors: Sequence[str] | None = None,
    metals=None, method: str = "vesta", distance_scale: float = 1.0, clash_scale: float = 1.0,
) -> dict:
    profile = resolve_profile(mode)
    names = _descriptor_names(profile, descriptors)
    result = _result_metadata(profile, mode, obj, names, metals, method, distance_scale, clash_scale)
    result.update({"error": message, "composite_status": "indeterminate", "valid": None})
    if profile.corrected_input:
        result["input_sanity"] = None
        result["errors"] = {"input_sanity": message}
        result["check_results"] = {
            field: {
                "status": "indeterminate", "severity": "error",
                "error": message if field == "input_sanity" else "not evaluated after structure failure",
                "rule": "required corrected-profile check", "atoms": [], "values": {},
            }
            for field in (*profile.presence_flags, *profile.problem_flags)
        }
    result.update(_profile_statuses(result, profile))
    result["composites"] = {
        profile.composite_name: {
            "status": "indeterminate", "value": None,
            "problem_flags": list(profile.problem_flags),
            "presence_flags": list(profile.presence_flags),
        }
    }
    return result


def check_structure(
    obj,
    *,
    descriptors: Sequence[str] | None = None,
    metals=None,
    method: str = "vesta",
    distance_scale: float = 1.0,
    clash_scale: float = 1.0,
    mode: str = DEFAULT_MODE,
) -> dict:
    """Run the selected diagnostics on one structure (graph built once).

    ``obj`` may be a pymatgen Structure, ASE Atoms, or a CIF path. Each
    descriptor is computed independently; a descriptor that errors is reported
    under ``errors`` rather than aborting the rest.

    ``distance_scale`` / ``clash_scale`` scale the bond-graph and atomic-overlap
    cutoffs respectively (see ``MOFChecker``); both default to ``1.0``.
    """
    profile = resolve_profile(mode)
    names = _descriptor_names(profile, descriptors)
    result = _result_metadata(
        profile, mode, obj, names, metals, method, distance_scale, clash_scale,
    )
    try:
        if mode == "0.9.6" and isinstance(obj, (str, Path)):
            checker = MOFChecker.from_cif(
                obj, mode=mode, metals=metals, method=method,
                distance_scale=distance_scale, clash_scale=clash_scale,
            )
            structure = checker.structure
        else:
            structure = normalize_structure(obj, corrected=profile.corrected_input)
            checker = MOFChecker(
                structure, mode=mode, metals=metals, method=method,
                distance_scale=distance_scale, clash_scale=clash_scale,
            )
    except Exception as exc:
        if not profile.corrected_input:
            raise
        return _failed_result(
            obj, mode=mode, message=f"{type(exc).__name__}: {exc}"[:200],
            descriptors=descriptors, metals=metals, method=method,
            distance_scale=distance_scale, clash_scale=clash_scale,
        )

    result["n_atoms"] = len(structure)
    if profile.corrected_input:
        result["input_sanity"] = True
    errors = {}
    for name in names:
        try:
            result[name] = getattr(checker, name)
        except Exception as exc:  # noqa: BLE001
            errors[name] = f"{type(exc).__name__}: {exc}"[:160]
    if errors:
        result["errors"] = errors
    if profile.corrected_input:
        result["check_results"] = _corrected_check_results(result, checker, profile)
    result.update(_profile_statuses(result, profile))
    validity = is_valid(result)
    result["composite_status"] = (
        "valid" if validity is True else "invalid" if validity is False else "indeterminate"
    )
    result["valid"] = validity  # compatibility alias for the named composite
    result["composites"] = {
        profile.composite_name: {
            "status": "pass" if validity is True else "fail" if validity is False else "indeterminate",
            "value": validity,
            "problem_flags": list(profile.problem_flags),
            "presence_flags": list(profile.presence_flags),
        }
    }
    return result


def _alarm(_signum, _frame):
    raise TimeoutError("structure check timed out")


def _worker(item, descriptors, metals, method, distance_scale, clash_scale, mode, on_error, timeout_s):
    index, obj = item
    profile = resolve_profile(mode)
    base = {
        "index": index,
        "id": _input_id(obj, index),
        "mode": mode,
        "resolved_profile": profile.id,
        "profile_provisional": profile.provisional,
    }
    old_handler = None
    try:
        if timeout_s:
            old_handler = signal.signal(signal.SIGALRM, _alarm)
            signal.setitimer(signal.ITIMER_REAL, float(timeout_s))
        return {**base, **check_structure(
            obj, descriptors=descriptors, metals=metals, method=method,
            distance_scale=distance_scale, clash_scale=clash_scale, mode=mode,
        )}
    except Exception as exc:  # noqa: BLE001
        if on_error == "raise":
            raise
        failed = _failed_result(
            obj, mode=mode, message=f"{type(exc).__name__}: {exc}"[:200],
            descriptors=descriptors, metals=metals, method=method,
            distance_scale=distance_scale, clash_scale=clash_scale,
        )
        return {**base, **failed}
    finally:
        if timeout_s:
            signal.setitimer(signal.ITIMER_REAL, 0.0)
            if old_handler is not None:
                signal.signal(signal.SIGALRM, old_handler)


def check_structures(
    inputs: Iterable,
    *,
    n_workers: int | None = None,
    descriptors: Sequence[str] | None = None,
    metals=None,
    method: str = "vesta",
    distance_scale: float = 1.0,
    clash_scale: float = 1.0,
    on_error: str = "record",
    chunksize: int = 1,
    progress: bool = False,
    timeout_s: float | None = None,
    mode: str = DEFAULT_MODE,
) -> list[dict]:
    """Validate many structures in parallel.

    Args:
        inputs: iterable of pymatgen Structures, ASE Atoms, and/or CIF paths.
        n_workers: process count (default: all CPUs). ``1`` runs serially.
        descriptors: subset of descriptor names (default: ``DEFAULT_DESCRIPTORS``,
            the diagnostics). Pass ``ALL_DESCRIPTORS`` to also get metadata,
            symmetry, and graph hashes.
        metals: metal symbol set (default: MOFChecker's METALS).
        distance_scale: multiplier on the bond-graph distance cutoffs (default
            ``1.0`` = MOFChecker behavior; >1 relaxes, recovering slightly-long
            bonds). Affects undercoordination, lone-molecule, and connectivity.
        clash_scale: multiplier on the atomic-overlap (clash) cutoffs (default
            ``1.0``). Affects ``has_atomic_overlaps``.
        on_error: ``"record"`` adds an ``error`` field to failed structures;
            ``"raise"`` propagates the first failure.
        progress: show a tqdm bar if tqdm is installed.
        timeout_s: optional per-structure wall-clock timeout; timed-out
            structures are recorded as errors when ``on_error='record'``.
        mode: ``"0.9.6"`` (default), ``"2.0"``, or provisional corrected ``"next"``.

    Returns:
        list of per-structure dicts (ordered to match ``inputs``), each with an
        ``index``, ``id``, ``n_atoms``, and the requested descriptors.
    """
    items = list(enumerate(inputs))
    if metals is not None:
        metals = frozenset(str(m) for m in metals)
    work = partial(
        _worker, descriptors=descriptors, metals=metals, method=method,
        distance_scale=distance_scale, clash_scale=clash_scale, on_error=on_error,
        timeout_s=timeout_s, mode=mode,
    )

    def _maybe_progress(iterator):
        if not progress:
            return iterator
        try:
            from tqdm import tqdm

            return tqdm(iterator, total=len(items))
        except ImportError:
            return iterator

    n_workers = n_workers or os.cpu_count() or 1
    if n_workers <= 1:
        results = list(_maybe_progress(map(work, items)))
    else:
        with Pool(n_workers) as pool:
            results = list(_maybe_progress(pool.imap_unordered(work, items, chunksize=chunksize)))

    results.sort(key=lambda r: r["index"])
    return results


def _cap_memory_from_env() -> None:
    gb = float(os.environ.get("MOFCHECKER_WORKER_MEM_GB", "8"))
    try:
        import resource

        nbytes = int(gb * 1024**3)
        soft, hard = resource.getrlimit(resource.RLIMIT_AS)
        cap = nbytes if hard == resource.RLIM_INFINITY else min(nbytes, hard)
        resource.setrlimit(resource.RLIMIT_AS, (cap, hard))
    except Exception:
        pass  # ponytail: best-effort OS guard; timeout still protects progress.


def _stream_worker_from_env() -> None:
    _cap_memory_from_env()
    descriptors = json.loads(os.environ.get("MOFCHECKER_DESCRIPTORS", "null"))
    mode = os.environ.get("MOFCHECKER_MODE", DEFAULT_MODE)
    print("READY", flush=True)
    for line in sys.stdin:
        line = line.strip()
        if not line:
            continue
        idx_s, raw_path = line.split("\t", 1)
        index = int(idx_s)
        path = Path(raw_path)
        try:
            r = check_structure(path, descriptors=descriptors, mode=mode)
        except Exception as exc:  # noqa: BLE001
            r = _failed_result(
                path, mode=mode, message=f"{type(exc).__name__}: {exc}"[:200],
                descriptors=descriptors,
            )
        r["index"] = index
        r["id"] = path.name
        print(json.dumps(r), flush=True)


def _run_stream(
    items: list[tuple[int, Path]], timeout_s: float,
    descriptors: Sequence[str] | None, mode: str,
) -> list[dict]:
    results: list[dict] = []
    i = 0

    def failed(index, path, message):
        return {
            "index": index, "id": path.name,
            **_failed_result(path, mode=mode, message=message, descriptors=descriptors),
        }

    env = os.environ | {
        "MOFCHECKER_DESCRIPTORS": json.dumps(list(descriptors) if descriptors is not None else None),
        "MOFCHECKER_MODE": mode,
    }
    while i < len(items):
        rest = items[i:]
        proc = subprocess.Popen(
            [sys.executable, "-m", "mofchecker_next.cli", "--_stream_worker"],
            stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.DEVNULL,
            text=True, env=env,
        )
        q: queue.Queue = queue.Queue()

        def _reader():
            assert proc.stdout is not None
            for line in proc.stdout:
                q.put(line)
            q.put(None)

        threading.Thread(target=_reader, daemon=True).start()
        assert proc.stdin is not None
        proc.stdin.write("".join(f"{idx}\t{p}\n" for idx, p in rest))
        proc.stdin.close()
        try:
            ready = q.get(timeout=60.0)
        except queue.Empty:
            ready = None
        if ready is None or ready.strip() != "READY":
            proc.kill()
            idx, path = rest[0]
            results.append(failed(idx, path, "worker failed to start"))
            i += 1
            continue
        for idx, path in rest:
            try:
                line = q.get(timeout=timeout_s)
            except queue.Empty:
                proc.kill()
                results.append(failed(idx, path, f"TimeoutExpired: >{timeout_s}s"))
                i += 1
                break
            if line is None:
                results.append(failed(idx, path, "worker died (child exited early)"))
                i += 1
                break
            try:
                results.append(json.loads(line))
            except Exception as exc:  # noqa: BLE001
                results.append(failed(idx, path, f"json parse failed: {type(exc).__name__}: {exc}"[:200]))
            i += 1
        else:
            proc.wait()
    return results


def check_cif_paths(
    paths: Iterable[str | Path],
    *,
    n_workers: int | None = None,
    descriptors: Sequence[str] | None = None,
    progress: bool = False,
    timeout_s: float | None = None,
    hard_timeout: bool = True,
    mode: str = DEFAULT_MODE,
) -> list[dict]:
    """Check CIF paths with the batch API; optional hard per-CIF timeout.

    ``check_structures(..., timeout_s=...)`` uses an in-process alarm. This path
    keeps one imported worker per shard but kills/respawns it if one CIF exceeds
    ``timeout_s``, so bad generated CIFs cannot stall a batch.
    """
    paths = [Path(p) for p in paths]
    if not hard_timeout or not timeout_s:
        return check_structures(
            paths, n_workers=n_workers, descriptors=descriptors, progress=progress,
            timeout_s=timeout_s, mode=mode,
        )
    n_workers = n_workers or os.cpu_count() or 1
    items = list(enumerate(paths))
    shards = [items[k::n_workers] for k in range(n_workers) if items[k::n_workers]]
    results_by_index = {}
    with cf.ThreadPoolExecutor(max_workers=n_workers) as ex:
        iterator = ex.map(lambda s: _run_stream(s, float(timeout_s), descriptors, mode), shards)
        for shard in iterator:
            for r in shard:
                results_by_index[r.get("index")] = r
    ordered = []
    for index, path in enumerate(paths):
        result = results_by_index.get(index)
        if result is None:
            result = {
                "index": index, "id": path.name,
                **_failed_result(path, mode=mode, message="missing result", descriptors=descriptors),
            }
        ordered.append(result)
    return ordered


def _profile_for_result(result: dict):
    return resolve_profile(result.get("mode", DEFAULT_MODE))


def _composite_fields(result: dict) -> tuple[tuple[str, ...], tuple[str, ...]] | None:
    profile = _profile_for_result(result)
    embedded = (
        tuple(result.get("composite_problem_flags", profile.problem_flags)),
        tuple(result.get("composite_presence_flags", profile.presence_flags)),
    )
    expected = (profile.problem_flags, profile.presence_flags)
    if result.get("resolved_profile", profile.id) != profile.id or embedded != expected:
        return None
    return embedded


def _composite_schema_error(result: dict) -> bool:
    return _composite_fields(result) is None


def _composite_error_fields(result: dict) -> set[str]:
    fields = _composite_fields(result)
    if fields is None:
        return set()
    problem_flags, presence_flags = fields
    required = set(problem_flags) | set(presence_flags)
    return required.intersection(result.get("errors", {}))


def is_valid(result: dict) -> bool | None:
    """Named-profile composite validity; missing/error fields are indeterminate."""
    fields = _composite_fields(result)
    if fields is None:
        return None
    problem_flags, presence_flags = fields
    required = set(problem_flags) | set(presence_flags)
    if result.get("error") or _composite_error_fields(result):
        return None
    if not required.issubset(result):
        return None
    if any(not isinstance(result[field], bool) for field in required):
        return None
    if any(result[flag] is True for flag in problem_flags):
        return False
    if any(result[flag] is not True for flag in presence_flags):
        return False
    return True


def summarize_results(results: Sequence[dict]) -> dict:
    n = len(results)
    valids = [is_valid(r) for r in results]
    actual_errors = [
        bool(r.get("error") or _composite_schema_error(r) or _composite_error_fields(r))
        for r in results
    ]
    n_errors = sum(actual_errors)
    n_indeterminate = sum(v is None for v in valids)
    n_valid = sum(v is True for v in valids)
    n_invalid = sum(v is False for v in valids)
    n_scored = n_valid + n_invalid
    descriptor_names = sorted({
        key for r in results for key, value in r.items() if isinstance(value, bool)
        and key not in {"valid", "profile_provisional"}
    })
    per_desc = {}
    for descriptor in descriptor_names:
        bools = [r[descriptor] for r in results if isinstance(r.get(descriptor), bool)]
        per_desc[descriptor] = sum(bools) / len(bools)

    numeric_descriptor_summary = {}
    for descriptor in (
        "max_abs_eqeq_charge", "eqeq_charge_sum", "eqeq_expected_total_charge",
        "eqeq_charge_residual", "eqeq_charge_threshold",
    ):
        values = [
            float(result[descriptor]) for result in results
            if isinstance(result.get(descriptor), (int, float))
            and not isinstance(result.get(descriptor), bool)
        ]
        if values:
            numeric_descriptor_summary[descriptor] = {
                "n": len(values),
                "min": min(values),
                "mean": sum(values) / len(values),
                "max": max(values),
            }

    descriptor_error_counts: dict[str, int] = {}
    error_categories: dict[str, int] = {}
    for result in results:
        messages = list(result.get("errors", {}).values())
        if result.get("error"):
            messages.append(result["error"])
        if _composite_schema_error(result):
            messages.append("CompositeSchemaError: result composite fields/profile do not match registry")
            descriptor_error_counts["__composite_schema__"] = descriptor_error_counts.get("__composite_schema__", 0) + 1
        for descriptor in result.get("errors", {}):
            descriptor_error_counts[descriptor] = descriptor_error_counts.get(descriptor, 0) + 1
        for message in messages:
            category = message.split(":", 1)[0]
            error_categories[category] = error_categories.get(category, 0) + 1

    modes = sorted({r.get("mode", DEFAULT_MODE) for r in results})
    profiles = sorted({
        r.get("resolved_profile", _profile_for_result(r).id) for r in results
    })
    provisional = sorted({
        r.get("profile_provisional", _profile_for_result(r).provisional) for r in results
    })
    status_names = sorted({
        key for result in results for key, value in result.items()
        if key.endswith("_status") and isinstance(value, str)
    })
    status_counts = {
        name: {
            status: sum(result.get(name) == status for result in results)
            for status in ("pass", "fail", "valid", "invalid", "indeterminate")
            if any(result.get(name) == status for result in results)
        }
        for name in status_names
    }
    error_records = {}
    for index, result in enumerate(results):
        messages = dict(result.get("errors", {}))
        if result.get("error"):
            messages["__structure__"] = result["error"]
        if _composite_schema_error(result):
            messages["__composite_schema__"] = "result composite fields/profile do not match registry"
        if messages:
            stable_index = str(result.get("index", index))
            error_records[stable_index] = {"id": result.get("id", stable_index), "errors": messages}
    summary = {
        "mode": modes[0] if len(modes) == 1 else modes,
        "resolved_profile": profiles[0] if len(profiles) == 1 else profiles,
        "profile_provisional": provisional[0] if len(provisional) == 1 else provisional,
        "n_structures": n,
        "n_requested": n,
        "n_errors": n_errors,
        "n_scored": n_scored,
        "n_valid": n_valid,
        "n_invalid": n_invalid,
        "n_indeterminate": n_indeterminate,
        "valid_rate": (n_valid / n_scored) if n_scored else 0.0,
        "valid_rate_incl_errors": (n_valid / n) if n else 0.0,
        "unconditional_valid_rate": (n_valid / n) if n else 0.0,
        "descriptor_true_rate": per_desc,
        "numeric_descriptor_summary": numeric_descriptor_summary,
        "descriptor_error_counts": descriptor_error_counts,
        "error_categories": error_categories,
        "n_descriptor_errors": sum(descriptor_error_counts.values()),
        "status_counts": status_counts,
        "errors": error_records,
    }
    if len(profiles) == 1:
        profile = resolve_profile(modes[0])
        summary["composite_name"] = profile.composite_name
        summary["composite_problem_flags"] = list(profile.problem_flags)
        summary["composite_presence_flags"] = list(profile.presence_flags)
        summary["manifest"] = profile_manifest(modes[0])
        protocols = {json.dumps(result.get("input_protocol", {}), sort_keys=True) for result in results}
        settings = {json.dumps(result.get("evaluation_settings", {}), sort_keys=True) for result in results}
        summary["input_protocols"] = [json.loads(value) for value in sorted(protocols)]
        summary["evaluation_settings"] = [json.loads(value) for value in sorted(settings)]
    return summary
