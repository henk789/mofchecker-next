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

import os
import signal
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Sequence

from mofchecker_next.core import DEFAULT_DESCRIPTORS as ALL_DESCRIPTORS  # full set (metadata + hashes + diagnostics)
from mofchecker_next.core import MOFChecker, normalize_structure

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

# Composite validity = ADiT / Mofasa "Validity rate (all passed)" (Mofasa paper
# Appendix G, Table 4): a structure is valid iff all 3 presence checks are True
# and all 12 problem checks are False. NB this differs from the older native
# composite -- it adds has_hydrogen + has_geometrically_exposed_metal and
# DROPS has_3d_connected_graph.
PROBLEM_FLAGS = (
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_stray_atom", "has_lone_molecule", "has_suspicious_terminal_oxo", "has_high_charges",
    "has_geometrically_exposed_metal",
)
PRESENCE_FLAGS = ("has_carbon", "has_hydrogen", "has_metal")


def _input_id(obj, index: int) -> str:
    if isinstance(obj, (str, Path)):
        return Path(obj).name
    return str(index)


def check_structure(
    obj,
    *,
    descriptors: Sequence[str] | None = None,
    metals=None,
    method: str = "vesta",
    distance_scale: float = 1.0,
    clash_scale: float = 1.0,
) -> dict:
    """Run the selected diagnostics on one structure (graph built once).

    ``obj`` may be a pymatgen Structure, ASE Atoms, or a CIF path. Each
    descriptor is computed independently; a descriptor that errors is reported
    under ``errors`` rather than aborting the rest.

    ``distance_scale`` / ``clash_scale`` scale the bond-graph and atomic-overlap
    cutoffs respectively (see ``MOFChecker``); both default to ``1.0``.
    """
    names = list(descriptors) if descriptors is not None else list(DEFAULT_DESCRIPTORS)
    structure = normalize_structure(obj)
    checker = MOFChecker(
        structure, metals=metals, method=method,
        distance_scale=distance_scale, clash_scale=clash_scale,
    )

    result: dict = {"n_atoms": len(structure)}
    errors = {}
    for name in names:
        try:
            result[name] = getattr(checker, name)
        except Exception as exc:  # noqa: BLE001
            errors[name] = f"{type(exc).__name__}: {exc}"[:160]
    if errors:
        result["errors"] = errors
    return result


def _alarm(_signum, _frame):
    raise TimeoutError("structure check timed out")


def _worker(item, descriptors, metals, method, distance_scale, clash_scale, on_error, timeout_s):
    index, obj = item
    base = {"index": index, "id": _input_id(obj, index)}
    old_handler = None
    try:
        if timeout_s:
            old_handler = signal.signal(signal.SIGALRM, _alarm)
            signal.setitimer(signal.ITIMER_REAL, float(timeout_s))
        return {**base, **check_structure(
            obj, descriptors=descriptors, metals=metals, method=method,
            distance_scale=distance_scale, clash_scale=clash_scale,
        )}
    except Exception as exc:  # noqa: BLE001
        if on_error == "raise":
            raise
        return {**base, "error": f"{type(exc).__name__}: {exc}"[:200]}
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
        timeout_s=timeout_s,
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


def is_valid(result: dict) -> bool | None:
    """ADiT / Mofasa "all 15 criteria passed" validity (see PROBLEM/PRESENCE_FLAGS).

    Returns ``None`` when the structure errored before scoring. In ADiT/Mofasa
    terms an errored structure is invalid but still counts in the denominator --
    that is exactly ``valid_rate_incl_errors`` in ``summarize_results`` (``None``
    is not counted as valid, denominator is all structures). ``valid_rate``
    excludes errors.
    """
    if result.get("error"):
        return None
    if any(result.get(flag) is True for flag in PROBLEM_FLAGS):
        return False
    if any(result.get(flag) is not True for flag in PRESENCE_FLAGS):
        return False
    return True


def summarize_results(results: Sequence[dict]) -> dict:
    n = len(results)
    n_err = sum(1 for r in results if r.get("error"))
    valids = [is_valid(r) for r in results]
    n_valid = sum(v is True for v in valids)
    n_scored = sum(v is not None for v in valids)
    per_desc = {}
    for d in DEFAULT_DESCRIPTORS:
        bools = [r.get(d) for r in results if not r.get("error") and isinstance(r.get(d), bool)]
        if bools:
            per_desc[d] = sum(bools) / len(bools)
    return {
        "n_structures": n,
        "n_errors": n_err,
        "n_scored": n_scored,
        "n_valid": n_valid,
        "valid_rate": (n_valid / n_scored) if n_scored else 0.0,
        "valid_rate_incl_errors": (n_valid / n) if n else 0.0,
        "composite_problem_flags": list(PROBLEM_FLAGS),
        "composite_presence_flags": list(PRESENCE_FLAGS),
        "descriptor_true_rate": per_desc,
        "errors": {r.get("id", str(i)): r.get("error") for i, r in enumerate(results) if r.get("error")},
    }
