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
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
    "has_high_charges",
)


def _input_id(obj, index: int) -> str:
    if isinstance(obj, (str, Path)):
        return Path(obj).name
    return str(index)


def check_structure(obj, *, descriptors: Sequence[str] | None = None, metals=None, method: str = "vesta") -> dict:
    """Run the selected diagnostics on one structure (graph built once).

    ``obj`` may be a pymatgen Structure, ASE Atoms, or a CIF path. Each
    descriptor is computed independently; a descriptor that errors is reported
    under ``errors`` rather than aborting the rest.
    """
    names = list(descriptors) if descriptors is not None else list(DEFAULT_DESCRIPTORS)
    structure = normalize_structure(obj)
    checker = MOFChecker(structure, metals=metals, method=method)

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


def _worker(item, descriptors, metals, method, on_error):
    index, obj = item
    base = {"index": index, "id": _input_id(obj, index)}
    try:
        return {**base, **check_structure(obj, descriptors=descriptors, metals=metals, method=method)}
    except Exception as exc:  # noqa: BLE001
        if on_error == "raise":
            raise
        return {**base, "error": f"{type(exc).__name__}: {exc}"[:200]}


def check_structures(
    inputs: Iterable,
    *,
    n_workers: int | None = None,
    descriptors: Sequence[str] | None = None,
    metals=None,
    method: str = "vesta",
    on_error: str = "record",
    chunksize: int = 1,
    progress: bool = False,
) -> list[dict]:
    """Validate many structures in parallel.

    Args:
        inputs: iterable of pymatgen Structures, ASE Atoms, and/or CIF paths.
        n_workers: process count (default: all CPUs). ``1`` runs serially.
        descriptors: subset of descriptor names (default: ``DEFAULT_DESCRIPTORS``,
            the diagnostics). Pass ``ALL_DESCRIPTORS`` to also get metadata,
            symmetry, and graph hashes.
        metals: metal symbol set (default: MOFChecker's METALS).
        on_error: ``"record"`` adds an ``error`` field to failed structures;
            ``"raise"`` propagates the first failure.
        progress: show a tqdm bar if tqdm is installed.

    Returns:
        list of per-structure dicts (ordered to match ``inputs``), each with an
        ``index``, ``id``, ``n_atoms``, and the requested descriptors.
    """
    items = list(enumerate(inputs))
    if metals is not None:
        metals = frozenset(str(m) for m in metals)
    work = partial(_worker, descriptors=descriptors, metals=metals, method=method, on_error=on_error)

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
            results = list(_maybe_progress(pool.imap(work, items, chunksize=chunksize)))

    results.sort(key=lambda r: r["index"])
    return results
