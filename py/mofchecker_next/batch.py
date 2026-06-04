"""Batch validation API.

Run the full MOFChecker-parity diagnostic set over many structures efficiently.

- Accepts pymatgen ``Structure``/``IStructure``, ASE ``Atoms``, or CIF path
  (``str``/``Path``) interchangeably -- or a mix.
- Builds each structure's graph **once** and reuses it across all graph-based
  checks (instead of rebuilding it per check).
- Parallelizes across structures with ``multiprocessing`` (the work is CPU-bound
  Python, so processes scale; structures are independent).
- Lets you pick a subset of descriptors (dropping ``has_oms`` alone removes the
  single most expensive check).

Example
-------
>>> from mofchecker_next.batch import check_structures
>>> results = check_structures(list_of_atoms_or_structures_or_paths, n_workers=16)
>>> bad = [r for r in results if r.get("has_atomic_overlaps")]

For a single structure use ``check_structure``.
"""

from __future__ import annotations

import os
from functools import partial
from multiprocessing import Pool
from pathlib import Path
from typing import Iterable, Sequence

# Radius constants (mirror MOFChecker 2.0's tables).
VDW_H_RADIUS = 1.1
COVALENT_MEDIAN = 1.49


def normalize_structure(obj):
    """Coerce a pymatgen Structure, ASE Atoms, or CIF path into a Structure."""
    from pymatgen.core import IStructure, Structure

    if isinstance(obj, (Structure, IStructure)):
        return obj
    if isinstance(obj, (str, Path)):
        return Structure.from_file(str(obj))
    # ASE Atoms (duck-typed to avoid importing ase eagerly).
    if hasattr(obj, "get_chemical_symbols") and hasattr(obj, "get_positions"):
        from pymatgen.io.ase import AseAtomsAdaptor

        return AseAtomsAdaptor.get_structure(obj)
    raise TypeError(
        f"Unsupported structure input {type(obj)!r}; expected pymatgen Structure, "
        "ASE Atoms, or a CIF path."
    )


def _input_id(obj, index: int) -> str:
    if isinstance(obj, (str, Path)):
        return Path(obj).name
    return str(index)


def _descriptor_table(structure, graph, metals):
    """Map descriptor name -> zero-arg callable computing it (graph reused)."""
    from mofchecker_next.checks import charge_oms as co
    from mofchecker_next.checks import composition as comp
    from mofchecker_next.checks import geometry as geo
    from mofchecker_next.checks import graph as g

    def atomic_overlaps():
        from pymatgen.core import Element

        atomic_numbers = [int(site.specie.Z) for site in structure]
        radii_by_z = {}
        for site in structure:
            sym = str(site.specie.symbol)
            radii_by_z[int(Element(sym).Z)] = co.COVALENT_RADII.get(sym, COVALENT_MEDIAN)
        matrix = geo.build_overlap_cutoff_matrix(atomic_numbers, radii_by_z, default_radius=COVALENT_MEDIAN)
        return len(geo.check_atomic_overlaps(structure, matrix)) > 0

    return {
        "has_carbon": lambda: comp.has_element(structure, "C"),
        "has_hydrogen": lambda: comp.has_element(structure, "H"),
        "has_nitrogen": lambda: comp.has_element(structure, "N"),
        "has_metal": lambda: comp.has_metal(structure, metals),
        "metal_number": lambda: comp.metal_number(structure, metals),
        "has_atomic_overlaps": atomic_overlaps,
        "has_overcoordinated_c": lambda: len(g.overcoordinated_carbon_indices_from_structure(structure, metals, graph=graph)) > 0,
        "has_overcoordinated_n": lambda: len(g.overcoordinated_nitrogen_indices_from_structure(structure, metals, graph=graph)) > 0,
        "has_overcoordinated_h": lambda: len(geo.overcoordinated_hydrogen_indices(structure, VDW_H_RADIUS)) > 0,
        "has_undercoordinated_c": lambda: len(g.undercoordinated_carbon_indices_from_structure(structure, metals, co.COVALENT_RADII, graph=graph)) > 0,
        "has_undercoordinated_n": lambda: len(g.undercoordinated_nitrogen_indices_from_structure(structure, metals, graph=graph)) > 0,
        "has_undercoordinated_rare_earth": lambda: len(g.undercoordinated_rare_earth_indices_from_structure(structure, graph=graph)) > 0,
        "has_undercoordinated_alkali_alkaline": lambda: len(g.undercoordinated_alkali_alkaline_indices_from_structure(structure, graph=graph)) > 0,
        "has_lone_molecule": lambda: len(g.floating_solvent_indices_from_structure(structure, graph=graph)) > 0,
        "has_3d_connected_graph": lambda: g.is_3d_connected_graph_from_structure(structure, graph=graph),
        "has_suspicious_terminal_oxo": lambda: len(g.false_oxo_indices_from_structure(structure, metals, graph=graph)) > 0,
        "has_geometrically_exposed_metal": lambda: len(g.geometrically_exposed_metal_indices_from_structure(structure, metals, graph=graph)) > 0,
        "possible_charged_fused_ring": lambda: len(co.fused_ring_indices(structure, graph)) > 0,
        "positive_charge_from_linkers": lambda: len(co.positive_charge_indices(structure, graph)),
        "negative_charge_from_linkers": lambda: len(co.negative_charge_indices(structure, graph)),
        "has_oms": lambda: len(co.oms_indices(structure, graph)) > 0,
        "oms_indices": lambda: sorted(co.oms_indices(structure, graph)),
        "has_high_charges": _has_high_charges_factory(structure),
    }


def _has_high_charges_factory(structure):
    def fn():
        from mofchecker_next.eqeq import has_high_charges

        return has_high_charges(structure)

    return fn


# Default descriptor set. ``oms_indices`` (the explicit index list) and
# ``has_high_charges`` (EQeq; GPL subpackage) are opt-in extras.
DEFAULT_DESCRIPTORS = (
    "has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
    "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
    "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
    "has_lone_molecule", "has_3d_connected_graph", "has_suspicious_terminal_oxo",
    "has_geometrically_exposed_metal", "possible_charged_fused_ring",
    "positive_charge_from_linkers", "negative_charge_from_linkers", "has_oms",
)

ALL_DESCRIPTORS = DEFAULT_DESCRIPTORS + ("oms_indices", "has_high_charges")

# Descriptors that do not require the StructureGraph (so it is only built when
# at least one graph-based descriptor is requested).
_NON_GRAPH_DESCRIPTORS = {
    "has_carbon", "has_hydrogen", "has_nitrogen", "has_metal", "metal_number",
    "has_atomic_overlaps", "has_overcoordinated_h", "has_high_charges",
}


def check_structure(obj, *, descriptors: Sequence[str] | None = None, metals=None, method: str = "vesta") -> dict:
    """Run the selected diagnostics on one structure (graph built once).

    ``obj`` may be a pymatgen Structure, ASE Atoms, or a CIF path. Each
    descriptor is computed independently; a descriptor that errors is reported
    under ``errors`` rather than aborting the rest.
    """
    from mofchecker_next.checks.charge_oms import METALS
    from mofchecker_next.checks.graph import build_structure_graph

    metals = METALS if metals is None else metals
    names = list(descriptors) if descriptors is not None else list(DEFAULT_DESCRIPTORS)

    structure = normalize_structure(obj)
    result: dict = {"n_atoms": len(structure)}
    graph = build_structure_graph(structure, method) if (set(names) - _NON_GRAPH_DESCRIPTORS) else None
    table = _descriptor_table(structure, graph, metals)

    errors = {}
    for name in names:
        fn = table.get(name)
        if fn is None:
            errors[name] = "unknown descriptor"
            continue
        try:
            result[name] = fn()
        except Exception as exc:  # noqa: BLE001
            errors[name] = f"{type(exc).__name__}: {exc}"[:160]
    if errors:
        result["errors"] = errors
    return result


def _worker(item, descriptors, metals, method, on_error):
    index, obj = item
    base = {"index": index, "id": _input_id(obj, index)}
    try:
        res = check_structure(obj, descriptors=descriptors, metals=metals, method=method)
        return {**base, **res}
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
        descriptors: subset of descriptor names (default: ``DEFAULT_DESCRIPTORS``).
            Use ``ALL_DESCRIPTORS`` to also get ``oms_indices`` and the EQeq
            ``has_high_charges``.
        metals: metal symbol set (default: MOFChecker's METALS).
        on_error: ``"record"`` puts an ``error`` field on failed structures;
            ``"raise"`` propagates the first failure.
        progress: show a tqdm bar if tqdm is installed.

    Returns:
        list of per-structure dicts (ordered to match ``inputs``), each with an
        ``index``, an ``id``, ``n_atoms``, and the requested descriptors.
    """
    items = list(enumerate(inputs))
    if metals is not None:
        metals = frozenset(str(m) for m in metals)
    work = partial(_worker, descriptors=descriptors, metals=metals, method=method, on_error=on_error)

    n_workers = n_workers or os.cpu_count() or 1
    if n_workers <= 1:
        iterator = map(work, items)
    else:
        pool = Pool(n_workers)
        iterator = pool.imap(work, items, chunksize=chunksize)

    if progress:
        try:
            from tqdm import tqdm

            iterator = tqdm(iterator, total=len(items))
        except ImportError:
            pass

    try:
        results = list(iterator)
    finally:
        if n_workers > 1:
            pool.close()
            pool.join()
    results.sort(key=lambda r: r["index"])
    return results
