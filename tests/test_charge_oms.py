"""Tests for linker-charge, fused-ring, and OMS diagnostics.

These assert our values on synthetic structures. The values were verified to
match the MOFChecker 2.0 oracle via ``scripts/compare_charge_oms_parity.py``
(full reference parity, 16/16 on the local reference CIFs plus these synthetic
positive-direction triggers).
"""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("structuregraph_helpers")
from pymatgen.core import Lattice, Structure  # noqa: E402

from mofchecker_next.checks.charge_oms import (  # noqa: E402
    _clean_cycles,
    negative_charge_from_linkers_from_structure,
    oms_indices_from_structure,
    positive_charge_from_linkers_from_structure,
    possible_charged_fused_ring_from_structure,
)


def _box(symbols, coords_angstrom, a=25.0):
    lattice = Lattice.cubic(a)
    frac = [lattice.get_fractional_coords([c[0] + 12, c[1] + 12, c[2] + 12]) for c in coords_angstrom]
    return Structure(lattice, symbols, frac)


def test_positive_charge_ge_and_sb():
    # Sb contributes 3 and Ge contributes 4 to the positive-charge count.
    structure = _box(["Ge", "Sb"], [[0, 0, 0], [8, 8, 8]])
    assert positive_charge_from_linkers_from_structure(structure) == 7


def test_positive_charge_quaternary_nitrogen():
    d = 1.47
    tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3) * d
    structure = _box(["N"] + ["C"] * 4, [[0, 0, 0]] + [list(v) for v in tet])
    assert positive_charge_from_linkers_from_structure(structure) == 1


def test_positive_charge_trisubstituted_oxygen():
    tri = np.array([[1, 0, 0], [-0.5, 0.866, 0], [-0.5, -0.866, 0]]) * 1.43
    structure = _box(["O"] + ["C"] * 3, [[0, 0, 0]] + [list(v) for v in tri])
    assert positive_charge_from_linkers_from_structure(structure) == 1


def test_negative_charge_isolated_oxygen():
    # An oxygen with no neighbors contributes 2 to the negative-charge count.
    structure = _box(["O", "O"], [[0, 0, 0], [10, 10, 10]])
    assert negative_charge_from_linkers_from_structure(structure) == 4


def _benzimidazole():
    atoms = {
        "C3a": (1.390, 0.000), "C4": (0.695, 1.204), "C5": (-0.695, 1.204),
        "C6": (-1.390, 0.000), "C7": (-0.695, -1.204), "C7a": (0.695, -1.204),
        "N3": (2.748, -0.289), "C2": (2.894, -1.671), "N1": (1.624, -2.236),
    }
    hydrogens = {
        "HN1": (1.416, -3.214), "HC2": (3.829, -2.211), "HC4": (1.235, 2.139),
        "HC5": (-1.235, 2.139), "HC6": (-2.470, 0.000), "HC7": (-1.235, -2.139),
    }
    symbols = ["N" if k.startswith("N") else "C" for k in atoms]
    coords = [[x, y, 0.0] for (x, y) in atoms.values()]
    symbols += ["H"] * len(hydrogens)
    coords += [[x, y, 0.0] for (x, y) in hydrogens.values()]
    return _box(symbols, coords)


def test_fused_ring_benzimidazole_flagged():
    assert possible_charged_fused_ring_from_structure(_benzimidazole()) is True


def test_fused_ring_negative_on_simple_molecule():
    # A lone quaternary-ammonium has no fused N-heterocycle.
    d = 1.47
    tet = np.array([[1, 1, 1], [1, -1, -1], [-1, 1, -1], [-1, -1, 1]]) / np.sqrt(3) * d
    structure = _box(["N"] + ["C"] * 4, [[0, 0, 0]] + [list(v) for v in tet])
    assert possible_charged_fused_ring_from_structure(structure) is False


def test_rust_cycles_match_networkx_reference():
    import networkx as nx
    from mofchecker_next.checks.graph import build_structure_graph
    from structuregraph_helpers.create import construct_clean_graph

    def canon(cycle):
        m = min(range(len(cycle)), key=cycle.__getitem__)
        rotated = cycle[m:] + cycle[:m]
        reversed_ = [rotated[0], *reversed(rotated[1:])]
        return tuple(min(rotated, reversed_))

    graph = build_structure_graph(_benzimidazole())
    reference = {canon(c) for c in nx.simple_cycles(construct_clean_graph(graph), length_bound=16)}
    assert {canon(c) for c in _clean_cycles(graph)} == reference


def test_oms_no_metal_returns_empty():
    # No metal -> no open metal site (we return empty rather than raising).
    structure = _box(["O", "O"], [[0, 0, 0], [10, 10, 10]])
    assert oms_indices_from_structure(structure) == []
