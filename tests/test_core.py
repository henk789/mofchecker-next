"""Tests for the drop-in MOFChecker class."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("structuregraph_helpers")
from pymatgen.core import Lattice, Structure  # noqa: E402

from mofchecker_next import MOFChecker  # noqa: E402
from mofchecker_next.core import DEFAULT_DESCRIPTORS  # noqa: E402


def _mof_like():
    lattice = Lattice.cubic(14.0)
    species = ["Zn", "O", "O", "C", "C", "H"]
    coords = [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1], [0.1, 0.2, 0.1],
              [0.25, 0.18, 0.1], [0.35, 0.25, 0.1], [0.4, 0.3, 0.1]]
    return Structure(lattice, species, coords)


def test_constructors(tmp_path):
    s = _mof_like()
    assert isinstance(MOFChecker(s).structure, Structure)
    p = tmp_path / "m.cif"
    s.to(filename=str(p), fmt="cif")
    c = MOFChecker.from_cif(p)
    assert c.name == "m" and c.path is not None and len(c.structure) == 6


def test_core_diagnostics_and_metadata():
    c = MOFChecker(_mof_like())
    assert c.has_carbon is True and c.has_metal is True
    assert c.metal_number == 1
    assert isinstance(c.has_atomic_overlaps, bool)
    assert isinstance(c.get_overlapping_indices(), list)
    assert isinstance(c.formula, str)
    assert c.spacegroup_number >= 1
    # all four graph hashes are strings
    for h in (c.graph_hash, c.undecorated_graph_hash, c.scaffold_hash, c.undecorated_scaffold_hash):
        assert isinstance(h, str) and h


def test_symmetry_hash_is_deterministic():
    s = _mof_like()
    assert MOFChecker(s).symmetry_hash == MOFChecker(s).symmetry_hash


def test_graph_built_once():
    c = MOFChecker(_mof_like())
    g1 = c.graph
    _ = c.has_oms, c.possible_charged_fused_ring, c.has_lone_molecule
    assert c.graph is g1  # cached, reused across checks


def test_floating_solvent_split_stray_vs_lone_molecule():
    c = MOFChecker(_mof_like())
    c.__dict__["floating_solvent_indices"] = [[1], [2, 3]]
    assert c.stray_atom_indices == [[1]]
    assert c.lone_molecule_indices == [[2, 3]]
    assert c.has_stray_atom is True
    assert c.has_lone_molecule is True


def test_get_mof_descriptors_default_and_subset():
    c = MOFChecker(_mof_like())
    full = c.get_mof_descriptors()
    assert set(full) == set(DEFAULT_DESCRIPTORS)
    assert full["has_metal"] is True
    sub = c.get_mof_descriptors(["has_carbon", "metal_number"])
    assert list(sub) == ["has_carbon", "metal_number"]


def test_healing_not_implemented():
    c = MOFChecker(_mof_like())
    with pytest.raises(NotImplementedError):
        _ = c.adding_hydrogen
    with pytest.raises(NotImplementedError):
        _ = c.adding_linker


def test_is_porous_none_without_zeopp():
    assert MOFChecker(_mof_like()).is_porous is None
