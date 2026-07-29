from types import SimpleNamespace

import pytest

pytest.importorskip("pymatgen")
from pymatgen.core import Lattice, Structure  # noqa: E402

from mofchecker_next import MOFChecker  # noqa: E402
from mofchecker_next.checks.graph import (  # noqa: E402
    undercoordinated_carbon_indices_from_structure,
    undercoordinated_carbon_indices_v096_from_structure,
    undercoordinated_nitrogen_indices_from_structure,
    undercoordinated_nitrogen_indices_v096_from_structure,
)
from mofchecker_next.core import LEGACY_DEFAULT_DESCRIPTORS  # noqa: E402


def _graph(structure, adjacency):
    return SimpleNamespace(
        get_connected_sites=lambda index: [
            SimpleNamespace(index=neighbor, site=structure[neighbor])
            for neighbor in adjacency.get(index, [])
        ]
    )


def test_v096_carbon_cn1_does_not_exempt_short_nitrile():
    structure = Structure(
        Lattice.cubic(20), ["C", "N"], [[10, 10, 10], [11.1, 10, 10]], coords_are_cartesian=True
    )
    graph = _graph(structure, {0: [1], 1: [0]})
    assert undercoordinated_carbon_indices_v096_from_structure(structure, graph=graph) == [0]
    assert undercoordinated_carbon_indices_from_structure(structure, set(), {}, graph=graph) == []


def test_v096_nitrogen_cn3_branch_is_active():
    structure = Structure(
        Lattice.cubic(20),
        ["N", "H", "H", "Fe"],
        [[10, 10, 10], [11, 10, 10], [9.5, 10.866, 10], [10, 10, 12]],
        coords_are_cartesian=True,
    )
    graph = _graph(structure, {0: [1, 2, 3], 1: [0], 2: [0], 3: [0]})
    assert undercoordinated_nitrogen_indices_v096_from_structure(structure, graph=graph) == [0]
    assert undercoordinated_nitrogen_indices_from_structure(structure, {"Fe"}, graph=graph) == []


def test_v096_mode_uses_legacy_floating_union_and_descriptor_names():
    checker = MOFChecker.__new__(MOFChecker)
    checker.mode = "0.9.6"
    checker.__dict__["floating_solvent_indices"] = [[1], [2, 3]]
    assert checker.stray_atom_indices == []
    assert checker.lone_molecule_indices == [[1], [2, 3]]
    assert "decorated_scaffold_hash" in LEGACY_DEFAULT_DESCRIPTORS
    assert "has_suspicicious_terminal_oxo" in LEGACY_DEFAULT_DESCRIPTORS


def test_v096_mode_uses_charge_threshold_three_and_reference_labels(monkeypatch):
    structure = Structure(Lattice.cubic(10), ["C"], [[0, 0, 0]])
    seen = []
    monkeypatch.setattr(
        "mofchecker_next.eqeq.has_high_charges",
        lambda _structure, threshold=4.0, reference_cif_labels=False: seen.append(
            (threshold, reference_cif_labels)
        )
        or False,
    )
    assert MOFChecker(structure, mode="0.9.6", symprec=None, angle_tolerance=None, primitive=False).has_high_charges is False
    assert MOFChecker(structure).has_high_charges is False
    assert MOFChecker(structure, mode="2.0").has_high_charges is False
    assert seen == [(3.0, True), (3.0, True), (4.0, False)]


def test_unknown_mode_rejected():
    structure = Structure(Lattice.cubic(10), ["C"], [[0, 0, 0]])
    with pytest.raises(ValueError):
        MOFChecker(structure, mode="legacy")


def test_reference_cif_labels_reproduce_eqeq_label_parsing():
    from mofchecker_next.eqeq.parameters import parameters_for

    # Two-letter symbols survive EQeq's 2-char label read.
    assert parameters_for("Ag", reference_cif_labels=True) == parameters_for("Ag")
    # One-letter symbols and elements outside the table share the Z=0 fallback,
    # which is hydrogen's ionization row without hydrogen's hI0 special case.
    fallback = parameters_for("C", reference_cif_labels=True)
    assert fallback == parameters_for("U", reference_cif_labels=True)
    assert fallback != parameters_for("C")
    assert fallback != parameters_for("H")


def test_v096_rare_earth_excludes_scandium_and_yttrium():
    from mofchecker_next.checks.graph import undercoordinated_rare_earth_indices_from_structure

    # Far-apart isolated atoms: both have CN 0, so only the element set matters.
    structure = Structure(Lattice.cubic(18.0), ["Y", "La"], [[0, 0, 0], [0.5, 0.5, 0.5]])
    assert undercoordinated_rare_earth_indices_from_structure(structure) == [0, 1]
    assert undercoordinated_rare_earth_indices_from_structure(structure, include_sc_y=False) == [1]
