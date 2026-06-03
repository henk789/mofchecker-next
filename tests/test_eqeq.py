"""Tests for the EQeq charge-equilibration subpackage.

The subpackage is a faithful translation of EQeq's C++ (GPLv2). These tests pin
the parameter derivation and the kernel's basic invariants; full bit-exact charge
parity against the pyeqeq oracle is exercised by
``scripts/compare_high_charges_parity.py`` (requires the .venv-ref oracle).
"""

from __future__ import annotations

import math

import numpy as np
import pytest

from mofchecker_next.eqeq import compute_charges, has_high_charges, high_charge_indices
from mofchecker_next.eqeq.parameters import parameters_for

pytest.importorskip("pymatgen")
from pymatgen.core import Lattice, Structure  # noqa: E402


def test_hydrogen_parameters_match_eqeq_formula():
    # X = 0.5 * (hI1 + hI0) = 0.5 * (13.598 + (-2.0)); J = hI1 - hI0.
    params = parameters_for("H")
    assert params.electronegativity == pytest.approx(0.5 * (13.598 - 2.0))
    assert params.hardness == pytest.approx(13.598 - (-2.0))


def test_carbon_parameters_neutral_charge_center():
    # Carbon charge center is 0: X = 0.5*(IE1 + EA), J = IE1 - EA.
    # EA(C)=1.26212, IE1(C)=11.26000 from the vendored ionization table.
    params = parameters_for("C")
    assert params.electronegativity == pytest.approx(0.5 * (11.26000 + 1.26212), abs=1e-5)
    assert params.hardness == pytest.approx(11.26000 - 1.26212, abs=1e-5)


def test_lanthanum_uses_charge_center_two():
    # La's charge center in EQeq's table is 2 (not 3): X uses IP[3],IP[2] and is
    # shifted by -cc*J. This is the value that makes AFIPAH match the reference.
    params = parameters_for("La")
    assert params.hardness > 0
    # Sanity: charge-center shift makes La strongly electropositive (negative X).
    assert params.electronegativity < 0


def _two_atom_cell(symbols, frac):
    lattice = Lattice.from_parameters(14.0, 14.0, 14.0, 90.0, 90.0, 90.0)
    return Structure(lattice, symbols, frac)


def test_charges_sum_to_zero():
    structure = _two_atom_cell(["C", "O"], [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1]])
    charges = compute_charges(structure)
    assert charges.sum() == pytest.approx(0.0, abs=1e-9)


def test_more_electronegative_atom_is_more_negative():
    structure = _two_atom_cell(["C", "O"], [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1]])
    charges = compute_charges(structure)
    # Oxygen (index 1) is more electronegative than carbon -> more negative.
    assert charges[1] < charges[0]


def test_charges_rounded_to_precision():
    structure = _two_atom_cell(["C", "O"], [[0.1, 0.1, 0.1], [0.25, 0.1, 0.1]])
    charges = compute_charges(structure, charge_precision=3)
    for q in charges:
        assert math.isclose(q, round(q, 3), abs_tol=1e-12)


def test_high_charge_helpers_agree():
    structure = _two_atom_cell(["C", "O"], [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1]])
    indices = high_charge_indices(structure, threshold=4.0)
    assert has_high_charges(structure, threshold=4.0) == bool(indices)
    # Physical charges here are small; nothing should exceed the threshold.
    assert indices == []


def test_single_atom_is_neutral():
    structure = Structure(
        Lattice.from_parameters(12.0, 12.0, 12.0, 90.0, 90.0, 90.0),
        ["Na"],
        [[0.0, 0.0, 0.0]],
    )
    charges = compute_charges(structure)
    assert charges[0] == pytest.approx(0.0, abs=1e-9)
