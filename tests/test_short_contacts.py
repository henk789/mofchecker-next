import random

import pytest
from pymatgen.core import Lattice, Structure

from mofchecker_next._rust import find_neighbor_candidates, find_short_contacts
from mofchecker_next.checks.geometry import build_overlap_cutoff_matrix, check_atomic_overlaps, find_short_contacts_py
from mofchecker_next.checks.geometry import find_neighbor_candidates_py
from mofchecker_next.checks.geometry import overcoordinated_hydrogen_indices


def assert_contacts_match(left, right):
    assert len(left) == len(right)
    for left_contact, right_contact in zip(left, right):
        assert left_contact["i"] == right_contact["i"]
        assert left_contact["j"] == right_contact["j"]
        assert left_contact["image_j"] == right_contact["image_j"]
        assert abs(left_contact["distance"] - right_contact["distance"]) <= 1e-10
        assert abs(left_contact["cutoff"] - right_contact["cutoff"]) <= 1e-12


def assert_neighbor_candidates_match(left, right):
    assert len(left) == len(right)
    for left_candidate, right_candidate in zip(left, right):
        assert left_candidate["i"] == right_candidate["i"]
        assert left_candidate["j"] == right_candidate["j"]
        assert left_candidate["image_j"] == right_candidate["image_j"]
        assert abs(left_candidate["distance"] - right_candidate["distance"]) <= 1e-10


def test_python_short_contacts_expected_values():
    frac_coords = [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.5, 0.5, 0.5]]
    atomic_numbers = [6, 6, 8]
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    cutoff_matrix = [[0.0 for _ in range(9)] for _ in range(9)]
    cutoff_matrix[6][6] = 0.76
    cutoff_matrix[6][8] = 0.66
    cutoff_matrix[8][6] = 0.66
    contacts = find_short_contacts_py(frac_coords, atomic_numbers, lattice, cutoff_matrix)
    assert contacts == [
        {"i": 0, "j": 1, "image_j": [0, 0, 0], "distance": 0.5, "cutoff": 0.76}
    ]


def test_overlap_cutoff_matrix_uses_minimum_radii():
    matrix = build_overlap_cutoff_matrix([1, 6, 8], {1: 0.31, 6: 0.76, 8: 0.66})
    assert matrix[6][6] == 0.76
    assert matrix[6][8] == 0.66
    assert matrix[8][6] == 0.66
    assert matrix[1][8] == 0.31


def test_overlap_cutoff_matrix_default_radius():
    matrix = build_overlap_cutoff_matrix([1, 118], {1: 0.31}, default_radius=0.75)
    assert matrix[118][118] == 0.75
    assert matrix[1][118] == 0.31


def test_overlap_cutoff_matrix_missing_radius_errors():
    try:
        build_overlap_cutoff_matrix([1, 118], {1: 0.31})
    except KeyError as exc:
        assert "118" in str(exc)
    else:
        raise AssertionError("Expected missing radius to raise KeyError")


def test_rust_short_contacts_match_python():
    frac_coords = [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.99, 0.0, 0.0]]
    atomic_numbers = [1, 1, 1]
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    cutoff_matrix = [[0.0, 0.0], [0.0, 0.7]]
    python = find_short_contacts_py(frac_coords, atomic_numbers, lattice, cutoff_matrix)
    rust = find_short_contacts(frac_coords, atomic_numbers, lattice, cutoff_matrix, 1.0)
    assert_contacts_match(rust, python)


def test_randomized_rust_vs_python_short_contacts():
    rng = random.Random(8675309)
    lattice = [[10.0, 0.0, 0.0], [0.5, 9.0, 0.0], [0.25, 0.5, 8.0]]
    cutoff_matrix = [[0.0 for _ in range(9)] for _ in range(9)]
    for zi in (1, 6, 8):
        for zj in (1, 6, 8):
            cutoff_matrix[zi][zj] = 1.0

    for _ in range(50):
        n = 8
        frac_coords = [[rng.uniform(-0.25, 1.25) for _ in range(3)] for _ in range(n)]
        atomic_numbers = [rng.choice([1, 6, 8]) for _ in range(n)]
        python = find_short_contacts_py(frac_coords, atomic_numbers, lattice, cutoff_matrix)
        rust = find_short_contacts(frac_coords, atomic_numbers, lattice, cutoff_matrix, 1.0)
        assert_contacts_match(rust, python)


def test_neighbor_candidates_match_python_expected_values():
    frac_coords = [[0.0, 0.0, 0.0], [0.04, 0.0, 0.0], [0.5, 0.5, 0.5]]
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    python = find_neighbor_candidates_py(frac_coords, lattice, 0.5)
    rust = find_neighbor_candidates(frac_coords, lattice, 0.5)
    assert_neighbor_candidates_match(rust, python)
    assert rust == [{"i": 0, "j": 1, "image_j": [0, 0, 0], "distance": 0.4}]


def test_randomized_neighbor_candidates_match_python():
    rng = random.Random(10101)
    lattice = [[9.0, 0.0, 0.0], [0.7, 8.0, 0.0], [0.4, 0.2, 7.0]]
    for _ in range(50):
        frac_coords = [[rng.uniform(-0.5, 1.5) for _ in range(3)] for _ in range(10)]
        cutoff = rng.uniform(0.5, 2.0)
        python = find_neighbor_candidates_py(frac_coords, lattice, cutoff)
        rust = find_neighbor_candidates(frac_coords, lattice, cutoff)
        assert_neighbor_candidates_match(rust, python)


def test_overcoordinated_hydrogen_indices_use_pymatgen_neighbors():
    structure = Structure(
        Lattice.cubic(10.0),
        ["H", "C", "C", "H"],
        [[0.0, 0.0, 0.0], [0.05, 0.0, 0.0], [0.0, 0.05, 0.0], [0.5, 0.5, 0.5]],
    )
    assert overcoordinated_hydrogen_indices(structure, 1.1) == [0]


def test_atomic_overlap_periodic_boundary_diagnostic_image():
    structure = Structure(Lattice.cubic(10.0), ["C", "C"], [[0.99, 0.0, 0.0], [0.01, 0.0, 0.0]])
    cutoff_matrix = build_overlap_cutoff_matrix([6, 6], {6: 0.76})
    diagnostics = check_atomic_overlaps(structure, cutoff_matrix)
    assert len(diagnostics) == 1
    assert diagnostics[0].atoms[1].image == (1, 0, 0)
    assert diagnostics[0].values["distance_angstrom"] == pytest.approx(0.2)


def test_atomic_overlap_skew_cell_diagnostic():
    structure = Structure(
        Lattice([[1.0, 0.0, 0.0], [0.99, 0.1, 0.0], [0.0, 0.0, 10.0]]),
        ["C", "C"],
        [[0.0, 0.0, 0.0], [0.49, 0.49, 0.0]],
    )
    cutoff_matrix = build_overlap_cutoff_matrix([6, 6], {6: 0.76})
    diagnostics = check_atomic_overlaps(structure, cutoff_matrix)
    assert len(diagnostics) == 1
    assert diagnostics[0].atoms[1].image == (0, -1, 0)


def test_atomic_overlap_hydrogen_and_metal_nonmetal_cases():
    h_structure = Structure(Lattice.cubic(10.0), ["H", "H"], [[0, 0, 0], [0.01, 0, 0]])
    h_cutoffs = build_overlap_cutoff_matrix([1, 1], {1: 0.31})
    assert len(check_atomic_overlaps(h_structure, h_cutoffs)) == 1

    metal_structure = Structure(Lattice.cubic(10.0), ["Na", "C"], [[0, 0, 0], [0.01, 0, 0]])
    metal_cutoffs = build_overlap_cutoff_matrix([11, 6], {11: 1.66, 6: 0.76})
    assert len(check_atomic_overlaps(metal_structure, metal_cutoffs)) == 1


def test_atomic_overlap_no_overlap_negative_case():
    structure = Structure(Lattice.cubic(10.0), ["C", "C"], [[0, 0, 0], [0.5, 0, 0]])
    cutoff_matrix = build_overlap_cutoff_matrix([6, 6], {6: 0.76})
    assert check_atomic_overlaps(structure, cutoff_matrix) == []
