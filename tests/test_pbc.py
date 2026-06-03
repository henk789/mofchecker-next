import random

from pymatgen.core import Lattice

from mofchecker_next._rust import minimum_image_distance
from mofchecker_next.checks.geometry import minimum_image_distance_py


def assert_close_sequence(left, right, tolerance):
    assert len(left) == len(right)
    for left_item, right_item in zip(left, right):
        assert abs(left_item - right_item) <= tolerance


def assert_mic_matches_python(frac_i, frac_j, lattice):
    rust = minimum_image_distance(frac_i, frac_j, lattice)
    python = minimum_image_distance_py(frac_i, frac_j, lattice)
    assert_close_sequence(rust["delta_frac"], python["delta_frac"], 1e-12)
    assert rust["image"] == python["image"]
    assert abs(rust["distance"] - python["distance"]) <= 1e-10


def test_python_mic_expected_values():
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    result = minimum_image_distance_py([0.0, 0.0, 0.0], [0.25, 0.0, 0.0], lattice)
    assert result["delta_frac"] == [0.25, 0.0, 0.0]
    assert result["image"] == [0, 0, 0]
    assert result["distance"] == 2.5


def test_python_mic_tie_behavior_at_half():
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    result = minimum_image_distance_py([0.0, 0.0, 0.0], [0.5, -0.5, 0.0], lattice)
    assert result["image"] == [0, 0, 0]
    assert result["delta_frac"] == [0.5, -0.5, 0.0]


def test_python_mic_skew_lattice_avoids_component_wrapping_pitfall():
    lattice = [[1.0, 0.0, 0.0], [0.99, 0.1, 0.0], [0.0, 0.0, 10.0]]
    result = minimum_image_distance_py([0.0, 0.0, 0.0], [0.49, 0.49, 0.0], lattice)
    assert result["image"] == [0, -1, 0]
    assert_close_sequence(result["delta_frac"], [0.49, -0.51, 0.0], 1e-12)
    assert result["distance"] < 0.06


def test_rust_mic_matches_python_cubic_boundary():
    lattice = [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]]
    assert_mic_matches_python([0.95, 0.2, 0.3], [0.05, 0.9, -0.4], lattice)


def test_rust_mic_matches_python_non_orthogonal():
    lattice = [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.25, 5.0]]
    assert_mic_matches_python([0.1, -0.2, 1.1], [0.8, 0.6, -0.3], lattice)


def test_randomized_rust_vs_python_mic():
    rng = random.Random(1729)
    for _ in range(200):
        frac_i = [rng.uniform(-2.0, 2.0) for _ in range(3)]
        frac_j = [rng.uniform(-2.0, 2.0) for _ in range(3)]
        lattice = [
            [rng.uniform(3.0, 8.0), 0.0, 0.0],
            [rng.uniform(-1.0, 1.0), rng.uniform(3.0, 8.0), 0.0],
            [rng.uniform(-1.0, 1.0), rng.uniform(-1.0, 1.0), rng.uniform(3.0, 8.0)],
        ]
        assert_mic_matches_python(frac_i, frac_j, lattice)


def assert_mic_distance_matches_pymatgen(frac_i, frac_j, lattice):
    rust = minimum_image_distance(frac_i, frac_j, lattice)
    pymatgen_distance, _ = Lattice(lattice).get_distance_and_image(frac_i, frac_j)
    assert abs(rust["distance"] - float(pymatgen_distance)) <= 1e-10


def test_rust_mic_matches_pymatgen_distance_and_image_for_representative_cases():
    cases = [
        (
            [0.95, 0.0, 0.0],
            [0.05, 0.0, 0.0],
            [[10.0, 0.0, 0.0], [0.0, 10.0, 0.0], [0.0, 0.0, 10.0]],
            [1, 0, 0],
        ),
        (
            [0.0, 0.0, 0.0],
            [0.49, 0.49, 0.0],
            [[1.0, 0.0, 0.0], [0.99, 0.1, 0.0], [0.0, 0.0, 10.0]],
            [0, -1, 0],
        ),
        (
            [-0.25, 1.2, 0.1],
            [1.3, -0.4, 0.9],
            [[4.0, 0.0, 0.0], [1.0, 3.0, 0.0], [0.5, 0.25, 5.0]],
            None,
        ),
    ]
    for frac_i, frac_j, lattice, expected_image in cases:
        rust = minimum_image_distance(frac_i, frac_j, lattice)
        pymatgen_distance, pymatgen_image = Lattice(lattice).get_distance_and_image(frac_i, frac_j)
        assert abs(rust["distance"] - float(pymatgen_distance)) <= 1e-10
        if expected_image is not None:
            assert rust["image"] == expected_image
            assert rust["image"] == [int(v) for v in pymatgen_image]


def test_randomized_rust_mic_distance_matches_pymatgen():
    rng = random.Random(20250220)
    for _ in range(200):
        frac_i = [rng.uniform(-2.0, 2.0) for _ in range(3)]
        frac_j = [rng.uniform(-2.0, 2.0) for _ in range(3)]
        lattice = [
            [rng.uniform(3.0, 10.0), 0.0, 0.0],
            [rng.uniform(-1.5, 1.5), rng.uniform(3.0, 10.0), 0.0],
            [rng.uniform(-1.5, 1.5), rng.uniform(-1.5, 1.5), rng.uniform(3.0, 10.0)],
        ]
        assert_mic_distance_matches_pymatgen(frac_i, frac_j, lattice)
