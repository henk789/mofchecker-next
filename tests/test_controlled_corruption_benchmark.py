import importlib.util
from pathlib import Path

import pytest

pytest.importorskip("pymatgen")
from pymatgen.core import Lattice, Structure

_SCRIPT = Path(__file__).parents[1] / "scripts" / "controlled_corruption_benchmark.py"
_spec = importlib.util.spec_from_file_location("controlled_corruption_benchmark", _SCRIPT)
benchmark = importlib.util.module_from_spec(_spec)
assert _spec.loader is not None
_spec.loader.exec_module(benchmark)


def test_wilson_and_confusion_are_fail_closed():
    assert benchmark.wilson_lower(100, 100) > 0.95
    metrics = benchmark.confusion(
        [False, True, True, False], [False, True, None, True]
    )
    assert metrics["tp"] == 1
    assert metrics["fp"] == 1
    assert metrics["tn"] == 1
    assert metrics["fn"] == 0
    assert metrics["indeterminate"] == 1
    assert metrics["precision"] == 0.5
    assert metrics["recall"] == 1.0


def test_overlap_corruption_is_deterministic_and_detected():
    structure = Structure(Lattice.cubic(12), ["C", "O"], [[0.1, 0.1, 0.1], [0.7, 0.7, 0.7]])
    corrupted = benchmark.overlap_corruption(structure)
    assert corrupted[0].frac_coords.tolist() == pytest.approx(corrupted[1].frac_coords.tolist())
    result = benchmark.check_structure(corrupted, mode="next", descriptors=["has_atomic_overlaps"])
    assert result["has_atomic_overlaps"] is True
    assert benchmark.split_for(b"same", 7) == benchmark.split_for(b"same", 7)
