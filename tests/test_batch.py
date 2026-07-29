"""Tests for the batch validation API."""

from __future__ import annotations

import numpy as np
import pytest

pytest.importorskip("pymatgen")
pytest.importorskip("structuregraph_helpers")
from pymatgen.core import Lattice, Structure  # noqa: E402

from mofchecker_next.batch import (  # noqa: E402
    DEFAULT_DESCRIPTORS,
    LEGACY_DEFAULT_DESCRIPTORS,
    PRESENCE_FLAGS,
    PROBLEM_FLAGS,
    check_structure,
    check_structures,
    is_valid,
    normalize_structure,
    summarize_results,
)


def _mof_like():
    # Small Zn + carboxylate-ish cell; enough to exercise the graph checks.
    lattice = Lattice.cubic(14.0)
    species = ["Zn", "O", "O", "C", "C", "H"]
    coords = [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1], [0.1, 0.2, 0.1],
              [0.25, 0.18, 0.1], [0.35, 0.25, 0.1], [0.4, 0.3, 0.1]]
    return Structure(lattice, species, coords)


def test_normalize_accepts_structure_and_path(tmp_path):
    s = _mof_like()
    assert normalize_structure(s) is s
    p = tmp_path / "x.cif"
    s.to(filename=str(p), fmt="cif")
    assert len(normalize_structure(str(p))) == len(s)
    assert len(normalize_structure(p)) == len(s)


def test_normalize_rejects_unknown():
    with pytest.raises(TypeError):
        normalize_structure(42)


def test_check_structure_defaults_to_legacy_descriptors():
    r = check_structure(_mof_like())
    assert "errors" not in r
    assert r["n_atoms"] == 6
    for name in LEGACY_DEFAULT_DESCRIPTORS:
        assert name in r
    assert r["has_metal"] is True
    assert "metal_number" not in r
    assert r["composite_status"] in {"valid", "invalid"}


def test_check_structure_2_0_descriptors_remain_available():
    r = check_structure(_mof_like(), mode="2.0")
    for name in DEFAULT_DESCRIPTORS:
        assert name in r
    assert isinstance(r["metal_number"], int)


def test_subset_descriptors_skip_graph():
    # Composition-only subset should not need the graph and must be exact.
    r = check_structure(_mof_like(), descriptors=["has_carbon", "has_metal", "metal_number"])
    assert {"n_atoms", "has_carbon", "has_metal", "metal_number", "mode",
            "resolved_profile", "profile_provisional", "composite_status", "valid",
            "composite_name", "composites", "implementation", "input_protocol",
            "evaluation_settings"} <= set(r)
    assert "has_atomic_overlaps" not in r
    assert r["has_carbon"] is True
    assert r["metal_number"] == 1


def test_check_structures_parallel_matches_serial():
    structures = [_mof_like() for _ in range(4)]
    serial = check_structures(structures, n_workers=1)
    parallel = check_structures(structures, n_workers=2)

    def strip(r):
        return {k: v for k, v in r.items() if k != "index"}

    assert [strip(a) for a in serial] == [strip(b) for b in parallel]
    assert [r["index"] for r in parallel] == [0, 1, 2, 3]


def test_check_structures_records_errors():
    results = check_structures([_mof_like(), 12345], n_workers=1, on_error="record")
    assert "error" in results[1]
    assert "errors" not in results[0] and results[0]["has_metal"] is True


def test_check_structures_mixed_inputs(tmp_path):
    s = _mof_like()
    p = tmp_path / "m.cif"
    s.to(filename=str(p), fmt="cif")
    results = check_structures([s, str(p)], n_workers=1)
    assert results[0]["n_atoms"] == results[1]["n_atoms"] == 6
    assert results[1]["id"] == "m.cif"


def _composite_result(**values):
    return {
        "mode": "2.0",
        **{field: True for field in PRESENCE_FLAGS},
        **{field: False for field in PROBLEM_FLAGS},
        **values,
    }


def test_summarize_results():
    good = _composite_result(id="good")
    bad = _composite_result(id="bad", has_stray_atom=True)
    err = {"id": "err", "error": "boom"}
    assert is_valid(good) is True
    assert is_valid(bad) is False
    assert is_valid(err) is None
    assert is_valid(good | {"errors": {"has_high_charges": "failed"}}) is None
    assert is_valid(good | {"errors": {"has_oms": "failed"}}) is True
    s = summarize_results([good, bad, err])
    assert s["n_structures"] == 3
    assert s["n_scored"] == 2
    assert s["n_valid"] == 1
    assert s["valid_rate"] == 0.5
    assert s["valid_rate_incl_errors"] == 1 / 3
    descriptor_error = summarize_results([good | {"errors": {"has_high_charges": "failed"}}])
    assert descriptor_error["n_errors"] == 1
    assert descriptor_error["n_scored"] == 0
    assert descriptor_error["errors"]["0"]["errors"] == {"has_high_charges": "failed"}


def test_missing_or_nonboolean_composite_field_is_indeterminate():
    assert is_valid({"has_carbon": True, "has_hydrogen": True, "has_metal": True}) is None
    assert is_valid(_composite_result(has_high_charges=None)) is None


def test_stray_atom_preserves_composite_validity():
    base = _composite_result()
    assert is_valid(base | {"has_stray_atom": True, "has_lone_molecule": False}) is False
    assert is_valid(base | {"has_stray_atom": False, "has_lone_molecule": True}) is False
