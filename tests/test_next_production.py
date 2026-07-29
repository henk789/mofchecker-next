import json
import subprocess
import sys
from dataclasses import asdict

import numpy as np
import pytest

pytest.importorskip("pymatgen")
from pymatgen.core import Lattice, Structure

from mofchecker_next.batch import check_structure, is_valid, summarize_results
from mofchecker_next.profiles import resolve_profile


def _probe():
    return Structure(
        Lattice.cubic(14), ["Zn", "O", "C", "H"],
        [[0.10, 0.10, 0.10], [0.18, 0.10, 0.10], [0.25, 0.10, 0.10], [0.31, 0.10, 0.10]],
    )


def _passing_composite():
    profile = resolve_profile("next")
    return {
        "mode": "next",
        "resolved_profile": profile.id,
        "composite_problem_flags": list(profile.problem_flags),
        "composite_presence_flags": list(profile.presence_flags),
        **{field: False for field in profile.problem_flags},
        **{field: True for field in profile.presence_flags},
    }


def test_compatibility_profiles_are_locked():
    legacy = resolve_profile("0.9.6")
    modern = resolve_profile("2.0")
    assert asdict(legacy) == {
        "id": "0.9.6", "descriptor_set": "legacy",
        "composite_name": "adit_mofasa_all_passed",
        "problem_flags": (
            "has_atomic_overlaps", "has_overcoordinated_c", "has_overcoordinated_n",
            "has_overcoordinated_h", "has_undercoordinated_c", "has_undercoordinated_n",
            "has_undercoordinated_rare_earth", "has_undercoordinated_alkali_alkaline",
            "has_lone_molecule", "has_suspicious_terminal_oxo", "has_high_charges",
            "has_geometrically_exposed_metal",
        ),
        "presence_flags": ("has_carbon", "has_hydrogen", "has_metal"),
        "status_groups": (), "report_only_fields": (), "legacy_input": True,
        "corrected_carbon_images": False, "corrected_input": False,
        "total_charge": None, "eqeq_threshold": 3.0,
        "eqeq_charge_sum_tolerance": 1e-9, "provisional": False,
    }
    assert modern.id == "2.0"
    assert modern.descriptor_set == "modern"
    assert modern.problem_flags[-1] == "has_geometrically_exposed_metal"
    assert modern.presence_flags == ("has_carbon", "has_hydrogen", "has_metal")
    assert modern.eqeq_threshold == 4.0
    assert modern.corrected_input is False


def test_corrected_composite_required_and_report_only_semantics():
    base = _passing_composite()
    assert is_valid(base) is True
    profile = resolve_profile("next")
    for field in profile.problem_flags:
        assert is_valid(base | {field: True}) is False
    for field in profile.presence_flags:
        assert is_valid(base | {field: False}) is False
    for field in (*profile.problem_flags, *profile.presence_flags):
        missing = dict(base)
        missing.pop(field)
        assert is_valid(missing) is None
        assert is_valid(base | {"errors": {field: "failed"}}) is None
    forged = base | {"composite_problem_flags": [], "composite_presence_flags": []}
    assert is_valid(forged) is None
    forged_summary = summarize_results([forged])
    assert forged_summary["n_errors"] == 1
    assert forged_summary["error_categories"] == {"CompositeSchemaError": 1}
    assert forged_summary["descriptor_error_counts"] == {"__composite_schema__": 1}
    assert is_valid(base | {"resolved_profile": "next-dev-1"}) is None
    assert is_valid(base | {"has_undercoordinated_c": True}) is True
    assert is_valid(base | {"errors": {"has_undercoordinated_c": "report-only failed"}}) is True


def test_corrected_input_rejects_nan_and_disorder_without_fabricated_boolean():
    nan_structure = Structure(Lattice.cubic(10), ["C"], [[0, 0, 0]])
    nan_structure.translate_sites([0], [np.nan, 0, 0], frac_coords=True, to_unit_cell=False)
    disordered = Structure(Lattice.cubic(10), [{"C": 0.5, "N": 0.5}], [[0, 0, 0]])
    for structure, text in ((nan_structure, "NaN or infinity"), (disordered, "partial occupancy")):
        result = check_structure(structure, mode="next")
        assert result["valid"] is None
        assert result["composite_status"] == "indeterminate"
        assert result["structure_status"] == "indeterminate"
        assert result["check_results"]["input_sanity"]["status"] == "indeterminate"
        assert text in result["errors"]["input_sanity"]
        assert "has_atomic_overlaps" not in result


def test_custom_metal_policy_is_recorded():
    result = check_structure(_probe(), mode="next", descriptors=["has_metal"], metals={"Zn", "Cu"})
    assert result["evaluation_settings"]["metals"] == ["Cu", "Zn"]


def test_corrected_check_evidence_is_json_safe_and_consistent(monkeypatch):
    monkeypatch.setattr(
        "mofchecker_next.eqeq.compute_charges",
        lambda structure, **kwargs: np.array([5.0, -2.0, -2.0, -1.0]),
    )
    result = check_structure(_probe(), mode="next")
    charge = result["check_results"]["has_high_charges"]
    assert result["has_high_charges"] is True
    assert charge["status"] == "fail"
    assert charge["atoms"] == [{"index": 0, "image": [0, 0, 0]}]
    assert charge["values"]["threshold"] == 4.0
    assert "true-element EQeq" in charge["rule"]
    assert result["composites"][result["composite_name"]]["value"] is False
    json.dumps(result)


def test_corrected_mechanical_fields_are_translation_permutation_and_cif_invariant(tmp_path):
    descriptors = [
        "has_carbon", "has_metal", "has_atomic_overlaps", "has_overcoordinated_c",
        "has_overcoordinated_n", "has_overcoordinated_h", "has_stray_atom",
        "has_lone_molecule",
    ]
    structure = _probe()
    translated = structure.copy()
    translated.translate_sites(range(len(translated)), [1, -2, 3], frac_coords=True, to_unit_cell=False)
    shifted = structure.copy()
    shifted.translate_sites(range(len(shifted)), [0.137, 0.211, 0.319], frac_coords=True, to_unit_cell=False)
    order = [2, 0, 3, 1]
    permuted = Structure(structure.lattice, [structure[i].species for i in order], [structure[i].frac_coords for i in order])
    supercell = structure.copy()
    supercell.make_supercell([2, 1, 1])
    path = tmp_path / "probe.cif"
    structure.to(filename=path, fmt="cif")
    rows = [check_structure(value, mode="next", descriptors=descriptors) for value in (structure, translated, shifted, permuted, supercell, path)]
    expected = tuple(rows[0][field] for field in descriptors)
    assert all(tuple(row[field] for field in descriptors) == expected for row in rows[1:])
    assert rows[-1]["input_protocol"]["coordinate_snapping"] is False


def test_cli_writes_reproducible_summary_and_ordered_results(tmp_path):
    path = tmp_path / "probe.cif"
    _probe().to(filename=path, fmt="cif")
    summary_path = tmp_path / "summary.json"
    results_path = tmp_path / "results.jsonl"
    subprocess.run([
        sys.executable, "-m", "mofchecker_next.cli", str(path), "--mode", "next",
        "--soft_timeout", "--no_progress", "--output_json", str(summary_path),
        "--output_results_jsonl", str(results_path),
    ], check=True, capture_output=True, text=True)
    summary = json.loads(summary_path.read_text())
    result = json.loads(results_path.read_text())
    assert summary["input_manifest"]["n_inputs"] == 1
    assert len(summary["input_manifest"]["ordered_paths_and_bytes_sha256"]) == 64
    assert summary["manifest"]["resolved_profile"] == "next-dev-2"
    assert result["index"] == 0
    assert result["composite_name"] == "generated_desolvated_organic_mof_v1"


def test_summary_contains_reproducible_manifest_and_stable_error_keys():
    good = check_structure(_probe(), mode="next")
    bad_a = good | {"index": 7, "id": "same.cif", "error": "TimeoutError: first", "valid": None}
    bad_b = good | {"index": 8, "id": "same.cif", "error": "TimeoutError: second", "valid": None}
    summary = summarize_results([bad_a, bad_b])
    assert set(summary["errors"]) == {"7", "8"}
    assert summary["composite_name"] == "generated_desolvated_organic_mof_v1"
    assert summary["manifest"]["resolved_profile"] == "next-dev-2"
    assert len(summary["manifest"]["source_tree_sha256"]) == 64
    assert len(summary["manifest"]["rust_extension_sha256"]) == 64
    assert summary["manifest"]["profile"]["total_charge"] == 0.0
    json.dumps(summary)
