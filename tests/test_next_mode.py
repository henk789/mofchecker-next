import json
from types import SimpleNamespace

import numpy as np
import pytest

pytest.importorskip("pymatgen")
from pymatgen.core import Lattice, Structure  # noqa: E402

from mofchecker_next import DEFAULT_MODE, MOFChecker  # noqa: E402
from mofchecker_next.batch import check_structure, summarize_results  # noqa: E402
from mofchecker_next.profiles import PROFILES, resolve_profile  # noqa: E402


def _co_cell():
    return Structure(Lattice.cubic(14), ["C", "O"], [[0.1, 0.1, 0.1], [0.2, 0.1, 0.1]])


def test_next_profile_resolution_and_default_are_explicit():
    assert DEFAULT_MODE == "0.9.6"
    assert resolve_profile("next").id == "next-dev-2"
    assert resolve_profile("next").provisional is True
    assert "has_stray_atom" not in resolve_profile("0.9.6").problem_flags
    assert "has_stray_atom" in resolve_profile("next").problem_flags
    assert resolve_profile("next").composite_name == "generated_desolvated_organic_mof_v1"
    assert "has_undercoordinated_c" not in resolve_profile("next").problem_flags
    assert "has_hydrogen" not in resolve_profile("next").presence_flags
    with pytest.raises(TypeError):
        PROFILES["next"] = resolve_profile("2.0")


def test_next_only_uses_connected_image_vectors_for_pbc_carbon_angle():
    structure = Structure(
        Lattice.cubic(10), ["C", "O", "O"],
        [[0.1, 0, 0], [0.8, 0, 0], [0.4, 0, 0]],
    )
    image_a = SimpleNamespace(index=1, site=SimpleNamespace(specie=structure[1].specie, coords=np.array([-2.0, 0, 0])))
    image_b = SimpleNamespace(index=2, site=SimpleNamespace(specie=structure[2].specie, coords=np.array([4.0, 0, 0])))
    graph = SimpleNamespace(get_connected_sites=lambda i: [image_a, image_b] if i == 0 else [])
    parity = MOFChecker(structure, mode="2.0")
    corrected = MOFChecker(structure, mode="next")
    parity.__dict__["graph"] = graph
    corrected.__dict__["graph"] = graph
    assert parity.undercoordinated_c_indices == [0]
    assert corrected.undercoordinated_c_indices == []


def test_next_eqeq_is_true_element_cached_and_exposes_metrics(monkeypatch):
    calls = []

    def fake(structure, **kwargs):
        calls.append(kwargs)
        return np.array([0.25, -0.25])

    monkeypatch.setattr("mofchecker_next.eqeq.compute_charges", fake)
    checker = MOFChecker(_co_cell(), mode="next")
    assert checker.max_abs_eqeq_charge == 0.25
    assert checker.eqeq_charge_sum == 0.0
    assert checker.eqeq_charge_threshold == 4.0
    assert checker.has_high_charges is False
    assert calls == [{"total_charge": 0.0, "reference_cif_labels": False}]


def test_next_default_descriptor_set_includes_compact_eqeq_metrics():
    result = check_structure(_co_cell(), mode="next")
    assert result["eqeq_charge_threshold"] == 4.0
    assert result["eqeq_charge_sum"] == pytest.approx(0.0, abs=1e-9)
    assert result["eqeq_expected_total_charge"] == 0.0
    assert result["eqeq_charge_residual"] == pytest.approx(0.0, abs=1e-9)
    assert isinstance(result["max_abs_eqeq_charge"], float)
    assert isinstance(result["has_high_charges"], bool)
    assert result["charge_status"] == "pass"
    assert result["scope_status"] == "fail"  # no metal in the two-atom probe
    assert {result[name] for name in (
        "structure_status", "local_chemistry_status", "component_status",
    )} <= {"pass", "fail"}


def test_explicit_per_atom_charges_are_json_safe():
    result = check_structure(_co_cell(), mode="next", descriptors=["eqeq_charges"])
    assert isinstance(result["eqeq_charges"], tuple)
    json.dumps(result)
    assert result["composite_status"] == "indeterminate"  # composite fields not requested


def test_next_eqeq_charge_sum_validation_fails_closed(monkeypatch):
    monkeypatch.setattr(
        "mofchecker_next.eqeq.compute_charges", lambda structure, **kwargs: np.array([0.2, 0.1])
    )
    result = check_structure(_co_cell(), mode="next", descriptors=["has_high_charges"])
    assert result["composite_status"] == "indeterminate"
    assert result["valid"] is None
    assert "charge sum" in result["errors"]["has_high_charges"]


def test_next_unsupported_eqeq_element_is_descriptor_error():
    structure = Structure(Lattice.cubic(12), ["Ac"], [[0, 0, 0]])
    result = check_structure(structure, mode="next", descriptors=["has_high_charges"])
    assert result["composite_status"] == "indeterminate"
    assert result["valid"] is None
    assert "outside EQeq's parameter table" in result["errors"]["has_high_charges"]
    assert result["charge_status"] == "indeterminate"


def test_next_result_metadata_and_summary_error_accounting():
    good = check_structure(_co_cell(), mode="next")
    assert good["mode"] == "next"
    assert good["resolved_profile"] == "next-dev-2"
    assert good["profile_provisional"] is True
    errored = good | {
        "errors": {"has_high_charges": "ValueError: unsupported"},
        "charge_status": "indeterminate",
        "composite_status": "indeterminate",
        "valid": None,
    }
    summary = summarize_results([good, errored])
    assert summary["mode"] == "next"
    assert summary["resolved_profile"] == "next-dev-2"
    assert summary["n_indeterminate"] == 1
    assert summary["n_invalid"] == 1
    assert summary["unconditional_valid_rate"] == 0.0
    assert summary["numeric_descriptor_summary"]["max_abs_eqeq_charge"]["n"] == 2
    assert summary["numeric_descriptor_summary"]["max_abs_eqeq_charge"]["max"] == good["max_abs_eqeq_charge"]
    assert summary["descriptor_error_counts"] == {"has_high_charges": 1}
    assert summary["error_categories"] == {"ValueError": 1}
    assert summary["status_counts"]["charge_status"] == {"pass": 1, "indeterminate": 1}
