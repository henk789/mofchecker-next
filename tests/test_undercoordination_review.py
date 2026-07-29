"""The adjudication harness' scoring math (see scripts/undercoordination_review.py)."""

from __future__ import annotations

import importlib.util
from pathlib import Path

import pytest

_PATH = Path(__file__).resolve().parents[1] / "scripts" / "undercoordination_review.py"
_SPEC = importlib.util.spec_from_file_location("undercoordination_review", _PATH)
review = importlib.util.module_from_spec(_SPEC)
_SPEC.loader.exec_module(review)


def _key():
    return [
        {"case": "case-001", "check": "undercoordinated_c", "flagged_by": "0.9.6"},
        {"case": "case-002", "check": "undercoordinated_c", "flagged_by": "0.9.6"},
        {"case": "case-003", "check": "undercoordinated_n", "flagged_by": "2.0"},
    ]


def test_selfcheck_runs():
    review.selfcheck(None)


def test_flagging_version_gets_the_true_positive_and_the_other_a_false_negative():
    out = review.score_review(
        _key(),
        {"case-001": "needs_bond", "case-002": "ok", "case-003": "needs_bond"},
        {"undercoordinated_c|v096_only": 100, "undercoordinated_n|v20_only": 4},
        {},
    )
    assert out["versions"]["0.9.6"]["tp"] == 1
    assert out["versions"]["0.9.6"]["fp"] == 1
    assert out["versions"]["0.9.6"]["fn"] == 1  # missed the atom 2.0 flagged
    assert out["versions"]["2.0"]["tn"] == 1  # correctly silent on the spurious one
    assert out["favored_calls"] == {"0.9.6": 1, "2.0": 2}


def test_unsure_and_blank_verdicts_are_excluded():
    out = review.score_review(_key(), {"case-001": "unsure", "case-002": ""}, {}, {})
    assert out["n_labeled"] == 0
    assert out["skipped"]["unlabeled_or_unsure"] == 3


def test_pool_weighting_extrapolates_the_real_flag_rate():
    out = review.score_review(
        _key(),
        {"case-001": "needs_bond", "case-002": "needs_bond", "case-003": "ok"},
        {"undercoordinated_c|v096_only": 1040},
        {},
    )
    weighted = out["pool_weighted"]["undercoordinated_c|v096_only"]
    assert weighted["real_rate"] == 1.0
    assert weighted["estimated_real_flags_in_pool"] == 1040.0
    low, high = weighted["real_rate_95ci"]
    assert 0.0 < low < 1.0 and high == 1.0


@pytest.mark.parametrize(
    "old,new,expected",
    [(0, 0, 1.0), (5, 5, 1.0), (10, 0, True), (9, 1, True)],
)
def test_sign_test_two_sided(old, new, expected):
    p = review._sign_test(old, new)
    assert p == expected if expected is not True else p < 0.05


def test_wilson_interval_brackets_the_point_estimate():
    low, high = review._wilson(8, 10)
    assert low < 0.8 < high
    assert review._wilson(0, 0) == (0.0, 1.0)
