from __future__ import annotations

import json

import pytest

from dark_siren_robustness.bootstrap_attempt_audit import canonical_inner_attempt_record


def test_canonical_inner_attempt_record_reconciles_successes_and_failures():
    result = {
        "n_attempted": 5,
        "n_successful": 3,
        "n_failed": 2,
        "success_mask": [True, False, True, False, True],
        "failure_codes": ["", "zero_variance", "", "zero_variance", ""],
    }
    record = canonical_inner_attempt_record(result)
    assert record["inner_attempts"] == 5
    assert record["inner_successes"] == 3
    assert record["inner_failures"] == 2
    assert json.loads(record["inner_failure_code_counts"]) == {"zero_variance": 2}


def test_canonical_inner_attempt_record_uses_stable_sorted_failure_codes():
    result = {
        "n_attempted": 3,
        "n_successful": 1,
        "n_failed": 2,
        "success_mask": [False, True, False],
        "failure_codes": ["z_failure", "", "a_failure"],
    }
    record = canonical_inner_attempt_record(result)
    assert record["inner_failure_code_counts"] == '{"a_failure":1,"z_failure":1}'


@pytest.mark.parametrize(
    "changes",
    [
        {"n_failed": 1},
        {"success_mask": [True, True, True, True, True]},
        {"failure_codes": ["", "", "", "", ""]},
        {"failure_codes": ["", "", "", "missing", ""]},
    ],
)
def test_canonical_inner_attempt_record_rejects_inconsistent_attempt_ledgers(changes):
    result = {
        "n_attempted": 5,
        "n_successful": 3,
        "n_failed": 2,
        "success_mask": [True, False, True, False, True],
        "failure_codes": ["", "zero_variance", "", "zero_variance", ""],
    }
    result.update(changes)
    with pytest.raises(ValueError):
        canonical_inner_attempt_record(result)
