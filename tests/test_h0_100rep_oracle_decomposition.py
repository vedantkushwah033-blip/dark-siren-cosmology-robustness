from pathlib import Path
import sys

import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_h0_100rep_oracle_decomposition import (
    assert_oracle_eligible_statuses, assert_planned_design, coverage_bounds, wilson_interval,
)


def test_non_evaluable_cell_is_retained_for_oracle_reconstruction():
    ledger = pd.DataFrame({"status": ["SUCCESS", "NON_EVALUABLE", "SUCCESS"]})
    assert_oracle_eligible_statuses(ledger)


def test_failed_or_interrupted_cell_is_not_silently_accepted():
    for status in ("FAILED", "RUNNING", "NOT_RUN"):
        ledger = pd.DataFrame({"status": ["SUCCESS", status]})
        with pytest.raises(AssertionError, match="unsupported statuses"):
            assert_oracle_eligible_statuses(ledger)


def test_all_attempt_coverage_bounds_keep_planned_denominator():
    result = coverage_bounds(success_covered=54, success=99, planned=100)
    assert result["success_only"] == pytest.approx(54/99)
    assert result["all_attempt_lower"] == pytest.approx(0.54)
    assert result["all_attempt_upper"] == pytest.approx(0.55)
    assert result["unresolved"] == 1


def test_coverage_bounds_validate_counts():
    with pytest.raises(ValueError):
        coverage_bounds(success_covered=3, success=2, planned=4)


def test_wilson_interval_is_ordered_and_bounded():
    low, high = wilson_interval(54, 99)
    assert 0 <= low < 54/99 < high <= 1



def _small_planned_design():
    rows = []
    for rep in range(2):
        for gamma in (-1.0, 0.0, 1.0):
            for nval in (100, 500):
                rows.append({
                    "replicate": rep, "gamma": gamma, "n_validation": nval,
                    "analysis_seed": 1000+rep, "validation_seed": 2000+rep,
                    "h0_true": 70.0,
                })
    return pd.DataFrame(rows)


def test_planned_design_checks_factorial_and_pairing():
    assert_planned_design(_small_planned_design(), expected_replicates=2)


def test_planned_design_rejects_unpaired_seeds():
    planned = _small_planned_design()
    planned.loc[(planned.replicate == 0) & (planned.gamma == 1.0), "analysis_seed"] += 1
    with pytest.raises(AssertionError, match="seeds are not paired"):
        assert_planned_design(planned, expected_replicates=2)
