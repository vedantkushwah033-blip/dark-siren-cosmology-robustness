"""Regression tests for multi-replicate lifecycle aggregation."""
from pathlib import Path
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from run_five_arm_replicates import ARMS, _aggregate, _seed_schedule


def test_seed_schedule_is_reproducible_and_collision_free():
    first = _seed_schedule(20261115, 20)
    second = _seed_schedule(20261115, 20)
    assert first == second
    flattened = [value for pair in first for value in pair]
    assert len(set(flattened)) == 40


def test_failure_inclusive_bounds_and_paired_contrasts():
    summary_rows = []
    ledger_rows = []
    errors = {
        "A_catalogue_only": 1.0,
        "B_parent_oracle": 0.5,
        "C_known_probability_ipw": 0.8,
        "D_map_plugin_ipw": 1.2,
        "E_uncertainty_integrated_ipw": 0.6,
    }
    for arm in ARMS:
        error = errors[arm]
        summary_rows.append({
            "replicate": 1,
            "validation_size": 100,
            "arm": arm,
            "status": "SUCCESS",
            "h0_true": 70.0,
            "median_error": error,
            "posterior_median": 70.0 + error,
            "q16": 65.0,
            "q84": 75.0,
        })
        ledger_rows.append({
            "replicate": 1, "validation_size": 100, "arm": arm,
            "status": "SUCCESS", "audit_status": "PASS",
        })
        ledger_rows.append({
            "replicate": 2, "validation_size": 100, "arm": arm,
            "status": "FAILED", "audit_status": "NOT_RUN",
        })
    summary = pd.DataFrame(summary_rows)
    ledger = pd.DataFrame(ledger_rows)
    aggregate, contrasts = _aggregate(summary, ledger, replicates=2, max_abs_error=20.0)
    row = aggregate.loc[aggregate.arm == "A_catalogue_only"].iloc[0]
    assert row.successful_replicates == 1
    assert row.failed_or_non_evaluable == 1
    assert row.failure_inclusive_coverage_lower == 0.5
    assert row.failure_inclusive_coverage_upper == 1.0
    assert row.failure_inclusive_rmse_lower == np.sqrt(0.5)
    assert row.failure_inclusive_rmse_upper == np.sqrt(200.5)
    assert row.failure_inclusive_mean_68pct_width_lower == 5.0
    assert row.failure_inclusive_mean_68pct_width_upper == 25.0
    assert len(contrasts) == 8
    assert set(contrasts.reference_arm) == {
        "A_catalogue_only", "C_known_probability_ipw"
    }


def test_all_failed_empty_summary_retains_aggregate_cells():
    ledger = pd.DataFrame([
        {"replicate": rep, "validation_size": 100, "arm": arm, "status": "FAILED"}
        for rep in (1, 2) for arm in ARMS
    ])
    summary = pd.DataFrame(columns=[
        "replicate", "validation_size", "arm", "status", "h0_true",
        "median_error", "q16", "q84",
    ])
    aggregate, contrasts = _aggregate(summary, ledger, replicates=2)
    assert len(aggregate) == len(ARMS)
    assert (aggregate.successful_replicates == 0).all()
    assert (aggregate.failed_or_non_evaluable == 2).all()
    assert contrasts.empty


def test_failure_bounds_expand_when_true_h0_is_off_centre():
    summary = pd.DataFrame([{
        "replicate": 1,
        "validation_size": 100,
        "arm": "A_catalogue_only",
        "status": "SUCCESS",
        "h0_true": 60.0,
        "median_error": 1.0,
        "posterior_median": 61.0,
        "q16": 55.0,
        "q84": 65.0,
    }])
    ledger = pd.DataFrame([
        {"replicate": 1, "validation_size": 100, "arm": "A_catalogue_only", "status": "SUCCESS"},
        {"replicate": 2, "validation_size": 100, "arm": "A_catalogue_only", "status": "FAILED"},
    ])
    # With H0_true=60 and prior [50, 90], the maximum absolute median error
    # is 30, not 20. The aggregation must use that actual support-derived bound.
    aggregate, _ = _aggregate(
        summary, ledger, replicates=2, max_abs_error=30.0
    )
    row = aggregate.iloc[0]
    assert row.failure_inclusive_rmse_upper == np.sqrt(450.5)
    assert row.failure_inclusive_abs_error_upper == 15.5
    assert row.failure_inclusive_bias_lower == -14.5
    assert row.failure_inclusive_bias_upper == 15.5
