from __future__ import annotations

import numpy as np
import pandas as pd

from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration


def test_procedure_and_inner_resample_denominators_are_separate_and_exclude_fpc_derivative():
    rows = []
    statuses = ["success", "non_evaluable", "failed", "not_run"]
    counts = [(2, 2, 0), (3, 1, 2), (None, None, None), (0, 0, 0)]
    for draw, (status, count) in enumerate(zip(statuses, counts)):
        outer_status = "failed" if status == "not_run" else "success"
        outer_mean = np.nan if outer_status == "failed" else 70.0 + draw
        attempted, successful, failed = count
        for method in ("paired_bootstrap", "design_aware_srswor", "paired_bootstrap_fpc"):
            method_status = status if method != "design_aware_srswor" else (
                "not_run" if outer_status == "failed" else "success"
            )
            method_counts = count if method != "design_aware_srswor" else (
                (0, 0, 0) if outer_status == "failed" else (2, 2, 0)
            )
            is_derived = method == "paired_bootstrap_fpc"
            metrics_ok = method_status == "success"
            rows.append({
                "scenario": 0,
                "validation_size": 20,
                "outer_draw": draw,
                "method": method,
                "outer_status": outer_status,
                "bootstrap_status": method_status,
                "outer_posterior_mean": outer_mean,
                "bootstrap_se": 1.0 if metrics_ok else np.nan,
                "interval_low": 68.0 if metrics_ok else np.nan,
                "interval_high": 72.0 if metrics_ok else np.nan,
                "finite_parent_target": 70.0,
                "bootstrap_is_derived": is_derived,
                "n_bootstrap_attempted": method_counts[0],
                "n_bootstrap_successful": method_counts[1],
                "n_bootstrap_failed": method_counts[2],
            })

    scenario, _ = summarize_bootstrap_calibration(pd.DataFrame(rows))
    paired = scenario.loc[scenario["method"].eq("paired_bootstrap")].iloc[0]
    derived = scenario.loc[scenario["method"].eq("paired_bootstrap_fpc")].iloc[0]

    assert paired["n_bootstrap_procedures_attempted"] == 3
    assert paired["n_bootstrap_procedures_successful"] == 1
    assert paired["n_bootstrap_procedures_non_evaluable"] == 1
    assert paired["n_bootstrap_procedures_failed"] == 1
    assert paired["n_bootstrap_procedures_not_run"] == 1
    assert paired["n_bootstrap_procedures_with_unknown_inner_count"] == 1
    assert paired["n_inner_resamples_attempted"] == 5
    assert paired["n_inner_resamples_successful"] == 3
    assert paired["n_inner_resamples_failed"] == 2
    assert np.isclose(paired["inner_resample_failure_fraction"], 0.4)
    assert np.isclose(paired["bootstrap_failure_fraction"], 2 / 3)

    assert derived["n_bootstrap_procedures_attempted"] == 0
    assert derived["n_inner_resamples_attempted"] == 0


def test_bootstrap_count_reconciliation_rejects_inconsistent_ledger():
    records = pd.DataFrame([{
        "scenario": 0,
        "validation_size": 20,
        "outer_draw": 0,
        "method": "paired_bootstrap",
        "outer_status": "success",
        "bootstrap_status": "success",
        "outer_posterior_mean": 70.0,
        "bootstrap_se": 1.0,
        "interval_low": 68.0,
        "interval_high": 72.0,
        "finite_parent_target": 70.0,
        "n_bootstrap_attempted": 3,
        "n_bootstrap_successful": 2,
        "n_bootstrap_failed": 0,
    }])
    try:
        summarize_bootstrap_calibration(records)
    except ValueError as exc:
        assert "do not reconcile" in str(exc)
    else:
        raise AssertionError("inconsistent bootstrap counts were accepted")
