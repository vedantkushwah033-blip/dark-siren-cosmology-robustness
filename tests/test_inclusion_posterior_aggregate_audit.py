import pandas as pd
import pytest

from scripts.audit_inclusion_posterior_aggregate import (
    _validate_all_attempt_bounds,
    _validate_paired_contrasts,
)


def _ledger():
    rows = []
    for size in (100, 300):
        rows.append({
            "cell_id": f"r1-n{size}", "replicate": 1,
            "experiment": "FIXED_TRUTH", "catalogue_floor": 0.0,
            "n_validation": size, "status": "SUCCESS",
            "h0_median": 70.0 + (1.0 if size == 100 else 0.5),
            "h0_true": 70.0,
            "h0_q16": 65.0, "h0_q84": 75.0 if size == 100 else 74.0,
            "h0_q05": 60.0, "h0_q95": 80.0 if size == 100 else 79.0,
            "h0_q025": 58.0, "h0_q975": 82.0 if size == 100 else 81.0,
            "covered_68": 1, "covered_90": 1, "covered_95": 1,
        })
    rows.append({
        "cell_id": "r2-n100", "replicate": 2,
        "experiment": "FIXED_TRUTH", "catalogue_floor": 0.0,
        "n_validation": 100, "status": "SUCCESS",
        "h0_median": 69.0, "h0_true": 70.0,
        "h0_q16": 64.0, "h0_q84": 76.0,
        "h0_q05": 59.0, "h0_q95": 81.0,
        "h0_q025": 57.0, "h0_q975": 83.0,
        "covered_68": 1, "covered_90": 1, "covered_95": 1,
    })
    rows.append({
        "cell_id": "r2-n300", "replicate": 2,
        "experiment": "FIXED_TRUTH", "catalogue_floor": 0.0,
        "n_validation": 300, "status": "FAILED",
    })
    return pd.DataFrame(rows)


def _stored_summary(frame):
    cells = {}
    for size, bounds in ((100, [1.0, 1.0]), (300, [0.5, 1.0])):
        subset = frame[frame.n_validation == size]
        item = {
            "experiment": "FIXED_TRUTH",
            "catalogue_floor": 0.0,
            "n_validation": size,
            "attempted": int((subset.status != "NOT_RUN").sum()),
            "success": int((subset.status == "SUCCESS").sum()),
        }
        for level in (68, 90, 95):
            item[f"coverage_{level}_all_attempt_bounds"] = bounds.copy()
        cells[f"FIXED_TRUTH|floor=0|n={size}"] = item

    # Hand-calculated expected results from the one common evaluable pair.
    metrics = {
        "posterior_median": -0.5,
        "signed_median_error": -0.5,
        "absolute_median_error": -0.5,
        "width_68": -1.0,
        "width_90": -1.0,
        "width_95": -1.0,
        "coverage_68": 0.0,
        "coverage_90": 0.0,
        "coverage_95": 0.0,
    }
    contrast = {
        "experiment": "FIXED_TRUTH",
        "catalogue_floor": 0.0,
        "smaller_validation_size": 100,
        "larger_validation_size": 300,
        "attempted_common_replicates": 2,
        "evaluable_pairs": 1,
        "excluded_pairs": 1,
        "difference_convention": "larger_validation_size_minus_smaller_validation_size",
        "metrics": {
            name: {
                "n_pairs": 1,
                "mean_delta": value,
                "paired_mc_se": None,
                "normal_approx_95_interval": None,
            }
            for name, value in metrics.items()
        },
    }
    return {
        "summary_by_cell": cells,
        "paired_validation_size_contrasts": {
            "FIXED_TRUTH|floor=0|100->300": contrast
        },
    }


def test_aggregate_auditor_reconstructs_failure_bounds_and_paired_contrasts():
    frame = _ledger()
    summary = _stored_summary(frame)
    assert _validate_all_attempt_bounds(frame, summary) == 6
    assert _validate_paired_contrasts(frame, summary) == 9


def test_aggregate_auditor_rejects_tampered_paired_summary():
    frame = _ledger()
    summary = _stored_summary(frame)
    summary["paired_validation_size_contrasts"][
        "FIXED_TRUTH|floor=0|100->300"
    ]["metrics"]["posterior_median"]["mean_delta"] += 0.25
    with pytest.raises(AssertionError, match="mean_delta"):
        _validate_paired_contrasts(frame, summary)


def test_aggregate_auditor_rejects_tampered_failure_inclusive_bound():
    frame = _ledger()
    summary = _stored_summary(frame)
    key = "FIXED_TRUTH|floor=0|n=300"
    summary["summary_by_cell"][key]["coverage_68_all_attempt_bounds"] = [0.0, 0.0]
    with pytest.raises(AssertionError, match="lower"):
        _validate_all_attempt_bounds(frame, summary)
