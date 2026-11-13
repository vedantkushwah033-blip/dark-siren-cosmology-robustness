import pandas as pd
import pytest

from scripts.audit_inclusion_posterior_aggregate import (
    _reconstruct_paired_contrasts,
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


def _summary_from(frame):
    cells = {}
    for size in (100, 300):
        subset = frame[(frame.n_validation == size)]
        cells[f"FIXED_TRUTH|floor=0|n={size}"] = {
            "experiment": "FIXED_TRUTH",
            "catalogue_floor": 0.0,
            "n_validation": size,
            "attempted": int((subset.status != "NOT_RUN").sum()),
            "success": int((subset.status == "SUCCESS").sum()),
        }
    summary = {"summary_by_cell": cells}
    summary["paired_validation_size_contrasts"] = _reconstruct_paired_contrasts(frame)
    for size in (100, 300):
        subset = frame[frame.n_validation == size]
        item = cells[f"FIXED_TRUTH|floor=0|n={size}"]
        for level in (68, 90, 95):
            from scripts.audit_inclusion_posterior_aggregate import _attempt_bounds
            item[f"coverage_{level}_all_attempt_bounds"] = _attempt_bounds(
                subset, f"covered_{level}"
            )
    return summary


def test_aggregate_auditor_reconstructs_failure_bounds_and_paired_contrasts():
    frame = _ledger()
    summary = _summary_from(frame)
    assert _validate_all_attempt_bounds(frame, summary) == 6
    assert _validate_paired_contrasts(frame, summary) == 9


def test_aggregate_auditor_rejects_tampered_paired_summary():
    frame = _ledger()
    summary = _summary_from(frame)
    key = next(iter(summary["paired_validation_size_contrasts"]))
    summary["paired_validation_size_contrasts"][key]["metrics"]["posterior_median"][
        "mean_delta"
    ] += 0.25
    with pytest.raises(AssertionError, match="mean_delta"):
        _validate_paired_contrasts(frame, summary)


def test_aggregate_auditor_rejects_tampered_failure_inclusive_bound():
    frame = _ledger()
    summary = _summary_from(frame)
    key = "FIXED_TRUTH|floor=0|n=300"
    summary["summary_by_cell"][key]["coverage_68_all_attempt_bounds"] = [0.0, 0.0]
    with pytest.raises(AssertionError, match="lower"):
        _validate_all_attempt_bounds(frame, summary)
