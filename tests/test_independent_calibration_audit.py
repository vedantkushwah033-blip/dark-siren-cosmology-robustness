import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration
from dark_siren_robustness.independent_calibration_audit import reconstruct_calibration_summaries


def fixture_records():
    rows = []
    for scenario, means in [(0, [68.0, 72.0]), (1, [69.0, 71.0, 69.0, 71.0])]:
        for method, se in [("paired_bootstrap", 1.0), ("design_aware_srswor", 0.5)]:
            for draw, mean in enumerate(means):
                rows.append({
                    "scenario": scenario, "validation_size": 100, "outer_draw": draw,
                    "method": method, "outer_status": "success",
                    "bootstrap_status": "success", "outer_posterior_mean": mean,
                    "bootstrap_se": se, "interval_low": 68.0, "interval_high": 72.0,
                    "finite_parent_target": 70.0,
                })
    return pd.DataFrame(rows)


def test_independent_summaries_match_reference_implementation():
    records = fixture_records()
    actual = reconstruct_calibration_summaries(records)
    expected = summarize_bootstrap_calibration(records)
    for left, right in zip(actual, expected):
        pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=False,
                                      rtol=1e-10, atol=1e-12)


def test_independent_summarizer_rejects_duplicate_design_rows():
    records = fixture_records()
    duplicate = pd.concat([records, records.iloc[[0]]], ignore_index=True)
    with pytest.raises(ValueError, match="unique"):
        reconstruct_calibration_summaries(duplicate)


def test_independent_summarizer_rejects_method_missing_from_outer_draw():
    records = fixture_records()
    records = records.loc[
        ~((records.scenario == 0) & (records.outer_draw == 0)
          & (records.method == "design_aware_srswor"))
    ]
    with pytest.raises(ValueError, match="method set"):
        reconstruct_calibration_summaries(records)


def test_independent_summarizer_rejects_disagreement_between_methods():
    records = fixture_records()
    records.loc[0, "outer_posterior_mean"] = 67.0
    with pytest.raises(ValueError, match="differs across methods"):
        reconstruct_calibration_summaries(records)


def test_stored_summary_comparison_rejects_changed_summary_value(tmp_path):
    from scripts.audit_nested_validation_bootstrap_replay import _require_summary_match

    scenario, _ = reconstruct_calibration_summaries(fixture_records())
    output = tmp_path / "scenario_summary.csv"
    scenario.to_csv(output, index=False)
    _require_summary_match(tmp_path, "scenario_summary.csv", scenario)

    altered = pd.read_csv(output)
    altered.loc[0, "finite_parent_target_interval_coverage"] = 0.123
    altered.to_csv(output, index=False)
    with pytest.raises(AssertionError, match="numeric column mismatch"):
        _require_summary_match(tmp_path, "scenario_summary.csv", scenario)
