from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration


def _fixture():
    rows = []
    # Scenario 0 has two successful outer draws with empirical SD 2.
    for draw, mean in enumerate([68.0, 72.0, np.nan]):
        status = "failed" if draw == 2 else "success"
        for method, se, lo, hi in [
            ("paired_bootstrap", 2.0, 66.0, 74.0),
            ("design_aware_srswor", 1.0, 69.0, 71.0),
        ]:
            rows.append({
                "scenario": 0, "validation_size": 100, "outer_draw": draw,
                "method": method, "status": status,
                "posterior_mean": mean if status == "success" else np.nan,
                "bootstrap_se": se if status == "success" else np.nan,
                "interval_low": lo if status == "success" else np.nan,
                "interval_high": hi if status == "success" else np.nan,
                "finite_parent_target": 70.0,
            })
    # Scenario 1 has many draws but remains one top-level scenario.
    for draw, mean in enumerate([69.0, 71.0, 69.0, 71.0]):
        for method, se, lo, hi in [
            ("paired_bootstrap", 1.0, 68.0, 72.0),
            ("design_aware_srswor", 0.5, 69.5, 70.5),
        ]:
            rows.append({
                "scenario": 1, "validation_size": 100, "outer_draw": draw,
                "method": method, "status": "success",
                "posterior_mean": mean, "bootstrap_se": se,
                "interval_low": lo, "interval_high": hi,
                "finite_parent_target": 70.0,
            })
    return pd.DataFrame(rows)


def test_scenario_balanced_calibration_uses_outer_srswor_sd_and_keeps_failures():
    scenario, aggregate = summarize_bootstrap_calibration(_fixture())
    first = scenario.query("scenario == 0 and method == 'paired_bootstrap'").iloc[0]
    assert first["n_outer_attempted"] == 3
    assert first["n_outer_successful"] == 2
    assert first["outer_failure_fraction"] == pytest.approx(1 / 3)
    assert first["empirical_srswor_sd"] == pytest.approx(np.std([68.0, 72.0], ddof=1))
    assert first["mean_se_to_empirical_sd_ratio"] == pytest.approx(1.0)

    # Scenario-level estimates are equally weighted, not pooled by outer draw.
    paired = aggregate.query("method == 'paired_bootstrap'").iloc[0]
    expected = np.mean([
        2.0 / np.std([68.0, 72.0], ddof=1),
        1.0 / np.std([69.0, 71.0, 69.0, 71.0], ddof=1),
    ])
    assert paired["scenario_balanced_mean_se_ratio"] == pytest.approx(expected)
    assert paired["n_scenarios"] == 2


def test_calibration_rejects_duplicate_grain_and_changing_finite_parent_target():
    frame = _fixture()
    with pytest.raises(ValueError, match="unique"):
        summarize_bootstrap_calibration(pd.concat([frame, frame.iloc[[0]]], ignore_index=True))
    frame.loc[0, "finite_parent_target"] = 71.0
    with pytest.raises(ValueError, match="fixed"):
        summarize_bootstrap_calibration(frame)


def test_calibration_rejects_invalid_success_metrics_and_reversed_intervals():
    frame = _fixture()
    frame.loc[0, "bootstrap_se"] = -1.0
    with pytest.raises(ValueError, match="nonnegative"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[0, "interval_low"] = 75.0
    with pytest.raises(ValueError, match="cannot exceed"):
        summarize_bootstrap_calibration(frame)


def test_calibration_requires_complete_attempt_records_and_known_status():
    frame = _fixture().drop(columns=["finite_parent_target"])
    with pytest.raises(ValueError, match="missing required"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[0, "status"] = "maybe"
    with pytest.raises(ValueError, match="status"):
        summarize_bootstrap_calibration(frame)
