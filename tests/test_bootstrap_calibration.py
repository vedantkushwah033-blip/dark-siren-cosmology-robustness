from __future__ import annotations

import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration


def _fixture():
    rows = []
    # Outer draw 2 fails for both methods; draw 1 has a design-aware inner failure.
    outer = [(0, 68.0, "success"), (1, 72.0, "success"), (2, np.nan, "failed")]
    outer_sd = float(np.std([68.0, 72.0], ddof=1))
    method_values = {
        "paired_bootstrap": [(outer_sd, 66.0, 74.0), (outer_sd, 66.0, 74.0), (np.nan, np.nan, np.nan)],
        "design_aware_srswor": [(1.0, 69.0, 71.0), (np.nan, np.nan, np.nan), (np.nan, np.nan, np.nan)],
    }
    for method, values in method_values.items():
        for (draw, mean, outer_status), (se, low, high) in zip(outer, values):
            bootstrap_status = "not_run" if outer_status == "failed" else (
                "failed" if np.isnan(se) else "success"
            )
            rows.append({
                "scenario": 0, "validation_size": 100, "outer_draw": draw,
                "method": method, "outer_status": outer_status,
                "bootstrap_status": bootstrap_status,
                "outer_posterior_mean": mean, "bootstrap_se": se,
                "interval_low": low, "interval_high": high,
                "finite_parent_target": 70.0,
            })
    # Scenario 1 has four outer draws, but remains one top-level scenario.
    for method, se, low, high in [
        ("paired_bootstrap", 1.0, 68.0, 72.0),
        ("design_aware_srswor", 0.5, 69.5, 70.5),
    ]:
        for draw, mean in enumerate([69.0, 71.0, 69.0, 71.0]):
            rows.append({
                "scenario": 1, "validation_size": 100, "outer_draw": draw,
                "method": method, "outer_status": "success",
                "bootstrap_status": "success", "outer_posterior_mean": mean,
                "bootstrap_se": se, "interval_low": low, "interval_high": high,
                "finite_parent_target": 70.0,
            })
    return pd.DataFrame(rows)


def test_scenario_balanced_calibration_separates_outer_and_inner_failures():
    scenario, aggregate = summarize_bootstrap_calibration(_fixture())
    first = scenario.query("scenario == 0 and method == 'paired_bootstrap'").iloc[0]
    assert first["n_outer_attempted"] == 3
    assert first["n_outer_successful"] == 2
    assert first["outer_failure_fraction"] == pytest.approx(1 / 3)
    assert first["bootstrap_failure_fraction"] == pytest.approx(0.0)
    assert first["empirical_srswor_sd"] == pytest.approx(np.std([68.0, 72.0], ddof=1))
    assert first["mean_se_to_empirical_sd_ratio"] == pytest.approx(1.0)

    design = scenario.query("scenario == 0 and method == 'design_aware_srswor'").iloc[0]
    assert design["empirical_srswor_sd"] == first["empirical_srswor_sd"]
    assert design["bootstrap_failure_fraction"] == pytest.approx(0.5)

    # Scenario-level estimates are equally weighted, not pooled by outer draw.
    paired = aggregate.query("method == 'paired_bootstrap'").iloc[0]
    expected = np.mean([
        1.0,
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


def test_calibration_requires_complete_attempt_records_and_consistent_outer_results():
    frame = _fixture().drop(columns=["finite_parent_target"])
    with pytest.raises(ValueError, match="missing required"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[1, "outer_posterior_mean"] = 71.0
    with pytest.raises(ValueError, match="identical across bootstrap methods"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[0, "outer_status"] = "maybe"
    with pytest.raises(ValueError, match="outer_status"):
        summarize_bootstrap_calibration(frame)


def test_calibration_rejects_invalid_design_keys_and_missing_finite_frame_target():
    frame = _fixture()
    frame.loc[0, "validation_size"] = 0
    with pytest.raises(ValueError, match="validation_size"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[0, "outer_draw"] = -1
    with pytest.raises(ValueError, match="outer_draw"):
        summarize_bootstrap_calibration(frame)
    frame = _fixture()
    frame.loc[0, "finite_parent_target"] = np.nan
    with pytest.raises(ValueError, match="finite_parent_target"):
        summarize_bootstrap_calibration(frame)


def test_calibration_rejects_method_specific_missing_outer_draw():
    frame = _fixture()
    frame = frame.loc[
        ~((frame["scenario"] == 0) & (frame["outer_draw"] == 1)
          & (frame["method"] == "design_aware_srswor"))
    ]
    with pytest.raises(ValueError, match="complete method set"):
        summarize_bootstrap_calibration(frame)
