from __future__ import annotations

import sys
from pathlib import Path

import pandas as pd
import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from summarize_selection_misspecification_sensitivity import summarize_paired_frames

ARMS = [
    "analytic_truth", "oracle", "catalogue_only", "correct_weighting",
    "shifted_transition", "shallow_slope",
]


def _frame():
    rows = []
    for replicate, seed in ((0, 20262000), (1, 20262001)):
        for arm in ARMS:
            mean = {
                "analytic_truth": 70.0, "oracle": 70.1, "catalogue_only": 60.0,
                "correct_weighting": 69.0, "shifted_transition": 69.8,
                "shallow_slope": 68.5,
            }[arm] + replicate * 0.1
            rows.append({
                "replicate": replicate, "seed": seed, "arm": arm,
                "h0_true": 70.0, "n_parent": 4000, "n_catalogue_observed": 2500,
                "n_detected": 30, "posterior_mean": mean,
                "ci68_low": mean - 3.0, "ci68_high": mean + 3.0,
                "ci68_width": 6.0, "density_ise": 0.1,
                "covered_68": bool(mean - 3.0 <= 70.0 <= mean + 3.0),
            })
    return pd.DataFrame(rows)


def test_identical_paired_configs_have_zero_change():
    frame = _frame()
    arms, primary = summarize_paired_frames(
        frame, frame.copy(), name="h0_321", seed=11, resamples=1000
    )
    assert len(arms) == 6
    assert all(row["max_abs_posterior_mean_change"] == pytest.approx(0.0) for row in arms)
    assert all(row["max_abs_ci_endpoint_change"] == pytest.approx(0.0) for row in arms)
    assert primary["contrast_change_vs_baseline"] == pytest.approx(0.0)
    assert primary["contrast_change_bootstrap_95_low"] == pytest.approx(0.0)
    assert primary["contrast_change_bootstrap_95_high"] == pytest.approx(0.0)


def test_primary_contrast_change_is_paired_and_directionally_correct():
    baseline = _frame()
    comparison = baseline.copy()
    mask = comparison["arm"] == "correct_weighting"
    comparison.loc[mask, "posterior_mean"] += 0.5
    _, primary = summarize_paired_frames(
        baseline, comparison, name="bandwidth_025", seed=12, resamples=1000
    )
    assert primary["contrast_change_vs_baseline"] == pytest.approx(0.5)


def test_pairing_mismatch_is_rejected():
    baseline = _frame()
    comparison = baseline.iloc[:-1].copy()
    with pytest.raises(ValueError, match="pairing differs"):
        summarize_paired_frames(baseline, comparison, name="h0_321", resamples=1000)
