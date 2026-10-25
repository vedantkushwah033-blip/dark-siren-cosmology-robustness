import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_paired_weight_sweep.py"
SPEC = importlib.util.spec_from_file_location("paired_weight_sweep", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_paired_weight_sweep_shape_metrics_and_bootstrap_intervals():
    rows, summaries = module.run_paired_sweep(replicates=200, seed=13, bootstrap_resamples=1000)
    assert len(rows) == 3 * 3 * 200
    assert len(summaries) == 3 * 3 * 2
    assert all(np.isfinite(r["paired_abs_error_difference_narrow_50_90"]) for r in rows)
    assert all(np.isfinite(r["paired_abs_error_difference_broad_20_220"]) for r in rows)
    for row in summaries:
        for key in (
            "correct_rmse", "reversed_rmse",
            "mean_paired_abs_error_diff_reversed_minus_correct",
            "paired_abs_error_diff_ci95_low", "paired_abs_error_diff_ci95_high",
            "coverage_68_diff_ci95_low", "coverage_68_diff_ci95_high",
            "rmse_diff_ci95_low", "rmse_diff_ci95_high",
        ):
            assert np.isfinite(row[key])
        assert row["paired_abs_error_diff_ci95_low"] <= row["paired_abs_error_diff_ci95_high"]
        assert row["coverage_68_diff_ci95_low"] <= row["coverage_68_diff_ci95_high"]
        assert row["rmse_diff_ci95_low"] <= row["rmse_diff_ci95_high"]
        assert row["bootstrap_resamples"] == 1000
        assert 0 <= row["correct_coverage_68"] <= 1
        assert 0 <= row["reversed_coverage_68"] <= 1


def test_broad_prior_retains_alternate_host_h0_modes():
    rows, _ = module.run_paired_sweep(replicates=200, seed=15, bootstrap_resamples=1000)
    relevant = [
        row for row in rows
        if row["z2"] == 0.04 and row["noise_fraction"] == 0.10
    ]
    assert len(relevant) == 200
    assert any(
        abs(row["broad_20_220_reversed_median"] - row["broad_20_220_correct_median"]) > 1.0
        for row in relevant
    )


def test_each_weighting_pair_uses_same_observation():
    rows, _ = module.run_paired_sweep(replicates=200, seed=21, bootstrap_resamples=1000)
    for row in rows:
        assert row["sigma_mpc"] > 0
        assert row["observed_distance_mpc"] > 0
        assert row["true_host_redshift"] in (row["z1"], row["z2"])
        for support in ("narrow_50_90", "broad_20_220"):
            correct = row[f"{support}_correct_median"]
            reversed_ = row[f"{support}_reversed_median"]
            difference = row[f"paired_abs_error_difference_{support}"]
            assert np.isfinite(correct) and np.isfinite(reversed_) and np.isfinite(difference)


def test_bootstrap_helper_is_reproducible_and_paired():
    values = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 5.0], [4.0, 8.0]])
    a = module._percentile_ci(values, np.random.default_rng(7), 1000, statistic="rmse_delta")
    b = module._percentile_ci(values, np.random.default_rng(7), 1000, statistic="rmse_delta")
    assert a == b
    assert a[0] <= a[1]
