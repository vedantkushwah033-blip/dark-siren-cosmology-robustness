import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_paired_weight_sweep.py"
SPEC = importlib.util.spec_from_file_location("paired_weight_sweep", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_paired_weight_sweep_shape_and_metrics():
    rows, summaries = module.run_paired_sweep(replicates=200, seed=13, bootstrap_resamples=1000)
    assert len(rows) == 3 * 3 * 200
    assert len(summaries) == 9
    assert all(np.isfinite(r["paired_abs_error_difference_reversed_minus_correct"]) for r in rows)
    assert all(0 <= r["correct_coverage_68"] <= 1 for r in summaries)
    assert all(0 <= r["reversed_coverage_68"] <= 1 for r in summaries)
    for row in summaries:
        for key in (
            "paired_abs_error_diff_ci95_low", "paired_abs_error_diff_ci95_high",
            "coverage_68_diff_ci95_low", "coverage_68_diff_ci95_high",
            "rmse_diff_ci95_low", "rmse_diff_ci95_high",
        ):
            assert np.isfinite(row[key])
        assert row["paired_abs_error_diff_ci95_low"] <= row["paired_abs_error_diff_ci95_high"]
        assert row["coverage_68_diff_ci95_low"] <= row["coverage_68_diff_ci95_high"]
        assert row["rmse_diff_ci95_low"] <= row["rmse_diff_ci95_high"]
        assert row["bootstrap_resamples"] == 1000


def test_paired_comparison_uses_same_observations():
    rows, _ = module.run_paired_sweep(replicates=200, seed=15, bootstrap_resamples=1000)
    for row in rows:
        assert row["correct_median"] != row["reversed_median"] or row["paired_abs_error_difference_reversed_minus_correct"] == 0
        assert row["sigma_mpc"] > 0
        assert row["observed_distance_mpc"] > 0


def test_bootstrap_helper_is_reproducible_and_paired():
    values = np.array([[1.0, 2.0], [2.0, 4.0], [3.0, 5.0], [4.0, 8.0]])
    a = module._percentile_ci(values, np.random.default_rng(7), 1000, statistic="rmse_delta")
    b = module._percentile_ci(values, np.random.default_rng(7), 1000, statistic="rmse_delta")
    assert a == b
    assert a[0] <= a[1]
