import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_paired_weight_sweep.py"
SPEC = importlib.util.spec_from_file_location("paired_weight_sweep", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_paired_weight_sweep_shape_and_metrics():
    rows, summaries = module.run_paired_sweep(replicates=200, seed=13)
    assert len(rows) == 3 * 3 * 200
    assert len(summaries) == 9
    assert all(np.isfinite(r["paired_abs_error_difference_reversed_minus_correct"]) for r in rows)
    assert all(0 <= r["correct_coverage_68"] <= 1 for r in summaries)
    assert all(0 <= r["reversed_coverage_68"] <= 1 for r in summaries)


def test_paired_comparison_uses_same_observations():
    rows, _ = module.run_paired_sweep(replicates=200, seed=15)
    for row in rows:
        assert row["correct_median"] != row["reversed_median"] or row["paired_abs_error_difference_reversed_minus_correct"] == 0
        assert row["sigma_mpc"] > 0
        assert row["observed_distance_mpc"] > 0
