import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "run_selection_paired_comparison.py"
SPEC = importlib.util.spec_from_file_location("selection_paired_comparison", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_generator_returns_only_detected_observations():
    rng = np.random.default_rng(27)
    for truth in module.TRUTHS:
        for _ in range(20):
            host, distance = module.draw_detected_observation(truth, rng)
            assert host in (0, 1)
            assert 0 < distance < module.DETECTION_DISTANCE_MAX


def test_paired_comparison_is_reproducible_and_finite():
    rows1, summaries1 = module.run_comparison(replicates_per_truth=100, seed=81, grid_size=1001)
    rows2, summaries2 = module.run_comparison(replicates_per_truth=100, seed=81, grid_size=1001)
    assert rows1 == rows2
    assert summaries1 == summaries2
    assert len(rows1) == 300
    assert len(summaries1) == 8
    for row in rows1:
        assert np.isfinite(row["uncorrected_median"])
        assert np.isfinite(row["corrected_median"])
        assert np.isfinite(row["paired_abs_error_change_corrected_minus_uncorrected"])
        assert 0 < row["observed_distance_mpc"] < module.DETECTION_DISTANCE_MAX


def test_comparison_rejects_too_few_replicates():
    with pytest.raises(ValueError, match="at least 100"):
        module.run_comparison(replicates_per_truth=20)
