import numpy as np
import pytest
from scripts.run_host_mixture_threshold_sensitivity import (
    cell_seed, expected_detected_host_fractions, run
)


def test_cell_seed_is_stable_and_design_specific():
    assert cell_seed(12, 70, (0.02, 0.04), (0.7, 0.3), 120) == cell_seed(12, 70, (0.02, 0.04), (0.7, 0.3), 120)
    assert cell_seed(12, 70, (0.02, 0.04), (0.7, 0.3), 120) != cell_seed(12, 70, (0.02, 0.05), (0.7, 0.3), 120)


def test_expected_detected_host_fractions_sum_to_one():
    p = expected_detected_host_fractions(70, (0.02, 0.04), (0.7, 0.3), 120)
    assert np.all(np.isfinite(p))
    assert np.all(p >= 0)
    assert np.isclose(p.sum(), 1.0)


def test_small_smoke_run_has_consistent_dimensions():
    rows, summaries, manifest = run(replicates=1, seed=17, grid_size=1001, bootstrap_draws=1000)
    assert len(rows) == 108
    assert manifest["cell_count"] == 108
    assert manifest["raw_rows"] == len(rows)
    assert all(0 < r["observed_distance_mpc"] < r["detection_max_mpc"] for r in rows)
    assert all(np.isfinite(r["observed_distance_mpc"]) for r in rows)
    assert len(summaries) == manifest["summary_rows"]
