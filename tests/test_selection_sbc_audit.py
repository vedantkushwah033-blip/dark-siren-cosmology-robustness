import importlib.util
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
MODULE_PATH = ROOT / "scripts" / "run_selection_sbc_audit.py"
SPEC = importlib.util.spec_from_file_location("selection_sbc_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_detected_observation_generator_respects_threshold():
    rng = np.random.default_rng(15)
    for h0 in (50.0, 70.0, 90.0):
        for _ in range(30):
            host, distance = module.draw_detected_observation(h0, rng)
            assert host in (0, 1)
            assert 0.0 < distance < module.DETECTION_DISTANCE_MAX


def test_selection_sbc_small_run_is_finite_and_reproducible():
    rows1, summary1 = module.run_audit(replicates=500, seed=33, grid_size=1001)
    rows2, summary2 = module.run_audit(replicates=500, seed=33, grid_size=1001)
    assert rows1 == rows2
    assert summary1 == summary2
    assert len(rows1) == 500
    ranks = np.array([row["posterior_cdf_rank"] for row in rows1])
    assert np.all(np.isfinite(ranks))
    assert np.all((ranks >= 0.0) & (ranks <= 1.0))
    assert np.all(np.isfinite(uncorrected_ranks))
    assert np.all((uncorrected_ranks >= 0.0) & (uncorrected_ranks <= 1.0))
    assert "uncorrected_ks_pvalue" in summary1
    assert 0 < summary1["beta_min_on_prior_grid"] <= summary1["beta_max_on_prior_grid"] <= 1


def test_selection_sbc_rejects_too_few_replicates():
    with pytest.raises(ValueError, match="at least 500"):
        module.run_audit(replicates=100)
