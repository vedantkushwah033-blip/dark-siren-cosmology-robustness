import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_multi_h0_coverage_audit.py"
SPEC = importlib.util.spec_from_file_location("multi_h0_coverage_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_multi_h0_coverage_audit_shape_and_bounds():
    rows, summaries = module.run_coverage_audit(replicates=200, seed=23)
    assert len(rows) == 7 * 3 * 2 * 200 * 2
    assert len(summaries) == 7 * 3 * 2 * 2
    assert all(np.isfinite(r["median_rmse"]) for r in summaries)
    assert all(0 <= r["coverage_68"] <= 1 for r in summaries)
    assert all(0 <= r["coverage_95"] <= 1 for r in summaries)
    assert all(r["coverage_68_mcse"] >= 0 for r in summaries)


def test_generator_and_inference_share_sigma_and_truth():
    rows, _ = module.run_coverage_audit(replicates=200, seed=29)
    for row in rows:
        assert row["sigma_mpc"] > 0
        assert module.PRIOR_MIN <= row["H0_true"] <= module.PRIOR_MAX
        assert row["observed_distance_mpc"] > 0
        assert row["inference_weighting"] in ("correct", "reversed")
