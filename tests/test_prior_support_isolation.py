import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_prior_support_isolation.py"
SPEC = importlib.util.spec_from_file_location("prior_support_isolation", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_prior_support_audit_has_expected_rows_and_finite_metrics():
    rows, summary = module.run_audit(replicates=100, seed=17, grid_size=1001)
    assert len(rows) == 7 * 2 * 100
    assert len(summary) == 7 * 2
    assert {r["prior"] for r in summary} == {"narrow_50_90", "broad_30_120"}
    for row in summary:
        assert row["replicates"] == 100
        for key in ("coverage_68", "coverage_95", "coverage_68_mc_se",
                    "coverage_95_mc_se", "rmse_median", "mean_width_68",
                    "mean_width_95", "boundary_touch_rate_68", "boundary_touch_rate_95"):
            assert np.isfinite(row[key])
        assert 0.0 <= row["coverage_68"] <= 1.0
        assert 0.0 <= row["coverage_95"] <= 1.0
        assert row["mean_width_95"] > row["mean_width_68"]


def test_same_simulated_distance_is_compared_under_both_priors():
    rows, _ = module.run_audit(replicates=100, seed=91, grid_size=1001)
    by_key = {}
    for row in rows:
        key = (row["true_H0"], row["replicate"])
        by_key.setdefault(key, {})[row["prior"]] = row
    assert len(by_key) == 7 * 100
    for pair in by_key.values():
        assert pair["narrow_50_90"]["observed_distance_mpc"] == pair["broad_30_120"]["observed_distance_mpc"]
