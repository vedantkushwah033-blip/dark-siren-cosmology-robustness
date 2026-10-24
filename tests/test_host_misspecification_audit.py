import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_host_misspecification_audit.py"
SPEC = importlib.util.spec_from_file_location("host_misspecification_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_generated_host_support_matches_declared_scenarios():
    rows = module.generate_dataset(200, 17)
    assert len(rows) == 800
    for row in rows:
        assert row["observed_distance_mpc"] > 0
        assert row["sigma_mpc"] > 0
        if row["scenario"] in ("correct_weights", "reversed_weights"):
            assert row["true_host_redshift"] in (0.02, 0.04)
        else:
            assert row["true_host_redshift"] in (0.02, 0.03, 0.04)


def test_stress_test_metrics_are_finite_and_bounded():
    rows, summaries = module.run_stress_test(replicates=200, seed=18)
    assert len(rows) == 800
    assert len(summaries) == 4
    for row in summaries:
        for key in ("mean_median_error", "rmse_median", "coverage_68", "coverage_95",
                    "mean_width_68", "mean_width_95"):
            assert np.isfinite(row[key])
        assert row["rmse_median"] >= 0
        assert 0 <= row["coverage_68"] <= 1
        assert 0 <= row["coverage_95"] <= 1
        assert row["mean_width_68"] > 0
        assert row["mean_width_95"] > row["mean_width_68"]
