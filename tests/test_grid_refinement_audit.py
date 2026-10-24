import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_grid_refinement_audit.py"
SPEC = importlib.util.spec_from_file_location("grid_refinement_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_grid_refinement_output_is_complete_and_finite():
    runs, cells = module.run_grid_audit(replicates=100, seed=4)
    assert len(runs) == 3 * 3 * 2 * 100 * 2 * 3
    assert len(cells) == 3 * 3 * 2 * 2
    summary_fields = ("posterior_median", "lower_68", "upper_68", "lower_95", "upper_95")
    assert all(np.isfinite(row[key]) for row in runs for key in summary_fields)
    assert all(row["lower_95"] <= row["lower_68"] <= row["posterior_median"]
               <= row["upper_68"] <= row["upper_95"] for row in runs)
    assert all(row["max_median_abs_median_diff_1001_vs_4001"] >= 0 for row in cells)

    # Each simulated observation must be held fixed across weighting and grid
    # resolution; otherwise differences would mix Monte Carlo variation with
    # numerical discretization. The observation is keyed by truth, z2, noise,
    # and replicate.
    observations = {}
    for row in runs:
        key = (row["H0_true"], row["z2"], row["noise_fraction"], row["replicate"])
        observations.setdefault(key, set()).add(row["observed_distance_mpc"])
    assert all(len(values) == 1 for values in observations.values())

    # A maximum absolute grid discrepancy cannot be smaller than its median
    # or 95th percentile within the same design cell.
    for row in cells:
        for metric in ("median_abs_median_diff_1001_vs_4001",
                       "abs_lower68_diff_1001_vs_4001", "abs_upper68_diff_1001_vs_4001",
                       "abs_lower95_diff_1001_vs_4001", "abs_upper95_diff_1001_vs_4001",
                       "median_abs_median_diff_2001_vs_4001"):
            assert row[f"max_{metric}"] >= row[f"median_{metric}"] >= 0
            assert row[f"max_{metric}"] >= row[f"p95_{metric}"] >= 0


def test_finer_grid_typically_changes_median_less():
    _, cells = module.run_grid_audit(replicates=100, seed=5)
    coarse = np.median([r["median_median_abs_median_diff_1001_vs_4001"] for r in cells])
    fine = np.median([r["median_median_abs_median_diff_2001_vs_4001"] for r in cells])
    assert fine <= coarse
