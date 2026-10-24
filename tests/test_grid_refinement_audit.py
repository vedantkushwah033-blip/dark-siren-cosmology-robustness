import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_grid_refinement_audit.py"
SPEC = importlib.util.spec_from_file_location("grid_refinement_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_grid_refinement_output_is_complete_and_finite():
    runs, cells = module.run_grid_audit(replicates=100, seed=4)
    assert len(runs) == 3 * 3 * 2 * 2 * 100 * 2 * 3
    assert len(cells) == 3 * 3 * 2 * 2
    assert all(np.isfinite(row["posterior_median"]) for row in runs)
    assert all(row["max_median_abs_median_diff_1001_vs_4001"] >= 0 for row in cells)


def test_finer_grid_typically_changes_median_less():
    _, cells = module.run_grid_audit(replicates=100, seed=5)
    coarse = np.median([r["median_median_abs_median_diff_1001_vs_4001"] for r in cells])
    fine = np.median([r["median_median_abs_median_diff_2001_vs_4001"] for r in cells])
    assert fine <= coarse
