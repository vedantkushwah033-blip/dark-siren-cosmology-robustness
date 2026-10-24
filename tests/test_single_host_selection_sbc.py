import importlib.util
from pathlib import Path
import numpy as np

PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_single_host_selection_sbc.py"
SPEC = importlib.util.spec_from_file_location("single_host_selection_sbc", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)

def test_single_host_sbc_is_reproducible_and_finite():
    rows1, summaries1 = module.run_audit(replicates=100, seed=17, grid_size=1001)
    rows2, summaries2 = module.run_audit(replicates=100, seed=17, grid_size=1001)
    assert rows1 == rows2
    assert summaries1 == summaries2
    assert len(rows1) == 100
    assert {r["method"] for r in summaries1} == {"corrected", "uncorrected"}
    assert all(0 < r["observed_distance_mpc"] < module.DETECTION_MAX for r in rows1)
    assert all(0 <= r["corrected_rank"] <= 1 and 0 <= r["uncorrected_rank"] <= 1 for r in rows1)
    numeric = ("mean_rank", "rank_variance", "ks_statistic", "ks_pvalue",
               "decile_chi_square", "decile_chi_square_pvalue",
               "beta_min_on_prior_grid", "beta_max_on_prior_grid")
    assert all(np.isfinite(float(row[key])) for row in summaries1 for key in numeric)
