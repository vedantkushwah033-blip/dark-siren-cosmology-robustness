import importlib.util
from pathlib import Path

import numpy as np

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_sbc_rank_audit.py"
SPEC = importlib.util.spec_from_file_location("sbc_rank_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_sbc_ranks_are_finite_and_in_unit_interval():
    rows, summary = module.run_sbc(replicates=200, seed=44, grid_size=1001)
    assert len(rows) == 400
    assert len(summary) == 2
    assert {row["prior"] for row in summary} == {"narrow_50_90", "broad_30_120"}
    ranks = np.array([row["posterior_cdf_rank"] for row in rows])
    assert np.all(np.isfinite(ranks))
    assert np.all((ranks >= 0.0) & (ranks <= 1.0))
    for row in summary:
        assert sum(row[f"decile_{i}_count"] for i in range(1, 11)) == 200
        assert np.isfinite(row["ks_statistic"])
        assert 0.0 <= row["ks_pvalue"] <= 1.0


def test_posterior_cdf_rank_for_symmetric_normal_is_near_half_at_mean():
    grid = np.linspace(-10.0, 10.0, 4001)
    density = np.exp(-0.5 * grid**2)
    assert abs(module.posterior_cdf_at(grid, density, 0.0) - 0.5) < 1e-4
