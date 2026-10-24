import importlib.util
from pathlib import Path

import numpy as np
import pytest

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


def test_posterior_cdf_excludes_zero_density_points_outside_prior_support():
    grid = np.linspace(0.0, 10.0, 101)
    density = np.where((grid >= 3.0) & (grid <= 8.0), 1.0, 0.0)
    # A constant density on [3, 8] has CDF(5) = (5 - 3) / (8 - 3) = 0.4.
    # Including exterior zero bins in trapezoidal integration creates artificial
    # boundary half-bins and shifts this value.
    assert module.posterior_cdf_at(grid, density, 5.0) == pytest.approx(0.4, abs=1e-12)


def test_posterior_cdf_rejects_malformed_inputs():
    with pytest.raises(ValueError, match="same-length"):
        module.posterior_cdf_at([0.0, 1.0], [1.0], 0.5)
    with pytest.raises(ValueError, match="strictly increasing"):
        module.posterior_cdf_at([0.0, 0.0], [1.0, 1.0], 0.0)
