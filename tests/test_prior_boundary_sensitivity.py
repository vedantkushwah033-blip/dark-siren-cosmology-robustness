import numpy as np
import pytest

from scripts.run_prior_boundary_sensitivity import cdf_from_density, interval


def test_cdf_normalizes_nonuniform_density():
    grid = np.linspace(30.0, 120.0, 1001)
    density = np.ones_like(grid)
    cdf = cdf_from_density(grid, density)
    assert np.isclose(cdf[0], 0.0)
    assert np.isclose(cdf[-1], 1.0)
    assert np.all(np.diff(cdf) >= 0.0)


def test_cdf_rejects_nonpositive_posterior():
    grid = np.linspace(30.0, 120.0, 11)
    with pytest.raises(ValueError, match="invalid posterior normalization"):
        cdf_from_density(grid, np.zeros_like(grid))


def test_central_interval_order_and_containment():
    grid = np.linspace(30.0, 120.0, 901)
    density = np.exp(-0.5 * ((grid - 70.0) / 5.0) ** 2)
    cdf = cdf_from_density(grid, density)
    lo68, hi68 = interval(grid, cdf, 0.68)
    lo95, hi95 = interval(grid, cdf, 0.95)
    assert lo95 < lo68 < 70.0 < hi68 < hi95
    assert hi95 - lo95 > hi68 - lo68


def test_interval_rejects_bad_cdf_by_construction():
    grid = np.array([0.0, 1.0, 2.0])
    cdf = np.array([0.0, 0.8, 1.0])
    lo, hi = interval(grid, cdf, 0.68)
    assert lo < hi
