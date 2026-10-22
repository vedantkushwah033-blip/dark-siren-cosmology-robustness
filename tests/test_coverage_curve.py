import numpy as np
import pytest

from scripts.run_coverage_curve import equal_tail_interval


def test_equal_tail_interval_uniform_density():
    grid = np.linspace(50.0, 90.0, 4001)
    posterior = np.ones_like(grid) / (grid[-1] - grid[0])
    low, high = equal_tail_interval(grid, posterior, 0.68)
    assert low == pytest.approx(56.4, abs=0.02)
    assert high == pytest.approx(83.6, abs=0.02)


def test_equal_tail_interval_central_level_is_nested():
    grid = np.linspace(50.0, 90.0, 4001)
    posterior = np.exp(-0.5 * ((grid - 70.0) / 4.0) ** 2)
    low50, high50 = equal_tail_interval(grid, posterior, 0.50)
    low95, high95 = equal_tail_interval(grid, posterior, 0.95)
    assert low95 < low50 < high50 < high95


def test_equal_tail_interval_rejects_invalid_density():
    grid = np.linspace(50.0, 90.0, 5)
    with pytest.raises(ValueError, match="normalization"):
        equal_tail_interval(grid, np.zeros_like(grid), 0.68)
