import importlib.util
from pathlib import Path

import numpy as np
import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_fixed_truth_coverage_audit.py"
SPEC = importlib.util.spec_from_file_location("fixed_truth_coverage_audit", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_normalized_cdf_is_monotone_and_ends_at_one():
    grid = np.linspace(50.0, 90.0, 401)
    posterior = np.exp(-0.5 * ((grid - 70.0) / 5.0) ** 2)
    cdf = module.normalized_cdf(grid, posterior)
    assert cdf[0] == pytest.approx(0.0)
    assert cdf[-1] == pytest.approx(1.0)
    assert np.all(np.diff(cdf) >= 0.0)


def test_central_interval_contains_median_and_is_ordered():
    grid = np.linspace(50.0, 90.0, 1601)
    posterior = np.exp(-0.5 * ((grid - 70.0) / 4.0) ** 2)
    cdf = module.normalized_cdf(grid, posterior)
    lo68, hi68 = module.central_interval(grid, cdf, 0.68)
    lo95, hi95 = module.central_interval(grid, cdf, 0.95)
    median = float(np.interp(0.5, cdf, grid))
    assert lo95 <= lo68 < median < hi68 <= hi95


def test_invalid_posterior_is_rejected():
    grid = np.array([0.0, 1.0, 2.0])
    with pytest.raises(ValueError):
        module.normalized_cdf(grid, np.array([1.0, -0.1, 1.0]))
    with pytest.raises(ValueError):
        module.normalized_cdf(grid, np.array([0.0, 0.0, 0.0]))


def test_interval_level_has_expected_width_order():
    grid = np.linspace(50.0, 90.0, 1601)
    posterior = np.exp(-0.5 * ((grid - 70.0) / 3.0) ** 2)
    cdf = module.normalized_cdf(grid, posterior)
    lo68, hi68 = module.central_interval(grid, cdf, 0.68)
    lo95, hi95 = module.central_interval(grid, cdf, 0.95)
    assert hi95 - lo95 > hi68 - lo68
