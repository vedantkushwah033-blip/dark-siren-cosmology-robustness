from __future__ import annotations

import sys
from pathlib import Path

import numpy as np
import pytest
from scipy.stats import norm, t

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from audit_expanded_prior_100rep_pair import posterior, t_interval


def test_uniform_posterior_quantiles_on_expanded_support():
    grid = np.linspace(30.0, 110.0, 561)
    result = posterior(grid, np.zeros_like(grid), 30.0, 110.0)
    assert result["median"] == pytest.approx(70.0, abs=1e-10)
    assert result["q025"] == pytest.approx(32.0, abs=1e-10)
    assert result["q16"] == pytest.approx(42.8, abs=1e-10)
    assert result["q84"] == pytest.approx(97.2, abs=1e-10)
    assert result["q975"] == pytest.approx(108.0, abs=1e-10)


def test_gaussian_posterior_quantiles_match_continuous_reference():
    grid = np.linspace(30.0, 110.0, 561)
    sigma = 5.0
    logl = -0.5 * ((grid - 70.0) / sigma) ** 2
    result = posterior(grid, logl, 30.0, 110.0)
    for key, probability in [("q025", 0.025), ("q16", 0.16), ("median", 0.5),
                             ("q84", 0.84), ("q975", 0.975)]:
        assert result[key] == pytest.approx(norm.ppf(probability, loc=70.0, scale=sigma), abs=0.02)


def test_posterior_rejects_support_without_exact_grid_endpoints():
    grid = np.linspace(30.0, 110.0, 561)
    with pytest.raises(AssertionError, match="support endpoints"):
        posterior(grid, np.zeros_like(grid), 31.0, 90.0)


def test_paired_t_interval_matches_reference():
    values = np.array([-1.0, 0.0, 1.0, 2.0, 3.0])
    mean, se, low, high = t_interval(values)
    expected_se = float(np.std(values, ddof=1) / np.sqrt(len(values)))
    half = float(t.ppf(0.975, len(values) - 1) * expected_se)
    assert mean == pytest.approx(float(values.mean()), abs=1e-14)
    assert se == pytest.approx(expected_se, abs=1e-14)
    assert low == pytest.approx(mean - half, abs=1e-14)
    assert high == pytest.approx(mean + half, abs=1e-14)


def test_paired_t_interval_constant_values_has_zero_width():
    mean, se, low, high = t_interval(np.full(20, 2.5))
    assert mean == pytest.approx(2.5)
    assert se == pytest.approx(0.0)
    assert low == pytest.approx(2.5)
    assert high == pytest.approx(2.5)
