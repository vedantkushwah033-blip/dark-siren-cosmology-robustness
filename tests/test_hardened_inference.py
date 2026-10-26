import numpy as np
import pytest
from scipy.stats import norm

from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary


def test_valid_posterior_is_finite_and_normalized():
    grid = np.linspace(50, 90, 1601)
    posterior = posterior_h0(grid, [0.02, 0.03, 0.04], 150.0, 20.0)
    assert np.all(np.isfinite(posterior))
    assert np.all(posterior >= 0)
    assert np.trapezoid(posterior, grid) == pytest.approx(1.0, abs=1e-10)


@pytest.mark.parametrize("grid", [[70.0], [70.0, 69.0], [70.0, np.nan], [[60.0, 70.0]]])
def test_invalid_h0_grid_rejected(grid):
    with pytest.raises(ValueError):
        posterior_h0(grid, [0.02], 100.0, 10.0)


@pytest.mark.parametrize("sigma", [0.0, -1.0, np.nan, np.inf])
def test_invalid_distance_uncertainty_rejected(sigma):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], [0.02], 100.0, sigma)


@pytest.mark.parametrize("z", [[np.nan], [np.inf], [-0.1], [[0.1]]])
def test_invalid_redshift_inputs_rejected(z):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], z, 100.0, 10.0)


@pytest.mark.parametrize("weights", [[1, np.nan], [1, np.inf], [-1, 2], [0, 0], [[1, 2]]])
def test_invalid_catalogue_weights_rejected(weights):
    with pytest.raises(ValueError):
        posterior_h0([50, 70, 90], [0.02, 0.03], 100.0, 10.0, weights=weights)


def test_summary_rejects_shape_mismatch_and_zero_area():
    with pytest.raises(ValueError, match="same length"):
        posterior_summary([50, 60, 70], [0.1, 0.2])
    with pytest.raises(ValueError, match="positive integrated area"):
        posterior_summary([50, 60, 70], [0, 0, 0])


def test_summary_normalizes_valid_density_and_returns_ordered_intervals():
    grid = np.linspace(50, 90, 401)
    density = np.exp(-0.5 * ((grid - 70) / 4) ** 2)
    summary = posterior_summary(grid, density)
    assert 50 < summary["lower_95"] < summary["lower_68"] < summary["median"]
    assert summary["median"] < summary["upper_68"] < summary["upper_95"] < 90


def test_missing_host_component_requires_support():
    with pytest.raises(ValueError, match="missing_redshifts"):
        posterior_h0([50, 70, 90], [0.02], 100.0, 10.0, missing_weight=0.2)


def test_posterior_normalization_uses_only_prior_supported_grid():
    # The prior interval is strictly inside the computational grid. The
    # trapezoidal rule must not count artificial half-bins that bridge the
    # discontinuous prior boundary to zero-valued exterior grid points.
    grid = np.linspace(25.0, 125.0, 1001)
    posterior = posterior_h0(
        grid, [0.043633490629581634, 0.06794687740863647, 0.07082449674972259],
        271.77310405878274, 7.639409914451161,
        weights=[0.28097647356639377, 0.3988538762771191, 0.32016965015648713],
        prior_min=50.0, prior_max=90.0,
    )
    support = (grid >= 50.0) & (grid <= 90.0)
    assert np.trapezoid(posterior[support], grid[support]) == pytest.approx(1.0, abs=1e-10)
    # The full-grid trapezoid would incorrectly add boundary half-bins because
    # the represented prior density jumps from nonzero to zero outside support.
    assert np.trapezoid(posterior, grid) > 1.0



def test_positive_truncated_likelihood_matches_independent_formula():
    grid = np.linspace(25.0, 125.0, 1001)
    redshifts = np.array([0.02, 0.04])
    distance, sigma = 100.0, 20.0
    weights = np.array([0.6, 0.4])
    means = 299792.458 * redshifts[None, :] / grid[:, None]
    component_density = norm.pdf(distance, loc=means, scale=sigma) / norm.cdf(means / sigma)
    expected = np.sum(weights[None, :] * component_density, axis=1)
    support = (grid >= 50.0) & (grid <= 90.0)
    expected[~support] = 0.0
    expected /= np.trapezoid(expected[support], grid[support])

    actual = posterior_h0(
        grid, redshifts, distance, sigma, weights=weights,
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    untruncated = posterior_h0(
        grid, redshifts, distance, sigma, weights=weights,
        prior_min=50.0, prior_max=90.0,
    )
    assert actual == pytest.approx(expected, abs=1e-12)
    assert not np.allclose(actual, untruncated, rtol=0.0, atol=1e-10)


def test_posterior_summary_ignores_zero_density_outside_prior_support():
    grid = np.linspace(25.0, 125.0, 1001)
    posterior = posterior_h0(
        grid, [0.043633490629581634, 0.06794687740863647, 0.07082449674972259],
        271.77310405878274, 7.639409914451161,
        weights=[0.28097647356639377, 0.3988538762771191, 0.32016965015648713],
        prior_min=50.0, prior_max=90.0,
    )
    summary = posterior_summary(grid, posterior)
    assert 50.0 <= summary["lower_95"] < summary["lower_68"] < summary["median"]
    assert summary["median"] < summary["upper_68"] < summary["upper_95"] <= 90.0
