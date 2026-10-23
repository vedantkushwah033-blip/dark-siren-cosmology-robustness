"""Independent reference checks for the simplified low-z H0 likelihood.

These tests deliberately calculate the Gaussian mixture directly, without using
logsumexp or the implementation's internal helper functions. They validate
numerical implementation only; they do not validate the astrophysical model.
"""
import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0


def _direct_reference(grid, redshifts, distance, sigma, weights, prior_min, prior_max):
    """Small, transparent reference implementation of the documented likelihood."""
    grid = np.asarray(grid, dtype=float)
    z = np.asarray(redshifts, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    likelihood = np.zeros(grid.shape, dtype=float)
    for i, h0 in enumerate(grid):
        host_terms = []
        for redshift, weight in zip(z, weights):
            predicted_distance = C_KM_S * redshift / h0
            host_terms.append(weight * norm.pdf(distance, loc=predicted_distance, scale=sigma))
        likelihood[i] = sum(host_terms)
    likelihood[(grid < prior_min) | (grid > prior_max)] = 0.0
    return likelihood / trapezoid(likelihood, grid)


def test_posterior_matches_independent_direct_gaussian_mixture():
    grid = np.linspace(40.0, 100.0, 1201)
    redshifts = np.array([0.012, 0.021, 0.034])
    weights = np.array([2.0, 5.0, 3.0])
    distance, sigma = 118.0, 16.0
    prior_min, prior_max = 50.0, 90.0

    actual = posterior_h0(
        grid, redshifts, distance, sigma, weights=weights,
        prior_min=prior_min, prior_max=prior_max,
    )
    expected = _direct_reference(
        grid, redshifts, distance, sigma, weights, prior_min, prior_max,
    )

    assert np.all(np.isfinite(actual))
    assert np.all(actual >= 0.0)
    assert trapezoid(actual, grid) == pytest.approx(1.0, abs=1e-10)
    # Same mathematical model, independently evaluated in ordinary probability space.
    assert np.allclose(actual, expected, rtol=2e-10, atol=1e-14)


def test_host_permutation_does_not_change_posterior():
    grid = np.linspace(50.0, 90.0, 1001)
    z = np.array([0.011, 0.019, 0.027, 0.041])
    w = np.array([1.0, 4.0, 2.0, 3.0])
    kwargs = dict(gw_distance=145.0, gw_sigma=18.0)
    first = posterior_h0(grid, z, weights=w, **kwargs)
    order = np.array([2, 0, 3, 1])
    second = posterior_h0(grid, z[order], weights=w[order], **kwargs)
    assert np.allclose(first, second, rtol=1e-12, atol=1e-14)


def test_common_weight_rescaling_does_not_change_posterior():
    grid = np.linspace(50.0, 90.0, 1001)
    z = [0.013, 0.024, 0.038]
    kwargs = dict(gw_distance=130.0, gw_sigma=20.0)
    first = posterior_h0(grid, z, weights=[1.0, 2.0, 7.0], **kwargs)
    second = posterior_h0(grid, z, weights=[11.0, 22.0, 77.0], **kwargs)
    assert np.allclose(first, second, rtol=1e-12, atol=1e-14)


def test_single_host_likelihood_mode_matches_distance_relation():
    grid = np.linspace(50.0, 90.0, 40001)
    redshift = 0.02
    observed_distance = 299792.458 * redshift / 70.0
    posterior = posterior_h0(
        grid, [redshift], observed_distance, 1.0,
        prior_min=50.0, prior_max=90.0,
    )
    # For one host and Gaussian distance error, the likelihood is maximized
    # where predicted distance equals observed distance: H0 = c*z/d.
    expected_mode = 299792.458 * redshift / observed_distance
    actual_mode = float(grid[np.argmax(posterior)])
    assert actual_mode == pytest.approx(expected_mode, abs=0.002)


def test_posterior_is_zero_outside_prior_even_with_wide_grid():
    grid = np.linspace(30.0, 120.0, 1801)
    posterior = posterior_h0(
        grid, [0.02, 0.03], 120.0, 15.0,
        prior_min=50.0, prior_max=90.0,
    )
    outside = (grid < 50.0) | (grid > 90.0)
    assert np.all(posterior[outside] == 0.0)
    assert np.any(posterior[~outside] > 0.0)
