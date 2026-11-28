"""Independent deterministic checks for the hidden-host posterior implementation.

The quadrature reference below is intentionally separate from the production
trapezoidal-grid calculation: it integrates the continuous likelihood with
scipy.integrate.quad and solves quantiles with scipy.optimize.brentq.
"""
import numpy as np
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import logsumexp
from scipy.stats import norm

from dark_siren_robustness.hardened_inference import (
    C_KM_S,
    posterior_h0,
    posterior_summary,
)


def _reference_log_likelihood(h0, z, weights, distance, sigma, *,
                              missing_fraction=0.0, missing_z=(),
                              missing_weights=(), positive_truncated=True):
    z = np.asarray(z, dtype=float)
    weights = np.asarray(weights, dtype=float)
    weights = weights / weights.sum()
    means = C_KM_S * z / h0
    component = norm.logpdf(distance, loc=means, scale=sigma)
    if positive_truncated:
        component = component - norm.logcdf(means / sigma)
    terms = np.log(weights) + component
    catalog = logsumexp(terms) + np.log1p(-missing_fraction)
    if missing_fraction <= 0:
        return float(catalog)
    mz = np.asarray(missing_z, dtype=float)
    mw = np.asarray(missing_weights, dtype=float)
    mw = mw / mw.sum()
    mmeans = C_KM_S * mz / h0
    mcomponent = norm.logpdf(distance, loc=mmeans, scale=sigma)
    if positive_truncated:
        mcomponent = mcomponent - norm.logcdf(mmeans / sigma)
    missing = logsumexp(np.log(mw) + mcomponent) + np.log(missing_fraction)
    return float(np.logaddexp(catalog, missing))


def _continuous_reference(z, weights, distance, sigma, prior_min, prior_max,
                          *, missing_fraction=0.0, missing_z=(),
                          missing_weights=()):
    midpoint = 0.5 * (prior_min + prior_max)
    shift = _reference_log_likelihood(
        midpoint, z, weights, distance, sigma,
        missing_fraction=missing_fraction, missing_z=missing_z,
        missing_weights=missing_weights,
    )

    def unnormalized(h0):
        return np.exp(_reference_log_likelihood(
            h0, z, weights, distance, sigma,
            missing_fraction=missing_fraction, missing_z=missing_z,
            missing_weights=missing_weights,
        ) - shift)

    normalizer = quad(
        unnormalized, prior_min, prior_max,
        epsabs=1e-11, epsrel=1e-10, limit=300,
    )[0]

    def cdf(h0):
        return quad(
            unnormalized, prior_min, h0,
            epsabs=1e-11, epsrel=1e-10, limit=300,
        )[0] / normalizer

    probabilities = (0.025, 0.16, 0.50, 0.84, 0.975)
    quantiles = [
        brentq(lambda x: cdf(x) - probability,
               prior_min, prior_max, xtol=1e-10)
        for probability in probabilities
    ]
    mean = quad(
        lambda h0: h0 * unnormalized(h0),
        prior_min, prior_max,
        epsabs=1e-9, epsrel=1e-10, limit=300,
    )[0] / normalizer
    return {
        "mean": mean,
        "lower_95": quantiles[0],
        "lower_68": quantiles[1],
        "median": quantiles[2],
        "upper_68": quantiles[3],
        "upper_95": quantiles[4],
    }


@pytest.mark.parametrize(
    "z,weights,distance,sigma,missing_fraction,missing_z,missing_weights",
    [
        ([0.03], [1.0], 130.0, 18.0, 0.0, (), ()),
        ([0.02, 0.04], [0.7, 0.3], 160.0, 20.0, 0.0, (), ()),
        ([0.02, 0.04], [0.7, 0.3], 180.0, 18.0, 0.30, [0.03], [1.0]),
        ([0.02, 0.04], [0.7, 0.3], 230.0, 9.0, 0.45, [0.035], [1.0]),
    ],
)
def test_grid_posterior_matches_independent_adaptive_quadrature(
    z, weights, distance, sigma, missing_fraction, missing_z, missing_weights
):
    grid = np.linspace(25.0, 125.0, 4001)
    kwargs = {}
    if missing_fraction > 0:
        kwargs = {
            "missing_weight": missing_fraction,
            "missing_redshifts": missing_z,
            "missing_weights": missing_weights,
        }
    density = posterior_h0(
        grid, z, distance, sigma, weights=weights,
        prior_min=50.0, prior_max=90.0,
        positive_truncated=True, **kwargs,
    )
    actual = posterior_summary(grid, density)
    expected = _continuous_reference(
        z, weights, distance, sigma, 50.0, 90.0,
        missing_fraction=missing_fraction,
        missing_z=missing_z, missing_weights=missing_weights,
    )
    for key in expected:
        assert actual[key] == pytest.approx(expected[key], abs=5e-4)


def test_zero_hidden_mass_is_exact_catalogue_only_reduction():
    grid = np.linspace(25.0, 125.0, 4001)
    args = (grid, [0.02, 0.04], 175.0, 18.0)
    base = posterior_h0(
        *args, weights=[0.7, 0.3], prior_min=50.0, prior_max=90.0,
        positive_truncated=True,
    )
    zero_mass = posterior_h0(
        *args, weights=[0.7, 0.3], missing_weight=0.0,
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    np.testing.assert_allclose(zero_mass, base, rtol=0.0, atol=1e-14)


def test_duplicate_component_splitting_does_not_change_posterior():
    grid = np.linspace(25.0, 125.0, 4001)
    base = posterior_h0(
        grid, [0.02, 0.04], 160.0, 20.0, weights=[0.7, 0.3],
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    duplicated = posterior_h0(
        grid, [0.02, 0.02, 0.04], 160.0, 20.0,
        weights=[0.2, 0.5, 0.3],
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    np.testing.assert_allclose(duplicated, base, rtol=2e-13, atol=2e-13)


def test_mixture_component_permutation_does_not_change_posterior():
    grid = np.linspace(25.0, 125.0, 4001)
    base = posterior_h0(
        grid, [0.02, 0.04], 160.0, 20.0, weights=[0.7, 0.3],
        missing_weight=0.25, missing_redshifts=[0.025, 0.035],
        missing_weights=[0.4, 0.6], prior_min=50.0, prior_max=90.0,
        positive_truncated=True,
    )
    permuted = posterior_h0(
        grid, [0.04, 0.02], 160.0, 20.0, weights=[0.3, 0.7],
        missing_weight=0.25, missing_redshifts=[0.035, 0.025],
        missing_weights=[0.6, 0.4], prior_min=50.0, prior_max=90.0,
        positive_truncated=True,
    )
    np.testing.assert_allclose(permuted, base, rtol=2e-13, atol=2e-13)


def test_one_host_is_invariant_to_splitting_identical_host_mass_between_components():
    grid = np.linspace(25.0, 125.0, 4001)
    single = posterior_h0(
        grid, [0.03], 130.0, 18.0, weights=[1.0],
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    same_host_as_catalogue_and_hidden = posterior_h0(
        grid, [0.03], 130.0, 18.0, weights=[1.0],
        missing_weight=0.4, missing_redshifts=[0.03],
        prior_min=50.0, prior_max=90.0, positive_truncated=True,
    )
    np.testing.assert_allclose(
        same_host_as_catalogue_and_hidden, single, rtol=2e-13, atol=2e-13
    )
