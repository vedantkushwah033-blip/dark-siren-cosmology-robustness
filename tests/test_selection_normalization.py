import numpy as np
from scipy.integrate import trapezoid
import pytest
from scipy.integrate import quad
from scipy.optimize import brentq
from scipy.special import ndtr
from scipy.stats import norm

from dark_siren_robustness.hardened_inference import C_KM_S, posterior_h0, posterior_summary
from dark_siren_robustness.selection import detection_efficiency_h0, posterior_h0_selected


def test_single_host_detection_efficiency_matches_analytic_cdf():
    grid = np.linspace(50.0, 90.0, 401)
    z, sigma, threshold = 0.03, 12.0, 150.0
    expected = ndtr((threshold - C_KM_S * z / grid) / sigma) - ndtr((0.0 - C_KM_S * z / grid) / sigma)
    actual = detection_efficiency_h0(grid, [z], sigma, threshold)
    assert actual == pytest.approx(expected, abs=1e-14)


def test_mixture_detection_efficiency_is_weighted_host_sum():
    grid = np.linspace(50.0, 90.0, 401)
    z = np.array([0.02, 0.04])
    weights = np.array([0.7, 0.3])
    sigma, threshold = 10.0, 180.0
    expected = (
        weights[0] * (ndtr((threshold - C_KM_S * z[0] / grid) / sigma) - ndtr(-C_KM_S * z[0] / grid / sigma))
        + weights[1] * (ndtr((threshold - C_KM_S * z[1] / grid) / sigma) - ndtr(-C_KM_S * z[1] / grid / sigma))
    )
    assert detection_efficiency_h0(grid, z, sigma, threshold, weights) == pytest.approx(expected, abs=1e-14)


def test_detection_efficiency_approaches_one_for_very_high_threshold():
    grid = np.linspace(50.0, 90.0, 401)
    beta = detection_efficiency_h0(grid, [0.02, 0.04], 10.0, 1e6, [0.7, 0.3])
    assert np.all(beta <= 1.0)
    assert np.all(beta > 1.0 - 1e-12)


def test_selected_posterior_reduces_to_baseline_when_selection_is_negligible():
    grid = np.linspace(25.0, 125.0, 4001)
    z = [0.02, 0.04]
    distance, sigma, threshold = 150.0, 10.0, 1e6
    baseline = posterior_h0(grid, z, distance, sigma, weights=[0.7, 0.3])
    selected = posterior_h0_selected(grid, z, distance, sigma, threshold, weights=[0.7, 0.3])
    support = (grid >= 50.0) & (grid <= 90.0)
    assert selected[support] == pytest.approx(baseline[support], abs=1e-10)


def test_selected_posterior_is_normalized_and_differs_when_selection_matters():
    grid = np.linspace(25.0, 125.0, 4001)
    z = [0.02, 0.08]
    distance, sigma, threshold = 150.0, 18.0, 160.0
    baseline = posterior_h0(grid, z, distance, sigma, weights=[0.5, 0.5])
    selected = posterior_h0_selected(grid, z, distance, sigma, threshold, weights=[0.5, 0.5])
    support = (grid >= 50.0) & (grid <= 90.0)
    assert trapezoid(selected[support], grid[support]) == pytest.approx(1.0, abs=1e-10)
    assert np.max(np.abs(selected[support] - baseline[support])) > 1e-5


@pytest.mark.parametrize("weights", [[0, 0], [-1, 2], [1, np.nan], [1]])
def test_invalid_host_weights_are_rejected(weights):
    with pytest.raises(ValueError):
        detection_efficiency_h0(np.linspace(50, 90, 11), [0.02, 0.03], 10.0, 150.0, weights)


def test_selected_inference_rejects_undetected_observation():
    with pytest.raises(ValueError, match="detection interval"):
        posterior_h0_selected(np.linspace(25, 125, 1001), [0.02], 200.0, 10.0, 150.0)


def test_selected_posterior_quantiles_match_independent_adaptive_quadrature():
    grid = np.linspace(25.0, 125.0, 16001)
    redshifts = np.array([0.02, 0.08])
    weights = np.array([0.5, 0.5])
    distance, sigma, threshold = 150.0, 18.0, 160.0
    density = posterior_h0_selected(
        grid, redshifts, distance, sigma, threshold, weights,
        prior_min=50.0, prior_max=90.0,
    )
    got = posterior_summary(grid, density)

    def unnormalized(h0):
        means = C_KM_S * redshifts / h0
        numerator = float(np.dot(weights, norm.pdf(distance, loc=means, scale=sigma)))
        beta = float(np.dot(weights, ndtr((threshold - means) / sigma) - ndtr((0.0 - means) / sigma)))
        return numerator / beta

    normalizer, _ = quad(unnormalized, 50.0, 90.0, epsabs=1e-12, epsrel=1e-10)
    def cdf(h0):
        if h0 <= 50.0:
            return 0.0
        if h0 >= 90.0:
            return 1.0
        area, _ = quad(unnormalized, 50.0, h0, epsabs=1e-12, epsrel=1e-10)
        return area / normalizer

    for name, probability in (("lower_68", 0.16), ("median", 0.5), ("upper_68", 0.84),
                              ("lower_95", 0.025), ("upper_95", 0.975)):
        reference = brentq(lambda h0: cdf(h0) - probability, 50.0, 90.0)
        assert got[name] == pytest.approx(reference, abs=0.003)


def test_detection_efficiency_matches_finite_positive_distance_interval():
    grid = np.linspace(50.0, 90.0, 101)
    z, sigma, lower, upper = 0.001, 100.0, 20.0, 180.0
    means = C_KM_S * z / grid
    expected = ndtr((upper - means) / sigma) - ndtr((lower - means) / sigma)
    actual = detection_efficiency_h0(grid, [z], sigma, upper, detection_distance_min=lower)
    assert actual == pytest.approx(expected, abs=1e-14)


def test_selected_inference_rejects_observation_below_lower_detection_bound():
    with pytest.raises(ValueError, match="detection interval"):
        posterior_h0_selected(
            np.linspace(25, 125, 1001), [0.02], 10.0, 10.0, 150.0,
            detection_distance_min=20.0,
        )
