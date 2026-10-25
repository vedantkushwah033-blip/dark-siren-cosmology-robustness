import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.special import ndtr
from scipy.stats import norm

from dark_siren_robustness.selection import (
    detection_efficiency_h0,
    posterior_h0_selected,
    posterior_h0_selected_with_missing_hosts,
)


H0 = np.linspace(40.0, 100.0, 6001)
DISTANCE = 299792.458 * 0.04 / 70.0


def test_zero_missing_mass_matches_existing_selected_posterior():
    baseline = posterior_h0_selected(
        H0, [0.02, 0.04], DISTANCE, 15.0, 300.0,
        weights=[0.7, 0.3], prior_min=40.0, prior_max=100.0,
    )
    combined = posterior_h0_selected_with_missing_hosts(
        H0, [0.02, 0.04], DISTANCE, 15.0, 300.0,
        catalog_weights=[0.7, 0.3], missing_weight=0.0,
        prior_min=40.0, prior_max=100.0,
    )
    np.testing.assert_allclose(combined, baseline, rtol=1e-12, atol=1e-12)


def test_complete_population_matches_explicit_full_mixture():
    z_catalog = [0.02, 0.04]
    z_missing = [0.06]
    catalog_weight = 0.65
    missing_weight = 0.35
    baseline = posterior_h0_selected(
        H0, z_catalog + z_missing, DISTANCE, 15.0, 300.0,
        weights=[0.65 * 0.7, 0.65 * 0.3, 0.35],
        prior_min=40.0, prior_max=100.0,
    )
    combined = posterior_h0_selected_with_missing_hosts(
        H0, z_catalog, DISTANCE, 15.0, 300.0,
        catalog_weights=[0.7, 0.3], missing_weight=missing_weight,
        missing_redshifts=z_missing, missing_weights=[1.0],
        prior_min=40.0, prior_max=100.0,
    )
    np.testing.assert_allclose(combined, baseline, rtol=1e-12, atol=1e-12)
    assert np.sum(np.r_[0.65 * np.array([0.7, 0.3]), 0.35]) == pytest.approx(1.0)


def test_selection_efficiency_uses_full_population_weights():
    beta = detection_efficiency_h0(
        H0, [0.02, 0.04, 0.06], 15.0, 300.0,
        weights=[0.455, 0.195, 0.35],
    )
    beta_manual = (
        0.455 * detection_efficiency_h0(H0, [0.02], 15.0, 300.0)
        + 0.195 * detection_efficiency_h0(H0, [0.04], 15.0, 300.0)
        + 0.35 * detection_efficiency_h0(H0, [0.06], 15.0, 300.0)
    )
    np.testing.assert_allclose(beta, beta_manual, rtol=1e-12, atol=1e-12)


@pytest.mark.parametrize("kwargs", [
    {"missing_weight": 0.2},
    {"missing_weight": -0.1, "missing_redshifts": [0.06]},
    {"missing_weight": 1.0, "missing_redshifts": [0.06]},
    {"missing_weight": 0.2, "missing_redshifts": [-0.06]},
    {"missing_weight": 0.2, "missing_redshifts": [0.06], "missing_weights": [0.0]},
])
def test_invalid_missing_component_rejected(kwargs):
    base = {"H0_grid": H0, "catalog_redshifts": [0.02, 0.04],
            "observed_distance": DISTANCE, "gw_sigma": 15.0,
            "detection_distance_max": 300.0, "prior_min": 40.0,
            "prior_max": 100.0}
    with pytest.raises(ValueError):
        posterior_h0_selected_with_missing_hosts(**base, **kwargs)


def test_no_event_selection_limit_matches_catalogue_only_posterior():
    from dark_siren_robustness.hardened_inference import posterior_h0

    baseline = posterior_h0(
        H0, [0.02, 0.04], DISTANCE, 15.0,
        weights=[0.7, 0.3], prior_min=40.0, prior_max=100.0,
    )
    selected = posterior_h0_selected(
        H0, [0.02, 0.04], DISTANCE, 15.0, 1.0e6,
        weights=[0.7, 0.3], prior_min=40.0, prior_max=100.0,
    )
    np.testing.assert_allclose(selected, baseline, rtol=1e-10, atol=1e-12)


def test_selected_posterior_matches_independent_direct_mixture_and_quadrature():
    # Independently evaluate p(d | H0) / P(detection | H0), then normalize
    # with trapezoidal quadrature over the prior support.
    z = np.array([0.02, 0.04])
    weights = np.array([0.7, 0.3])
    sigma = 15.0
    lower, upper = 50.0, 260.0
    grid = H0
    means = 299792.458 * z[None, :] / grid[:, None]
    likelihood = np.sum(
        weights[None, :] * norm.pdf(DISTANCE, loc=means, scale=sigma), axis=1
    )
    beta_hosts = ndtr((upper - means) / sigma) - ndtr((lower - means) / sigma)
    beta = beta_hosts @ weights
    support = (grid >= 40.0) & (grid <= 100.0)
    expected = np.zeros_like(grid)
    expected[support] = likelihood[support] / beta[support]
    expected[support] /= trapezoid(expected[support], grid[support])

    actual = posterior_h0_selected(
        grid, z, DISTANCE, sigma, upper, weights=weights,
        prior_min=40.0, prior_max=100.0, detection_distance_min=lower,
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=1e-12)
