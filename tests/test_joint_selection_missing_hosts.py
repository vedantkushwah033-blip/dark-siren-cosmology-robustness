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



def test_missing_host_zero_weight_rejects_unused_missing_weights():
    with pytest.raises(ValueError, match="missing_weights were supplied but missing_weight is zero"):
        posterior_h0_selected_with_missing_hosts(
            H0, [0.02, 0.04], DISTANCE, 15.0, 300.0,
            catalog_weights=[0.7, 0.3], missing_weight=0.0,
            missing_weights=[1.0], prior_min=40.0, prior_max=100.0,
        )


def test_missing_host_mixture_is_invariant_to_rescaling_component_weights():
    kwargs = dict(
        H0_grid=H0, catalog_redshifts=[0.02, 0.04],
        observed_distance=DISTANCE, gw_sigma=15.0,
        detection_distance_max=300.0, missing_weight=0.35,
        missing_redshifts=[0.06, 0.08], prior_min=40.0, prior_max=100.0,
    )
    a = posterior_h0_selected_with_missing_hosts(
        **kwargs, catalog_weights=[7.0, 3.0], missing_weights=[2.0, 1.0]
    )
    b = posterior_h0_selected_with_missing_hosts(
        **kwargs, catalog_weights=[70.0, 30.0], missing_weights=[20.0, 10.0]
    )
    np.testing.assert_allclose(a, b, rtol=1e-12, atol=1e-12)


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


def test_detection_efficiency_matches_cdf_difference_in_negative_tail():
    # For H0=40 and z=0.001, the mean distance is above [5, 6], so both
    # standardized interval bounds are negative. Check the non-survival branch.
    # Use a two-point grid because the public helper requires a proper grid.
    h0 = np.array([40.0, 41.0])
    z = np.array([0.001])
    means = 299792.458 * z[0] / h0
    lower_z = (5.0 - means) / 1.0
    upper_z = (6.0 - means) / 1.0
    expected = ndtr(upper_z) - ndtr(lower_z)
    actual = detection_efficiency_h0(
        h0, z, 1.0, 6.0, detection_distance_min=5.0,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-12, atol=0.0)


def test_detection_efficiency_preserves_representable_positive_tail_probability():
    # Direct CDF subtraction rounds both Phi(6) and Phi(5) to values near one;
    # the survival-function identity preserves the small but representable mass.
    beta = detection_efficiency_h0(
        np.array([40.0, 100.0]), [0.0], 1.0, 6.0,
        detection_distance_min=5.0,
    )
    expected = ndtr(-5.0) - ndtr(-6.0)
    assert np.all(beta > 0.0)
    np.testing.assert_allclose(beta, expected, rtol=1e-12, atol=0.0)


def test_detection_efficiency_fails_closed_when_probability_underflows():
    # At these deliberately extreme values, both Gaussian CDF terms underflow
    # to the same value. Returning machine-tiny would fabricate nonzero support.
    with pytest.raises(ValueError, match="not representable as a positive probability"):
        detection_efficiency_h0(
            np.array([40.0, 100.0]), [100.0], 1.0, 1.0,
            detection_distance_min=0.0,
        )


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
    np.testing.assert_allclose(selected, baseline, rtol=1e-6, atol=1e-8)


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


def test_prior_predictive_sbc_is_approximately_uniform_for_matched_selection_model():
    """Prior-predictive PIT check with H0 drawn once per replicate.

    Conditional on each sampled H0, rejection sampling continues over host and
    distance draws until the event passes the same detection window used by the
    likelihood. Redrawing H0 after rejected events would instead weight the
    simulated truth distribution by detection efficiency beta(H0), breaking
    the intended prior-predictive calibration design.
    """
    rng = np.random.default_rng(20261032)
    grid = np.linspace(40.0, 100.0, 1201)
    z = np.array([0.02, 0.04])
    weights = np.array([0.7, 0.3])
    sigma = 15.0
    lower, upper = 50.0, 260.0
    pits = []
    for _ in range(300):
        # Draw the parameter once from the declared prior. Selection conditions
        # the event data given this truth; it must not resample the truth itself.
        truth = rng.uniform(40.0, 100.0)
        while True:
            host = rng.choice(len(z), p=weights)
            distance = rng.normal(299792.458 * z[host] / truth, sigma)
            if lower < distance < upper:
                break
        posterior = posterior_h0_selected(
            grid, z, distance, sigma, upper, weights=weights,
            prior_min=40.0, prior_max=100.0, detection_distance_min=lower,
        )
        cdf = np.concatenate(([0.0], np.cumsum(
            0.5 * (posterior[1:] + posterior[:-1]) * np.diff(grid)
        )))
        cdf /= cdf[-1]
        pits.append(float(np.interp(truth, grid, cdf)))
    pits = np.asarray(pits)
    assert np.all((pits >= 0.0) & (pits <= 1.0))
    # Broad deterministic tolerances guard gross generator/likelihood mismatch
    # without pretending 300 draws can establish exact rank uniformity.
    assert abs(float(np.mean(pits)) - 0.5) < 0.06
    counts, _ = np.histogram(pits, bins=np.linspace(0.0, 1.0, 11))
    assert np.max(np.abs(counts / len(pits) - 0.1)) < 0.07


def test_missing_host_selected_posterior_matches_independent_full_mixture_quadrature():
    """Independently reconstruct both numerator and selection denominator.

    This deliberately avoids calling either production posterior helper for the
    expected result, so a catalogue-only normalization or omitted missing-host
    contribution cannot cancel out unnoticed.
    """
    from scipy.stats import norm

    grid = np.linspace(40.0, 100.0, 6001)
    catalog_z = np.array([0.02, 0.04])
    missing_z = np.array([0.06, 0.08])
    catalog_weights = np.array([0.8, 0.2])
    missing_weights = np.array([0.25, 0.75])
    missing_mass = 0.35
    sigma = 18.0
    lower, upper = 45.0, 270.0
    observed = 180.0

    full_z = np.concatenate((catalog_z, missing_z))
    full_weights = np.concatenate((
        (1.0 - missing_mass) * catalog_weights,
        missing_mass * missing_weights,
    ))
    means = 299792.458 * full_z[None, :] / grid[:, None]
    numerator = np.sum(
        full_weights[None, :] * norm.pdf(observed, loc=means, scale=sigma),
        axis=1,
    )
    host_efficiency = ndtr((upper - means) / sigma) - ndtr((lower - means) / sigma)
    beta = host_efficiency @ full_weights
    expected = numerator / beta
    expected /= trapezoid(expected, grid)

    actual = posterior_h0_selected_with_missing_hosts(
        grid, catalog_z, observed, sigma, upper,
        catalog_weights=catalog_weights, missing_weight=missing_mass,
        missing_redshifts=missing_z, missing_weights=missing_weights,
        prior_min=40.0, prior_max=100.0, detection_distance_min=lower,
    )
    np.testing.assert_allclose(actual, expected, rtol=2e-10, atol=1e-12)


@pytest.mark.parametrize("distance", [0.0, 45.0, 270.0, 300.0, np.nan, np.inf])
def test_selected_posterior_rejects_observations_outside_open_detection_interval(distance):
    with pytest.raises(ValueError):
        posterior_h0_selected(
            H0, [0.02, 0.04], distance, 15.0, 270.0,
            weights=[0.7, 0.3], prior_min=40.0, prior_max=100.0,
            detection_distance_min=45.0,
        )


def test_detection_efficiency_and_selected_posterior_are_invariant_to_weight_scale():
    """Only relative host weights may affect a finite-mixture posterior."""
    z = np.array([0.015, 0.035, 0.07])
    weights = np.array([2.0, 3.0, 5.0])
    scaled_weights = 17.0 * weights
    lower, upper = 40.0, 250.0
    observed = 160.0

    beta = detection_efficiency_h0(
        H0, z, 20.0, upper, weights=weights,
        detection_distance_min=lower,
    )
    beta_scaled = detection_efficiency_h0(
        H0, z, 20.0, upper, weights=scaled_weights,
        detection_distance_min=lower,
    )
    np.testing.assert_allclose(beta_scaled, beta, rtol=2e-14, atol=2e-14)

    posterior = posterior_h0_selected(
        H0, z, observed, 20.0, upper, weights=weights,
        prior_min=40.0, prior_max=100.0, detection_distance_min=lower,
    )
    posterior_scaled = posterior_h0_selected(
        H0, z, observed, 20.0, upper, weights=scaled_weights,
        prior_min=40.0, prior_max=100.0, detection_distance_min=lower,
    )
    np.testing.assert_allclose(posterior_scaled, posterior, rtol=2e-13, atol=2e-14)


def test_zero_missing_mass_reduces_to_catalogue_only_even_if_missing_support_is_supplied():
    """With zero missing mass, latent missing-host support contributes no mass."""
    baseline = posterior_h0_selected(
        H0, [0.02, 0.04], DISTANCE, 15.0, 300.0,
        weights=[0.7, 0.3], prior_min=40.0, prior_max=100.0,
    )
    with_zero_mass_hosts = posterior_h0_selected_with_missing_hosts(
        H0, [0.02, 0.04], DISTANCE, 15.0, 300.0,
        catalog_weights=[0.7, 0.3], missing_weight=0.0,
        missing_redshifts=[0.06, 0.08], prior_min=40.0, prior_max=100.0,
    )
    np.testing.assert_allclose(with_zero_mass_hosts, baseline, rtol=1e-12, atol=1e-12)



def test_normal_interval_probability_preserves_representable_negative_tail_mass():
    # Direct subtraction of two tiny CDF values is poorly conditioned here.
    # Compare with a log-CDF/expm1 reference that retains the relative gap.
    from scipy.special import log_ndtr
    from dark_siren_robustness.selection import _normal_interval_probability

    lower = np.array([-37.500001])
    upper = np.array([-37.5])
    log_lower = log_ndtr(lower)
    log_upper = log_ndtr(upper)
    expected = np.exp(log_upper) * (-np.expm1(log_lower - log_upper))
    actual = _normal_interval_probability(lower, upper)
    naive = ndtr(-lower) - ndtr(-upper)
    assert expected[0] > 0.0
    assert actual[0] > 0.0
    # Prove this case distinguishes the stable implementation from the old
    # direct survival-function subtraction, not merely that both are positive.
    assert abs((naive[0] - expected[0]) / expected[0]) > 1e-11
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=0.0)


def test_normal_interval_probability_rejects_nonfinite_bounds():
    from dark_siren_robustness.selection import _normal_interval_probability

    with pytest.raises(ValueError, match="bounds must be finite"):
        _normal_interval_probability(np.array([-np.inf]), np.array([0.0]))


def test_out_of_prior_underflow_does_not_invalidate_selected_posterior():
    # H0=1 is far outside the prior and makes the Gaussian detection
    # probability underflow, while the prior-supported H0 values have valid
    # positive efficiencies. The posterior must ignore the irrelevant point.
    grid = np.linspace(1.0, 100.0, 100)
    observed_distance = 299792.458 * 0.04 / 70.0
    posterior = posterior_h0_selected(
        grid, [0.04], observed_distance, 15.0, 260.0,
        prior_min=50.0, prior_max=90.0,
    )
    support = (grid >= 50.0) & (grid <= 90.0)
    assert np.all(np.isfinite(posterior))
    assert np.all(posterior[~support] == 0.0)
    assert trapezoid(posterior[support], grid[support]) == pytest.approx(1.0)


def test_normal_interval_probability_preserves_representable_positive_tail_mass():
    # Symmetric counterpart to the negative-tail regression: direct survival
    # subtraction loses relative precision for close bounds far above zero.
    from scipy.special import log_ndtr
    from dark_siren_robustness.selection import _normal_interval_probability

    lower = np.array([37.5])
    upper = np.array([37.500001])
    log_lower_survival = log_ndtr(-lower)
    log_upper_survival = log_ndtr(-upper)
    expected = np.exp(log_lower_survival) * (
        -np.expm1(log_upper_survival - log_lower_survival)
    )
    actual = _normal_interval_probability(lower, upper)
    assert expected[0] > 0.0
    assert actual[0] > 0.0
    np.testing.assert_allclose(actual, expected, rtol=2e-12, atol=0.0)
