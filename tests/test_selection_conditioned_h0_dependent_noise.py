from __future__ import annotations

import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
    C_KM_S,
    dataset_posterior_h0_selection_conditioned,
    event_log_likelihood_h0_dependent_noise,
    population_detection_efficiency_h0_dependent_noise,
    selection_conditioned_dataset_log_likelihood,
)


def test_single_host_likelihood_matches_direct_heteroscedastic_normal():
    h0 = np.array([60.0, 70.0, 80.0])
    z, distance, frac, floor = 0.035, 410.0, 0.12, 8.0
    mu = C_KM_S * z / h0
    sigma = np.hypot(floor, frac * mu)
    actual = event_log_likelihood_h0_dependent_noise(
        h0, redshifts=[z], gw_distance=distance, fractional_sigma=frac,
        sigma_floor_mpc=floor,
    )
    np.testing.assert_allclose(actual, norm.logpdf(distance, loc=mu, scale=sigma), rtol=1e-13, atol=1e-13)


def test_host_mixture_matches_independent_weighted_logsumexp():
    h0, z = np.array([60.0, 70.0, 80.0]), np.array([0.02, 0.035, 0.06])
    weights = np.array([2.0, 3.0, 5.0])
    distance, frac, floor = 500.0, 0.08, 15.0
    mu = C_KM_S * z[None, :] / h0[:, None]
    sigma = np.hypot(floor, frac * mu)
    terms = norm.logpdf(distance, loc=mu, scale=sigma) + np.log(weights / weights.sum())[None, :]
    expected = np.logaddexp.reduce(terms, axis=1)
    actual = event_log_likelihood_h0_dependent_noise(
        h0, redshifts=z, gw_distance=distance, fractional_sigma=frac,
        sigma_floor_mpc=floor, weights=weights,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)


def test_positive_truncation_uses_per_host_h0_dependent_normalizer():
    h0 = np.array([60.0, 70.0, 80.0])
    z, distance, frac, floor = 0.03, 350.0, 0.20, 5.0
    mu = C_KM_S * z / h0
    sigma = np.hypot(floor, frac * mu)
    expected = norm.logpdf(distance, loc=mu, scale=sigma) - norm.logcdf(mu / sigma)
    actual = event_log_likelihood_h0_dependent_noise(
        h0, redshifts=[z], gw_distance=distance, fractional_sigma=frac,
        sigma_floor_mpc=floor, positive_truncated=True,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)


def test_population_efficiency_matches_analytic_single_host_survival():
    h0 = np.array([60.0, 70.0, 80.0])
    z, frac, floor, threshold = 0.04, 0.1, 10.0, 300.0
    mu = C_KM_S * z / h0
    sigma = np.hypot(floor, frac * mu)
    expected = norm.sf((threshold - mu) / sigma)
    actual = population_detection_efficiency_h0_dependent_noise(
        h0, population_redshifts=[z], fractional_sigma=frac,
        sigma_floor_mpc=floor, detection_distance_min=threshold,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)


def test_truncated_population_efficiency_is_conditional_survival_ratio():
    h0 = np.array([60.0, 70.0, 80.0])
    z, frac, floor, threshold = 0.04, 0.15, 10.0, 300.0
    mu = C_KM_S * z / h0
    sigma = np.hypot(floor, frac * mu)
    expected = norm.sf((threshold - mu) / sigma) / norm.sf(-mu / sigma)
    actual = population_detection_efficiency_h0_dependent_noise(
        h0, population_redshifts=[z], fractional_sigma=frac,
        sigma_floor_mpc=floor, detection_distance_min=threshold,
        positive_truncated=True,
    )
    np.testing.assert_allclose(actual, expected, rtol=1e-13, atol=1e-13)


def test_selected_dataset_uses_h0_dependent_detection_normalization():
    h0 = np.array([60.0, 70.0, 80.0])
    events = [
        {"redshifts": [0.02, 0.04], "weights": [0.6, 0.4], "gw_distance": 410.0},
        {"redshifts": [0.025, 0.05], "weights": [0.25, 0.75], "gw_distance": 450.0},
    ]
    result = selection_conditioned_dataset_log_likelihood(
        h0, events, population_redshifts=[0.015, 0.03, 0.05, 0.08],
        population_weights=[0.1, 0.3, 0.4, 0.2], detection_distance_min=300.0,
        fractional_sigma=0.15, sigma_floor_mpc=5.0,
    )
    direct_joint = sum(event_log_likelihood_h0_dependent_noise(
        h0, redshifts=e["redshifts"], weights=e["weights"], gw_distance=e["gw_distance"],
        fractional_sigma=0.15, sigma_floor_mpc=5.0,
    ) for e in events)
    assert result["unconditional_log_likelihood"] == pytest.approx(direct_joint)
    assert result["selected_log_likelihood"] == pytest.approx(
        direct_joint - len(events) * np.log(result["detection_efficiency"])
    )


def test_posterior_normalizes_including_non_grid_prior_endpoints():
    h0 = np.linspace(40.0, 100.0, 301)
    events = [{"redshifts": [0.02, 0.04], "gw_distance": 410.0},
              {"redshifts": [0.03, 0.05], "gw_distance": 470.0}]
    result = dataset_posterior_h0_selection_conditioned(
        h0, events, prior_min=50.07, prior_max=89.91,
        population_redshifts=np.linspace(0.01, 0.1, 200),
        detection_distance_min=300.0, fractional_sigma=0.12, sigma_floor_mpc=8.0,
    )
    assert result["normalization_grid"][0] == pytest.approx(50.07)
    assert result["normalization_grid"][-1] == pytest.approx(89.91)
    assert trapezoid(result["normalization_density"], result["normalization_grid"]) == pytest.approx(1.0, abs=2e-12)
    assert result["posterior_integral"] == pytest.approx(1.0, abs=2e-12)
    assert np.all(result["posterior_density"][h0 < 50.07] == 0.0)
    assert np.all(result["posterior_density"][h0 > 89.91] == 0.0)


def test_detection_efficiency_is_monte_carlo_calibrated():
    rng = np.random.default_rng(20261029)
    h0_true, frac, floor, threshold = 70.0, 0.12, 12.0, 350.0
    z = np.array([0.02, 0.035, 0.05, 0.075])
    weights = np.array([0.15, 0.35, 0.35, 0.15])
    predicted = population_detection_efficiency_h0_dependent_noise(
        [h0_true, 75.0], population_redshifts=z, population_weights=weights,
        fractional_sigma=frac, sigma_floor_mpc=floor, detection_distance_min=threshold,
    )[0]
    n = 200_000
    host_z = rng.choice(z, size=n, p=weights)
    mu = C_KM_S * host_z / h0_true
    sigma = np.hypot(floor, frac * mu)
    empirical = np.mean(rng.normal(mu, sigma) >= threshold)
    se = np.sqrt(predicted * (1.0 - predicted) / n)
    assert abs(empirical - predicted) < 5.0 * se + 0.001


def test_invalid_inputs_and_below_threshold_events_fail_closed():
    with pytest.raises(ValueError, match="H0_grid"):
        event_log_likelihood_h0_dependent_noise(
            [70.0, 60.0], redshifts=[0.03], gw_distance=400.0, fractional_sigma=0.1,
        )
    for frac in (0.0, -0.1, np.nan, np.inf):
        with pytest.raises(ValueError, match="fractional_sigma"):
            event_log_likelihood_h0_dependent_noise(
                [60.0, 70.0], redshifts=[0.03], gw_distance=400.0, fractional_sigma=frac,
            )
    for floor in (-1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="sigma_floor_mpc"):
            event_log_likelihood_h0_dependent_noise(
                [60.0, 70.0], redshifts=[0.03], gw_distance=400.0,
                fractional_sigma=0.1, sigma_floor_mpc=floor,
            )
    with pytest.raises(ValueError, match="detection_distance_min"):
        selection_conditioned_dataset_log_likelihood(
            [60.0, 70.0], [{"redshifts": [0.03], "gw_distance": 299.0}],
            population_redshifts=[0.03], detection_distance_min=300.0, fractional_sigma=0.1,
        )
    for threshold in (0.0, -1.0, np.nan, np.inf):
        with pytest.raises(ValueError, match="detection_distance_min"):
            population_detection_efficiency_h0_dependent_noise(
                [60.0, 70.0], population_redshifts=[0.03],
                detection_distance_min=threshold, fractional_sigma=0.1,
            )


def test_zero_weight_host_is_equivalent_to_omitting_host():
    h0 = np.array([60.0, 70.0, 80.0])
    all_hosts = event_log_likelihood_h0_dependent_noise(
        h0, redshifts=[0.02, 0.04], weights=[1.0, 0.0], gw_distance=400.0, fractional_sigma=0.1,
    )
    one_host = event_log_likelihood_h0_dependent_noise(
        h0, redshifts=[0.02], weights=[1.0], gw_distance=400.0, fractional_sigma=0.1,
    )
    np.testing.assert_allclose(all_hosts, one_host, rtol=0, atol=1e-14)

def _prior_predictive_rank_calibration(n_replicates=160, *, include_selection=True):
    """Small deterministic SBC diagnostic for the H0-dependent-noise toy."""
    from scipy.integrate import cumulative_trapezoid
    from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
        simulate_detected_distance_events,
    )

    redshifts = np.linspace(0.025, 0.16, 12)
    population_weights = np.linspace(1.0, 2.5, redshifts.size)
    population_weights = population_weights / population_weights.sum()
    h0_grid = np.linspace(50.0, 90.0, 161)
    ranks, cover68, cover90 = [], [], []
    for replicate in range(n_replicates):
        rng = np.random.default_rng(6172026 + replicate)
        h0_true = float(rng.uniform(50.0, 90.0))
        events = simulate_detected_distance_events(
            rng=rng, n_events=3, population_redshifts=redshifts,
            population_weights=population_weights, H0_true=h0_true,
            detection_distance_min=300.0, fractional_sigma=0.12,
            sigma_floor_mpc=8.0,
        )
        result = selection_conditioned_dataset_log_likelihood(
            h0_grid, events, population_redshifts=redshifts,
            population_weights=population_weights, detection_distance_min=300.0,
            fractional_sigma=0.12, sigma_floor_mpc=8.0,
        )
        loglike = result["selected_log_likelihood"] if include_selection else result["unconditional_log_likelihood"]
        density = np.exp(loglike - np.max(loglike))
        density /= trapezoid(density, h0_grid)
        cdf = np.r_[0.0, cumulative_trapezoid(density, h0_grid)]
        cdf /= cdf[-1]
        ranks.append(float(np.interp(h0_true, h0_grid, cdf)))
        q16, q84, q05, q95 = np.interp([0.16, 0.84, 0.05, 0.95], cdf, h0_grid)
        cover68.append(bool(q16 <= h0_true <= q84))
        cover90.append(bool(q05 <= h0_true <= q95))
    return {
        "mean_rank": float(np.mean(ranks)),
        "coverage68": float(np.mean(cover68)),
        "coverage90": float(np.mean(cover90)),
        "ranks": np.asarray(ranks),
    }


def test_prior_predictive_calibration_with_h0_dependent_noise_and_selection():
    result = _prior_predictive_rank_calibration()
    assert 0.43 < result["mean_rank"] < 0.57
    assert 0.54 < result["coverage68"] < 0.82
    assert 0.80 < result["coverage90"] < 0.99


def test_h0_dependent_selection_normalization_is_a_detectable_negative_control():
    correct = _prior_predictive_rank_calibration(include_selection=True)
    wrong = _prior_predictive_rank_calibration(include_selection=False)
    assert abs(wrong["mean_rank"] - 0.5) > 0.035 or abs(wrong["coverage68"] - 0.68) > 0.12
    assert np.mean(np.abs(wrong["ranks"] - 0.5)) > np.mean(np.abs(correct["ranks"] - 0.5))


def test_simulator_returns_only_detected_distance_events_and_is_reproducible():
    from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
        simulate_detected_distance_events,
    )
    z = np.array([0.03, 0.06, 0.10, 0.14])
    weights = np.array([1.0, 2.0, 3.0, 4.0])
    kwargs = dict(
        n_events=5, population_redshifts=z, population_weights=weights,
        H0_true=70.0, detection_distance_min=300.0,
        fractional_sigma=0.12, sigma_floor_mpc=8.0,
    )
    first = simulate_detected_distance_events(rng=np.random.default_rng(912), **kwargs)
    second = simulate_detected_distance_events(rng=np.random.default_rng(912), **kwargs)
    assert [e["gw_distance"] for e in first] == [e["gw_distance"] for e in second]
    assert len(first) == 5
    assert all(e["gw_distance"] >= 300.0 for e in first)
    np.testing.assert_allclose(first[0]["weights"].sum(), 1.0)
    np.testing.assert_array_equal(first[0]["redshifts"], z)


def test_simulator_rejects_bad_count_and_rng():
    from dark_siren_robustness.selection_conditioned_h0_dependent_noise import (
        simulate_detected_distance_events,
    )
    kwargs = dict(population_redshifts=[0.03, 0.08], H0_true=70.0,
                  detection_distance_min=300.0, fractional_sigma=0.12)
    with pytest.raises(ValueError, match="rng"):
        simulate_detected_distance_events(rng=123, n_events=3, **kwargs)
    with pytest.raises(ValueError, match="n_events"):
        simulate_detected_distance_events(rng=np.random.default_rng(1), n_events=0, **kwargs)

