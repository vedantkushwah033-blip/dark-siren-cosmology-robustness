import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.selection_toy import (
    C_KM_S,
    normalized_host_density,
    selection_aware_toy_likelihood,
)


def test_uniform_host_density_normalizes_on_nonuniform_grid():
    z = np.array([0.01, 0.02, 0.05, 0.08, 0.12])
    density = normalized_host_density(z, np.ones_like(z))
    assert trapezoid(density, z) == pytest.approx(1.0)
    assert density == pytest.approx(np.full(z.size, 1.0 / (z[-1] - z[0])))


def test_host_density_includes_declared_merger_weight():
    z = np.linspace(0.01, 0.12, 101)
    intensity = 1.0 + 4.0 * np.exp(-0.5 * ((z - 0.08) / 0.01) ** 2)
    rate = 1.0 + z
    density = normalized_host_density(z, intensity, rate)
    assert np.all(density > 0)
    assert trapezoid(density, z) == pytest.approx(1.0, abs=1e-12)
    assert density[np.argmax(z)] > density[0]


def test_detection_efficiency_matches_independent_quadrature():
    z = np.linspace(0.01, 0.15, 2001)
    density = normalized_host_density(z, 1.0 + 2.0 * z)
    h0 = np.array([60.0, 70.0, 80.0])
    threshold, sigma = 500.0, 100.0
    result = selection_aware_toy_likelihood(
        z, density, h0, [700.0, 850.0], sigma=sigma,
        detection_threshold=threshold,
    )
    expected = np.array([
        trapezoid(density * norm.cdf((C_KM_S * z / value - threshold) / sigma), z)
        for value in h0
    ])
    assert result["detection_efficiency"] == pytest.approx(expected, rel=1e-12, abs=1e-12)
    assert np.all((result["detection_efficiency"] > 0) & (result["detection_efficiency"] <= 1))


def test_selected_likelihood_uses_minus_n_log_efficiency():
    z = np.linspace(0.01, 0.15, 1001)
    density = normalized_host_density(z, np.ones_like(z))
    h0 = np.array([60.0, 70.0, 80.0])
    distances = np.array([700.0, 850.0, 950.0])
    result = selection_aware_toy_likelihood(
        z, density, h0, distances, sigma=100.0, detection_threshold=500.0
    )
    expected = (
        result["unconditional_log_likelihood"]
        - len(distances) * np.log(result["detection_efficiency"])
    )
    assert result["selected_log_likelihood"] == pytest.approx(expected)
    assert result["n_detected"] == len(distances)


def test_selected_events_below_threshold_are_rejected():
    z = np.linspace(0.01, 0.15, 101)
    density = normalized_host_density(z, np.ones_like(z))
    with pytest.raises(ValueError, match="detection_threshold"):
        selection_aware_toy_likelihood(
            z, density, [70.0, 75.0], [400.0], sigma=100.0,
            detection_threshold=500.0,
        )


def test_invalid_density_and_noise_scale_are_rejected():
    z = np.linspace(0.01, 0.15, 101)
    density = normalized_host_density(z, np.ones_like(z))
    with pytest.raises(ValueError, match="integrate to one"):
        selection_aware_toy_likelihood(
            z, density * 2.0, [70.0, 75.0], [700.0],
            sigma=100.0, detection_threshold=500.0,
        )
    with pytest.raises(ValueError, match="sigma"):
        selection_aware_toy_likelihood(
            z, density, [70.0, 75.0], [700.0],
            sigma=0.0, detection_threshold=500.0,
        )

def test_full_likelihood_matches_direct_quadrature_reference():
    """Cross-check event mixtures and selection correction without implementation helpers."""
    z = np.linspace(0.015, 0.14, 1703)
    raw = 1.0 + 3.0 * np.exp(-0.5 * ((z - 0.09) / 0.018) ** 2)
    density = normalized_host_density(z, raw, 1.0 + 0.5 * z)
    h0 = np.array([58.0, 71.0, 83.0])
    distances = np.array([610.0, 780.0, 910.0])
    sigma, threshold = 115.0, 480.0

    result = selection_aware_toy_likelihood(
        z, density, h0, distances, sigma=sigma, detection_threshold=threshold
    )

    # Deliberately use scipy's direct quadrature and probability density, not
    # the implementation's precomputed trapezoid weights or log-sum-exp path.
    reference_joint = []
    reference_beta = []
    for candidate_h0 in h0:
        event_integrals = [
            trapezoid(
                density * norm.pdf(
                    observed_distance, loc=C_KM_S * z / candidate_h0, scale=sigma
                ),
                z,
            )
            for observed_distance in distances
        ]
        reference_joint.append(np.sum(np.log(event_integrals)))
        detection_curve = norm.cdf(
            (C_KM_S * z / candidate_h0 - threshold) / sigma
        )
        reference_beta.append(trapezoid(density * detection_curve, z))

    reference_joint = np.asarray(reference_joint)
    reference_beta = np.asarray(reference_beta)
    assert result["unconditional_log_likelihood"] == pytest.approx(
        reference_joint, rel=1e-11, abs=1e-11
    )
    assert result["detection_efficiency"] == pytest.approx(
        reference_beta, rel=1e-11, abs=1e-11
    )
    assert result["selected_log_likelihood"] == pytest.approx(
        reference_joint - distances.size * np.log(reference_beta),
        rel=1e-11,
        abs=1e-11,
    )


def test_detection_efficiency_predicts_empirical_threshold_rate():
    """Monte Carlo calibration of beta against the declared generative model."""
    z = np.linspace(0.01, 0.15, 1201)
    density = normalized_host_density(z, 1.0 + 2.0 * z)
    h0, sigma, threshold = 70.0, 100.0, 500.0
    result = selection_aware_toy_likelihood(
        z, density, [h0, 75.0], [700.0], sigma=sigma,
        detection_threshold=threshold,
    )
    # Draw from the trapezoidal discrete approximation to the host density;
    # its node masses are density(z) times the quadrature weights. The
    # independently simulated measurement noise is Gaussian as specified.
    delta = np.diff(z)
    weights = np.empty_like(z)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    node_probability = density * weights
    node_probability /= node_probability.sum()

    rng = np.random.default_rng(20261028)
    n_draws = 80_000
    host_z = rng.choice(z, size=n_draws, p=node_probability)
    measured_distance = C_KM_S * host_z / h0 + rng.normal(0.0, sigma, n_draws)
    empirical_rate = np.mean(measured_distance >= threshold)
    predicted_rate = result["detection_efficiency"][0]
    monte_carlo_se = np.sqrt(predicted_rate * (1.0 - predicted_rate) / n_draws)
    assert abs(empirical_rate - predicted_rate) < 5.0 * monte_carlo_se + 0.001

