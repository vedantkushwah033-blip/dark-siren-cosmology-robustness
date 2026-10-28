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
