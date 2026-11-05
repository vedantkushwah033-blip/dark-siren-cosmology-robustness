import numpy as np
from scipy.integrate import trapezoid

from dark_siren_robustness.catalogue_density import estimate_host_density


def _adaptive_kde_directional_derivative(grid, sample, base_weights, delta, rate):
    """Independent derivative of the reflected, adaptive production KDE."""
    a = base_weights / np.sum(base_weights)
    delta_bar = float(np.sum(a * delta))
    a_prime = -a * (delta - delta_bar)

    mean = float(np.sum(a * sample))
    centered2 = (sample - mean) ** 2
    variance = float(np.sum(a * centered2))
    ess = float(1.0 / np.sum(a**2))
    variance_prime = -float(np.sum(a * (centered2 - variance) * delta))
    ess_log_prime = 2.0 * (
        float(np.sum(a**2 * delta) / np.sum(a**2)) - delta_bar
    )
    bandwidth = np.sqrt(variance) * ess ** (-1.0 / 5.0)
    bandwidth_log_prime = 0.5 * variance_prime / variance - 0.2 * ess_log_prime
    bandwidth_prime = bandwidth * bandwidth_log_prime

    lower, upper = float(grid[0]), float(grid[-1])
    reflected = np.stack([sample, 2.0 * lower - sample, 2.0 * upper - sample])
    offsets = grid[:, None, None] - reflected[None, :, :]
    kernels = np.exp(-0.5 * (offsets / bandwidth) ** 2) / (
        np.sqrt(2.0 * np.pi) * bandwidth
    )
    # Axis 1 indexes the three reflected copies; axis 2 indexes observations.
    weighted_kernel = kernels * (a[None, None, :] / 3.0)
    raw = np.sum(weighted_kernel)
    del raw  # The grid-dependent raw density is calculated below.

    raw_density = np.sum(weighted_kernel, axis=(1, 2))
    log_kernel_bandwidth_derivative = -1.0 / bandwidth + offsets**2 / bandwidth**3
    raw_derivative = np.sum(
        kernels * (a_prime[None, None, :] / 3.0)
        + weighted_kernel * log_kernel_bandwidth_derivative * bandwidth_prime,
        axis=(1, 2),
    )

    unnormalized = raw_density * rate
    area = float(trapezoid(unnormalized, grid))
    density = unnormalized / area
    relative_raw_derivative = raw_derivative / raw_density
    centering = float(trapezoid(density * relative_raw_derivative, grid))
    return density * (relative_raw_derivative - centering), bandwidth, bandwidth_prime


def test_adaptive_reflected_ipw_kde_derivative_matches_production_finite_difference():
    rng = np.random.default_rng(202611056)
    grid = np.linspace(0.01, 0.40, 401)
    sample = np.sort(rng.uniform(0.025, 0.385, 90))
    inclusion = 0.18 + 0.72 * (1.0 / (1.0 + np.exp(-2.0 * (sample - 0.20) / 0.195)))
    base_weights = 1.0 / inclusion
    delta = 0.55 * np.sin(2.7 * (sample - 0.01) / 0.39) + 0.2 * sample / 0.40
    rate = 1.0 + 0.8 * (grid - grid[0]) / (grid[-1] - grid[0])

    derivative, bandwidth, bandwidth_prime = _adaptive_kde_directional_derivative(
        grid, sample, base_weights, delta, rate
    )
    assert bandwidth > 0.0
    assert abs(bandwidth_prime) > 1e-5

    epsilon = 2e-5
    plus = estimate_host_density(
        grid, sample, inclusion_probability=inclusion * (1.0 + epsilon * delta),
        merger_weight=rate,
    )
    minus = estimate_host_density(
        grid, sample, inclusion_probability=inclusion * (1.0 - epsilon * delta),
        merger_weight=rate,
    )
    finite_difference = (plus - minus) / (2.0 * epsilon)

    assert np.all(np.isfinite(derivative))
    assert np.all(np.isfinite(finite_difference))
    assert np.max(np.abs(derivative)) > 0.0
    assert np.allclose(finite_difference, derivative, rtol=2e-4, atol=2e-6)
    assert abs(float(trapezoid(derivative, grid))) < 1e-10
