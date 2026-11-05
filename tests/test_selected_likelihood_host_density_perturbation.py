"""First-order host-density perturbation checks for the selected toy likelihood.

The derivative is conditional on a fixed redshift grid, fixed quadrature rule,
fixed event distances, and a normalized host-density perturbation. It includes
both event-integral and detection-efficiency normalization terms.
"""
import numpy as np
import pytest
from scipy.integrate import trapezoid
from scipy.stats import norm

from dark_siren_robustness.selection_toy import (
    C_KM_S,
    selection_aware_toy_likelihood,
)


def _trap_weights(grid):
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _first_order_loglikelihood_derivative(
    z, density, h0, distances, sigma, threshold, centered_delta
):
    """Derivative d ell(H0; eps)/d eps at zero for p_eps=p*(1-eps*g)."""
    quad = _trap_weights(z)
    means = C_KM_S * z / h0
    event_terms = []
    for distance in distances:
        kernel = norm.pdf(distance, loc=means, scale=sigma)
        mass = density * kernel * quad
        event_terms.append(float(np.sum(mass * centered_delta) / np.sum(mass)))

    selection = norm.cdf((means - threshold) / sigma)
    selected_mass = density * selection * quad
    selection_expectation = float(
        np.sum(selected_mass * centered_delta) / np.sum(selected_mass)
    )
    # d log I_j / d eps = -E_event[g]; d(-n log beta)/d eps
    # = +n E_selected[g]. Both pieces are essential.
    return -sum(event_terms) + len(distances) * selection_expectation


def test_selected_loglikelihood_first_order_derivative_includes_beta_term():
    z = np.linspace(0.01, 0.40, 1601)
    raw = 1.0 + 2.1 * np.exp(-0.5 * ((z - 0.11) / 0.035) ** 2)
    raw += 1.3 * np.exp(-0.5 * ((z - 0.29) / 0.055) ** 2)
    density = raw / trapezoid(raw, z)
    h0_grid = np.array([62.0, 70.0, 79.0])
    distances = np.array([520.0, 680.0, 830.0, 1010.0])
    sigma, threshold = 135.0, 430.0

    # Smooth, nonconstant perturbation. Centering under p ensures the
    # first-order density perturbation integrates to zero.
    raw_delta = 0.65 * np.sin(2.0 * np.pi * (z - z[0]) / (z[-1] - z[0]))
    raw_delta += 0.25 * (z - trapezoid(z * density, z))
    delta_mean = trapezoid(density * raw_delta, z)
    centered_delta = raw_delta - delta_mean
    assert abs(trapezoid(density * centered_delta, z)) < 1e-12

    baseline = selection_aware_toy_likelihood(
        z, density, h0_grid, distances,
        sigma=sigma, detection_threshold=threshold,
    )
    derivative = np.array([
        _first_order_loglikelihood_derivative(
            z, density, h0, distances, sigma, threshold, centered_delta
        )
        for h0 in h0_grid
    ])

    errors = []
    for epsilon in (2e-3, 1e-3, 5e-4):
        # This is an exactly normalized, positive density family whose
        # derivative at epsilon=0 is -p(z)*centered_delta(z).
        perturbed_raw = density * (1.0 - epsilon * centered_delta)
        perturbed_density = perturbed_raw / trapezoid(perturbed_raw, z)
        assert np.all(perturbed_density > 0)
        perturbed = selection_aware_toy_likelihood(
            z, perturbed_density, h0_grid, distances,
            sigma=sigma, detection_threshold=threshold,
        )
        finite_difference = (
            perturbed["selected_log_likelihood"]
            - baseline["selected_log_likelihood"]
        ) / epsilon
        errors.append(float(np.max(np.abs(finite_difference - derivative))))

    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 2e-3


def test_selection_normalizer_term_is_not_optional_in_perturbation():
    z = np.linspace(0.01, 0.40, 1201)
    density = (1.0 + 1.8 * z)
    density /= trapezoid(density, z)
    h0 = 70.0
    distances = np.array([510.0, 640.0, 780.0])
    sigma, threshold = 120.0, 460.0
    centered_delta = (z - trapezoid(z * density, z)) / 0.39

    full = _first_order_loglikelihood_derivative(
        z, density, h0, distances, sigma, threshold, centered_delta
    )
    quad = _trap_weights(z)
    means = C_KM_S * z / h0
    event_expectations = []
    for distance in distances:
        mass = density * norm.pdf(distance, loc=means, scale=sigma) * quad
        event_expectations.append(np.sum(mass * centered_delta) / np.sum(mass))
    event_only = -float(np.sum(event_expectations))
    selection = norm.cdf((means - threshold) / sigma)
    selected_mass = density * selection * quad
    beta_term = len(distances) * float(
        np.sum(selected_mass * centered_delta) / np.sum(selected_mass)
    )
    assert full == pytest.approx(event_only + beta_term, abs=1e-12)
    assert abs(beta_term) > 1e-4
