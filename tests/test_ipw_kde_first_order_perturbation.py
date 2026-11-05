"""Deterministic checks for the fixed-bandwidth IPW-KDE perturbation.

This isolates the mathematical operator with a fixed bandwidth. It does not
claim to validate the production estimator's adaptive bandwidth derivative.
"""
import numpy as np
from scipy.integrate import trapezoid


def _reflected_kernel_operator(grid, sample, weights, bandwidth):
    left = 2.0 * grid[0] - sample
    right = 2.0 * grid[-1] - sample
    augmented = np.concatenate((sample, left, right))
    augmented_weights = np.tile(weights / 3.0, 3)
    standardized = (grid[:, None] - augmented[None, :]) / bandwidth
    kernel = np.exp(-0.5 * standardized**2) / (
        np.sqrt(2.0 * np.pi) * bandwidth
    )
    return kernel @ augmented_weights


def _normalize(values, grid):
    return values / trapezoid(values, grid)


def test_fixed_bandwidth_reflected_ipw_kde_has_kernel_averaged_first_order_error():
    # Fixed synthetic catalogue support, deliberately including boundary-near
    # observations so the reflected operator is exercised.
    grid = np.linspace(0.01, 0.40, 801)
    sample = np.array([0.012, 0.025, 0.055, 0.09, 0.13, 0.18,
                       0.23, 0.29, 0.34, 0.385, 0.398])
    base = np.array([0.7, 1.0, 1.2, 1.1, 0.9, 1.3,
                     1.0, 0.8, 1.1, 0.9, 0.6])
    base = base / base.sum()
    delta = np.array([-0.8, -0.5, -0.2, 0.1, 0.35, 0.7,
                      0.9, 0.55, 0.2, -0.25, -0.6])
    bandwidth = 0.035
    merger_weight = 1.0 + 0.8 * (grid - grid[0]) / (grid[-1] - grid[0])

    # Baseline smoothed density after merger weighting and finite-grid
    # normalization, using exactly the same fixed kernel operator below.
    baseline_galaxy = _reflected_kernel_operator(
        grid, sample, base, bandwidth
    )
    baseline = _normalize(baseline_galaxy * merger_weight, grid)

    # Kernel-local average of delta under the baseline kernel contributions.
    # The three reflected copies carry the same one-third mass as production.
    left = 2.0 * grid[0] - sample
    right = 2.0 * grid[-1] - sample
    augmented = np.concatenate((sample, left, right))
    augmented_delta = np.tile(delta, 3)
    standardized = (grid[:, None] - augmented[None, :]) / bandwidth
    kernel = np.exp(-0.5 * standardized**2) / (
        np.sqrt(2.0 * np.pi) * bandwidth
    )
    contribution = kernel * np.tile(base / 3.0, 3)[None, :]
    local_delta = (contribution @ augmented_delta) / contribution.sum(axis=1)

    # After merger weighting and final normalization, the local perturbation
    # must be recentered under the baseline merger-weighted KDE.
    delta_bar = trapezoid(baseline * local_delta, grid)
    predicted_relative = -(local_delta - delta_bar)

    errors = []
    for epsilon in (2.0e-2, 1.0e-2, 5.0e-3):
        perturbed_weights = base / (1.0 + epsilon * delta)
        perturbed_weights /= perturbed_weights.sum()
        perturbed_galaxy = _reflected_kernel_operator(
            grid, sample, perturbed_weights, bandwidth
        )
        perturbed = _normalize(perturbed_galaxy * merger_weight, grid)
        observed_relative = (perturbed - baseline) / baseline
        errors.append(
            np.sqrt(trapezoid(
                baseline * (observed_relative / epsilon - predicted_relative)**2,
                grid,
            ))
        )

    # First-order remainder should shrink approximately linearly in epsilon
    # after dividing the density difference by epsilon (O(epsilon) error).
    assert errors[2] < errors[1] < errors[0]
    assert errors[2] < 0.01
    assert np.isclose(trapezoid(baseline, grid), 1.0, atol=1e-12)
