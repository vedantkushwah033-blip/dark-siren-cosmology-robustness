"""Downstream H0 convergence across validation-posterior quadrature orders.

Coefficient-posterior moment convergence is not sufficient: this runner
re-evaluates the complete shared-node selected likelihood at each quadrature
order and compares both its relative log-likelihood curve and normalized H0
posterior. It is a numerical diagnostic, not a calibration experiment.
"""
from __future__ import annotations

from dataclasses import dataclass
from collections.abc import Mapping, Sequence

import numpy as np
from numpy.typing import ArrayLike

from .h0_posterior import H0GridPosterior, compare_h0_grid_posteriors, normalize_h0_posterior
from .inclusion_posterior import (
    InclusionPosteriorQuadrature,
    fit_inclusion_posterior_quadrature,
)
from .shared_inclusion_likelihood import (
    SharedInclusionLikelihoodResult,
    shared_inclusion_selected_log_likelihood,
)


@dataclass(frozen=True)
class InclusionH0QuadratureConvergence:
    """Full per-order outputs and downstream convergence diagnostics."""

    orders: tuple[int, ...]
    fits: tuple[InclusionPosteriorQuadrature, ...]
    likelihoods: tuple[SharedInclusionLikelihoodResult, ...]
    posteriors: tuple[H0GridPosterior, ...]
    centered_log_likelihood_max_abs_change: tuple[float, ...]
    posterior_density_max_abs_change: tuple[float, ...]
    posterior_quantile_max_abs_change: tuple[float, ...]
    status: str


def run_inclusion_h0_quadrature_convergence(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    H0_grid: ArrayLike,
    catalogue_redshift: ArrayLike,
    events: Sequence[Mapping[str, float]],
    *,
    detection_threshold_mpc: float,
    population_base_weights: ArrayLike | None = None,
    degree: int = 1,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    coordinate_center: float | None = None,
    coordinate_scale: float | None = None,
    orders: tuple[int, ...] = (7, 11, 15, 19, 23),
    centered_log_likelihood_tolerance: float = 1e-3,
    posterior_density_tolerance: float = 1e-3,
    posterior_quantile_tolerance: float = 0.05,
    quantile_probabilities: tuple[float, ...] = (0.05, 0.16, 0.50, 0.84, 0.95),
    consecutive_passes: int = 2,
) -> InclusionH0QuadratureConvergence:
    """Recompute the selected H0 posterior at each coefficient-GH order.

    The comparison uses the maximum absolute difference between log-likelihood
    curves after subtracting each curve's own maximum. This removes an
    irrelevant H0-independent numerical offset while retaining shape changes.
    It also compares normalized posterior density on the common H0 grid and
    selected posterior quantiles. A convergence pass requires all three
    criteria for the requested number of consecutive order transitions.
    """
    order_tuple = tuple(orders)
    if len(order_tuple) < 2 or any(
        isinstance(order, bool) or not isinstance(order, int) or order < 3
        for order in order_tuple
    ) or any(right <= left for left, right in zip(order_tuple, order_tuple[1:])):
        raise ValueError("orders must be a strictly increasing tuple of integers >= 3")
    for value, label in (
        (centered_log_likelihood_tolerance, "centered_log_likelihood_tolerance"),
        (posterior_density_tolerance, "posterior_density_tolerance"),
        (posterior_quantile_tolerance, "posterior_quantile_tolerance"),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be positive and finite")
    if isinstance(consecutive_passes, bool) or not isinstance(consecutive_passes, int) or consecutive_passes < 1:
        raise ValueError("consecutive_passes must be a positive integer")
    probabilities = tuple(float(q) for q in quantile_probabilities)
    if not probabilities or any(not np.isfinite(q) or not 0.0 <= q <= 1.0 for q in probabilities):
        raise ValueError("quantile_probabilities must be a non-empty sequence in [0, 1]")

    h0 = np.asarray(H0_grid, dtype=float)
    if h0.ndim != 1 or h0.size < 2 or not np.all(np.isfinite(h0)) or np.any(np.diff(h0) <= 0):
        raise ValueError("H0_grid must be finite, one-dimensional, and strictly increasing")

    fits = []
    likelihoods = []
    posteriors = []
    log_changes = []
    density_changes = []
    quantile_changes = []
    passing_streak = 0
    converged = False

    for order in order_tuple:
        fit = fit_inclusion_posterior_quadrature(
            validation_redshift, retained, degree=degree,
            prior_mean=prior_mean, prior_sd=prior_sd,
            coordinate_center=coordinate_center,
            coordinate_scale=coordinate_scale, order=order,
        )
        likelihood = shared_inclusion_selected_log_likelihood(
            h0, catalogue_redshift, events, fit,
            detection_threshold_mpc=detection_threshold_mpc,
            population_base_weights=population_base_weights,
        )
        posterior = normalize_h0_posterior(
            h0, likelihood.marginal_log_likelihood,
            prior_min=float(h0[0]), prior_max=float(h0[-1]),
        )
        fits.append(fit)
        likelihoods.append(likelihood)
        posteriors.append(posterior)

        if len(fits) == 1:
            continue
        previous_log = likelihoods[-2].marginal_log_likelihood
        current_log = likelihood.marginal_log_likelihood
        previous_centered = previous_log - np.max(previous_log)
        current_centered = current_log - np.max(current_log)
        log_change = float(np.max(np.abs(current_centered - previous_centered)))
        posterior_comparison = compare_h0_grid_posteriors(
            posteriors[-2], posterior,
            density_tolerance=posterior_density_tolerance,
            quantile_tolerance=posterior_quantile_tolerance,
            quantile_probabilities=probabilities,
        )
        log_changes.append(log_change)
        density_changes.append(posterior_comparison.max_abs_density_difference)
        quantile_changes.append(posterior_comparison.max_abs_quantile_difference)
        passed = (
            log_change <= centered_log_likelihood_tolerance
            and posterior_comparison.passed
        )
        passing_streak = passing_streak + 1 if passed else 0
        if passing_streak >= consecutive_passes:
            converged = True
            break

    return InclusionH0QuadratureConvergence(
        orders=tuple(fit.order for fit in fits),
        fits=tuple(fits),
        likelihoods=tuple(likelihoods),
        posteriors=tuple(posteriors),
        centered_log_likelihood_max_abs_change=tuple(log_changes),
        posterior_density_max_abs_change=tuple(density_changes),
        posterior_quantile_max_abs_change=tuple(quantile_changes),
        status="CONVERGED" if converged else "NON_CONVERGED",
    )
