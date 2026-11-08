"""Order-ladder convergence checks for inclusion-posterior quadrature.

Convergence is assessed on posterior moments, predictive inclusion probabilities,
and the quadrature log normalizer. ESS and maximum node mass are reported as
diagnostics, not used as substitutes for convergence.
"""
from __future__ import annotations

from dataclasses import dataclass
from typing import Literal

import numpy as np
from numpy.typing import ArrayLike

from .inclusion_posterior import (
    InclusionPosteriorQuadrature,
    fit_inclusion_posterior_quadrature,
    posterior_predictive_inclusion,
)


@dataclass(frozen=True)
class QuadratureOrderCheck:
    """Comparison between two consecutive quadrature orders."""

    previous_order: int
    current_order: int
    max_abs_mean_difference: float
    max_abs_sd_difference: float
    max_abs_predictive_difference: float
    abs_log_normalizer_difference: float
    quadrature_ess: float
    maximum_normalized_weight: float
    passed: bool


@dataclass(frozen=True)
class QuadratureConvergenceResult:
    """Final order-ladder result; non-convergence is an explicit outcome."""

    status: Literal["CONVERGED", "NON_CONVERGED"]
    posterior: InclusionPosteriorQuadrature
    checks: tuple[QuadratureOrderCheck, ...]
    attempted_orders: tuple[int, ...]


def fit_inclusion_posterior_with_order_ladder(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    *,
    degree: int = 1,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    orders: tuple[int, ...] = (7, 11, 15, 19, 23),
    predictive_redshift: ArrayLike | None = None,
    mean_tolerance: float = 1e-3,
    sd_tolerance: float = 1e-3,
    predictive_tolerance: float = 1e-3,
    log_normalizer_tolerance: float = 1e-3,
) -> QuadratureConvergenceResult:
    """Fit increasing GH orders until all consecutive-order checks pass.

    Predictive comparisons use predictive_redshift, or validation redshifts
    when omitted. Tolerances are absolute. The full check history is retained.
    """
    order_tuple = tuple(orders)
    if len(order_tuple) < 2:
        raise ValueError("orders must contain at least two increasing values")
    if any(isinstance(o, bool) or not isinstance(o, int) or o < 3 for o in order_tuple):
        raise ValueError("every quadrature order must be an integer >= 3")
    if any(b <= a for a, b in zip(order_tuple, order_tuple[1:])):
        raise ValueError("orders must be strictly increasing")

    tolerances = np.asarray(
        [mean_tolerance, sd_tolerance, predictive_tolerance, log_normalizer_tolerance],
        dtype=float,
    )
    if not np.all(np.isfinite(tolerances)) or np.any(tolerances <= 0):
        raise ValueError("all convergence tolerances must be finite and positive")

    z = np.asarray(validation_redshift, dtype=float)
    probe = z if predictive_redshift is None else np.asarray(predictive_redshift, dtype=float)
    if probe.ndim != 1 or probe.size == 0 or not np.all(np.isfinite(probe)):
        raise ValueError("predictive_redshift must be a non-empty finite 1D array")

    previous = None
    previous_predictive = None
    checks: list[QuadratureOrderCheck] = []
    attempted: list[int] = []
    latest = None

    for order in order_tuple:
        latest = fit_inclusion_posterior_quadrature(
            z, retained, degree=degree, prior_mean=prior_mean,
            prior_sd=prior_sd, order=order,
        )
        attempted.append(order)
        current_predictive = posterior_predictive_inclusion(latest, probe)
        if previous is not None:
            previous_sd = np.sqrt(np.diag(previous.posterior_covariance))
            current_sd = np.sqrt(np.diag(latest.posterior_covariance))
            mean_diff = float(np.max(np.abs(latest.posterior_mean - previous.posterior_mean)))
            sd_diff = float(np.max(np.abs(current_sd - previous_sd)))
            predictive_diff = float(np.max(np.abs(current_predictive - previous_predictive)))
            normalizer_diff = float(abs(
                latest.log_evidence_unnormalized - previous.log_evidence_unnormalized
            ))
            passed = (
                mean_diff <= mean_tolerance
                and sd_diff <= sd_tolerance
                and predictive_diff <= predictive_tolerance
                and normalizer_diff <= log_normalizer_tolerance
            )
            checks.append(QuadratureOrderCheck(
                previous_order=previous.order, current_order=latest.order,
                max_abs_mean_difference=mean_diff,
                max_abs_sd_difference=sd_diff,
                max_abs_predictive_difference=predictive_diff,
                abs_log_normalizer_difference=normalizer_diff,
                quadrature_ess=latest.quadrature_ess,
                maximum_normalized_weight=latest.maximum_normalized_weight,
                passed=passed,
            ))
            if passed:
                return QuadratureConvergenceResult(
                    status="CONVERGED", posterior=latest,
                    checks=tuple(checks), attempted_orders=tuple(attempted),
                )
        previous = latest
        previous_predictive = current_predictive

    assert latest is not None
    return QuadratureConvergenceResult(
        status="NON_CONVERGED", posterior=latest,
        checks=tuple(checks), attempted_orders=tuple(attempted),
    )
