"""Separately versioned damped-Newton logistic inclusion estimators.

V2 is intentionally opt-in and does not replace the v1 L-BFGS-B estimator.
It solves the box-constrained convex logistic negative log-likelihood with
analytic gradient/Hessian, projected-KKT stopping, and monotone backtracking.
A finite coefficient box constrains separation numerically; boundary-limited
fits must be reported and are not evidence of positivity or model adequacy.
"""
from __future__ import annotations

from collections.abc import Mapping

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import expit, log_expit


def _validate_design(
    redshift: ArrayLike, included: ArrayLike, *, z_min: float, z_max: float,
    degree: int,
) -> tuple[NDArray[np.float64], NDArray[np.float64], NDArray[np.float64]]:
    z = np.asarray(redshift, dtype=float)
    raw_y = np.asarray(included)
    if z.ndim != 1 or z.size < 10 or not np.all(np.isfinite(z)):
        raise ValueError("validation_redshift must contain at least 10 finite values")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or not z_min >= 0 or not z_min < z_max:
        raise ValueError("z_min and z_max must be finite with 0 <= z_min < z_max")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("validation_redshift must lie inside the declared domain")
    if raw_y.shape != z.shape or raw_y.ndim != 1:
        raise ValueError("included must be a one-dimensional array matching validation_redshift")
    if not np.all(np.isin(raw_y, [0, 1, False, True])):
        raise ValueError("included must contain only binary 0/1 values")
    y = raw_y.astype(float)
    if np.unique(y).size != 2:
        raise ValueError("validation sample must contain both included and not-included outcomes")
    midpoint = (float(z_min) + float(z_max)) / 2.0
    half_range = (float(z_max) - float(z_min)) / 2.0
    x = (z - midpoint) / half_range
    columns = [np.ones(z.size), x]
    if degree == 2:
        columns.append(x * x)
    return z, y, np.column_stack(columns)


def _logistic_nll_gradient_hessian(
    design: NDArray[np.float64], y: NDArray[np.float64], theta: NDArray[np.float64],
) -> tuple[float, NDArray[np.float64], NDArray[np.float64]]:
    """Return stable logistic NLL, analytic gradient, and observed Hessian."""
    eta = design @ theta
    probability = expit(eta)
    residual = probability - y
    variance = probability * (1.0 - probability)
    nll = float(-np.sum(y * log_expit(eta) + (1.0 - y) * log_expit(-eta)))
    gradient = design.T @ residual
    hessian = design.T @ (design * variance[:, None])
    return nll, gradient, hessian


def _fit_box_logistic(
    design: NDArray[np.float64], y: NDArray[np.float64], bounds: tuple[tuple[float, float], ...],
    *, max_iter: int = 500, kkt_tolerance: float = 1e-7,
) -> dict[str, object]:
    lower = np.asarray([b[0] for b in bounds], dtype=float)
    upper = np.asarray([b[1] for b in bounds], dtype=float)
    if design.shape[1] != lower.size or np.any(lower >= upper):
        raise ValueError("coefficient bounds do not match the design")
    theta = np.zeros(design.shape[1], dtype=float)
    def objective(value: NDArray[np.float64]) -> float:
        return _logistic_nll_gradient_hessian(design, y, value)[0]

    def derivatives(value: NDArray[np.float64]) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
        _, gradient, hessian = _logistic_nll_gradient_hessian(design, y, value)
        return gradient, hessian

    value = objective(theta)
    accepted_steps = 0
    rejected_steps = 0
    converged = False
    for iteration in range(1, max_iter + 1):
        gradient, hessian = derivatives(theta)
        projected = theta - np.clip(theta - gradient, lower, upper)
        kkt = float(np.max(np.abs(projected)))
        if not np.isfinite(kkt):
            raise RuntimeError("v2 logistic fit produced a non-finite projected KKT residual")
        if kkt <= kkt_tolerance:
            converged = True
            break

        # Identify coordinates whose active box constraints already satisfy
        # the KKT sign condition. A full-space Newton step can point outside
        # such a face; clipping that step may destroy descent in the remaining
        # coordinates. Solve the Newton system on the free subspace instead.
        active_lower = (theta <= lower + 1e-10) & (gradient > 0.0)
        active_upper = (theta >= upper - 1e-10) & (gradient < 0.0)
        free = ~(active_lower | active_upper)

        direction = None
        if np.any(free):
            free_hessian = hessian[np.ix_(free, free)]
            free_gradient = gradient[free]
            ridge = 0.0
            for _ in range(12):
                try:
                    free_direction = np.linalg.solve(
                        free_hessian + ridge * np.eye(int(np.sum(free))),
                        -free_gradient,
                    )
                except np.linalg.LinAlgError:
                    free_direction = np.full_like(free_gradient, np.nan)
                if (
                    np.all(np.isfinite(free_direction))
                    and float(free_gradient @ free_direction) < 0.0
                ):
                    direction = np.zeros_like(theta)
                    direction[free] = free_direction
                    break
                ridge = 1e-10 if ridge == 0 else ridge * 10.0

        # If the reduced Hessian is singular/ill-conditioned, use a
        # projected-gradient direction. The line search still enforces descent.
        if direction is None:
            direction = -gradient / max(1.0, float(np.linalg.norm(gradient, ord=np.inf)))

        step = 1.0
        accepted = False
        for _ in range(60):
            candidate = np.clip(theta + step * direction, lower, upper)
            delta = candidate - theta
            if not np.any(delta):
                step *= 0.5
                rejected_steps += 1
                continue
            candidate_value = objective(candidate)
            directional_change = float(gradient @ delta)
            if (directional_change < 0.0 and np.isfinite(candidate_value)
                    and candidate_value <= value + 1e-4 * directional_change):
                theta = candidate
                value = candidate_value
                accepted = True
                accepted_steps += 1
                break
            step *= 0.5
            rejected_steps += 1
        if not accepted:
            # A projected-gradient fallback is attempted once, without silently
            # declaring convergence if it cannot make a valid step.
            direction = -gradient / max(1.0, float(np.linalg.norm(gradient, ord=np.inf)))
            step = 1.0
            for _ in range(60):
                candidate = np.clip(theta + step * direction, lower, upper)
                delta = candidate - theta
                candidate_value = objective(candidate)
                if np.any(delta) and np.isfinite(candidate_value) and candidate_value < value:
                    theta, value = candidate, candidate_value
                    accepted = True
                    accepted_steps += 1
                    break
                step *= 0.5
                rejected_steps += 1
        if not accepted:
            break

    _, gradient, final_hessian = _logistic_nll_gradient_hessian(design, y, theta)
    projected = theta - np.clip(theta - gradient, lower, upper)
    kkt = float(np.max(np.abs(projected)))
    eigenvalues = np.linalg.eigvalsh(final_hessian)
    largest_eigenvalue = float(max(0.0, eigenvalues[-1]))
    rank_tolerance = np.finfo(float).eps * max(final_hessian.shape) * largest_eigenvalue
    hessian_rank = int(np.count_nonzero(eigenvalues > rank_tolerance))
    positive_eigenvalues = eigenvalues[eigenvalues > rank_tolerance]
    hessian_condition = (
        float(largest_eigenvalue / positive_eigenvalues[0])
        if hessian_rank == final_hessian.shape[0] else float("inf")
    )
    if not converged and kkt <= kkt_tolerance:
        converged = True
    if not converged:
        raise RuntimeError(
            f"v2 logistic fit did not converge: iterations={iteration}, "
            f"projected_gradient_inf_norm={kkt:.6g}, objective={value:.12g}"
        )
    bound_hit = bool(np.any(np.isclose(theta, lower, atol=1e-6, rtol=0.0))
                     or np.any(np.isclose(theta, upper, atol=1e-6, rtol=0.0)))
    return {
        "coefficients": theta,
        "negative_log_likelihood": value,
        "iterations": iteration,
        "accepted_steps": accepted_steps,
        "rejected_steps": rejected_steps,
        "projected_gradient_inf_norm": kkt,
        "hessian_rank": hessian_rank,
        "parameter_dimension": int(final_hessian.shape[0]),
        "hessian_condition_number": hessian_condition,
        "hessian_min_eigenvalue": float(eigenvalues[0]),
        "hessian_max_eigenvalue": largest_eigenvalue,
        "coefficient_bound_hit": bound_hit,
        "fit_status": "boundary_limited" if bound_hit else "converged",
    }


def fit_logistic_inclusion_model_v2(
    validation_redshift: ArrayLike, included: ArrayLike, *, z_min: float, z_max: float,
) -> dict[str, object]:
    """Fit the linear inclusion model with the V2 damped-Newton solver."""
    z, y, design = _validate_design(
        validation_redshift, included, z_min=z_min, z_max=z_max, degree=1
    )
    fit = _fit_box_logistic(design, y, ((-25.0, 25.0), (-50.0, 50.0)))
    intercept, slope = fit.pop("coefficients")
    return {
        "model": "logistic_inclusion_v2", "solver_version": "v2",
        "intercept": float(intercept), "slope": float(slope),
        "z_min": float(z_min), "z_max": float(z_max),
        "n_validation": int(z.size), "n_included": int(y.sum()),
        "n_not_included": int(y.size - y.sum()), **fit,
    }


def fit_quadratic_logistic_inclusion_model_v2(
    validation_redshift: ArrayLike, included: ArrayLike, *, z_min: float, z_max: float,
) -> dict[str, object]:
    """Fit the quadratic inclusion model with the V2 damped-Newton solver."""
    z, y, design = _validate_design(
        validation_redshift, included, z_min=z_min, z_max=z_max, degree=2
    )
    fit = _fit_box_logistic(
        design, y, ((-25.0, 25.0), (-50.0, 50.0), (-50.0, 50.0))
    )
    intercept, slope, quadratic = fit.pop("coefficients")
    return {
        "model": "quadratic_logistic_inclusion_v2", "solver_version": "v2",
        "intercept": float(intercept), "slope": float(slope),
        "quadratic": float(quadratic), "z_min": float(z_min), "z_max": float(z_max),
        "n_validation": int(z.size), "n_included": int(y.sum()),
        "n_not_included": int(y.size - y.sum()), **fit,
    }


def predict_inclusion_probability_v2(
    redshift: ArrayLike, model: Mapping[str, object],
) -> NDArray[np.float64]:
    """Predict probabilities from a V2 linear or quadratic fit result."""
    kind = model.get("model")
    if kind not in {"logistic_inclusion_v2", "quadratic_logistic_inclusion_v2"}:
        raise ValueError("model must be a supported V2 logistic inclusion result")
    keys = ("intercept", "slope", "z_min", "z_max")
    try:
        intercept, slope, z_min, z_max = (float(model[key]) for key in keys)
        quadratic = float(model.get("quadratic", 0.0))
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError("model is missing valid coefficients or domain bounds") from exc
    if not np.all(np.isfinite([intercept, slope, quadratic, z_min, z_max])) or z_min >= z_max:
        raise ValueError("model coefficients and domain bounds must be finite and valid")
    if abs(intercept) > 25 or abs(slope) > 50 or abs(quadratic) > 50:
        raise ValueError("model coefficients exceed the V2 fit bounds")
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)) or np.any((z < z_min) | (z > z_max)):
        raise ValueError("redshift must be finite, one-dimensional, and inside the fitted domain")
    x = (z - (z_min + z_max) / 2) / ((z_max - z_min) / 2)
    return np.asarray(expit(intercept + slope * x + quadratic * x**2), dtype=float)
