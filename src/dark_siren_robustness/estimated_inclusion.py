"""Finite-validation logistic models for catalogue inclusion probabilities.

This module estimates P(retained | redshift) from validation labels only.
It is a synthetic-methodology component, not a real-survey completeness model.
The fitter reports separation/conditioning failures rather than silently
regularizing them into apparently successful probability estimates.
"""
from __future__ import annotations
from dataclasses import dataclass
import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.special import expit

@dataclass(frozen=True)
class InclusionFit:
    """Fit result plus diagnostics; coefficients use standardized redshift."""
    status: str
    message: str
    degree: int
    center: float
    scale: float
    coefficients: NDArray[np.float64] | None
    covariance: NDArray[np.float64] | None
    n_validation: int
    n_retained: int
    n_not_retained: int
    design_rank: int
    design_condition_number: float
    gradient_inf_norm: float
    hessian_condition_number: float
    max_abs_coefficient: float
    iterations: int

    def predict(self, redshift: ArrayLike) -> NDArray[np.float64]:
        """Predict probabilities; unavailable for failed fits."""
        if self.status != "SUCCESS" or self.coefficients is None:
            raise RuntimeError(f"Cannot predict from fit with status={self.status}")
        z = np.asarray(redshift, dtype=float)
        if z.ndim != 1 or not np.all(np.isfinite(z)):
            raise ValueError("redshift must be a finite one-dimensional array")
        x = (z - self.center) / self.scale
        design = np.column_stack([x**power for power in range(self.degree + 1)])
        return np.asarray(expit(design @ self.coefficients), dtype=float)


def _design(redshift: NDArray[np.float64], degree: int):
    center = float(np.mean(redshift))
    scale = float(np.std(redshift, ddof=0))
    if not np.isfinite(scale) or scale <= 0:
        scale = 1.0
    x = (redshift - center) / scale
    return np.column_stack([x**power for power in range(degree + 1)]), center, scale


def fit_inclusion_logistic(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    *,
    degree: int = 1,
    max_iterations: int = 2000,
    separation_coefficient_limit: float = 25.0,
) -> InclusionFit:
    """Fit an unpenalized logistic inclusion model to validation labels.

    Degree 1 is linear in standardized redshift; degree 2 adds curvature.
    Validation sampling must be independent of the retained label.
    Failed or non-evaluable fits retain explicit status and diagnostics.
    """
    z = np.asarray(validation_redshift, dtype=float)
    y_raw = np.asarray(retained)
    if z.ndim != 1 or y_raw.ndim != 1 or z.shape != y_raw.shape:
        raise ValueError("validation_redshift and retained must be matching 1D arrays")
    if z.size < 5 or not np.all(np.isfinite(z)):
        raise ValueError("at least five finite validation redshifts are required")
    if np.any((y_raw != 0) & (y_raw != 1)):
        raise ValueError("retained labels must be binary (0 or 1)")
    if degree not in (1, 2):
        raise ValueError("degree must be 1 or 2")
    if not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    if not np.isfinite(separation_coefficient_limit) or separation_coefficient_limit <= 0:
        raise ValueError("separation_coefficient_limit must be finite and positive")

    y = y_raw.astype(float)
    n1, n0 = int(np.sum(y)), int(y.size - np.sum(y))
    xmat, center, scale = _design(z, degree)
    rank = int(np.linalg.matrix_rank(xmat))
    design_condition = float(np.linalg.cond(xmat))
    empty = dict(
        degree=degree, center=center, scale=scale, coefficients=None,
        covariance=None, n_validation=int(y.size), n_retained=n1,
        n_not_retained=n0, design_rank=rank,
        design_condition_number=design_condition,
        gradient_inf_norm=float("inf"), hessian_condition_number=float("inf"),
        max_abs_coefficient=float("inf"), iterations=0,
    )

    def failed(message: str, status: str = "FAILED", **updates) -> InclusionFit:
        return InclusionFit(status=status, message=message, **(empty | updates))

    if n0 == 0 or n1 == 0:
        return failed("Validation sample contains only one outcome class.", "NON_EVALUABLE")
    if rank != degree + 1 or not np.isfinite(design_condition):
        return failed("Validation design is rank deficient or ill-conditioned.")
    if design_condition > 1e10:
        return failed("Validation design condition number exceeds 1e10.")

    def objective(beta):
        eta = xmat @ beta
        return float(np.sum(np.logaddexp(0.0, eta) - y * eta))

    def gradient(beta):
        return xmat.T @ (expit(xmat @ beta) - y)

    fit = minimize(
        objective, np.zeros(degree + 1), jac=gradient, method="BFGS",
        options={"maxiter": max_iterations, "gtol": 1e-8},
    )
    beta = np.asarray(fit.x, dtype=float)
    grad_norm = (
        float(np.linalg.norm(gradient(beta), ord=np.inf))
        if np.all(np.isfinite(beta)) else float("inf")
    )

    # BFGS can report precision loss when it has reached a nearly converged
    # point.  A short damped-Newton polish uses the analytic score/information
    # to distinguish that numerical status from a genuinely unresolved fit.
    polish_iterations = 0
    if np.all(np.isfinite(beta)) and (not fit.success or grad_norm > 1e-10):
        for _ in range(50):
            score = gradient(beta)
            if not np.all(np.isfinite(score)) or np.linalg.norm(score, ord=np.inf) <= 1e-10:
                break
            probability = expit(xmat @ beta)
            information = xmat.T @ (
                (probability * (1.0 - probability))[:, None] * xmat
            )
            eigenvalues = np.linalg.eigvalsh(information)
            if not np.all(np.isfinite(eigenvalues)) or np.min(eigenvalues) <= 0:
                break
            direction = np.linalg.solve(information, score)
            decrement = float(score @ direction)
            if not np.isfinite(decrement) or decrement <= 0:
                break
            current_objective = objective(beta)
            step_length = 1.0
            accepted = False
            while step_length >= 2.0 ** -30:
                candidate = beta - step_length * direction
                if objective(candidate) <= (
                    current_objective - 1e-4 * step_length * decrement
                ):
                    beta = candidate
                    accepted = True
                    break
                step_length *= 0.5
            if not accepted:
                break
            polish_iterations += 1
        grad_norm = float(np.linalg.norm(gradient(beta), ord=np.inf))
    probability = expit(xmat @ beta)
    hessian = xmat.T @ ((probability * (1.0 - probability))[:, None] * xmat)
    eigenvalues = np.linalg.eigvalsh(hessian)
    hessian_condition = (
        float(np.max(eigenvalues) / np.min(eigenvalues))
        if np.min(eigenvalues) > 0 else float("inf")
    )
    max_beta = float(np.max(np.abs(beta)))
    common = dict(
        gradient_inf_norm=grad_norm,
        hessian_condition_number=hessian_condition,
        max_abs_coefficient=max_beta,
        iterations=int(getattr(fit, "nit", 0)) + polish_iterations,
    )

    if max_beta >= separation_coefficient_limit or hessian_condition > 1e12:
        return failed(
            "Possible complete/quasi separation or weak identification; "
            "unpenalized covariance is not reliable.", **common,
        )
    if not np.all(np.isfinite(beta)) or grad_norm > 1e-6:
        return failed(f"Score convergence criterion not met: {fit.message}", **common)
    if np.min(eigenvalues) <= 0 or not np.all(np.isfinite(eigenvalues)):
        return failed("Observed information matrix is not positive definite.", **common)

    covariance = np.linalg.inv(hessian)
    convergence_message = (
        "Converged unpenalized logistic MLE."
        if fit.success
        else f"Accepted by score tolerance despite optimizer message: {fit.message}"
    )
    return InclusionFit(
        status="SUCCESS", message=convergence_message,
        **(empty | {"coefficients": beta.copy(), "covariance": covariance, **common}),
    )
