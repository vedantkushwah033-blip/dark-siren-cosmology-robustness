"""First-order SRSWOR sandwich covariance for finite-frame logistic roots.

This is a design-based linearization primitive. It is deliberately separate
from the inverse-Hessian/model-based covariance of a Bayesian or likelihood
fit. The caller must supply the exact fixed-coordinate design matrix used by
the fitted coefficient vector; sample-estimated centering/scaling is not
silently treated as fixed.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.special import expit


@dataclass(frozen=True)
class SRSWORLogisticVariance:
    """Plug-in linearization covariance conditional on the realized frame."""

    covariance: NDArray[np.float64]
    sensitivity: NDArray[np.float64]
    score_covariance: NDArray[np.float64]
    fpc: float
    n_validation: int
    population_size: int
    design_condition_number: float


def srswor_logistic_linearization_covariance(
    design_matrix: ArrayLike,
    retained: ArrayLike,
    coefficients: ArrayLike,
    *,
    population_size: int,
) -> SRSWORLogisticVariance:
    """Estimate coefficient design variance for an SRSWOR validation sample.

    Let u_i(beta)=x_i * (y_i-expit(x_i' beta)). For a regular finite-frame
    pseudo-true score root, the first-order variance is estimated by

        A_hat^-1 * ((1 - n/N) / n * S_hat_u) * A_hat^-T,

    where A_hat is the sample mean of x_i x_i' p_i(1-p_i), and S_hat_u
    is the usual sample covariance of the individual score contributions.

    The supplied design matrix must use fixed, predeclared coordinates (or
    coordinates whose estimation uncertainty is separately propagated).
    This function does not turn an arbitrary inverse-Hessian covariance into
    a design variance, and does not include parent-frame or retention-process
    variability. It is not valid for separation, active coefficient bounds,
    a singular sensitivity matrix, or a non-regular score root.
    """
    x = np.asarray(design_matrix, dtype=float)
    y_raw = np.asarray(retained)
    beta = np.asarray(coefficients, dtype=float)

    if x.ndim != 2 or x.shape[0] < 2 or x.shape[1] < 1:
        raise ValueError("design_matrix must be a 2D array with at least two rows and one column")
    n, p = x.shape
    if y_raw.shape != (n,) or np.any((y_raw != 0) & (y_raw != 1)):
        raise ValueError("retained must be a binary vector matching design_matrix rows")
    if beta.shape != (p,) or not np.all(np.isfinite(beta)):
        raise ValueError("coefficients must be a finite vector matching design_matrix columns")
    if not np.all(np.isfinite(x)):
        raise ValueError("design_matrix must contain only finite values")
    if isinstance(population_size, bool) or not isinstance(population_size, int):
        raise ValueError("population_size must be an integer")
    if population_size < n:
        raise ValueError("population_size must be at least the validation sample size")

    y = y_raw.astype(float)
    if np.linalg.matrix_rank(x) != p:
        raise ValueError("design_matrix must have full column rank")
    eta = x @ beta
    probability = expit(eta)
    variance = probability * (1.0 - probability)
    sensitivity = (x.T * variance) @ x / n
    condition = float(np.linalg.cond(sensitivity))
    if not np.isfinite(condition) or condition > 1e12:
        raise ValueError("logistic sensitivity is singular or ill-conditioned")

    score = x * (y - probability)[:, None]
    centered = score - np.mean(score, axis=0, keepdims=True)
    score_covariance = centered.T @ centered / (n - 1)
    fpc = 1.0 - n / population_size
    left_solved = np.linalg.solve(sensitivity, score_covariance)
    covariance = (
        np.linalg.solve(sensitivity, left_solved.T).T
        * (fpc / n)
    )
    covariance = 0.5 * (covariance + covariance.T)
    eigenvalues = np.linalg.eigvalsh(covariance)
    if not np.all(np.isfinite(covariance)) or np.min(eigenvalues) < -1e-10:
        raise ArithmeticError("design covariance is non-finite or materially non-positive-semidefinite")

    return SRSWORLogisticVariance(
        covariance=covariance,
        sensitivity=sensitivity,
        score_covariance=score_covariance,
        fpc=float(fpc),
        n_validation=int(n),
        population_size=int(population_size),
        design_condition_number=condition,
    )
