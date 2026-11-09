"""Independent reference reconstruction of the validation logistic posterior.

This module intentionally does not import the production inclusion-posterior
fitter. It uses a separate BFGS mode search and independently reconstructs
the Gaussian proposal, tensor Hermite nodes, and target/proposal weights.
It is an audit implementation, not a second production inference API.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize
from scipy.special import expit, logsumexp, roots_hermitenorm


@dataclass(frozen=True)
class ReferencePosterior:
    nodes: np.ndarray
    weights: np.ndarray
    degree: int
    center: float
    scale: float
    mode: np.ndarray
    covariance: np.ndarray
    posterior_mean: np.ndarray
    posterior_covariance: np.ndarray
    log_normalizer: float
    mode_score_inf: float


def reconstruct_logistic_posterior(
    redshift: ArrayLike,
    retained: ArrayLike,
    *,
    degree: int,
    prior_mean: ArrayLike,
    prior_sd: ArrayLike,
    center: float,
    scale: float,
    order: int,
) -> ReferencePosterior:
    """Rebuild the proper Gaussian-prior logistic posterior from raw labels.

    The validation likelihood is Bernoulli-logistic and the coefficient prior
    is independent Normal. Coordinates are fixed externally by the analysis
    contract. The returned quadrature is a Laplace-proposal importance rule;
    normalized weights target the full logistic posterior, not the Gaussian
    approximation.
    """
    z = np.asarray(redshift, dtype=np.float64)
    y_raw = np.asarray(retained)
    mu = np.asarray(prior_mean, dtype=np.float64)
    sd = np.asarray(prior_sd, dtype=np.float64)

    if z.ndim != 1 or y_raw.ndim != 1 or z.shape != y_raw.shape or z.size < 5:
        raise ValueError("redshift and retained must be matching vectors with n >= 5")
    if not np.all(np.isfinite(z)) or not np.all(np.isfinite(y_raw)):
        raise ValueError("validation inputs must be finite")
    if np.any((y_raw != 0) & (y_raw != 1)):
        raise ValueError("retained labels must be binary")
    if degree not in (1, 2):
        raise ValueError("degree must be one or two")
    dimension = degree + 1
    if mu.shape != (dimension,) or sd.shape != (dimension,):
        raise ValueError("prior vectors must match the polynomial dimension")
    if not np.all(np.isfinite(mu)) or not np.all(np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("prior must be finite with strictly positive scales")
    if not np.isfinite(center) or not np.isfinite(scale) or scale <= 0:
        raise ValueError("coordinate center/scale are invalid")
    if isinstance(order, bool) or not isinstance(order, int) or order < 3:
        raise ValueError("quadrature order must be an integer >= 3")

    y = y_raw.astype(np.float64, copy=False)
    x = (z - float(center)) / float(scale)
    design = np.vander(x, N=dimension, increasing=True)
    if np.linalg.matrix_rank(design) != dimension:
        raise ValueError("validation design is rank deficient")
    precision = np.reciprocal(sd * sd)

    # Independent objective and score, expressed directly in matrix form.
    def negative_log_posterior(beta: np.ndarray) -> float:
        eta = design @ beta
        likelihood = np.sum(y * eta - np.logaddexp(0.0, eta))
        delta = (beta - mu) / sd
        prior = -0.5 * np.dot(delta, delta)
        return float(-(likelihood + prior))

    def score(beta: np.ndarray) -> np.ndarray:
        p = expit(design @ beta)
        return design.T @ (p - y) + precision * (beta - mu)

    # BFGS is deliberately distinct from the production trust-exact /
    # L-BFGS-B sequence. Its success flag is not trusted without a score check.
    result = minimize(
        negative_log_posterior,
        mu.copy(),
        jac=score,
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 5000},
    )
    mode = np.asarray(result.x, dtype=np.float64)
    score_norm = float(np.linalg.norm(score(mode), ord=np.inf))
    if not np.all(np.isfinite(mode)) or not np.isfinite(score_norm) or score_norm > 2e-7:
        raise ArithmeticError(
            f"independent BFGS mode failed score check: success={result.success}, "
            f"message={result.message}, score_inf={score_norm:.6g}"
        )

    probability = expit(design @ mode)
    curvature = probability * (1.0 - probability)
    observed_information = (
        design.T @ (curvature[:, None] * design)
        + np.diag(precision)
    )
    eigenvalues = np.linalg.eigvalsh(observed_information)
    if not np.all(np.isfinite(eigenvalues)) or np.min(eigenvalues) <= 0:
        raise ArithmeticError("reference observed information is not positive definite")
    covariance = np.linalg.solve(observed_information, np.eye(dimension))
    covariance = 0.5 * (covariance + covariance.T)
    lower = np.linalg.cholesky(covariance)

    # Construct the tensor rule independently from Cartesian multi-indices.
    abscissae, univariate_weights = roots_hermitenorm(order)
    multi_indices = tuple(product(range(order), repeat=dimension))
    standard = np.asarray(
        [[abscissae[j] for j in index] for index in multi_indices],
        dtype=np.float64,
    )
    log_rule_weights = np.asarray(
        [
            sum(math.log(float(univariate_weights[j])) for j in index)
            - 0.5 * dimension * math.log(2.0 * math.pi)
            for index in multi_indices
        ],
        dtype=np.float64,
    )
    nodes = mode[None, :] + standard @ lower.T

    log_abs_det_lower = float(np.log(np.diag(lower)).sum())
    log_standard_normal = (
        -0.5 * dimension * math.log(2.0 * math.pi)
        -0.5 * np.einsum("ij,ij->i", standard, standard)
    )
    log_proposal_density = log_standard_normal - log_abs_det_lower

    eta_nodes = design @ nodes.T
    log_likelihood_nodes = np.sum(
        y[:, None] * eta_nodes - np.logaddexp(0.0, eta_nodes),
        axis=0,
    )
    standardized = (nodes - mu[None, :]) / sd[None, :]
    log_prior_nodes = (
        -0.5 * np.einsum("ij,ij->i", standardized, standardized)
        -float(np.log(sd).sum())
        -0.5 * dimension * math.log(2.0 * math.pi)
    )
    log_target = log_likelihood_nodes + log_prior_nodes
    log_terms = log_rule_weights + log_target - log_proposal_density
    log_normalizer = float(logsumexp(log_terms))
    weights = np.exp(log_terms - log_normalizer)
    total = float(weights.sum())
    if not np.isfinite(total) or total <= 0 or not np.all(np.isfinite(weights)):
        raise ArithmeticError("reference posterior weights are invalid")
    weights /= total

    mean = weights @ nodes
    deviations = nodes - mean
    posterior_covariance = (deviations * weights[:, None]).T @ deviations
    return ReferencePosterior(
        nodes=nodes,
        weights=weights,
        degree=degree,
        center=float(center),
        scale=float(scale),
        mode=mode,
        covariance=covariance,
        posterior_mean=mean,
        posterior_covariance=posterior_covariance,
        log_normalizer=log_normalizer,
        mode_score_inf=score_norm,
    )
