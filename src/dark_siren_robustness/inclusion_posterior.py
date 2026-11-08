"""Proper Bayesian logistic inclusion posterior via adaptive Gauss-Hermite quadrature.

The posterior target is the validation Bernoulli likelihood times a proper
Gaussian prior. The Laplace Gaussian is only an importance proposal; normalized
quadrature weights correct it to the declared target posterior.
"""
from __future__ import annotations

from dataclasses import dataclass
from itertools import product
import math

import numpy as np
from numpy.typing import ArrayLike
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.special import roots_hermitenorm


@dataclass(frozen=True)
class InclusionPosteriorQuadrature:
    """Posterior quadrature nodes and normalized target-posterior weights."""

    nodes: np.ndarray
    weights: np.ndarray
    degree: int
    center: float
    scale: float
    order: int
    mode: np.ndarray
    laplace_covariance: np.ndarray
    log_evidence_unnormalized: float
    quadrature_ess: float
    maximum_normalized_weight: float
    posterior_mean: np.ndarray
    posterior_covariance: np.ndarray
    min_predicted_probability: float
    max_predicted_probability: float


def _validation_inputs(redshift: ArrayLike, retained: ArrayLike, degree: int):
    z = np.asarray(redshift, dtype=float)
    raw = np.asarray(retained)
    if z.ndim != 1 or raw.ndim != 1 or z.shape != raw.shape or z.size < 5:
        raise ValueError("redshift and retained must be matching 1D arrays with n >= 5")
    if not np.all(np.isfinite(z)):
        raise ValueError("redshift must be finite")
    if np.any((raw != 0) & (raw != 1)):
        raise ValueError("retained labels must be binary")
    if degree not in (1, 2):
        raise ValueError("degree must be 1 or 2")
    center = float(np.mean(z))
    scale = float(np.std(z, ddof=0))
    if not np.isfinite(scale) or scale <= 0:
        raise ValueError("validation redshifts must have positive finite spread")
    x = (z - center) / scale
    design = np.column_stack([x ** power for power in range(degree + 1)])
    if np.linalg.matrix_rank(design) != degree + 1:
        raise ValueError("validation design matrix is rank deficient")
    return z, raw.astype(float), design, center, scale


def _prior_vectors(dimension: int, prior_mean, prior_sd):
    mean = np.zeros(dimension) if prior_mean is None else np.asarray(prior_mean, float)
    if prior_sd is None:
        sd = np.asarray([5.0] + [2.5] * (dimension - 1), dtype=float)
    else:
        sd = np.asarray(prior_sd, dtype=float)
    if mean.shape != (dimension,) or sd.shape != (dimension,):
        raise ValueError("prior_mean and prior_sd must match coefficient dimension")
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("prior parameters must be finite and prior_sd strictly positive")
    return mean, sd


def _log_target(theta, design, retained, prior_mean, prior_sd):
    eta = design @ theta
    log_likelihood = np.sum(retained * eta - np.logaddexp(0.0, eta))
    standardized = (theta - prior_mean) / prior_sd
    log_prior = -0.5 * np.dot(standardized, standardized)
    log_prior -= np.sum(np.log(prior_sd)) + 0.5 * theta.size * math.log(2.0 * math.pi)
    return float(log_likelihood + log_prior)


def fit_inclusion_posterior_quadrature(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    *,
    degree: int = 1,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    order: int = 15,
    max_iterations: int = 2000,
) -> InclusionPosteriorQuadrature:
    """Approximate the proper logistic posterior with corrected Laplace-GH weights.

    Tensor Gauss-Hermite quadrature is used in the standardized coefficient
    space. This is intended for the project's 2- or 3-parameter models.
    """
    z, y, design, center, scale = _validation_inputs(
        validation_redshift, retained, degree
    )
    if isinstance(order, bool) or not isinstance(order, int) or order < 3:
        raise ValueError("order must be an integer >= 3")
    if isinstance(max_iterations, bool) or not isinstance(max_iterations, int) or max_iterations < 1:
        raise ValueError("max_iterations must be a positive integer")
    dimension = degree + 1
    mu, sd = _prior_vectors(dimension, prior_mean, prior_sd)
    prior_precision = 1.0 / (sd * sd)

    def objective(beta):
        return -_log_target(beta, design, y, mu, sd)

    def gradient(beta):
        probability = expit(design @ beta)
        return design.T @ (probability - y) + (beta - mu) * prior_precision

    def hessian(beta):
        probability = expit(design @ beta)
        curvature = probability * (1.0 - probability)
        return design.T @ (curvature[:, None] * design) + np.diag(prior_precision)

    result = minimize(
        objective, mu.copy(), jac=gradient, hess=hessian,
        method="trust-exact", options={"maxiter": max_iterations, "gtol": 1e-10},
    )
    mode = np.asarray(result.x, dtype=float)
    score_norm = float(np.linalg.norm(gradient(mode), ord=np.inf))
    if not np.all(np.isfinite(mode)) or not np.isfinite(score_norm) or score_norm > 1e-6:
        raise ArithmeticError(f"posterior mode did not converge: {result.message}; score={score_norm}")
    information = hessian(mode)
    eigenvalues = np.linalg.eigvalsh(information)
    if not np.all(np.isfinite(eigenvalues)) or np.min(eigenvalues) <= 0:
        raise ArithmeticError("posterior negative-log Hessian is not positive definite")
    covariance = np.linalg.inv(information)
    chol = np.linalg.cholesky(covariance)
    log_det_chol = float(np.sum(np.log(np.diag(chol))))

    gh_nodes, gh_weights = roots_hermitenorm(order)
    index_tuples = list(product(range(order), repeat=dimension))
    standard_nodes = np.asarray(
        [[gh_nodes[index] for index in indices] for indices in index_tuples],
        dtype=float,
    )
    base_log_weights = np.asarray(
        [sum(math.log(gh_weights[index]) for index in indices)
         - 0.5 * dimension * math.log(2.0 * math.pi)
         for indices in index_tuples],
        dtype=float,
    )
    nodes = mode[None, :] + standard_nodes @ chol.T

    # q(theta) is N(mode, covariance); the Jacobian is |chol|.
    log_proposal = (
        -0.5 * dimension * math.log(2.0 * math.pi)
        -0.5 * np.sum(standard_nodes * standard_nodes, axis=1)
        -log_det_chol
    )
    eta = design @ nodes.T
    log_likelihood = np.sum(
        y[:, None] * eta - np.logaddexp(0.0, eta), axis=0
    )
    standardized_prior = (nodes - mu[None, :]) / sd[None, :]
    log_prior = (
        -0.5 * np.sum(standardized_prior * standardized_prior, axis=1)
        -np.sum(np.log(sd))
        -0.5 * dimension * math.log(2.0 * math.pi)
    )
    log_target = log_likelihood + log_prior
    log_importance = base_log_weights + log_target - log_proposal
    log_evidence = float(logsumexp(log_importance))
    weights = np.exp(log_importance - log_evidence)
    if not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ArithmeticError("posterior quadrature weights are not finite and non-negative")
    weight_sum = float(weights.sum())
    if not np.isfinite(weight_sum) or weight_sum <= 0:
        raise ArithmeticError("posterior quadrature weights have invalid total mass")
    weights /= weight_sum

    posterior_mean = weights @ nodes
    centered = nodes - posterior_mean[None, :]
    posterior_covariance = (centered * weights[:, None]).T @ centered
    predicted = expit(design @ nodes.T)
    return InclusionPosteriorQuadrature(
        nodes=nodes,
        weights=weights,
        degree=degree,
        center=center,
        scale=scale,
        order=order,
        mode=mode,
        laplace_covariance=covariance,
        log_evidence_unnormalized=log_evidence,
        quadrature_ess=float(1.0 / np.sum(weights * weights)),
        maximum_normalized_weight=float(np.max(weights)),
        posterior_mean=posterior_mean,
        posterior_covariance=posterior_covariance,
        min_predicted_probability=float(np.min(predicted)),
        max_predicted_probability=float(np.max(predicted)),
    )


def posterior_predictive_inclusion(
    quadrature: InclusionPosteriorQuadrature, redshift: ArrayLike
) -> np.ndarray:
    """Return validation-posterior mean inclusion probability at supplied z."""
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or not np.all(np.isfinite(z)):
        raise ValueError("redshift must be a finite one-dimensional array")
    x = (z - quadrature.center) / quadrature.scale
    design = np.column_stack([x ** power for power in range(quadrature.degree + 1)])
    probabilities = expit(design @ quadrature.nodes.T)
    return probabilities @ quadrature.weights


def log_marginalize_shared_parameter(
    log_likelihood_by_node: ArrayLike,
    posterior_weights: ArrayLike,
) -> np.ndarray:
    """Marginalize a joint likelihood over one shared validation-posterior node.

    The first likelihood axis must index the same coefficient node across all
    events/other likelihood dimensions. This prevents event-wise independent
    marginalization of a shared catalogue-inclusion model.
    """
    log_likelihood = np.asarray(log_likelihood_by_node, dtype=float)
    weights = np.asarray(posterior_weights, dtype=float)
    if log_likelihood.ndim < 1 or log_likelihood.shape[0] != weights.size:
        raise ValueError("first likelihood dimension must match posterior weights")
    if weights.ndim != 1 or not np.all(np.isfinite(weights)) or np.any(weights < 0):
        raise ValueError("posterior_weights must be finite and non-negative")
    if not np.isclose(weights.sum(), 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError("posterior_weights must sum to one")
    if np.any(np.isnan(log_likelihood)) or np.any(np.isposinf(log_likelihood)):
        raise ValueError("log likelihoods may be finite or -inf, but not NaN/+inf")
    log_weights = np.full(weights.shape, -np.inf)
    positive = weights > 0
    log_weights[positive] = np.log(weights[positive])
    return logsumexp(log_likelihood + log_weights.reshape(
        (weights.size,) + (1,) * (log_likelihood.ndim - 1)
    ), axis=0)
