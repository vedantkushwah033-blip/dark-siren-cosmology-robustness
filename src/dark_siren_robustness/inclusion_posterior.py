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


@dataclass(frozen=True)
class InclusionPosteriorQuadratureLadder:
    """Order-refinement record; CONVERGED is a numerical, not scientific, claim."""

    fits: tuple[InclusionPosteriorQuadrature, ...]
    status: str
    mean_max_abs_change: tuple[float, ...]
    covariance_max_abs_change: tuple[float, ...]
    predictive_max_abs_change: tuple[float, ...]


def fit_inclusion_posterior_quadrature_ladder(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    *,
    orders: tuple[int, ...] = (7, 11, 15, 19, 23),
    degree: int = 1,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | None = None,
    predictive_redshift: ArrayLike | None = None,
    mean_tolerance: float = 1e-3,
    covariance_tolerance: float = 1e-3,
    predictive_tolerance: float = 1e-3,
    consecutive_passes: int = 2,
) -> InclusionPosteriorQuadratureLadder:
    """Fit increasing tensor-GH orders and require repeated summary stability.

    The convergence decision compares posterior coefficient moments and
    posterior-mean inclusion probabilities on a fixed prediction grid. It does
    not establish convergence of a downstream H0 likelihood; that must be
    checked after the shared-event likelihood is implemented.
    """
    order_tuple = tuple(orders)
    if len(order_tuple) < 2 or any(
        isinstance(o, bool) or not isinstance(o, int) or o < 3 for o in order_tuple
    ) or any(b <= a for a, b in zip(order_tuple, order_tuple[1:])):
        raise ValueError("orders must be a strictly increasing tuple of integers >= 3")
    for value, label in (
        (mean_tolerance, "mean_tolerance"),
        (covariance_tolerance, "covariance_tolerance"),
        (predictive_tolerance, "predictive_tolerance"),
    ):
        if not np.isfinite(value) or value <= 0:
            raise ValueError(f"{label} must be positive and finite")
    if isinstance(consecutive_passes, bool) or not isinstance(consecutive_passes, int) or consecutive_passes < 1:
        raise ValueError("consecutive_passes must be a positive integer")

    z = np.asarray(validation_redshift, dtype=float)
    prediction_grid = (
        np.linspace(float(np.min(z)), float(np.max(z)), 101)
        if predictive_redshift is None else np.asarray(predictive_redshift, dtype=float)
    )
    if prediction_grid.ndim != 1 or not np.all(np.isfinite(prediction_grid)):
        raise ValueError("predictive_redshift must be a finite one-dimensional array")

    fits = []
    mean_changes = []
    covariance_changes = []
    predictive_changes = []
    passing_streak = 0
    converged = False
    for order_value in order_tuple:
        fit = fit_inclusion_posterior_quadrature(
            validation_redshift, retained, degree=degree, prior_mean=prior_mean,
            prior_sd=prior_sd, order=order_value,
        )
        fits.append(fit)
        if len(fits) == 1:
            continue
        previous, current = fits[-2], fits[-1]
        dm = float(np.max(np.abs(current.posterior_mean - previous.posterior_mean)))
        dc = float(np.max(np.abs(current.posterior_covariance - previous.posterior_covariance)))
        dp = float(np.max(np.abs(
            posterior_predictive_inclusion(current, prediction_grid)
            - posterior_predictive_inclusion(previous, prediction_grid)
        )))
        mean_changes.append(dm)
        covariance_changes.append(dc)
        predictive_changes.append(dp)
        passed = dm <= mean_tolerance and dc <= covariance_tolerance and dp <= predictive_tolerance
        passing_streak = passing_streak + 1 if passed else 0
        if passing_streak >= consecutive_passes:
            converged = True
            break

    return InclusionPosteriorQuadratureLadder(
        fits=tuple(fits),
        status="CONVERGED" if converged else "NON_CONVERGED",
        mean_max_abs_change=tuple(mean_changes),
        covariance_max_abs_change=tuple(covariance_changes),
        predictive_max_abs_change=tuple(predictive_changes),
    )


@dataclass(frozen=True)
class NodeSpecificHostWeights:
    """Stable inverse-inclusion weights and diagnostics for each posterior node."""

    inclusion_probability: np.ndarray
    raw_inverse_weight: np.ndarray
    normalized_weight: np.ndarray
    kish_ess: np.ndarray
    entropy_ess: np.ndarray
    maximum_normalized_weight: np.ndarray
    top_one_percent_mass: np.ndarray
    top_five_percent_mass: np.ndarray
    minimum_inclusion_probability: np.ndarray
    maximum_raw_inverse_weight: np.ndarray


def build_node_specific_host_weights(
    quadrature: InclusionPosteriorQuadrature,
    catalogue_redshift: ArrayLike,
) -> NodeSpecificHostWeights:
    """Construct untruncated IPW host weights for every validation-posterior node.

    Normalization is performed in log space. Raw inverse weights are retained
    as diagnostics; if they cannot be represented in float64, the routine
    fails explicitly rather than clipping probabilities or weights.
    """
    z = np.asarray(catalogue_redshift, dtype=float)
    if z.ndim != 1 or z.size == 0 or not np.all(np.isfinite(z)):
        raise ValueError("catalogue_redshift must be a non-empty finite 1D array")
    if quadrature.nodes.ndim != 2 or quadrature.weights.shape != (quadrature.nodes.shape[0],):
        raise ValueError("quadrature nodes and weights have inconsistent shapes")
    if quadrature.nodes.shape[1] != quadrature.degree + 1:
        raise ValueError("quadrature coefficient dimension does not match degree")
    if np.any(~np.isfinite(quadrature.nodes)) or np.any(~np.isfinite(quadrature.weights)):
        raise ValueError("quadrature nodes and weights must be finite")
    if np.any(quadrature.weights < 0) or not np.isclose(quadrature.weights.sum(), 1.0):
        raise ValueError("quadrature weights must be non-negative and sum to one")

    x = (z - quadrature.center) / quadrature.scale
    design = np.column_stack([x ** power for power in range(quadrature.degree + 1)])
    eta = quadrature.nodes @ design.T
    log_probability = -np.logaddexp(0.0, -eta)
    log_raw = -log_probability
    log_float_max = math.log(np.finfo(float).max)
    log_float_tiny = math.log(np.nextafter(0.0, 1.0))
    if np.any(log_probability < log_float_tiny):
        raise ArithmeticError("inclusion probability underflows float64; refusing to clip")
    if np.any(log_raw > log_float_max):
        raise ArithmeticError("raw inverse-inclusion weight overflows float64; refusing to clip")

    probability = np.exp(log_probability)
    raw = np.exp(log_raw)
    log_normalizer = logsumexp(log_raw, axis=1, keepdims=True)
    normalized = np.exp(log_raw - log_normalizer)
    kish = 1.0 / np.sum(normalized * normalized, axis=1)
    entropy = -np.sum(
        np.where(normalized > 0, normalized * np.log(np.maximum(normalized, np.finfo(float).tiny)), 0.0),
        axis=1,
    )
    entropy_ess = np.exp(entropy)
    maximum = np.max(normalized, axis=1)
    n = z.size
    k1 = max(1, int(math.ceil(0.01 * n)))
    k5 = max(1, int(math.ceil(0.05 * n)))
    top1 = np.sum(np.partition(normalized, n-k1, axis=1)[:, -k1:], axis=1)
    top5 = np.sum(np.partition(normalized, n-k5, axis=1)[:, -k5:], axis=1)
    return NodeSpecificHostWeights(
        inclusion_probability=probability,
        raw_inverse_weight=raw,
        normalized_weight=normalized,
        kish_ess=kish,
        entropy_ess=entropy_ess,
        maximum_normalized_weight=maximum,
        top_one_percent_mass=top1,
        top_five_percent_mass=top5,
        minimum_inclusion_probability=np.min(probability, axis=1),
        maximum_raw_inverse_weight=np.max(raw, axis=1),
    )
