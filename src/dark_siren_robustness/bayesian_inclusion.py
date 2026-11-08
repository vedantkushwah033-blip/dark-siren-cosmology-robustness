"""Bayesian logistic inclusion posterior via fixed-basis t importance sampling.

The target is the validation-only modular posterior. Importance samples retain
their normalized importance weights; they are not mislabeled as equally
weighted posterior draws. This module is synthetic-methodology infrastructure.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.optimize import minimize
from scipy.special import expit, logsumexp
from scipy.stats import multivariate_t


@dataclass(frozen=True)
class InclusionPosteriorApproximation:
    status: str
    message: str
    degree: int
    z_min: float
    z_max: float
    prior_mean: NDArray[np.float64]
    prior_sd: NDArray[np.float64]
    map_coefficients: NDArray[np.float64] | None
    laplace_covariance: NDArray[np.float64] | None
    coefficient_proposals: NDArray[np.float64] | None
    normalized_importance_weights: NDArray[np.float64] | None
    log_importance_weights: NDArray[np.float64] | None
    posterior_mean: NDArray[np.float64] | None
    posterior_covariance: NDArray[np.float64] | None
    importance_ess: float
    max_normalized_weight: float
    log_evidence_estimate: float
    log_evidence_relative_mcse: float
    n_draws: int
    newton_iterations: int

    def predict(self, redshift: ArrayLike) -> NDArray[np.float64]:
        """Return posterior-mean inclusion probabilities at supplied redshifts."""
        if self.status != "SUCCESS" or self.coefficient_proposals is None:
            raise RuntimeError(f"Cannot predict from posterior approximation: {self.status}")
        z = np.asarray(redshift, dtype=float)
        if z.ndim != 1 or not np.all(np.isfinite(z)):
            raise ValueError("redshift must be a finite one-dimensional array")
        if np.any((z < self.z_min) | (z > self.z_max)):
            raise ValueError("redshift lies outside the fixed prior domain")
        u = (z - self.z_min) / (self.z_max - self.z_min)
        design = np.column_stack([u**power for power in range(self.degree + 1)])
        probabilities = expit(self.coefficient_proposals @ design.T)
        return np.asarray(self.normalized_importance_weights @ probabilities, dtype=float)


def fixed_redshift_design(
    redshift: ArrayLike, *, degree: int, z_min: float = 0.01, z_max: float = 0.40
) -> NDArray[np.float64]:
    """Construct the fixed domain basis [1,u,...,u^degree], u in [0,1]."""
    z = np.asarray(redshift, dtype=float)
    if z.ndim != 1 or z.size < 1 or not np.all(np.isfinite(z)):
        raise ValueError("redshift must be a nonempty finite one-dimensional array")
    if degree not in (1, 2):
        raise ValueError("degree must be 1 or 2")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or z_min >= z_max:
        raise ValueError("require finite z_min < z_max")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("redshift must lie inside the fixed domain")
    u = (z - z_min) / (z_max - z_min)
    return np.column_stack([u**power for power in range(degree + 1)])


def _log_posterior_and_derivatives(beta, design, labels, prior_mean, prior_sd):
    eta = design @ beta
    log_likelihood = float(np.sum(labels * eta - np.logaddexp(0.0, eta)))
    standardized = (beta - prior_mean) / prior_sd
    log_prior = float(
        -0.5 * np.sum(standardized**2)
        - np.sum(np.log(prior_sd))
        - 0.5 * beta.size * np.log(2.0 * np.pi)
    )
    probability = expit(eta)
    gradient = design.T @ (labels - probability) - (beta - prior_mean) / prior_sd**2
    information = design.T @ (
        (probability * (1.0 - probability))[:, None] * design
    ) + np.diag(1.0 / prior_sd**2)
    return log_likelihood + log_prior, gradient, information


def _find_map(design, labels, prior_mean, prior_sd, max_iterations):
    beta = prior_mean.copy()
    iterations = 0
    for iterations in range(1, max_iterations + 1):
        current, score, information = _log_posterior_and_derivatives(
            beta, design, labels, prior_mean, prior_sd
        )
        if not np.isfinite(current) or not np.all(np.isfinite(score)):
            raise ArithmeticError("Non-finite posterior objective or score")
        if np.linalg.norm(score, ord=np.inf) <= 1e-9:
            return beta, information, iterations
        eigenvalues = np.linalg.eigvalsh(information)
        if not np.all(np.isfinite(eigenvalues)) or np.min(eigenvalues) <= 0:
            raise ArithmeticError("Posterior information is not positive definite")
        direction = np.linalg.solve(information, score)
        decrement = float(score @ direction)
        if not np.isfinite(decrement) or decrement <= 0:
            raise ArithmeticError("Newton decrement is not positive")
        step = 1.0
        accepted = False
        while step >= 2.0**-30:
            candidate = beta + step * direction
            candidate_logp = _log_posterior_and_derivatives(
                candidate, design, labels, prior_mean, prior_sd
            )[0]
            if candidate_logp >= current + 1e-4 * step * decrement:
                beta = candidate
                accepted = True
                break
            step *= 0.5
        if not accepted:
            if np.linalg.norm(score, ord=np.inf) <= 1e-7:
                return beta, information, iterations
            raise ArithmeticError("Damped Newton line search failed")
    raise ArithmeticError("Posterior mode did not converge within max_iterations")


def importance_sample_inclusion_posterior(
    validation_redshift: ArrayLike,
    retained: ArrayLike,
    *,
    degree: int = 1,
    z_min: float = 0.01,
    z_max: float = 0.40,
    prior_mean: ArrayLike | None = None,
    prior_sd: ArrayLike | float = 1.0,
    n_draws: int = 20000,
    seed: int = 20261108,
    proposal_df: float = 7.0,
    proposal_scale: float = 1.5,
    min_ess: float = 200.0,
    min_ess_fraction: float = 0.05,
    max_newton_iterations: int = 200,
) -> InclusionPosteriorApproximation:
    """Approximate p(beta | validation) with self-normalized t importance sampling.

    The prior is defined in the fixed domain basis, not the validation-sample
    mean/SD basis used by the separate frequentist MLE fitter. The returned
    proposals must be integrated using their normalized importance weights.
    """
    z = np.asarray(validation_redshift, dtype=float)
    raw_labels = np.asarray(retained)
    if z.ndim != 1 or raw_labels.ndim != 1 or z.shape != raw_labels.shape:
        raise ValueError("validation_redshift and retained must be matching 1D arrays")
    if z.size < 5 or not np.all(np.isfinite(z)):
        raise ValueError("at least five finite validation redshifts are required")
    if np.any((raw_labels != 0) & (raw_labels != 1)):
        raise ValueError("retained labels must be binary (0 or 1)")
    if not isinstance(n_draws, int) or n_draws < 1000:
        raise ValueError("n_draws must be an integer >= 1000")
    if not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if not np.isfinite(proposal_df) or proposal_df <= 2:
        raise ValueError("proposal_df must be finite and > 2")
    if not np.isfinite(proposal_scale) or proposal_scale <= 0:
        raise ValueError("proposal_scale must be finite and positive")
    if not np.isfinite(min_ess) or min_ess < 1:
        raise ValueError("min_ess must be finite and >= 1")
    if not np.isfinite(min_ess_fraction) or not 0 < min_ess_fraction <= 1:
        raise ValueError("min_ess_fraction must lie in (0,1]")
    if not isinstance(max_newton_iterations, int) or max_newton_iterations < 1:
        raise ValueError("max_newton_iterations must be a positive integer")

    design = fixed_redshift_design(z, degree=degree, z_min=z_min, z_max=z_max)
    labels = raw_labels.astype(float)
    dimension = degree + 1
    mean = (
        np.zeros(dimension, dtype=float)
        if prior_mean is None else np.asarray(prior_mean, dtype=float)
    )
    sd_raw = np.asarray(prior_sd, dtype=float)
    sd = np.full(dimension, float(sd_raw)) if sd_raw.ndim == 0 else sd_raw.copy()
    if mean.shape != (dimension,) or sd.shape != (dimension,):
        raise ValueError("prior_mean and prior_sd must match the coefficient dimension")
    if not np.all(np.isfinite(mean)) or not np.all(np.isfinite(sd)) or np.any(sd <= 0):
        raise ValueError("prior parameters must be finite and prior_sd strictly positive")
    if np.unique(labels).size < 2:
        return InclusionPosteriorApproximation(
            status="NON_EVALUABLE", message="Validation labels contain only one outcome class.",
            degree=degree, z_min=float(z_min), z_max=float(z_max),
            prior_mean=mean.copy(), prior_sd=sd.copy(), map_coefficients=None,
            laplace_covariance=None, coefficient_proposals=None,
            normalized_importance_weights=None, log_importance_weights=None,
            posterior_mean=None, posterior_covariance=None, importance_ess=0.0,
            max_normalized_weight=float("inf"), log_evidence_estimate=float("nan"),
            log_evidence_relative_mcse=float("inf"), n_draws=n_draws, newton_iterations=0,
        )

    try:
        mode, information, iterations = _find_map(
            design, labels, mean, sd, max_newton_iterations
        )
        covariance = np.linalg.inv(information)
        covariance = 0.5 * (covariance + covariance.T)
        shape = proposal_scale * covariance
        # Validate positive definiteness before drawing from the proposal.
        np.linalg.cholesky(shape)
        rng = np.random.default_rng(seed)
        proposals = np.asarray(
            multivariate_t.rvs(
                loc=mode, shape=shape, df=proposal_df, size=n_draws, random_state=rng
            ),
            dtype=float,
        ).reshape(n_draws, dimension)

        eta = proposals @ design.T
        log_likelihood = np.sum(
            labels[None, :] * eta - np.logaddexp(0.0, eta), axis=1
        )
        standardized = (proposals - mean[None, :]) / sd[None, :]
        log_prior = (
            -0.5 * np.sum(standardized**2, axis=1)
            - np.sum(np.log(sd))
            - 0.5 * dimension * np.log(2.0 * np.pi)
        )
        log_target = log_likelihood + log_prior
        log_proposal = np.asarray(
            multivariate_t.logpdf(
                proposals, loc=mode, shape=shape, df=proposal_df
            ),
            dtype=float,
        )
        log_weights = log_target - log_proposal
        log_weight_sum = float(logsumexp(log_weights))
        normalized_weights = np.exp(log_weights - log_weight_sum)
        ess = float(1.0 / np.sum(normalized_weights**2))
        max_weight = float(np.max(normalized_weights))
        posterior_mean = normalized_weights @ proposals
        centered = proposals - posterior_mean
        posterior_covariance = centered.T @ (
            centered * normalized_weights[:, None]
        )
        log_evidence = log_weight_sum - np.log(n_draws)
        scaled = np.exp(log_weights - np.max(log_weights))
        relative_mcse = float(
            np.std(scaled, ddof=1) / np.sqrt(n_draws) / np.mean(scaled)
        )
        status = "SUCCESS"
        message = "Fixed-basis validation posterior approximated by weighted t importance sampling."
        if ess < max(min_ess, min_ess_fraction * n_draws):
            status = "NON_EVALUABLE"
            message = (
                "Importance proposal produced insufficient effective sample size; "
                "increase proposal draws or revise the proposal."
            )
        return InclusionPosteriorApproximation(
            status=status, message=message, degree=degree, z_min=float(z_min),
            z_max=float(z_max), prior_mean=mean.copy(), prior_sd=sd.copy(),
            map_coefficients=mode.copy(), laplace_covariance=covariance.copy(),
            coefficient_proposals=proposals,
            normalized_importance_weights=normalized_weights,
            log_importance_weights=log_weights,
            posterior_mean=posterior_mean,
            posterior_covariance=posterior_covariance,
            importance_ess=ess, max_normalized_weight=max_weight,
            log_evidence_estimate=float(log_evidence),
            log_evidence_relative_mcse=relative_mcse,
            n_draws=n_draws, newton_iterations=iterations,
        )
    except (ArithmeticError, FloatingPointError, np.linalg.LinAlgError, ValueError) as exc:
        return InclusionPosteriorApproximation(
            status="FAILED", message=str(exc), degree=degree, z_min=float(z_min),
            z_max=float(z_max), prior_mean=mean.copy(), prior_sd=sd.copy(),
            map_coefficients=None, laplace_covariance=None,
            coefficient_proposals=None, normalized_importance_weights=None,
            log_importance_weights=None, posterior_mean=None,
            posterior_covariance=None, importance_ess=0.0,
            max_normalized_weight=float("inf"), log_evidence_estimate=float("nan"),
            log_evidence_relative_mcse=float("inf"), n_draws=n_draws,
            newton_iterations=0,
        )


def weighted_log_marginal_likelihood(
    log_likelihood_draws: ArrayLike, log_importance_weights: ArrayLike
) -> NDArray[np.float64]:
    """Compute log E_p(theta|V)[L(H0;theta)] from weighted proposal draws.

    Inputs have shape (draw,) or (draw, H0). Importance weights may be
    unnormalized; self-normalization is performed stably in log space.
    """
    log_likelihood = np.asarray(log_likelihood_draws, dtype=float)
    log_weights = np.asarray(log_importance_weights, dtype=float)
    if log_likelihood.ndim not in (1, 2) or log_likelihood.shape[0] < 1:
        raise ValueError("log_likelihood_draws must have shape (draw,) or (draw,H0)")
    if log_weights.ndim != 1 or log_weights.shape[0] != log_likelihood.shape[0]:
        raise ValueError("log_importance_weights must match the draw dimension")
    if not np.all(np.isfinite(log_likelihood)) or not np.all(np.isfinite(log_weights)):
        raise ValueError("log likelihoods and log importance weights must be finite")
    denominator = logsumexp(log_weights)
    if log_likelihood.ndim == 1:
        return np.asarray(logsumexp(log_likelihood + log_weights) - denominator)
    return np.asarray(
        logsumexp(log_likelihood + log_weights[:, None], axis=0) - denominator,
        dtype=float,
    )
