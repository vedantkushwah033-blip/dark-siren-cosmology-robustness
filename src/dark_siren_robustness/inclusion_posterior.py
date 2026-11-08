"""Proper-prior logistic validation posterior with multi-chain random-walk MH.

This experimental sampler does not reuse the MLE covariance as a posterior.
The Gaussian prior is part of the target distribution; the inverse Hessian at
the MAP is used only to precondition proposals.
"""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.optimize import minimize
from scipy.special import expit


@dataclass(frozen=True)
class InclusionPosterior:
    status: str
    message: str
    degree: int
    z_min: float
    z_max: float
    prior_sd: float
    map_coefficients: np.ndarray | None
    map_covariance: np.ndarray | None
    draws: np.ndarray | None
    acceptance_rates: np.ndarray | None
    rhat: np.ndarray | None
    ess: np.ndarray | None
    design_rank: int
    design_condition_number: float
    warmup: int
    draws_per_chain: int
    seed: int

    def predict(self, redshift):
        """Predict inclusion probabilities for every retained posterior draw.

        Returns shape (chain, draw, redshift), using the same fixed-domain
        coordinate for validation and analysis-catalogue redshifts.
        """
        if self.draws is None:
            raise RuntimeError(f"Cannot predict from posterior status={self.status}")
        z = np.asarray(redshift, dtype=float)
        if z.ndim != 1 or not np.all(np.isfinite(z)):
            raise ValueError("redshift must be a finite one-dimensional array")
        x = (z - self.z_min) / (self.z_max - self.z_min)
        design = np.column_stack([x**k for k in range(self.degree + 1)])
        return expit(np.einsum("cdp,np->cdn", self.draws, design))


def _design(redshift, degree, z_min, z_max):
    """Fixed-coordinate polynomial basis; independent of validation moments."""
    x = (redshift - z_min) / (z_max - z_min)
    return np.column_stack([x**k for k in range(degree + 1)])


def _split_rhat(chains):
    """Classical split-Rhat; chains have shape (chain, draw, parameter)."""
    _, n_draw, _ = chains.shape
    half = n_draw // 2
    split = np.concatenate([chains[:, :half, :], chains[:, -half:, :]], axis=0)
    _, n_half, _ = split.shape
    means = split.mean(axis=1)
    variances = split.var(axis=1, ddof=1)
    within = variances.mean(axis=0)
    between = n_half * means.var(axis=0, ddof=1)
    variance_hat = ((n_half - 1) / n_half) * within + between / n_half
    with np.errstate(divide="ignore", invalid="ignore"):
        result = np.sqrt(variance_hat / within)
    zero = within == 0
    result[zero] = np.where(between[zero] == 0, 1.0, np.inf)
    return result


def _ess_initial_positive(chains):
    """Approximate ESS using the initial-positive adjacent-autocorrelation sum."""
    n_chain, n_draw, n_parameter = chains.shape
    result = np.empty(n_parameter)
    fft_size = 1 << (2 * n_draw - 1).bit_length()
    for parameter in range(n_parameter):
        autocorrelations = []
        for chain in range(n_chain):
            values = chains[chain, :, parameter]
            values = values - values.mean()
            if float(values @ values) <= 0:
                autocorrelations.append(np.r_[1.0, np.zeros(n_draw - 1)])
                continue
            spectrum = np.fft.rfft(values, n=fft_size)
            autocovariance = np.fft.irfft(
                spectrum * np.conjugate(spectrum), n=fft_size
            )[:n_draw]
            autocovariance /= np.arange(n_draw, 0, -1)
            autocorrelations.append(autocovariance / autocovariance[0])
        rho = np.mean(autocorrelations, axis=0)
        positive_sum = 0.0
        for lag in range(1, n_draw - 1, 2):
            pair = float(rho[lag] + rho[lag + 1])
            if not np.isfinite(pair) or pair <= 0:
                break
            positive_sum += pair
        integrated_time = max(1.0, 1.0 + 2.0 * positive_sum)
        result[parameter] = n_chain * n_draw / integrated_time
    return result


def sample_inclusion_posterior(
    validation_redshift,
    retained,
    *,
    degree=1,
    prior_sd=1.0,
    z_min=0.01,
    z_max=0.40,
    n_chains=4,
    warmup=2000,
    draws_per_chain=2000,
    seed=20261108,
    target_acceptance=0.234,
    max_rhat=1.05,
    min_ess=100.0,
):
    """Sample p(theta | validation) under independent N(0, prior_sd^2) priors.

    The polynomial basis uses the fixed coordinate u=(z-z_min)/(z_max-z_min),
    not validation-sample standardization. The proper Gaussian prior is thus
    invariant to the realized validation sample. Returned draws are correlated
    MCMC draws, not independent samples; classical split-Rhat, approximate
    initial-positive-sequence ESS, and post-warmup acceptance are reported.
    """
    z = np.asarray(validation_redshift, dtype=float)
    labels = np.asarray(retained)
    if z.ndim != 1 or labels.ndim != 1 or z.shape != labels.shape or z.size < 5:
        raise ValueError("validation inputs must be matching 1D arrays with at least five rows")
    if not np.all(np.isfinite(z)) or np.any((labels != 0) & (labels != 1)):
        raise ValueError("redshifts must be finite and labels binary")
    if degree not in (1, 2):
        raise ValueError("degree must be 1 or 2")
    if not np.isfinite(prior_sd) or prior_sd <= 0:
        raise ValueError("prior_sd must be finite and positive")
    if not np.isfinite(z_min) or not np.isfinite(z_max) or z_min >= z_max:
        raise ValueError("fixed redshift domain must be finite with z_min < z_max")
    if np.any((z < z_min) | (z > z_max)):
        raise ValueError("validation redshifts must lie inside the fixed prior domain")
    if n_chains < 4 or warmup < 500 or draws_per_chain < 500:
        raise ValueError("use at least four chains and 500 warmup/draws per chain")
    if not 0.05 < target_acceptance < 0.8:
        raise ValueError("target_acceptance must lie in (0.05, 0.8)")
    if not np.isfinite(max_rhat) or max_rhat <= 1:
        raise ValueError("max_rhat must be finite and greater than one")
    if not np.isfinite(min_ess) or min_ess <= 0:
        raise ValueError("min_ess must be finite and positive")

    y = labels.astype(float)
    design = _design(z, degree, z_min, z_max)
    rank = int(np.linalg.matrix_rank(design))
    condition = float(np.linalg.cond(design))
    if rank != degree + 1 or not np.isfinite(condition) or condition > 1e10:
        return InclusionPosterior(
            "NON_EVALUABLE", "Validation design is rank deficient or ill-conditioned.",
            degree, z_min, z_max, prior_sd, None, None, None, None, None, None,
            rank, condition, warmup, draws_per_chain, seed,
        )
    # The proper Gaussian prior regularizes complete/quasi separation,
    # including all-zero or all-one validation labels. Posterior diagnostics
    # below determine evaluability; a single observed class is not itself a
    # mathematical failure of this Bayesian model.

    prior_precision = 1.0 / prior_sd**2

    def log_posterior(beta):
        eta = design @ beta
        log_likelihood = -np.sum(np.logaddexp(0.0, eta) - y * eta)
        log_prior = -0.5 * prior_precision * float(beta @ beta)
        return float(log_likelihood + log_prior)

    def gradient(beta):
        return -(design.T @ (expit(design @ beta) - y)) - prior_precision * beta

    optimum = minimize(
        lambda beta: -log_posterior(beta),
        np.zeros(degree + 1),
        jac=lambda beta: -gradient(beta),
        method="BFGS",
        options={"gtol": 1e-10, "maxiter": 3000},
    )
    mode = np.asarray(optimum.x, dtype=float)
    score_norm = float(np.linalg.norm(gradient(mode), ord=np.inf))
    probability = expit(design @ mode)
    hessian = (
        design.T @ ((probability * (1.0 - probability))[:, None] * design)
        + prior_precision * np.eye(degree + 1)
    )
    eigenvalues = np.linalg.eigvalsh(hessian)
    if (
        not np.all(np.isfinite(mode))
        or score_norm > 1e-6
        or not np.all(np.isfinite(eigenvalues))
        or np.min(eigenvalues) <= 0
    ):
        return InclusionPosterior(
            "FAILED", f"MAP optimization/Hessian check failed (score={score_norm:g}).",
            degree, z_min, z_max, prior_sd, None, None, None, None, None, None,
            rank, condition, warmup, draws_per_chain, seed,
        )

    map_covariance = np.linalg.inv(hessian)
    cholesky = np.linalg.cholesky(map_covariance)
    dimension = degree + 1
    initial_scale = 2.38 / np.sqrt(dimension)
    seed_sequences = np.random.SeedSequence(int(seed)).spawn(n_chains)
    chains = np.empty((n_chains, draws_per_chain, dimension), dtype=float)
    acceptance_rates = np.zeros(n_chains, dtype=float)

    for chain_index, sequence in enumerate(seed_sequences):
        rng = np.random.default_rng(sequence)
        beta = mode + 0.1 * (cholesky @ rng.normal(size=dimension))
        current_log_posterior = log_posterior(beta)
        log_scale = float(np.log(initial_scale))
        accepted_warmup_window = 0
        accepted_sampling = 0

        for iteration in range(warmup + draws_per_chain):
            proposal = beta + np.exp(log_scale) * (cholesky @ rng.normal(size=dimension))
            proposal_log_posterior = log_posterior(proposal)
            if np.log(rng.random()) < proposal_log_posterior - current_log_posterior:
                beta = proposal
                current_log_posterior = proposal_log_posterior
                if iteration < warmup:
                    accepted_warmup_window += 1
                else:
                    accepted_sampling += 1

            if iteration < warmup and (iteration + 1) % 50 == 0:
                batch = (iteration + 1) // 50
                observed_rate = accepted_warmup_window / 50.0
                gain = min(0.5, 1.0 / np.sqrt(batch))
                log_scale += gain * (observed_rate - target_acceptance)
                log_scale = float(np.clip(log_scale, np.log(0.03), np.log(10.0)))
                accepted_warmup_window = 0
            elif iteration >= warmup:
                chains[chain_index, iteration - warmup, :] = beta

        acceptance_rates[chain_index] = accepted_sampling / draws_per_chain

    rhat = _split_rhat(chains)
    ess = _ess_initial_positive(chains)
    if (
        not np.all(np.isfinite(rhat))
        or not np.all(np.isfinite(ess))
        or np.max(rhat) > max_rhat
        or np.min(ess) < min_ess
    ):
        status = "NON_EVALUABLE"
        message = (
            f"MCMC diagnostics failed: max split-Rhat={np.max(rhat):.4f}, "
            f"min approximate ESS={np.min(ess):.1f}."
        )
    else:
        status = "SUCCESS"
        message = "Proper-prior logistic posterior sampled; declared MCMC diagnostics passed."

    return InclusionPosterior(
        status, message, degree, z_min, z_max, prior_sd, mode, map_covariance,
        chains, acceptance_rates, rhat, ess, rank, condition, warmup,
        draws_per_chain, int(seed),
    )
