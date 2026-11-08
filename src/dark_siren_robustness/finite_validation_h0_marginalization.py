"""Selection-conditioned likelihood propagation over finite-validation draws."""
from __future__ import annotations

from dataclasses import dataclass

import numpy as np
from scipy.special import expit

from dark_siren_robustness.bayesian_inclusion import (
    InclusionPosteriorApproximation,
    fixed_redshift_design,
    weighted_log_marginal_likelihood,
)
from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.finite_validation_h0 import selected_log_likelihood_by_draw


@dataclass(frozen=True)
class FiniteValidationH0Marginalization:
    host_density_draws: np.ndarray
    unconditional_log_likelihood_by_draw: np.ndarray
    detection_efficiency_by_draw: np.ndarray
    selected_log_likelihood_by_draw: np.ndarray
    normalized_importance_weights: np.ndarray
    marginalized_selected_log_likelihood: np.ndarray
    importance_ess: float
    max_normalized_importance_weight: float
    marginalization_batch_log_likelihoods: np.ndarray | None
    marginalization_log_likelihood_mcse: np.ndarray | None
    marginalization_batch_ess: np.ndarray | None
    marginalization_mcse_status: str


def selection_correct_each_draw(
    unconditional_log_likelihood, detection_efficiency, *, n_detected
):
    """Apply beta_GW(H0; theta_b) correction within each coefficient draw."""
    log_likelihood = np.asarray(unconditional_log_likelihood, dtype=float)
    beta = np.asarray(detection_efficiency, dtype=float)
    if (
        log_likelihood.ndim != 2
        or log_likelihood.shape[0] < 1
        or log_likelihood.shape[1] < 2
    ):
        raise ValueError(
            "unconditional log likelihood must have shape (draws>=1, H0 grid>=2)"
        )
    if beta.shape != log_likelihood.shape:
        raise ValueError("draw-specific detection efficiencies must match likelihood shape")
    if not np.all(np.isfinite(log_likelihood)):
        raise ValueError("unconditional log likelihoods must be finite; invalid draws are not dropped")
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1.0 + 1e-12):
        raise ValueError("every draw-specific beta_GW must lie in (0,1]")
    if (
        isinstance(n_detected, (bool, np.bool_))
        or not isinstance(n_detected, (int, np.integer))
        or n_detected < 1
    ):
        raise ValueError("n_detected must be a positive integer")
    return log_likelihood - int(n_detected) * np.log(beta)


def marginalize_selected_draws(selected_log_likelihood, log_importance_weights):
    """Self-normalized importance expectation of selected likelihood over theta."""
    log_likelihood = np.asarray(selected_log_likelihood, dtype=float)
    log_weights = np.asarray(log_importance_weights, dtype=float)
    if (
        log_likelihood.ndim != 2
        or log_likelihood.shape[0] < 1
        or log_likelihood.shape[1] < 2
    ):
        raise ValueError(
            "selected log likelihood must have shape (draws>=1, H0 grid>=2)"
        )
    if log_weights.ndim != 1 or log_weights.size != log_likelihood.shape[0]:
        raise ValueError("log importance weights must match coefficient draws")
    if not np.all(np.isfinite(log_likelihood)) or not np.all(np.isfinite(log_weights)):
        raise ValueError("selected log likelihoods and log weights must be finite")
    return weighted_log_marginal_likelihood(log_likelihood, log_weights)


def marginalize_selected_draws_with_mcse(
    selected_log_likelihood,
    log_importance_weights,
    *,
    n_batches=10,
    min_batch_ess=10.0,
):
    """Return the full SNIS estimate and a disjoint-batch log-scale MCSE.

    Importance proposals are iid. Disjoint batches provide a Monte Carlo
    stability diagnostic; the full-sample estimate remains the reported
    marginalized likelihood.
    """
    log_likelihood = np.asarray(selected_log_likelihood, dtype=float)
    log_weights = np.asarray(log_importance_weights, dtype=float)
    if (
        log_likelihood.ndim != 2
        or log_likelihood.shape[0] < 1
        or log_likelihood.shape[1] < 2
    ):
        raise ValueError("selected log likelihood must have shape (draws>=1, H0 grid>=2)")
    if log_weights.ndim != 1 or log_weights.size != log_likelihood.shape[0]:
        raise ValueError("log importance weights must match coefficient draws")
    if not np.all(np.isfinite(log_likelihood)) or not np.all(np.isfinite(log_weights)):
        raise ValueError("selected log likelihoods and log weights must be finite")
    if (
        isinstance(n_batches, (bool, np.bool_))
        or not isinstance(n_batches, (int, np.integer))
        or n_batches < 2
    ):
        raise ValueError("n_batches must be an integer >= 2")
    if not np.isfinite(min_batch_ess) or min_batch_ess < 1:
        raise ValueError("min_batch_ess must be finite and >= 1")

    full_estimate = marginalize_selected_draws(log_likelihood, log_weights)
    if log_likelihood.shape[0] < 2 * n_batches:
        return {
            "marginalized_selected_log_likelihood": full_estimate,
            "batch_log_likelihoods": None,
            "log_likelihood_mcse": None,
            "batch_ess": None,
            "status": "INSUFFICIENT_DRAWS",
        }

    batches = np.array_split(np.arange(log_likelihood.shape[0]), n_batches)
    batch_values = []
    batch_ess = []
    for indices in batches:
        batch_log_weights = log_weights[indices]
        scaled = np.exp(batch_log_weights - np.max(batch_log_weights))
        normalized = scaled / scaled.sum()
        batch_ess.append(float(1.0 / np.sum(normalized**2)))
        batch_values.append(
            marginalize_selected_draws(log_likelihood[indices], batch_log_weights)
        )
    batch_values = np.asarray(batch_values, dtype=float)
    batch_ess = np.asarray(batch_ess, dtype=float)
    mcse = batch_values.std(axis=0, ddof=1) / np.sqrt(n_batches)
    status = "PASS" if np.min(batch_ess) >= min_batch_ess else "LOW_BATCH_ESS"
    return {
        "marginalized_selected_log_likelihood": full_estimate,
        "batch_log_likelihoods": batch_values,
        "log_likelihood_mcse": mcse,
        "batch_ess": batch_ess,
        "status": status,
    }


def propagate_inclusion_posterior_to_h0(
    posterior: InclusionPosteriorApproximation,
    analysis_redshift,
    redshift_grid,
    H0_grid,
    distances,
    *,
    merger_weight,
    scott_multiplier,
    sigma,
    detection_threshold,
    max_output_bytes=1_000_000_000,
    marginalization_batches=10,
    min_batch_ess=10.0,
):
    """Build one IPW density and selected likelihood per weighted theta proposal.

    The same coefficient proposal determines the host density, event-mixture
    likelihood, and beta_GW vector across the full H0 grid. Importance weights
    are retained and used in the final log-space marginalization.
    """
    if posterior.status != "SUCCESS":
        raise ValueError(
            f"Cannot propagate posterior approximation with status={posterior.status}"
        )
    proposals = np.asarray(posterior.coefficient_proposals, dtype=float)
    weights = np.asarray(posterior.normalized_importance_weights, dtype=float)
    log_weights = np.asarray(posterior.log_importance_weights, dtype=float)
    if (
        proposals.ndim != 2
        or weights.shape != (proposals.shape[0],)
        or log_weights.shape != weights.shape
    ):
        raise ValueError("posterior proposal and importance-weight shapes are inconsistent")
    if (
        not np.all(np.isfinite(proposals))
        or not np.all(np.isfinite(weights))
        or np.any(weights < 0)
        or not np.any(weights > 0)
    ):
        raise ValueError(
            "posterior proposals and normalized importance weights must be finite and nonnegative"
        )
    if not np.isclose(weights.sum(), 1.0, rtol=1e-10, atol=1e-12):
        raise ValueError("posterior importance weights must sum to one")
    normalized_from_logs = np.exp(log_weights - np.max(log_weights))
    normalized_from_logs /= normalized_from_logs.sum()
    if not np.allclose(normalized_from_logs, weights, rtol=1e-9, atol=1e-14):
        raise ValueError("normalized and log importance weights are inconsistent")

    catalogue = np.asarray(analysis_redshift, dtype=float)
    z_grid = np.asarray(redshift_grid, dtype=float)
    h0_grid = np.asarray(H0_grid, dtype=float)
    distances = np.asarray(distances, dtype=float)
    merger = np.asarray(merger_weight, dtype=float)
    if catalogue.ndim != 1 or catalogue.size < 2 or not np.all(np.isfinite(catalogue)):
        raise ValueError("analysis_redshift must be a finite 1D catalogue with at least two objects")
    if (
        z_grid.ndim != 1
        or z_grid.size < 3
        or not np.all(np.isfinite(z_grid))
        or np.any(np.diff(z_grid) <= 0)
    ):
        raise ValueError("redshift_grid must be finite and strictly increasing")
    if (
        h0_grid.ndim != 1
        or h0_grid.size < 2
        or not np.all(np.isfinite(h0_grid))
        or np.any(np.diff(h0_grid) <= 0)
    ):
        raise ValueError("H0_grid must be finite and strictly increasing")
    if merger.shape != z_grid.shape or not np.all(np.isfinite(merger)) or np.any(merger <= 0):
        raise ValueError("merger_weight must be finite, positive and match redshift_grid")
    if distances.ndim != 1 or distances.size < 1 or not np.all(np.isfinite(distances)):
        raise ValueError("distances must be a nonempty finite vector")
    if proposals.shape[1] != posterior.degree + 1:
        raise ValueError("proposal coefficient dimension does not match posterior degree")
    if not np.isclose(z_grid[0], posterior.z_min, rtol=0.0, atol=1e-12) or not np.isclose(
        z_grid[-1], posterior.z_max, rtol=0.0, atol=1e-12
    ):
        raise ValueError("redshift_grid must span the full fixed inclusion-prior domain")
    if (
        isinstance(max_output_bytes, (bool, np.bool_))
        or not isinstance(max_output_bytes, (int, np.integer))
        or max_output_bytes < 1
    ):
        raise ValueError("max_output_bytes must be a positive integer")
    estimated_bytes = proposals.shape[0] * (
        z_grid.size + 3 * h0_grid.size + 1
    ) * np.dtype(float).itemsize
    if proposals.shape[0] >= 2 * marginalization_batches:
        estimated_bytes += marginalization_batches * h0_grid.size * np.dtype(float).itemsize
    if estimated_bytes > max_output_bytes:
        raise MemoryError(
            f"Requested draw-level outputs need about {estimated_bytes} bytes, "
            f"above max_output_bytes={max_output_bytes}; reduce proposal draws or grid size."
        )

    design = fixed_redshift_design(
        catalogue,
        degree=posterior.degree,
        z_min=posterior.z_min,
        z_max=posterior.z_max,
    )
    density_draws = []
    for coefficients in proposals:
        inclusion_probability = expit(design @ coefficients)
        density = estimate_host_density(
            z_grid,
            catalogue,
            inclusion_probability=inclusion_probability,
            merger_weight=merger,
            scott_multiplier=scott_multiplier,
        )
        density_draws.append(density)

    density_draws = np.asarray(density_draws)
    likelihood_draws = selected_log_likelihood_by_draw(
        z_grid,
        density_draws,
        h0_grid,
        distances,
        sigma=sigma,
        detection_threshold=detection_threshold,
    )
    unconditional = likelihood_draws["unconditional_log_likelihoods"]
    efficiencies = likelihood_draws["detection_efficiencies"]
    selected = selection_correct_each_draw(
        unconditional, efficiencies, n_detected=distances.size
    )
    marginalization = marginalize_selected_draws_with_mcse(
        selected, log_weights, n_batches=marginalization_batches,
        min_batch_ess=min_batch_ess,
    )
    marginalized = marginalization["marginalized_selected_log_likelihood"]
    return FiniteValidationH0Marginalization(
        host_density_draws=density_draws,
        unconditional_log_likelihood_by_draw=unconditional,
        detection_efficiency_by_draw=efficiencies,
        selected_log_likelihood_by_draw=selected,
        normalized_importance_weights=weights.copy(),
        marginalized_selected_log_likelihood=marginalized,
        importance_ess=float(1.0 / np.sum(weights**2)),
        max_normalized_importance_weight=float(np.max(weights)),
        marginalization_batch_log_likelihoods=marginalization["batch_log_likelihoods"],
        marginalization_log_likelihood_mcse=marginalization["log_likelihood_mcse"],
        marginalization_batch_ess=marginalization["batch_ess"],
        marginalization_mcse_status=marginalization["status"],
    )
