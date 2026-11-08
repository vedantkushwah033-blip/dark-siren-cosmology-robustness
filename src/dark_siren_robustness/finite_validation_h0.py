"""Finite-validation host-density draws and selection-conditioned H0 marginalization.

This module is a synthetic low-redshift toy-inference component. Catalogue
inclusion uncertainty is represented by weighted coefficient proposals; GW
detection selection is evaluated separately within every coefficient draw.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike, NDArray
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import expit
from scipy.stats import norm

from dark_siren_robustness.bayesian_inclusion import (
    fixed_redshift_design,
    weighted_log_marginal_likelihood,
)
from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.selection_toy import C_KM_S


def _strict_grid(values: ArrayLike, name: str, *, positive: bool = False) -> NDArray[np.float64]:
    result = np.asarray(values, dtype=float)
    if result.ndim != 1 or result.size < 3 or not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must be a finite 1D grid with at least three points")
    if np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be strictly increasing")
    if positive and np.any(result <= 0):
        raise ValueError(f"{name} must be strictly positive")
    return result


def _trapezoid_weights(grid: NDArray[np.float64]) -> NDArray[np.float64]:
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def reconstruct_host_density_draws(
    redshift_grid: ArrayLike,
    analysis_redshift: ArrayLike,
    coefficient_proposals: ArrayLike,
    *,
    degree: int,
    merger_weight: ArrayLike,
    z_min: float = 0.01,
    z_max: float = 0.40,
    scott_multiplier: float = 0.15,
) -> NDArray[np.float64]:
    """Reconstruct one weighted-KDE host density for every coefficient proposal.

    No oracle inclusion probabilities are accepted. Every proposal must yield
    a finite, strictly positive inclusion probability for every analysis object;
    invalid draws raise rather than being silently dropped.
    """
    grid = _strict_grid(redshift_grid, "redshift_grid", positive=True)
    sample = np.asarray(analysis_redshift, dtype=float)
    coefficients = np.asarray(coefficient_proposals, dtype=float)
    rate = np.asarray(merger_weight, dtype=float)
    if sample.ndim != 1 or sample.size < 5 or not np.all(np.isfinite(sample)):
        raise ValueError("analysis_redshift must contain at least five finite values")
    if np.any((sample < grid[0]) | (sample > grid[-1])):
        raise ValueError("analysis redshifts must lie inside redshift_grid bounds")
    if coefficients.ndim != 2 or coefficients.shape[0] < 1:
        raise ValueError("coefficient_proposals must have shape (draw, coefficient)")
    if coefficients.shape[1] != degree + 1 or degree not in (1, 2):
        raise ValueError("coefficient dimension must match degree 1 or 2")
    if not np.all(np.isfinite(coefficients)):
        raise ValueError("coefficient proposals must be finite")
    if rate.shape != grid.shape or not np.all(np.isfinite(rate)) or np.any(rate <= 0):
        raise ValueError("merger_weight must be finite, positive and match redshift_grid")
    if not np.isfinite(scott_multiplier) or scott_multiplier <= 0:
        raise ValueError("scott_multiplier must be finite and positive")

    analysis_design = fixed_redshift_design(
        sample, degree=degree, z_min=z_min, z_max=z_max
    )
    inclusion = expit(coefficients @ analysis_design.T)
    if not np.all(np.isfinite(inclusion)) or np.any(inclusion <= 0) or np.any(inclusion > 1):
        raise ArithmeticError("At least one coefficient proposal gives invalid inclusion probabilities")

    densities = np.empty((coefficients.shape[0], grid.size), dtype=float)
    for index, probabilities in enumerate(inclusion):
        try:
            densities[index] = estimate_host_density(
                grid,
                sample,
                inclusion_probability=probabilities,
                merger_weight=rate,
                scott_multiplier=scott_multiplier,
            )
        except (ValueError, ArithmeticError, FloatingPointError) as exc:
            raise ArithmeticError(
                f"Host-density reconstruction failed for coefficient proposal {index}: {exc}"
            ) from exc
    if not np.all(np.isfinite(densities)) or np.any(densities < 0):
        raise ArithmeticError("Reconstructed host densities contain invalid values")
    integrals = trapezoid(densities, grid, axis=1)
    if not np.allclose(integrals, 1.0, rtol=1e-8, atol=1e-10):
        raise ArithmeticError("At least one reconstructed host density is not normalized")
    return densities


def selected_log_likelihood_by_draw(
    redshift_grid: ArrayLike,
    host_density_draws: ArrayLike,
    H0_grid: ArrayLike,
    distances: ArrayLike,
    *,
    sigma: float,
    detection_threshold: float,
    h0_batch_size: int = 16,
) -> dict:
    """Evaluate unconditional and selection-conditioned log likelihood per draw.

    The Gaussian distance kernel is evaluated in batches over H0. For each
    coefficient draw b, beta_GW(H0; b) is calculated from that same density
    before subtracting N_detected * log(beta_GW). No averaging of beta occurs.
    """
    z = _strict_grid(redshift_grid, "redshift_grid", positive=True)
    h0 = _strict_grid(H0_grid, "H0_grid", positive=True)
    densities = np.asarray(host_density_draws, dtype=float)
    d = np.asarray(distances, dtype=float)
    sigma = float(sigma)
    threshold = float(detection_threshold)
    if densities.ndim != 2 or densities.shape[0] < 1 or densities.shape[1] != z.size:
        raise ValueError("host_density_draws must have shape (draw, redshift_grid)")
    if not np.all(np.isfinite(densities)) or np.any(densities < 0):
        raise ValueError("host densities must be finite and nonnegative")
    if not np.allclose(trapezoid(densities, z, axis=1), 1.0, rtol=1e-7, atol=1e-9):
        raise ValueError("every host density must integrate to one")
    if d.ndim != 1 or d.size < 1 or not np.all(np.isfinite(d)) or np.any(d <= 0):
        raise ValueError("distances must be a nonempty finite positive vector")
    if not np.isfinite(sigma) or sigma <= 0 or not np.isfinite(threshold):
        raise ValueError("sigma must be positive and detection_threshold finite")
    if np.any(d < threshold):
        raise ValueError("every supplied event must satisfy the detection threshold")
    if not isinstance(h0_batch_size, int) or h0_batch_size < 1:
        raise ValueError("h0_batch_size must be a positive integer")

    quadrature = _trapezoid_weights(z)
    mass = densities * quadrature[None, :]
    if not np.allclose(mass.sum(axis=1), 1.0, rtol=1e-7, atol=1e-9):
        raise ArithmeticError("Discrete host-redshift masses do not sum to one")

    n_draws = densities.shape[0]
    unconditional = np.empty((n_draws, h0.size), dtype=float)
    beta_all = np.empty((n_draws, h0.size), dtype=float)
    selected = np.empty((n_draws, h0.size), dtype=float)
    mass_transpose = mass.T

    for start in range(0, h0.size, h0_batch_size):
        stop = min(start + h0_batch_size, h0.size)
        batch_h0 = h0[start:stop]
        means = C_KM_S * z[None, :] / batch_h0[:, None]
        log_kernel = norm.logpdf(
            d[None, :, None], loc=means[:, None, :], scale=sigma
        )
        max_log_kernel = np.max(log_kernel, axis=2)
        scaled_kernel = np.exp(log_kernel - max_log_kernel[:, :, None])
        detection_probability = norm.cdf((means - threshold) / sigma)

        # Shapes: (H0 batch, event, draw) and (H0 batch, draw).
        mixture = np.einsum("hez,bz->heb", scaled_kernel, mass, optimize=True)
        beta = np.einsum("hz,bz->hb", detection_probability, mass, optimize=True)
        if not np.all(np.isfinite(mixture)) or np.any(mixture <= 0):
            raise ArithmeticError("An event/redshift mixture likelihood is numerically invalid")
        if (
            not np.all(np.isfinite(beta))
            or np.any(beta <= 0)
            or np.any(beta > 1.0 + 1e-12)
        ):
            raise ArithmeticError("A per-draw GW detection efficiency is outside (0,1]")

        log_unconditional = np.sum(
            np.log(mixture) + max_log_kernel[:, :, None], axis=1
        )
        log_selected = log_unconditional - d.size * np.log(beta)
        unconditional[:, start:stop] = log_unconditional.T
        beta_all[:, start:stop] = beta.T
        selected[:, start:stop] = log_selected.T

    return {
        "H0_grid": h0.copy(),
        "unconditional_log_likelihoods": unconditional,
        "detection_efficiencies": beta_all,
        "selected_log_likelihoods": selected,
        "n_detected": int(d.size),
    }


def marginalize_selection_conditioned_likelihood(
    selected_log_likelihoods: ArrayLike, log_importance_weights: ArrayLike
) -> NDArray[np.float64]:
    """Integrate selected likelihood over the validation-only coefficient posterior."""
    return weighted_log_marginal_likelihood(
        selected_log_likelihoods, log_importance_weights
    )


def normalize_h0_posterior(
    H0_grid: ArrayLike,
    log_marginal_likelihood: ArrayLike,
    *,
    prior_density: ArrayLike,
) -> dict:
    """Normalize a declared-prior H0 posterior and calculate 68% summaries."""
    h0 = _strict_grid(H0_grid, "H0_grid", positive=True)
    log_like = np.asarray(log_marginal_likelihood, dtype=float)
    prior = np.asarray(prior_density, dtype=float)
    if log_like.shape != h0.shape or not np.all(np.isfinite(log_like)):
        raise ValueError("log_marginal_likelihood must be finite and match H0_grid")
    if prior.shape != h0.shape or not np.all(np.isfinite(prior)) or np.any(prior < 0):
        raise ValueError("prior_density must be finite, nonnegative and match H0_grid")
    if not np.any(prior > 0):
        raise ValueError("prior_density must be positive somewhere on H0_grid")
    log_prior = np.full_like(prior, -np.inf)
    positive = prior > 0
    log_prior[positive] = np.log(prior[positive])
    log_posterior = log_like + log_prior
    shift = float(np.max(log_posterior))
    unnormalized = np.exp(log_posterior - shift)
    area = float(trapezoid(unnormalized, h0))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("H0 posterior normalization failed")
    density = unnormalized / area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h0)]
    cdf /= cdf[-1]
    q16, median, q84 = (
        float(np.interp(q, cdf, h0)) for q in (0.16, 0.50, 0.84)
    )
    mean = float(trapezoid(h0 * density, h0))
    five_percent = 0.05 * (h0[-1] - h0[0])
    left_cdf = float(np.interp(h0[0] + five_percent, h0, cdf))
    right_cdf = float(np.interp(h0[-1] - five_percent, h0, cdf))
    return {
        "posterior_density": density,
        "posterior_integral": float(trapezoid(density, h0)),
        "posterior_mean": mean,
        "posterior_median": median,
        "ci68_low": q16,
        "ci68_high": q84,
        "ci68_width": q84 - q16,
        "boundary_mass_fraction": left_cdf + (1.0 - right_cdf),
    }
