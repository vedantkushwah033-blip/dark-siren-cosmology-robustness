"""Selection-conditioned dataset inference with H0-dependent distance noise.

A transparent synthetic discrete-host model, not a gravitational-wave detector
likelihood. Event likelihoods and detection efficiency share one noise law.
"""
from __future__ import annotations

from collections.abc import Mapping, Sequence
import numpy as np
from scipy.integrate import trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _vector(value, name, *, min_size=1, positive=False, increasing=False):
    result = np.asarray(value, dtype=float)
    if result.ndim != 1 or result.size < min_size:
        raise ValueError(f"{name} must be a one-dimensional array with at least {min_size} value(s).")
    if not np.all(np.isfinite(result)):
        raise ValueError(f"{name} must contain only finite values.")
    if positive and np.any(result <= 0):
        raise ValueError(f"{name} must be strictly positive.")
    if increasing and np.any(np.diff(result) <= 0):
        raise ValueError(f"{name} must be strictly increasing.")
    return result


def _normalized_weights(value, size, name):
    if value is None:
        raw = np.ones(size, dtype=float)
    else:
        raw = _vector(value, name, min_size=1)
        if raw.size != size:
            raise ValueError(f"{name} must match the number of redshifts.")
        if np.any(raw < 0):
            raise ValueError(f"{name} must be non-negative.")
    maximum = float(np.max(raw))
    if not np.isfinite(maximum) or maximum <= 0:
        raise ValueError(f"{name} must have positive finite total mass.")
    scaled = raw / maximum
    total = float(np.sum(scaled))
    if not np.isfinite(total) or total <= 0:
        raise ValueError(f"{name} must have positive finite total mass.")
    return scaled / total


def _noise_parameters(fractional_sigma, sigma_floor_mpc):
    fractional_sigma = float(fractional_sigma)
    sigma_floor_mpc = float(sigma_floor_mpc)
    if not np.isfinite(fractional_sigma) or fractional_sigma <= 0:
        raise ValueError("fractional_sigma must be finite and positive.")
    if not np.isfinite(sigma_floor_mpc) or sigma_floor_mpc < 0:
        raise ValueError("sigma_floor_mpc must be finite and non-negative.")
    return fractional_sigma, sigma_floor_mpc


def _mean_and_sigma(H0_grid, redshifts, fractional_sigma, sigma_floor_mpc):
    means = C_KM_S * redshifts[None, :] / H0_grid[:, None]
    sigmas = np.hypot(sigma_floor_mpc, fractional_sigma * means)
    if not np.all(np.isfinite(means)) or not np.all(np.isfinite(sigmas)) or np.any(sigmas <= 0):
        raise ValueError("computed means and distance-noise scales must be finite, with positive scales.")
    return means, sigmas


def event_log_likelihood_h0_dependent_noise(
    H0_grid, *, redshifts, gw_distance, fractional_sigma, sigma_floor_mpc=0.0,
    weights=None, positive_truncated=False,
):
    """Return one event's host-mixture log likelihood on an H0 grid.

    For a host at redshift z, mu=c*z/H0 and
    sigma=sqrt(sigma_floor_mpc**2 + (fractional_sigma*mu)**2). If
    positive_truncated=True, condition each host's distance density on d>0.
    Host weights are discrete probability masses, not a density per redshift.
    """
    h0 = _vector(H0_grid, "H0_grid", min_size=2, positive=True, increasing=True)
    z = _vector(redshifts, "redshifts", min_size=1, positive=True)
    distance = float(gw_distance)
    if not np.isfinite(distance) or distance <= 0:
        raise ValueError("gw_distance must be finite and positive.")
    frac, floor = _noise_parameters(fractional_sigma, sigma_floor_mpc)
    host_weights = _normalized_weights(weights, z.size, "weights")
    means, sigmas = _mean_and_sigma(h0, z, frac, floor)
    terms = norm.logpdf(distance, loc=means, scale=sigmas)
    if positive_truncated:
        terms = terms - norm.logcdf(means / sigmas)
    log_weights = np.full(host_weights.shape, -np.inf)
    positive = host_weights > 0
    log_weights[positive] = np.log(host_weights[positive])
    result = logsumexp(terms + log_weights[None, :], axis=1)
    if not np.all(np.isfinite(result)):
        raise ArithmeticError("event log likelihood is not finite on the requested H0 grid.")
    return result


def population_detection_efficiency_h0_dependent_noise(
    H0_grid, *, population_redshifts, fractional_sigma, detection_distance_min,
    sigma_floor_mpc=0.0, population_weights=None, positive_truncated=False,
):
    """Compute P(detected|H0) under the same heteroscedastic discrete-host model.

    Detection is d>=detection_distance_min. Population weights are discrete
    masses; for a redshift grid approximating a continuous density, provide
    density times quadrature weights. With zero-truncated noise, divide the
    threshold survival probability by P(d>0).
    """
    h0 = _vector(H0_grid, "H0_grid", min_size=2, positive=True, increasing=True)
    z = _vector(population_redshifts, "population_redshifts", min_size=1, positive=True)
    frac, floor = _noise_parameters(fractional_sigma, sigma_floor_mpc)
    threshold = float(detection_distance_min)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_distance_min must be finite and positive.")
    weights = _normalized_weights(population_weights, z.size, "population_weights")
    means, sigmas = _mean_and_sigma(h0, z, frac, floor)
    survival = norm.sf((threshold - means) / sigmas)
    if positive_truncated:
        denominator = norm.sf(-means / sigmas)
        with np.errstate(divide="ignore", invalid="ignore"):
            survival = np.divide(survival, denominator, out=np.zeros_like(survival), where=denominator > 0)
    efficiency = survival @ weights
    if not np.all(np.isfinite(efficiency)) or np.any(efficiency <= 0) or np.any(efficiency > 1.0 + 1e-12):
        raise ArithmeticError("detection efficiency must be finite and in (0, 1]; inspect model inputs.")
    return efficiency


def selection_conditioned_dataset_log_likelihood(
    H0_grid, events: Sequence[Mapping], *, population_redshifts, detection_distance_min,
    fractional_sigma, sigma_floor_mpc=0.0, population_weights=None, positive_truncated=False,
):
    """Condition a fixed-size sample on detection with beta(H0) from the same noise law.

    Returns logL_selected = sum_i logL_i - N_detected*log(beta(H0)).
    The supplied events must pass the threshold. This does not infer event rates
    or counts and assumes the supplied population mixture is appropriate for
    selection normalization.
    """
    h0 = _vector(H0_grid, "H0_grid", min_size=2, positive=True, increasing=True)
    if not isinstance(events, Sequence) or isinstance(events, (str, bytes)) or len(events) < 1:
        raise ValueError("events must be a nonempty sequence of event mappings.")
    threshold = float(detection_distance_min)
    if not np.isfinite(threshold) or threshold <= 0:
        raise ValueError("detection_distance_min must be finite and positive.")
    frac, floor = _noise_parameters(fractional_sigma, sigma_floor_mpc)
    joint = np.zeros(h0.size, dtype=float)
    for index, event in enumerate(events):
        if not isinstance(event, Mapping):
            raise ValueError(f"events[{index}] must be a mapping.")
        distance = float(event.get("gw_distance", np.nan))
        if not np.isfinite(distance) or distance < threshold:
            raise ValueError(f"events[{index}] must have finite gw_distance >= detection_distance_min.")
        joint += event_log_likelihood_h0_dependent_noise(
            h0, redshifts=event.get("redshifts"), gw_distance=distance,
            fractional_sigma=frac, sigma_floor_mpc=floor, weights=event.get("weights"),
            positive_truncated=positive_truncated,
        )
    efficiency = population_detection_efficiency_h0_dependent_noise(
        h0, population_redshifts=population_redshifts, fractional_sigma=frac,
        sigma_floor_mpc=floor, detection_distance_min=threshold,
        population_weights=population_weights, positive_truncated=positive_truncated,
    )
    selected = joint - len(events) * np.log(efficiency)
    if not np.all(np.isfinite(selected)):
        raise ArithmeticError("selection-conditioned dataset log likelihood is not finite.")
    return {"H0_grid": h0.copy(), "unconditional_log_likelihood": joint,
            "detection_efficiency": efficiency, "selected_log_likelihood": selected,
            "n_detected": int(len(events))}


def dataset_posterior_h0_selection_conditioned(
    H0_grid, events: Sequence[Mapping], *, prior_min=50.0, prior_max=90.0,
    population_redshifts, detection_distance_min, fractional_sigma,
    sigma_floor_mpc=0.0, population_weights=None, positive_truncated=False,
):
    """Return an H0-grid posterior and endpoint-exact normalization arrays.

    The prior is uniform on [prior_min, prior_max]. normalization_grid inserts
    exact prior endpoints; integrate normalization_density on that grid to get
    exactly the declared prior-support normalization even when endpoints fall
    between nodes in H0_grid. posterior_density is evaluated on H0_grid and is
    zero outside the prior support.
    """
    h0 = _vector(H0_grid, "H0_grid", min_size=2, positive=True, increasing=True)
    lo, hi = float(prior_min), float(prior_max)
    if not np.isfinite(lo) or not np.isfinite(hi) or lo >= hi:
        raise ValueError("prior_min and prior_max must be finite with prior_min < prior_max.")
    if h0[0] > lo or h0[-1] < hi:
        raise ValueError("H0_grid must span the full prior interval.")
    selected = selection_conditioned_dataset_log_likelihood(
        h0, events, population_redshifts=population_redshifts,
        detection_distance_min=detection_distance_min, fractional_sigma=fractional_sigma,
        sigma_floor_mpc=sigma_floor_mpc, population_weights=population_weights,
        positive_truncated=positive_truncated,
    )["selected_log_likelihood"]
    in_prior = (h0 >= lo) & (h0 <= hi)
    if np.count_nonzero(in_prior) < 2:
        raise ValueError("H0_grid must contain at least two points within prior support.")
    norm_grid = np.unique(np.concatenate(([lo], h0[in_prior], [hi])))
    norm_loglike = np.interp(norm_grid, h0, selected)
    shifted = norm_loglike - np.max(norm_loglike)
    density = np.exp(shifted)
    area = float(trapezoid(density, norm_grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("selection-conditioned H0 posterior normalization failed.")
    norm_density = density / area
    posterior = np.zeros_like(h0)
    posterior[in_prior] = np.exp(selected[in_prior] - np.max(norm_loglike)) / area
    return {"H0_grid": h0.copy(), "posterior_density": posterior,
            "normalization_grid": norm_grid, "normalization_density": norm_density,
            "selected_log_likelihood": selected, "prior_min": lo, "prior_max": hi,
            "posterior_integral": float(trapezoid(norm_density, norm_grid))}
