"""Paired synthetic sensitivity study of distance-uncertainty likelihood assumptions."""
from __future__ import annotations

import math
from collections import defaultdict

import numpy as np

from dark_siren_robustness.uncertainty_likelihoods import (
    C_KM_S,
    log_likelihood_fixed_reference_sigma,
    log_likelihood_fixed_sigma,
    log_likelihood_fractional_sigma,
)

TRUTHS = (50., 55., 60., 65., 70., 75., 80., 85., 90.)
FRACTIONS = (.05, .10, .15, .20)
ARMS = ("A_fixed_realized_sigma", "B_fractional_sigma_H0", "C_fixed_reference_sigma")


def posterior_summary(grid, log_likelihood):
    """Summarize a normalized grid posterior with trapezoidal integration."""
    grid = np.asarray(grid, dtype=float)
    log_likelihood = np.asarray(log_likelihood, dtype=float)
    if grid.ndim != 1 or grid.size < 3 or log_likelihood.shape != grid.shape:
        raise ValueError("grid and log likelihood must be matching 1D arrays of length >= 3")
    if not np.all(np.isfinite(grid)) or not np.all(np.isfinite(log_likelihood)):
        raise ValueError("grid and log likelihood must be finite")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("grid must be strictly increasing")
    density = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(np.trapezoid(density, grid))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("invalid posterior normalization")
    density /= area
    cdf = np.r_[0., np.cumsum(.5 * (density[1:] + density[:-1]) * np.diff(grid))]
    cdf /= cdf[-1]
    mean = float(np.trapezoid(grid * density, grid))
    variance = float(np.trapezoid((grid - mean) ** 2 * density, grid))
    if not math.isfinite(variance) or variance <= 0:
        raise ValueError("posterior variance must be positive and finite")
    q = lambda p: float(np.interp(p, cdf, grid))
    return {"median": q(.5), "mean": mean, "posterior_sd": math.sqrt(variance),
            "lower_68": q(.16), "upper_68": q(.84),
            "lower_95": q(.025), "upper_95": q(.975)}


def run_experiment(replicates=300, seed=20261026, grid_points=1601,
                   truths=TRUTHS, fractional_errors=FRACTIONS):
    """Return raw rows and cell summaries; Arm A/B share observations exactly."""
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates < 2:
        raise ValueError("replicates must be an integer >= 2")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if isinstance(grid_points, bool) or not isinstance(grid_points, int) or grid_points < 101:
        raise ValueError("grid_points must be an integer >= 101")
    if not truths or any(not math.isfinite(t) or not 50 <= t <= 90 for t in truths):
        raise ValueError("truths must be within [50, 90]")
    if not fractional_errors or any(not math.isfinite(f) or f <= 0 for f in fractional_errors):
        raise ValueError("fractional errors must be positive and finite")

    rng, grid, raw = np.random.default_rng(seed), np.linspace(50., 90., grid_points), []
    for truth in truths:
        for frac in fractional_errors:
            for rep in range(replicates):
                z, eps = float(rng.uniform(.005, .08)), float(rng.normal())
                distance = C_KM_S * z / truth
                sigma_a = frac * distance
                sigma_c = frac * C_KM_S * z / 70.
                obs_ab, obs_c = distance + sigma_a * eps, distance + sigma_c * eps
                for arm in ARMS:
                    if arm == ARMS[0]:
                        obs, sigma = obs_ab, sigma_a
                        ll = log_likelihood_fixed_sigma(obs, grid, z, sigma)
                    elif arm == ARMS[1]:
                        obs, sigma = obs_ab, sigma_a
                        ll = log_likelihood_fractional_sigma(obs, grid, z, frac)
                    else:
                        obs, sigma = obs_c, sigma_c
                        ll = log_likelihood_fixed_reference_sigma(obs, grid, z, frac, 70.)
                    raw.append({"true_h0": truth, "fractional_error": frac, "replicate": rep,
                                "arm": arm, "host_redshift": z, "standard_normal": eps,
                                "true_distance_mpc": distance, "observed_distance_mpc": float(obs),
                                "generated_sigma_mpc": float(sigma), **posterior_summary(grid, ll)})

    groups = defaultdict(list)
    for row in raw:
        groups[(row["true_h0"], row["fractional_error"], row["arm"])].append(row)
    summaries = []
    for (truth, frac, arm), rows in sorted(groups.items()):
        med_err = np.array([r["median"] - truth for r in rows])
        mean_err = np.array([r["mean"] - truth for r in rows])
        n = len(rows)
        c68 = sum(r["lower_68"] <= truth <= r["upper_68"] for r in rows)
        c95 = sum(r["lower_95"] <= truth <= r["upper_95"] for r in rows)
        summaries.append({"true_h0": truth, "fractional_error": frac, "arm": arm, "n": n,
                          "coverage_68": c68/n, "coverage_68_successes": c68,
                          "coverage_95": c95/n, "coverage_95_successes": c95,
                          "mean_median_error": float(med_err.mean()),
                          "rmse_median": float(np.sqrt(np.mean(med_err**2))),
                          "mean_posterior_mean_error": float(mean_err.mean()),
                          "rmse_mean": float(np.sqrt(np.mean(mean_err**2))),
                          "mean_width_68": float(np.mean([r["upper_68"]-r["lower_68"] for r in rows])),
                          "mean_width_95": float(np.mean([r["upper_95"]-r["lower_95"] for r in rows])),
                          "truth_on_prior_boundary": truth in (50., 90.)})
    return raw, summaries
