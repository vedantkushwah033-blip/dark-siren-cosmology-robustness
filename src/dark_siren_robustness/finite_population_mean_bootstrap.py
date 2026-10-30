"""Studentized pseudo-population bootstrap for a finite-population mean.

This is a linear benchmark for checking finite-population resampling mechanics.
It is not an uncertainty procedure for the nonlinear dark-siren posterior
functional; that requires separate empirical validation.
"""
from __future__ import annotations

import numpy as np
from numpy.typing import ArrayLike


def studentized_pseudopopulation_srswor_bootstrap(
    sample: ArrayLike,
    population_size: int,
    *,
    n_bootstrap: int = 999,
    seed: int = 20261031,
    confidence: float = 0.95,
) -> dict[str, object]:
    """Construct a studentized pseudo-population SRSWOR interval for a mean.

    For each replicate, generate a size-N multinomial pseudo-population from
    the observed n records, then draw n units without replacement. The
    bootstrap t statistic is centered on that replicate's pseudo-population
    mean and studentized with its SRSWOR variance estimate. The observed
    interval inverts the bootstrap t quantiles around the observed sample mean.

    This follows the single-stage SRSWOR construction of Wang, Peng & Kim
    (2022). Their asymptotic result concerns a finite-population mean under
    stated regularity conditions; it does not establish validity for nonlinear
    fitted-model or posterior functionals.
    """
    x = np.asarray(sample, dtype=float)
    if x.ndim != 1 or x.size < 2 or not np.all(np.isfinite(x)):
        raise ValueError("sample must be a finite one-dimensional array of length >= 2")
    for name, value, minimum in (
        ("population_size", population_size, int(x.size)),
        ("n_bootstrap", n_bootstrap, 2),
        ("seed", seed, 0),
    ):
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
        if int(value) < minimum:
            raise ValueError(f"{name} must be >= {minimum}")
    N, n = int(population_size), int(x.size)
    if N <= n:
        raise ValueError("population_size must exceed sample size for studentized SRSWOR")
    if not np.isfinite(confidence) or not 0.0 < confidence < 1.0:
        raise ValueError("confidence must lie strictly between zero and one")

    sample_variance = float(np.var(x, ddof=1))
    fpc_variance = 1.0 - n / N
    observed_se = float(np.sqrt(fpc_variance * sample_variance / n))
    if not np.isfinite(observed_se) or observed_se <= 0.0:
        raise ValueError("observed sample must have positive finite variance")

    rng = np.random.default_rng(int(seed))
    t_statistics = np.full(int(n_bootstrap), np.nan, dtype=float)
    pseudo_means = np.full(int(n_bootstrap), np.nan, dtype=float)
    replicate_means = np.full(int(n_bootstrap), np.nan, dtype=float)
    replicate_ses = np.full(int(n_bootstrap), np.nan, dtype=float)
    success = np.zeros(int(n_bootstrap), dtype=bool)

    probabilities = np.full(n, 1.0 / n, dtype=float)
    for b in range(int(n_bootstrap)):
        counts = rng.multinomial(N, probabilities)
        pseudo_population = np.repeat(x, counts)
        positions = rng.choice(N, size=n, replace=False)
        bootstrap_sample = pseudo_population[positions]
        pseudo_mean = float(np.mean(pseudo_population))
        bootstrap_mean = float(np.mean(bootstrap_sample))
        bootstrap_variance = float(np.var(bootstrap_sample, ddof=1))
        bootstrap_se = float(np.sqrt(fpc_variance * bootstrap_variance / n))
        if not np.isfinite(bootstrap_se) or bootstrap_se <= 0.0:
            continue
        t_statistics[b] = (bootstrap_mean - pseudo_mean) / bootstrap_se
        pseudo_means[b] = pseudo_mean
        replicate_means[b] = bootstrap_mean
        replicate_ses[b] = bootstrap_se
        success[b] = True

    valid_t = t_statistics[success]
    alpha = 1.0 - float(confidence)
    if valid_t.size >= 2:
        q_low, q_high = np.quantile(valid_t, [alpha / 2.0, 1.0 - alpha / 2.0])
        estimate = float(np.mean(x))
        interval = [
            float(estimate - q_high * observed_se),
            float(estimate - q_low * observed_se),
        ]
        status = "ok"
    else:
        q_low = q_high = float("nan")
        interval = [float("nan"), float("nan")]
        status = "insufficient_valid_bootstrap_replicates"

    return {
        "schema": "studentized_finite_population_mean_srswor_v1",
        "status": status,
        "seed": int(seed),
        "population_size": N,
        "sample_size": n,
        "n_attempted": int(n_bootstrap),
        "n_successful": int(np.sum(success)),
        "n_failed": int(n_bootstrap - np.sum(success)),
        "confidence": float(confidence),
        "estimate": float(np.mean(x)),
        "sample_variance": sample_variance,
        "finite_population_correction_variance": fpc_variance,
        "estimated_design_se": observed_se,
        "bootstrap_t_quantiles": [float(q_low), float(q_high)],
        "interval": interval,
        "bootstrap_t_statistics": t_statistics,
        "bootstrap_pseudopopulation_means": pseudo_means,
        "bootstrap_sample_means": replicate_means,
        "bootstrap_design_ses": replicate_ses,
        "success_mask": success,
    }
