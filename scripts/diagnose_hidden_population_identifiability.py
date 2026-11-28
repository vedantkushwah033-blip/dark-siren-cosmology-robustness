"""Monte Carlo observed-Hessian approximation to local hidden-population information.

This diagnostic is local and model-conditional. It is not a proof of global
identifiability and does not replace finite-sample posterior recovery.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np

from dark_siren_robustness.hidden_population import log_likelihood
from scripts.run_continuous_hidden_population_recovery import (
    SIGMA_U, TRUTH, generate_dataset, LISTED_Z, LISTED_W, Z_BOUNDS,
)

PARAMETER_NAMES = ("H0", "fraction", "mu_h", "sigma_h")
STEPS = np.array([1e-2, 1e-4, 1e-5, 1e-5])
SCALES = np.array([10.0, 0.1, 0.01, 0.01])


def _loglik(u, theta, quadrature_order):
    return log_likelihood(
        u, theta[0], theta[1], theta[2], theta[3], SIGMA_U,
        listed_redshifts=LISTED_Z, listed_weights=LISTED_W,
        z_bounds=Z_BOUNDS, quadrature_order=quadrature_order,
    ) / len(u)


def numerical_hessian(u, theta=None, steps=STEPS, quadrature_order=48):
    """Central finite-difference Hessian of mean log likelihood."""
    x = np.asarray(
        [TRUTH[k] for k in ("H0", "fraction", "mu_h", "sigma_h")]
        if theta is None else theta, dtype=float
    )
    h = np.asarray(steps, dtype=float)
    if x.shape != (4,) or h.shape != (4,) or np.any(h <= 0):
        raise ValueError("theta and steps must be positive four-vectors")
    base = _loglik(u, x, quadrature_order)
    out = np.zeros((4, 4), dtype=float)
    for i in range(4):
        plus, minus = x.copy(), x.copy()
        plus[i] += h[i]
        minus[i] -= h[i]
        out[i, i] = (_loglik(u, plus, quadrature_order)-2*base+_loglik(u, minus, quadrature_order))/(h[i]**2)
        for j in range(i):
            pp, pm, mp, mm = (x.copy() for _ in range(4))
            pp[i] += h[i]; pp[j] += h[j]
            pm[i] += h[i]; pm[j] -= h[j]
            mp[i] -= h[i]; mp[j] += h[j]
            mm[i] -= h[i]; mm[j] -= h[j]
            value = (
                _loglik(u, pp, quadrature_order)-_loglik(u, pm, quadrature_order)
                -_loglik(u, mp, quadrature_order)+_loglik(u, mm, quadrature_order)
            )/(4*h[i]*h[j])
            out[i, j] = out[j, i] = value
    return out


def run(replicates=2, n_events=500, seed=20261128, quadrature_order=48):
    if replicates < 1 or n_events < 10:
        raise ValueError("replicates must be positive and n_events at least 10")
    rng = np.random.default_rng(seed)
    hessians = []
    for _ in range(replicates):
        u, _, _ = generate_dataset(n_events, rng)
        hessians.append(numerical_hessian(u, quadrature_order=quadrature_order))
    observed_information = -np.mean(hessians, axis=0)
    observed_information = (observed_information + observed_information.T)/2
    scale = np.diag(SCALES)
    standardized = scale @ observed_information @ scale
    eigenvalues = np.linalg.eigvalsh(standardized)
    positive = eigenvalues[eigenvalues > 0]
    condition = float(positive.max()/positive.min()) if len(positive)==4 else float("inf")

    result = {
        "status": "FINITE" if len(positive)==4 else "NON_POSITIVE_CURVATURE",
        "replicates": replicates, "events_per_replicate": n_events,
        "seed": seed, "quadrature_order": quadrature_order,
        "parameters": list(PARAMETER_NAMES),
        "finite_difference_steps": STEPS.tolist(),
        "standardization_scales": SCALES.tolist(),
        "mean_negative_hessian_per_event": observed_information.tolist(),
        "standardized_information": standardized.tolist(),
        "standardized_eigenvalues": eigenvalues.tolist(),
        "standardized_condition_number": condition,
    }
    if len(positive)==4:
        covariance_standardized = np.linalg.inv(standardized)
        covariance = scale @ covariance_standardized @ scale
        sd = np.sqrt(np.diag(covariance))
        correlation = covariance / np.outer(sd, sd)
        h0_known_nuisance_variance = 1.0/observed_information[0, 0]
        result.update({
            "local_asymptotic_sd_per_sqrt_N": sd.tolist(),
            "local_parameter_correlation": correlation.tolist(),
            "H0_variance_inflation_all_nuisance_free_vs_fixed": float(
                covariance[0, 0]/h0_known_nuisance_variance
            ),
        })
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--n-events", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261128)
    parser.add_argument("--quadrature-order", type=int, default=48)
    parser.add_argument("--output", default="results/continuous-hidden-population-recovery/local_information.json")
    args = parser.parse_args()
    result = run(args.replicates, args.n_events, args.seed, args.quadrature_order)
    path = Path(args.output)
    path.parent.mkdir(parents=True, exist_ok=True)
    path.write_text(json.dumps(result, indent=2)+"\n", encoding="utf-8")
    print(json.dumps({k:v for k,v in result.items() if k not in {
        "mean_negative_hessian_per_event", "standardized_information",
        "local_parameter_correlation"
    }}, indent=2))


if __name__ == "__main__":
    main()
