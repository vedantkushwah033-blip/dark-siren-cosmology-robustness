"""Quantify degree-two coefficient-prior mismatch under centered/raw quadratic bases."""
from __future__ import annotations

import argparse
import json

import numpy as np


def audit_prior_mapping(
    quadratic_center: float = 0.2845005579191241,
    prior_sd_intercept: float = 5.0,
    prior_sd_slope: float = 2.5,
    prior_sd_curvature: float = 2.5,
) -> dict:
    values = np.asarray([
        quadratic_center, prior_sd_intercept, prior_sd_slope, prior_sd_curvature
    ], dtype=float)
    if not np.all(np.isfinite(values)) or min(values[1:]) <= 0:
        raise ValueError("center must be finite and prior scales finite and positive")
    # theta_raw = (theta0, theta1, theta2) and phi_centered = (alpha,beta,gamma)
    # satisfy theta0 = alpha - c*gamma, theta1=beta, theta2=gamma.
    transform = np.array([
        [1.0, 0.0, -quadratic_center],
        [0.0, 1.0, 0.0],
        [0.0, 0.0, 1.0],
    ])
    inverse = np.linalg.inv(transform)
    centered_covariance = np.diag([
        prior_sd_intercept**2, prior_sd_slope**2, prior_sd_curvature**2
    ])
    raw_covariance_induced_by_centered_prior = (
        transform @ centered_covariance @ transform.T
    )
    raw_independent_covariance = np.diag(np.diag(centered_covariance))
    centered_covariance_induced_by_raw_prior = (
        inverse @ raw_independent_covariance @ inverse.T
    )
    alpha_gamma_cov = float(centered_covariance_induced_by_raw_prior[0, 2])
    alpha_sd = float(np.sqrt(centered_covariance_induced_by_raw_prior[0, 0]))
    gamma_sd = float(np.sqrt(centered_covariance_induced_by_raw_prior[2, 2]))
    alpha_gamma_corr = alpha_gamma_cov / (alpha_sd * gamma_sd)
    target_precision = np.linalg.inv(centered_covariance)
    sign_target, logdet_target = np.linalg.slogdet(centered_covariance)
    sign_current, logdet_current = np.linalg.slogdet(centered_covariance_induced_by_raw_prior)
    if sign_target <= 0 or sign_current <= 0:
        raise ArithmeticError("prior covariance is not positive definite")
    # KL(current raw-basis-induced prior || independent centered-basis prior).
    dimension = 3
    kl = 0.5 * (
        np.trace(target_precision @ centered_covariance_induced_by_raw_prior)
        - dimension + logdet_target - logdet_current
    )
    return {
        "status": "PASS",
        "quadratic_center": float(quadratic_center),
        "coefficient_mapping": "theta0=alpha-c*gamma; theta1=beta; theta2=gamma",
        "degree1_prior_sd_alpha_beta": [prior_sd_intercept, prior_sd_slope],
        "degree2_raw_independent_prior_sd": [
            prior_sd_intercept, prior_sd_slope, prior_sd_curvature
        ],
        "centered_prior_covariance_in_raw_coefficients":
            raw_covariance_induced_by_centered_prior.tolist(),
        "raw_independent_prior_covariance_in_centered_coefficients":
            centered_covariance_induced_by_raw_prior.tolist(),
        "raw_prior_induced_alpha_sd": alpha_sd,
        "raw_prior_induced_beta_sd": float(np.sqrt(centered_covariance_induced_by_raw_prior[1, 1])),
        "raw_prior_induced_gamma_sd": gamma_sd,
        "raw_prior_induced_alpha_gamma_correlation": float(alpha_gamma_corr),
        "kl_nats_raw_induced_vs_independent_centered_prior": float(kl),
        "interpretation": (
            "Independent Gaussian priors on raw [1,x,x^2] coefficients induce a "
            "correlated prior on centered [1,x,x^2-c] coefficients. A D1/D2 "
            "comparison with these defaults compares complete procedures and "
            "does not isolate polynomial degree."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--quadratic-center", type=float, default=0.2845005579191241)
    p.add_argument("--prior-sd-intercept", type=float, default=5.0)
    p.add_argument("--prior-sd-slope", type=float, default=2.5)
    p.add_argument("--prior-sd-curvature", type=float, default=2.5)
    args = p.parse_args()
    print(json.dumps(audit_prior_mapping(
        args.quadratic_center, args.prior_sd_intercept,
        args.prior_sd_slope, args.prior_sd_curvature,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
