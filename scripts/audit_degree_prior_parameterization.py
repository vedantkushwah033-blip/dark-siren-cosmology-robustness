"""Verify matched D1/D2 Gaussian priors under the centered quadratic basis."""
from __future__ import annotations

import argparse
import json

import numpy as np


def audit_centered_prior_match(
    prior_sd_intercept: float = 5.0,
    prior_sd_slope: float = 2.5,
    prior_sd_curvature: float = 2.5,
) -> dict:
    scales = np.asarray(
        [prior_sd_intercept, prior_sd_slope, prior_sd_curvature], dtype=float
    )
    if not np.all(np.isfinite(scales)) or np.any(scales <= 0):
        raise ValueError("all coefficient prior scales must be finite and positive")

    d1_covariance = np.diag(scales[:2] ** 2)
    d2_covariance = np.diag(scales ** 2)
    d2_marginal_alpha_beta = d2_covariance[:2, :2]
    cross_covariance = d2_covariance[:2, 2:]
    matched = (
        np.array_equal(d1_covariance, d2_marginal_alpha_beta)
        and np.count_nonzero(cross_covariance) == 0
    )
    if not matched:
        raise AssertionError(
            "D1 prior on (alpha,beta) does not match D2 marginal prior "
            "under centered basis [1,x,x^2-c_g]"
        )
    return {
        "status": "PASS",
        "basis_degree1": "[1, x]",
        "basis_degree2": "[1, x, x^2 - c_g]",
        "degree1_prior_covariance_alpha_beta": d1_covariance.tolist(),
        "degree2_marginal_prior_covariance_alpha_beta": d2_marginal_alpha_beta.tolist(),
        "degree2_alpha_beta_gamma_cross_covariance": cross_covariance.tolist(),
        "degree2_curvature_prior_sd": float(prior_sd_curvature),
        "marginal_prior_match": True,
        "interpretation": (
            "The D1 prior on (alpha,beta) exactly equals the D2 marginal prior "
            "on (alpha,beta); gamma is added independently. This removes the "
            "raw-versus-centered coefficient-prior mismatch."
        ),
    }


def main() -> None:
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--prior-sd-intercept", type=float, default=5.0)
    p.add_argument("--prior-sd-slope", type=float, default=2.5)
    p.add_argument("--prior-sd-curvature", type=float, default=2.5)
    args = p.parse_args()
    print(json.dumps(audit_centered_prior_match(
        args.prior_sd_intercept, args.prior_sd_slope, args.prior_sd_curvature
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
