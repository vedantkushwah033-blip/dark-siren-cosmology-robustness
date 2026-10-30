"""Independent replay audit for posterior-level validation bootstrap artifacts.

This module intentionally reconstructs the logistic fit, reflected weighted KDE,
selection-aware likelihood, and posterior summaries without calling the producer
pipeline's estimator or likelihood functions.
"""
from __future__ import annotations

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import minimize
from scipy.special import expit, log_expit, logsumexp
from scipy.stats import norm


def _independent_posterior(z, h0, observed, validation_z, validation_y, distances,
                           sigma, threshold, rate, bandwidth_factor):
    midpoint = 0.5 * (z[0] + z[-1])
    half = 0.5 * (z[-1] - z[0])
    x = (validation_z - midpoint) / half
    y = validation_y.astype(float)

    def objective(theta):
        eta = theta[0] + theta[1] * x
        return float(-np.sum(y * log_expit(eta) + (1.0 - y) * log_expit(-eta)))

    fit = minimize(
        objective, np.zeros(2), method="L-BFGS-B",
        bounds=((-25.0, 25.0), (-50.0, 50.0)),
        options={"maxiter": 2000, "ftol": 1e-12, "gtol": 1e-8},
    )
    if not fit.success or not np.all(np.isfinite(fit.x)):
        raise RuntimeError("independent logistic refit failed")
    intercept, slope = map(float, fit.x)
    probability = expit(intercept + slope * ((observed - midpoint) / half))
    inverse = 1.0 / probability
    scaled = inverse / np.max(inverse)
    ess = float(np.sum(scaled) ** 2 / np.sum(scaled ** 2))
    mean = float(np.average(observed, weights=scaled))
    variance = float(np.average((observed - mean) ** 2, weights=scaled))
    if not np.isfinite(variance) or variance <= 0:
        raise ArithmeticError("independent weighted variance invalid")
    factor = ess ** (-0.2) if bandwidth_factor is None else float(bandwidth_factor)
    bandwidth = factor * np.sqrt(variance)
    reflected = np.concatenate((observed, 2.0 * z[0] - observed, 2.0 * z[-1] - observed))
    kernel_weights = np.tile((scaled / np.sum(scaled)) / 3.0, 3)
    kernels = np.exp(-0.5 * ((z[:, None] - reflected[None, :]) / bandwidth) ** 2)
    density = (kernels @ kernel_weights) / (np.sqrt(2.0 * np.pi) * bandwidth)
    density *= rate
    density /= trapezoid(density, z)
    q = np.empty_like(z)
    delta = np.diff(z)
    q[0], q[-1] = delta[0] / 2.0, delta[-1] / 2.0
    q[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    valid = (density > 0) & (q > 0)
    log_mass = np.log(density[valid]) + np.log(q[valid])
    means = 299792.458 * z[valid][None, :] / h0[:, None]
    terms = norm.logpdf(
        distances[None, :, None], loc=means[:, None, :], scale=sigma
    ) + log_mass[None, None, :]
    joint = np.sum(logsumexp(terms, axis=2), axis=1)
    detect = norm.cdf((299792.458 * z[None, :] / h0[:, None] - threshold) / sigma)
    beta = np.sum(density[None, :] * q[None, :] * detect, axis=1)
    selected = joint - distances.size * np.log(beta)
    posterior = np.exp(selected - np.max(selected))
    posterior /= trapezoid(posterior, h0)
    cdf = np.r_[0.0, cumulative_trapezoid(posterior, h0)]
    cdf /= cdf[-1]
    summary = {
        "posterior_mean": float(trapezoid(h0 * posterior, h0)),
        "posterior_median": float(np.interp(0.50, cdf, h0)),
        "ci68_low": float(np.interp(0.16, cdf, h0)),
        "ci68_high": float(np.interp(0.84, cdf, h0)),
    }
    diagnostics = {
        "ess": ess,
        "max_inverse": float(np.max(inverse)),
        "boundary": bool(
            np.isclose(abs(intercept), 25.0, atol=1e-5, rtol=0.0)
            or np.isclose(abs(slope), 50.0, atol=1e-5, rtol=0.0)
        ),
    }
    return posterior, summary, diagnostics


def _assert_close(actual, expected, label, *, atol=5e-7, rtol=5e-6):
    if not np.isclose(float(actual), float(expected), atol=atol, rtol=rtol):
        raise AssertionError(f"{label} mismatch: stored={actual}, rebuilt={expected}")


def audit_posterior_validation_bootstrap(
    result: dict[str, object],
    redshift_grid,
    h0_grid,
    observed_redshift,
    validation_redshift,
    validation_included,
    observed_distances_mpc,
    *,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    parent_population_size: int,
    merger_weight=None,
    bandwidth_factor=None,
) -> dict[str, object]:
    """Replay RNG indices and independently reconstruct every posterior attempt."""
    z = np.asarray(redshift_grid, dtype=float)
    h0 = np.asarray(h0_grid, dtype=float)
    observed = np.asarray(observed_redshift, dtype=float)
    vz = np.asarray(validation_redshift, dtype=float)
    vy = np.asarray(validation_included, dtype=int)
    distances = np.asarray(observed_distances_mpc, dtype=float)
    rate = np.ones_like(z) if merger_weight is None else np.asarray(merger_weight, dtype=float)
    method = result.get("method")
    if result.get("schema") != "validation_h0_posterior_resampling_v1":
        raise AssertionError("unsupported posterior-bootstrap artifact schema")
    if method not in {"paired_bootstrap", "design_aware_srswor"}:
        raise AssertionError("unknown resampling method")
    attempted = int(result["n_attempted"])
    success = np.asarray(result["success_mask"], dtype=bool)
    if success.shape != (attempted,):
        raise AssertionError("success mask has incorrect shape")
    if int(np.sum(success)) != int(result["n_successful"]) or int(np.sum(~success)) != int(result["n_failed"]):
        raise AssertionError("attempt/success/failure accounting does not reconcile")
    means = np.asarray(result["attempt_posterior_mean"], dtype=float)
    references = np.asarray(result["reference_posterior_mean_by_attempt"], dtype=float)
    if means.shape != (attempted,) or references.shape != (attempted,):
        raise AssertionError("attempt-aligned posterior arrays have incorrect shape")
    samples = result["sample_indices"]
    pseudo = result["pseudo_population_indices"]
    positions = result["selected_pseudo_population_positions"]
    if len(samples) != attempted:
        raise AssertionError("sample index list must contain every attempt")
    if method == "design_aware_srswor" and not (
        len(pseudo) == attempted and len(positions) == attempted
    ):
        raise AssertionError("design-aware arrays must contain every attempt")

    rng = np.random.default_rng(int(result["seed"]))
    if method == "paired_bootstrap":
        _, full_ref, _ = _independent_posterior(
            z, h0, observed, vz, vy, distances, sigma_mpc,
            detection_threshold_mpc, rate, bandwidth_factor,
        )
        fixed_reference = full_ref["posterior_mean"]
        _assert_close(result["fixed_full_validation_reference_mean"], fixed_reference,
                      "full-validation reference")
    else:
        fixed_reference = float("nan")

    posterior_draws = np.asarray(result["posterior_density_draws"], dtype=float)
    if posterior_draws.shape != (int(np.sum(success)), h0.size):
        raise AssertionError("posterior-density draw matrix has incorrect shape")
    success_cursor = 0
    failure_indices = []
    for b in range(attempted):
        if method == "paired_bootstrap":
            expected = rng.integers(0, vz.size, size=vz.size, dtype=np.int64)
            if not np.array_equal(np.asarray(samples[b]), expected):
                raise AssertionError(f"paired-bootstrap indices do not replay at attempt {b}")
            sample_z, sample_y = vz[expected], vy[expected]
            reference = fixed_reference
        else:
            counts = rng.multinomial(parent_population_size, np.full(vz.size, 1.0 / vz.size))
            expected_pseudo = np.repeat(np.arange(vz.size, dtype=np.int64), counts)
            expected_positions = rng.choice(parent_population_size, size=vz.size, replace=False)
            expected_sample = expected_pseudo[expected_positions]
            if not np.array_equal(np.asarray(pseudo[b]), expected_pseudo):
                raise AssertionError(f"pseudo-population indices do not replay at attempt {b}")
            if not np.array_equal(np.asarray(positions[b]), expected_positions):
                raise AssertionError(f"selected positions do not replay at attempt {b}")
            if not np.array_equal(np.asarray(samples[b]), expected_sample):
                raise AssertionError(f"mapped SRSWOR indices do not replay at attempt {b}")
            pseudo_z, pseudo_y = vz[expected_pseudo], vy[expected_pseudo]
            sample_z, sample_y = vz[expected_sample], vy[expected_sample]
            try:
                _, pseudo_summary, _ = _independent_posterior(
                    z, h0, observed, pseudo_z, pseudo_y, distances, sigma_mpc,
                    detection_threshold_mpc, rate, bandwidth_factor,
                )
                reference = pseudo_summary["posterior_mean"]
            except (ValueError, RuntimeError, ArithmeticError, FloatingPointError):
                reference = float("nan")

        try:
            posterior, summary, diagnostics = _independent_posterior(
                z, h0, observed, sample_z, sample_y, distances, sigma_mpc,
                detection_threshold_mpc, rate, bandwidth_factor,
            )
            independently_successful = True
        except (ValueError, RuntimeError, ArithmeticError, FloatingPointError):
            independently_successful = False
            posterior = None
            summary = {}
            diagnostics = {}

        if success[b] != independently_successful:
            raise AssertionError(f"success/failure label does not replay at attempt {b}")
        if success[b]:
            _assert_close(means[b], summary["posterior_mean"], f"posterior mean {b}")
            _assert_close(references[b], reference, f"reference posterior mean {b}")
            _assert_close(result["posterior_median_draws_by_attempt"][b],
                          summary["posterior_median"], f"posterior median {b}")
            _assert_close(result["ci68_low_draws_by_attempt"][b],
                          summary["ci68_low"], f"credible lower endpoint {b}")
            _assert_close(result["ci68_high_draws_by_attempt"][b],
                          summary["ci68_high"], f"credible upper endpoint {b}")
            _assert_close(result["inverse_weight_ess_by_attempt"][b],
                          diagnostics["ess"], f"inverse-weight ESS {b}")
            _assert_close(result["inverse_weight_maximum_by_attempt"][b],
                          diagnostics["max_inverse"], f"maximum inverse weight {b}")
            if bool(result["coefficient_boundary_limited_by_attempt"][b]) != diagnostics["boundary"]:
                raise AssertionError(f"boundary flag mismatch at attempt {b}")
            if not np.allclose(posterior_draws[success_cursor], posterior, atol=5e-7, rtol=5e-6):
                raise AssertionError(f"posterior density mismatch at attempt {b}")
            success_cursor += 1
        else:
            failure_indices.append(b)

    records = result.get("failure_records", [])
    if [int(row["bootstrap_index"]) for row in records] != failure_indices:
        raise AssertionError("failure records do not exactly match failed attempt indices")
    if int(result["n_attempted"]) != int(result["n_successful"]) + len(records):
        raise AssertionError("failure records do not reconcile with attempt counts")
    if int(result["n_successful"]) >= 2:
        successful_means = means[success]
        expected_se = float(np.std(successful_means, ddof=1))
        _assert_close(result["raw_bootstrap_standard_error"], expected_se, "bootstrap SE")
        expected_interval = np.quantile(successful_means, [0.025, 0.975])
        if not np.allclose(result["posterior_mean_percentile_interval_95"],
                           expected_interval, atol=5e-7, rtol=5e-6):
            raise AssertionError("posterior-mean percentile interval mismatch")
    return {
        "status": "PASS",
        "n_attempted": attempted,
        "n_successful": int(np.sum(success)),
        "n_failed": int(np.sum(~success)),
        "independently_reconstructed_posteriors": success_cursor,
        "method": str(method),
    }
