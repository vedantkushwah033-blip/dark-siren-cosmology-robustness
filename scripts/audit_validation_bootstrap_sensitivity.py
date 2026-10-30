"""Independently audit persisted validation-bootstrap posterior propagation."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

C_KM_S = 299_792.458


def _quad_weights(grid: np.ndarray) -> np.ndarray:
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _posterior_from_density(
    z: np.ndarray, density: np.ndarray, h0: np.ndarray,
    distances: np.ndarray, sigma: float, threshold: float,
) -> tuple[float, float, float, float]:
    q = _quad_weights(z)
    valid = (density > 0) & (q > 0)
    log_mass = np.log(density[valid]) + np.log(q[valid])
    means = C_KM_S * z[valid][None, :] / h0[:, None]
    event_terms = norm.logpdf(
        distances[None, :, None],
        loc=means[:, None, :],
        scale=sigma,
    ) + log_mass[None, None, :]
    joint = np.sum(logsumexp(event_terms, axis=2), axis=1)
    detection = norm.cdf(
        (C_KM_S * z[None, :] / h0[:, None] - threshold) / sigma
    )
    beta = np.sum(density[None, :] * q[None, :] * detection, axis=1)
    if np.any(~np.isfinite(beta)) or np.any(beta <= 0) or np.any(beta > 1 + 1e-12):
        raise AssertionError("independent reconstruction found invalid detection efficiency")
    selected_loglike = joint - distances.size * np.log(beta)
    posterior = np.exp(selected_loglike - np.max(selected_loglike))
    posterior /= trapezoid(posterior, h0)
    cdf = np.r_[0.0, cumulative_trapezoid(posterior, h0)]
    cdf /= cdf[-1]
    low, median, high = (float(np.interp(p, cdf, h0)) for p in (0.16, 0.50, 0.84))
    mean = float(trapezoid(h0 * posterior, h0))
    return mean, median, low, high


def _assert_close(actual: np.ndarray, expected: np.ndarray, label: str) -> None:
    if actual.shape != expected.shape or not np.allclose(
        actual, expected, rtol=2e-7, atol=2e-8, equal_nan=False
    ):
        raise AssertionError(f"{label} does not match independent reconstruction")


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir", type=Path, required=True)
    parser.add_argument("--pilot-csv", type=Path, required=True)
    args = parser.parse_args()

    summary = json.loads((args.artifact_dir / "summary.json").read_text(encoding="utf-8"))
    inputs = json.loads((args.artifact_dir / "inputs.json").read_text(encoding="utf-8"))
    if summary.get("schema") != "validation_bootstrap_artifact_v1":
        raise AssertionError("unsupported bootstrap artifact schema")
    if inputs.get("schema") != "validation_bootstrap_inputs_v1":
        raise AssertionError("unsupported bootstrap input schema")
    if summary.get("finite_population_correction_applied") is not False:
        raise AssertionError("artifact must explicitly disclose absent finite-population correction")
    payload = inputs.get("inference_payload")
    if not isinstance(payload, dict):
        raise AssertionError("missing inference-only input payload")
    forbidden = {"parent_redshift", "true_inclusion_probability", "estimated_inclusion_probability"}
    if forbidden.intersection(payload):
        raise AssertionError("hidden-parent/scoring truth leaked into bootstrap inference inputs")

    pilot = pd.read_csv(args.pilot_csv)
    rows = pilot.loc[
        (pilot["replicate"] == int(inputs["replicate"]))
        & (pilot["arm"] == str(inputs["arm"]))
        & (pilot["status"] == "success")
    ]
    if len(rows) != 1:
        raise AssertionError("pilot CSV does not contain exactly one matching successful row")
    pilot_payload = json.loads(rows.iloc[0]["inference_audit_payload_json"])
    if pilot_payload != payload:
        raise AssertionError("persisted bootstrap inputs differ from the selected pilot inference payload")

    with np.load(args.artifact_dir / "bootstrap_draws.npz", allow_pickle=False) as data:
        arrays = {name: np.asarray(data[name]) for name in data.files}
    required = {
        "density_draws", "intercept_draws", "slope_draws",
        "coefficient_boundary_limited_draws", "inverse_weight_ess_draws",
        "inverse_weight_maximum_draws", "inverse_weight_p99_draws",
        "posterior_mean_draws", "posterior_median_draws",
        "posterior_ci68_low_draws", "posterior_ci68_high_draws",
    }
    if set(arrays) != required:
        raise AssertionError("bootstrap NPZ has missing or unexpected arrays")

    n_attempted = int(summary["n_attempted"])
    n_successful = int(summary["n_successful"])
    n_failed = int(summary["n_failed"])
    if n_attempted < 2 or n_successful + n_failed != n_attempted:
        raise AssertionError("attempted/successful/failed bootstrap counts do not reconcile")
    if len(summary.get("failure_reasons", [])) != n_failed:
        raise AssertionError("failure reason count does not match failed resamples")
    if summary.get("seed") is None or not isinstance(summary["seed"], int):
        raise AssertionError("bootstrap seed is missing")
    z = np.asarray(payload["redshift_grid"], dtype=float)
    h0 = np.asarray(payload["h0_grid"], dtype=float)
    density_draws = arrays["density_draws"]
    if density_draws.shape != (n_successful, z.size):
        raise AssertionError("density draw dimensions disagree with declared success count")
    if not np.all(np.isfinite(density_draws)) or np.any(density_draws < 0):
        raise AssertionError("density draws contain invalid values")
    if n_successful and not np.allclose(
        trapezoid(density_draws, z, axis=1), 1.0, rtol=1e-7, atol=1e-8
    ):
        raise AssertionError("one or more bootstrap densities fail unit normalization")

    scalar_draw_names = (
        "intercept_draws", "slope_draws", "coefficient_boundary_limited_draws",
        "inverse_weight_ess_draws", "inverse_weight_maximum_draws",
        "inverse_weight_p99_draws", "posterior_mean_draws",
        "posterior_median_draws", "posterior_ci68_low_draws",
        "posterior_ci68_high_draws",
    )
    for name in scalar_draw_names:
        if arrays[name].shape != (n_successful,):
            raise AssertionError(f"{name} length disagrees with successful draw count")
    for name in scalar_draw_names:
        if name != "coefficient_boundary_limited_draws" and not np.all(np.isfinite(arrays[name])):
            raise AssertionError(f"{name} contains non-finite values")
    if n_successful and (
        np.any(arrays["inverse_weight_ess_draws"] <= 0)
        or np.any(arrays["inverse_weight_maximum_draws"] <= 0)
        or np.any(arrays["inverse_weight_p99_draws"] <= 0)
        or np.any(arrays["posterior_ci68_low_draws"] >= arrays["posterior_ci68_high_draws"])
    ):
        raise AssertionError("positive-weight or interval-order invariant failed")

    reconstructed = [
        _posterior_from_density(
            z, density, h0,
            np.asarray(payload["observed_distances_mpc"], dtype=float),
            float(payload["sigma_mpc"]),
            float(payload["detection_threshold_mpc"]),
        )
        for density in density_draws
    ]
    if n_successful:
        values = np.asarray(reconstructed, dtype=float)
        _assert_close(arrays["posterior_mean_draws"], values[:, 0], "posterior means")
        _assert_close(arrays["posterior_median_draws"], values[:, 1], "posterior medians")
        _assert_close(arrays["posterior_ci68_low_draws"], values[:, 2], "posterior lower limits")
        _assert_close(arrays["posterior_ci68_high_draws"], values[:, 3], "posterior upper limits")
    expected_percentiles = (
        np.quantile(arrays["posterior_mean_draws"], [0.025, 0.975]).tolist()
        if n_successful >= 2 else None
    )
    stored_percentiles = summary.get("posterior_mean_bootstrap_percentile_95")
    if expected_percentiles is None:
        if stored_percentiles is not None:
            raise AssertionError("percentile interval should be null with fewer than two successful draws")
    elif not np.allclose(stored_percentiles, expected_percentiles, rtol=1e-10, atol=1e-10):
        raise AssertionError("stored posterior-mean bootstrap percentiles do not match draws")

    print(json.dumps({
        "status": "PASS",
        "replicate": int(inputs["replicate"]),
        "arm": inputs["arm"],
        "n_attempted": n_attempted,
        "n_successful": n_successful,
        "n_failed": n_failed,
        "posterior_draws_independently_reconstructed": n_successful,
        "scope": "posterior propagation audit only; bootstrap refit generation is not independently rerun",
    }, indent=2))


if __name__ == "__main__":
    main()
