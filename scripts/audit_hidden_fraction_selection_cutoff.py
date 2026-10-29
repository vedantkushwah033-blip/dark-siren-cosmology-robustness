"""Independently reconstruct hidden-fraction/cutoff pilot artifacts.

This audit intentionally does not import the experiment runner or its likelihood
function. It recomputes component densities, selection efficiencies, likelihood
grids, posterior summaries, coverage and aggregate metrics from saved event data.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _density(z: np.ndarray, mean: float, sd: float) -> np.ndarray:
    value = np.exp(-0.5 * ((z - mean) / sd) ** 2)
    return value / trapezoid(value, z)


def _posterior(grid: np.ndarray, loglike: np.ndarray) -> dict:
    weight = np.exp(loglike - np.max(loglike))
    area = float(trapezoid(weight, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("independent posterior normalization failed")
    density = weight / area
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    q = {p: float(np.interp(p, cdf, grid)) for p in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "density": density,
        "mean": float(trapezoid(grid * density, grid)),
        "ci95_low": q[0.025], "ci68_low": q[0.16], "median": q[0.5],
        "ci68_high": q[0.84], "ci95_high": q[0.975],
        "integral": float(trapezoid(density, grid)),
    }


def _independent_loglike(z, g_cat, g_miss, h0, distances, sigma, f, d_min, d_max):
    means = 299_792.458 * z[None, :] / h0[:, None]
    accept = norm.cdf((d_max - means) / sigma) - norm.cdf((d_min - means) / sigma)
    alpha_c = trapezoid(accept * g_cat[None, :], z, axis=1)
    alpha_m = trapezoid(accept * g_miss[None, :], z, axis=1)
    beta = (1.0 - f) * alpha_c + f * alpha_m
    if np.any(beta <= 0) or not np.all(np.isfinite(beta)):
        raise ValueError("independent mixture selection efficiency invalid")
    total = np.zeros_like(h0)
    for distance in distances:
        q_c = trapezoid(norm.pdf(distance, loc=means, scale=sigma) * g_cat[None, :], z, axis=1)
        q_m = trapezoid(norm.pdf(distance, loc=means, scale=sigma) * g_miss[None, :], z, axis=1)
        numerator = (1.0 - f) * q_c + f * q_m
        if np.any(numerator <= 0) or not np.all(np.isfinite(numerator)):
            raise ValueError("independent mixture numerator invalid")
        total += np.log(numerator)
    return total - len(distances) * np.log(beta), beta


def audit(output_dir: Path) -> dict:
    events_path = output_dir / "events.csv"
    fits_path = output_dir / "fits.csv"
    summary_path = output_dir / "summary.json"
    manifest_path = output_dir / "manifest.json"
    for path in (events_path, fits_path, summary_path, manifest_path):
        if not path.is_file():
            raise FileNotFoundError(f"missing required artifact: {path.name}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for name in ("events.csv", "fits.csv", "summary.json"):
        record = manifest.get("files", {}).get(name)
        if not isinstance(record, dict) or record.get("sha256") != _sha(output_dir / name):
            raise ValueError(f"manifest checksum mismatch for {name}")
    events = pd.read_csv(events_path)
    fits = pd.read_csv(fits_path)
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    required_event_cols = {"replicate", "data_seed", "d_max_true_mpc", "event_index",
                           "z_true_for_audit_only", "component_for_audit_only", "distance_mpc",
                           "d_min_true_mpc", "d_max_true_mpc"}
    required_fit_cols = {"replicate", "data_seed", "d_max_true_mpc", "d_max_model_mpc",
                         "f_model_preselection", "h0_true", "posterior_mean", "posterior_median",
                         "ci68_low", "ci68_high", "ci95_low", "ci95_high", "ci68_covered",
                         "ci95_covered", "posterior_density_grid_json", "log_likelihood_grid_json",
                         "beta_grid_json", "h0_grid_json"}
    if not required_event_cols.issubset(events.columns):
        raise ValueError(f"event columns missing: {sorted(required_event_cols - set(events.columns))}")
    if not required_fit_cols.issubset(fits.columns):
        raise ValueError(f"fit columns missing: {sorted(required_fit_cols - set(fits.columns))}")
    if events.duplicated(["replicate", "d_max_true_mpc", "event_index"]).any():
        raise ValueError("duplicate event keys")
    if fits.duplicated(["replicate", "d_max_true_mpc", "f_model_preselection"]).any():
        raise ValueError("duplicate fit keys")
    nrep = int(summary["replicates_per_cutoff"])
    nevents = int(summary["n_events_per_dataset"])
    cutoffs = [float(x) for x in summary["d_max_levels_mpc"]]
    f_models = [float(x) for x in summary["f_model_levels"]]
    expected_events = nrep * nevents * len(cutoffs)
    expected_fits = nrep * len(cutoffs) * len(f_models)
    if len(events) != expected_events or len(fits) != expected_fits:
        raise ValueError(f"row count mismatch: events {len(events)}/{expected_events}, fits {len(fits)}/{expected_fits}")
    if not np.all(np.isfinite(events["distance_mpc"])) or not np.all(np.isfinite(events["z_true_for_audit_only"])):
        raise ValueError("non-finite event data")
    if not events["component_for_audit_only"].isin(["catalogue", "missing"]).all():
        raise ValueError("unknown component label")
    if (events["distance_mpc"] < events["d_min_true_mpc"]).any() or (events["distance_mpc"] > events["d_max_true_mpc"]).any():
        raise ValueError("event distance violates the true selection region")
    if not np.allclose(fits["d_max_true_mpc"], fits["d_max_model_mpc"], rtol=0, atol=1e-12):
        raise ValueError("this pilot requires correctly specified selection cutoff in the model")

    z_spec = summary["redshift_support"]
    z = np.linspace(float(z_spec[0]), float(z_spec[1]), int(z_spec[2]))
    h_spec = summary["h0_grid"]
    h0 = np.linspace(float(h_spec[0]), float(h_spec[1]), int(h_spec[2]))
    g_cat = _density(z, 0.12, 0.05)
    g_miss = _density(z, 0.25, 0.07)
    sigma = float(summary["sigma_mpc"])
    h0_true = float(summary["h0_true"])
    posterior_reconstructions = 0
    max_loglike_abs_diff = 0.0
    max_density_abs_diff = 0.0

    for (rep, cutoff), group in events.groupby(["replicate", "d_max_true_mpc"], sort=True):
        group = group.sort_values("event_index")
        if len(group) != nevents or not np.array_equal(group["event_index"].to_numpy(), np.arange(nevents)):
            raise ValueError("event count or event index sequence mismatch")
        distances = group["distance_mpc"].to_numpy(dtype=float)
        for _, fit in fits[(fits.replicate == rep) & (fits.d_max_true_mpc == cutoff)].iterrows():
            f = float(fit["f_model_preselection"])
            loglike, beta = _independent_loglike(z, g_cat, g_miss, h0, distances, sigma, f, 0.0, float(cutoff))
            saved_h0 = np.asarray(json.loads(fit["h0_grid_json"]), dtype=float)
            saved_loglike = np.asarray(json.loads(fit["log_likelihood_grid_json"]), dtype=float)
            saved_density = np.asarray(json.loads(fit["posterior_density_grid_json"]), dtype=float)
            saved_beta = np.asarray(json.loads(fit["beta_grid_json"]), dtype=float)
            if not np.array_equal(saved_h0, h0):
                raise ValueError("saved H0 grid differs from configuration")
            if saved_loglike.shape != loglike.shape or saved_density.shape != h0.shape or saved_beta.shape != h0.shape:
                raise ValueError("saved grid arrays have invalid shapes")
            log_diff = float(np.max(np.abs(saved_loglike - loglike)))
            max_loglike_abs_diff = max(max_loglike_abs_diff, log_diff)
            if not np.allclose(saved_loglike, loglike, rtol=2e-7, atol=2e-7):
                raise ValueError(f"log-likelihood reconstruction mismatch at rep={rep}, cutoff={cutoff}, f={f}")
            if not np.allclose(saved_beta, beta, rtol=2e-7, atol=1e-10):
                raise ValueError("selection efficiency reconstruction mismatch")
            post = _posterior(h0, loglike)
            dens_diff = float(np.max(np.abs(saved_density - post["density"])))
            max_density_abs_diff = max(max_density_abs_diff, dens_diff)
            if not np.allclose(saved_density, post["density"], rtol=2e-7, atol=2e-8):
                raise ValueError("posterior density reconstruction mismatch")
            for key, expected in (
                ("posterior_mean", post["mean"]), ("posterior_median", post["median"]),
                ("ci68_low", post["ci68_low"]), ("ci68_high", post["ci68_high"]),
                ("ci95_low", post["ci95_low"]), ("ci95_high", post["ci95_high"]),
            ):
                if not np.isclose(float(fit[key]), expected, rtol=2e-7, atol=2e-6):
                    raise ValueError(f"posterior summary mismatch: {key}")
            if bool(fit["ci68_covered"]) != (post["ci68_low"] <= h0_true <= post["ci68_high"]):
                raise ValueError("68% coverage flag mismatch")
            if bool(fit["ci95_covered"]) != (post["ci95_low"] <= h0_true <= post["ci95_high"]):
                raise ValueError("95% coverage flag mismatch")
            posterior_reconstructions += 1

    # Independently reconstruct the primary grouped metrics and paired contrasts.
    for item in summary["aggregate_metrics"]:
        sub = fits[(fits.d_max_true_mpc == float(item["d_max_true_mpc"])) &
                   (fits.f_model_preselection == float(item["f_model_preselection"]))]
        err = sub["posterior_median"].to_numpy(dtype=float) - h0_true
        if len(sub) != nrep:
            raise ValueError("aggregate group has wrong replicate count")
        expected = {
            "bias_median": float(np.mean(err)),
            "rmse_median": float(np.sqrt(np.mean(err**2))),
            "coverage68": float(sub["ci68_covered"].astype(bool).mean()),
            "coverage95": float(sub["ci95_covered"].astype(bool).mean()),
            "mean_ci68_width": float((sub["ci68_high"] - sub["ci68_low"]).mean()),
        }
        for key, value in expected.items():
            if not np.isclose(float(item[key]), value, rtol=2e-7, atol=2e-6):
                raise ValueError(f"aggregate mismatch: {key}")
        if int(item["coverage68_count"]) != int(sub["ci68_covered"].astype(bool).sum()):
            raise ValueError("coverage count mismatch")
        if int(item["coverage95_count"]) != int(sub["ci95_covered"].astype(bool).sum()):
            raise ValueError("95% coverage count mismatch")

    return {
        "status": "PASS",
        "events_reconstructed_structure": len(events),
        "posterior_grids_independently_reconstructed": posterior_reconstructions,
        "expected_posterior_grids": expected_fits,
        "max_abs_log_likelihood_difference": max_loglike_abs_diff,
        "max_abs_posterior_density_difference": max_density_abs_diff,
        "manifest_checksums_verified": 3,
        "scope": "synthetic artifact and likelihood reconstruction; not external or physical validation",
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
