"""Independent reconstruction audit for smooth-selection pilot artifacts.

This module does not import the experiment runner or production likelihood.
"""
from __future__ import annotations

import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import expit, roots_hermitenorm
from scipy.stats import norm


def _sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _density(z, mean, sd):
    value = np.exp(-0.5 * ((z - mean) / sd) ** 2)
    return value / trapezoid(value, z)


def _selection(distance, d50, width):
    return expit(-(np.asarray(distance, dtype=float) - d50) / width)


def _alpha(z, density, h0, sigma, d50, width, order):
    nodes, weights = roots_hermitenorm(order)
    weights = weights / np.sqrt(2.0 * np.pi)
    means = 299_792.458 * z[None, :] / h0[:, None]
    trial_d = means[:, :, None] + sigma * nodes[None, None, :]
    accept = np.sum(_selection(trial_d, d50, width) * weights[None, None, :], axis=2)
    return trapezoid(accept * density[None, :], z, axis=1)


def _loglike(z, gc, gm, h0, distances, sigma, f, d50, width, order):
    ac = _alpha(z, gc, h0, sigma, d50, width, order)
    am = _alpha(z, gm, h0, sigma, d50, width, order)
    beta = (1.0 - f) * ac + f * am
    if np.any(beta <= 0) or not np.all(np.isfinite(beta)):
        raise ValueError("invalid independently reconstructed selection efficiency")
    means = 299_792.458 * z[None, :] / h0[:, None]
    total = np.zeros_like(h0)
    for d in distances:
        qc = trapezoid(norm.pdf(d, loc=means, scale=sigma) * gc[None, :], z, axis=1)
        qm = trapezoid(norm.pdf(d, loc=means, scale=sigma) * gm[None, :], z, axis=1)
        numerator = (1.0 - f) * qc + f * qm
        if np.any(numerator <= 0) or not np.all(np.isfinite(numerator)):
            raise ValueError("invalid independently reconstructed event density")
        total += np.log(numerator) + np.log(_selection(d, d50, width))
    return total - len(distances) * np.log(beta), beta


def _posterior(grid, loglike):
    w = np.exp(loglike - np.max(loglike))
    area = float(trapezoid(w, grid))
    density = w / area
    cdf = np.r_[0.0, cumulative_trapezoid(density, grid)]
    cdf /= cdf[-1]
    quantile = {p: float(np.interp(p, cdf, grid)) for p in (0.025, 0.16, 0.5, 0.84, 0.975)}
    return {
        "density": density,
        "mean": float(trapezoid(grid * density, grid)),
        "ci95_low": quantile[0.025], "ci68_low": quantile[0.16],
        "median": quantile[0.5], "ci68_high": quantile[0.84],
        "ci95_high": quantile[0.975],
    }


def audit(output_dir: Path):
    required = {name: output_dir / name for name in
                ("events.csv", "fits.csv", "summary.json", "manifest.json")}
    for name, path in required.items():
        if not path.is_file():
            raise FileNotFoundError(f"missing artifact: {name}")
    manifest = json.loads(required["manifest.json"].read_text(encoding="utf-8"))
    for name in ("events.csv", "fits.csv", "summary.json"):
        record = manifest.get("files", {}).get(name, {})
        if record.get("sha256") != _sha(required[name]):
            raise ValueError(f"checksum mismatch: {name}")
    events = pd.read_csv(required["events.csv"])
    fits = pd.read_csv(required["fits.csv"])
    summary = json.loads(required["summary.json"].read_text(encoding="utf-8"))
    if events.duplicated(["replicate", "d50_true_mpc", "event_index"]).any():
        raise ValueError("duplicate event keys")
    if fits.duplicated(["replicate", "d50_true_mpc", "d50_model_mpc", "f_model_preselection"]).any():
        raise ValueError("duplicate fit keys")
    reps = int(summary["replicates_per_true_curve"])
    nevents = int(summary["n_events_per_dataset"])
    true_levels = [float(x) for x in summary["d50_true_levels_mpc"]]
    model_levels = [float(x) for x in summary["d50_model_levels_mpc"]]
    f_levels = [float(x) for x in summary["f_model_levels"]]
    expected_events = reps * nevents * len(true_levels)
    expected_fits = reps * len(true_levels) * len(model_levels) * len(f_levels)
    if len(events) != expected_events or len(fits) != expected_fits:
        raise ValueError("artifact row counts do not match declared configuration")
    if not np.all(np.isfinite(events["distance_mpc"])):
        raise ValueError("non-finite event distance")
    if not events["component_for_audit_only"].isin(["catalogue", "missing"]).all():
        raise ValueError("invalid latent component label")
    zspec, hspec = summary["redshift_support"], summary["h0_grid"]
    z = np.linspace(float(zspec[0]), float(zspec[1]), int(zspec[2]))
    h0 = np.linspace(float(hspec[0]), float(hspec[1]), int(hspec[2]))
    gc, gm = _density(z, 0.12, 0.05), _density(z, 0.25, 0.07)
    sigma, width = float(summary["sigma_mpc"]), float(summary["width_mpc"])
    truth = float(summary["h0_true"])
    n_reconstructed = 0
    max_ll_diff = 0.0
    max_density_diff = 0.0

    for (rep, dtrue), group in events.groupby(["replicate", "d50_true_mpc"], sort=True):
        group = group.sort_values("event_index")
        if len(group) != nevents or not np.array_equal(group["event_index"].to_numpy(), np.arange(nevents)):
            raise ValueError("event indices/counts inconsistent")
        distances = group["distance_mpc"].to_numpy(dtype=float)
        cellfits = fits[(fits.replicate == rep) & (fits.d50_true_mpc == dtrue)]
        for _, row in cellfits.iterrows():
            dmodel, fmodel = float(row["d50_model_mpc"]), float(row["f_model_preselection"])
            ll, beta = _loglike(z, gc, gm, h0, distances, sigma, fmodel, dmodel, width, 40)
            saved_grid = np.asarray(json.loads(row["h0_grid_json"]), dtype=float)
            saved_ll = np.asarray(json.loads(row["log_likelihood_grid_json"]), dtype=float)
            saved_density = np.asarray(json.loads(row["posterior_density_grid_json"]), dtype=float)
            saved_beta = np.asarray(json.loads(row["beta_grid_json"]), dtype=float)
            if not np.array_equal(saved_grid, h0):
                raise ValueError("H0 grid mismatch")
            if not np.allclose(saved_ll, ll, rtol=2e-7, atol=2e-7):
                raise ValueError("likelihood grid reconstruction mismatch")
            if not np.allclose(saved_beta, beta, rtol=2e-7, atol=2e-9):
                raise ValueError("selection efficiency reconstruction mismatch")
            max_ll_diff = max(max_ll_diff, float(np.max(np.abs(saved_ll - ll))))
            post = _posterior(h0, ll)
            max_density_diff = max(max_density_diff, float(np.max(np.abs(saved_density - post["density"]))))
            if not np.allclose(saved_density, post["density"], rtol=2e-7, atol=2e-8):
                raise ValueError("posterior density reconstruction mismatch")
            for key in ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci95_low", "ci95_high"):
                mapped = {"posterior_mean": "mean", "posterior_median": "median"}.get(key, key)
                if not np.isclose(float(row[key]), post[mapped], rtol=2e-7, atol=2e-6):
                    raise ValueError(f"posterior summary mismatch: {key}")
            if bool(row["ci68_covered"]) != (post["ci68_low"] <= truth <= post["ci68_high"]):
                raise ValueError("68% coverage flag mismatch")
            if bool(row["ci95_covered"]) != (post["ci95_low"] <= truth <= post["ci95_high"]):
                raise ValueError("95% coverage flag mismatch")
            n_reconstructed += 1

    for item in summary["aggregate_metrics"]:
        group = fits[(fits.d50_true_mpc == float(item["d50_true_mpc"])) &
                     (fits.d50_model_mpc == float(item["d50_model_mpc"])) &
                     (fits.f_model_preselection == float(item["f_model_preselection"]))]
        err = group["posterior_median"].to_numpy(dtype=float) - truth
        checks = {
            "bias_median": float(np.mean(err)),
            "rmse_median": float(np.sqrt(np.mean(err ** 2))),
            "coverage68": float(group["ci68_covered"].astype(bool).mean()),
            "coverage95": float(group["ci95_covered"].astype(bool).mean()),
            "mean_ci68_width": float(group["ci68_width"].mean()),
        }
        for key, expected in checks.items():
            if not np.isclose(float(item[key]), expected, rtol=2e-7, atol=2e-6):
                raise ValueError(f"aggregate mismatch: {key}")
    return {
        "status": "PASS",
        "events_structurally_checked": len(events),
        "posterior_fits_independently_reconstructed": n_reconstructed,
        "expected_posterior_fits": expected_fits,
        "max_abs_log_likelihood_difference": max_ll_diff,
        "max_abs_posterior_density_difference": max_density_diff,
        "manifest_checksums_verified": 3,
        "scope": "synthetic artifact reconstruction; not physical or external-data validation",
    }


def main():
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
