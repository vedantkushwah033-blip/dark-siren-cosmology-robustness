"""Run a synthetic smooth-selection misspecification experiment.

This experiment tests a two-sided shift in a logistic selection curve. It is a
toy methodological study, not a model of a real GW detector or an LVK analysis.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
import time
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import roots_hermitenorm
from scipy.stats import norm

from dark_siren_robustness.selection_conditioned_hidden_host import C_KM_S
from dark_siren_robustness.smooth_selection_hidden_host import (
    _component_selection_efficiency,
    selection_conditioned_smooth_hidden_host_log_likelihood,
    smooth_selection_probability,
)
from dark_siren_robustness.host_density_benchmark import posterior_summary

TRUE_D50_LEVELS = (700.0, 1000.0, 1300.0)
MODEL_D50_LEVELS = (700.0, 1000.0, 1300.0)
F_TRUE = 0.30
F_MODEL_LEVELS = (0.10, 0.30, 0.50)
H0_TRUE = 70.0
SIGMA_MPC = 140.0
WIDTH_MPC = 120.0
QUADRATURE_ORDER = 40


def sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def stable_seed(seed: int, config: dict) -> int:
    encoded = json.dumps([int(seed), config], sort_keys=True, separators=(",", ":"))
    return int.from_bytes(hashlib.sha256(encoded.encode()).digest()[:8], "big") % (2**32 - 1)


def normalized_density(z: np.ndarray, mean: float, sd: float) -> np.ndarray:
    raw = np.exp(-0.5 * ((z - mean) / sd) ** 2)
    return raw / float(trapezoid(raw, z))


def selected_z_sample(z, density, h0, sigma, d50, width, rng):
    nodes, weights = roots_hermitenorm(QUADRATURE_ORDER)
    weights = weights / np.sqrt(2.0 * np.pi)
    mu = C_KM_S * z / h0
    d = mu[:, None] + sigma * nodes[None, :]
    acceptance = np.sum(
        smooth_selection_probability(d, d50_mpc=d50, width_mpc=width)
        * weights[None, :], axis=1
    )
    weighted = density * acceptance
    cdf = np.r_[0.0, cumulative_trapezoid(weighted, z)]
    cdf /= cdf[-1]
    return float(np.interp(rng.random(), cdf, z))


def generate_selected_dataset(z, gc, gm, *, h0, f, sigma, d50, width, n_events, seed):
    rng = np.random.default_rng(seed)
    alpha_c = float(_component_selection_efficiency(
        z, gc, np.array([h0]), sigma, d50, width, QUADRATURE_ORDER
    )[0])
    alpha_m = float(_component_selection_efficiency(
        z, gm, np.array([h0]), sigma, d50, width, QUADRATURE_ORDER
    )[0])
    beta = (1.0 - f) * alpha_c + f * alpha_m
    p_missing_detected = f * alpha_m / beta
    rows = []
    for i in range(n_events):
        missing = bool(rng.random() < p_missing_detected)
        density = gm if missing else gc
        # Draw z from its selected marginal, then draw d conditional on
        # acceptance at that fixed z. Redrawing z after a failed distance draw
        # would apply the selection weight twice and distort the population.
        z_true = selected_z_sample(z, density, h0, sigma, d50, width, rng)
        while True:
            distance = float(rng.normal(C_KM_S * z_true / h0, sigma))
            p_accept = float(smooth_selection_probability(
                distance, d50_mpc=d50, width_mpc=width
            ))
            if rng.random() < p_accept:
                break
        rows.append({
            "event_index": i,
            "distance_mpc": distance,
            "z_true_for_audit_only": z_true,
            "component_for_audit_only": "missing" if missing else "catalogue",
            "d50_true_mpc": float(d50),
        })
    return rows


def run_experiment(*, replicates, n_events, seed, h0_points, z_points):
    if replicates < 1 or n_events < 1:
        raise ValueError("replicates and n_events must be positive")
    if h0_points < 31 or z_points < 101:
        raise ValueError("h0_points must be >=31 and z_points >=101")
    z = np.linspace(0.01, 0.40, z_points)
    h0_grid = np.linspace(50.0, 90.0, h0_points)
    gc = normalized_density(z, 0.12, 0.05)
    gm = normalized_density(z, 0.25, 0.07)
    event_rows, fit_rows = [], []
    started = time.perf_counter()

    for d50_true in TRUE_D50_LEVELS:
        for rep in range(replicates):
            config = {"d50_true_mpc": d50_true, "replicate": rep}
            data_seed = stable_seed(seed, config)
            events = generate_selected_dataset(
                z, gc, gm, h0=H0_TRUE, f=F_TRUE, sigma=SIGMA_MPC,
                d50=d50_true, width=WIDTH_MPC, n_events=n_events, seed=data_seed,
            )
            for row in events:
                event_rows.append({
                    "replicate": rep, "data_seed": data_seed,
                    "d50_true_mpc": d50_true, **row,
                })
            distances = np.array([x["distance_mpc"] for x in events])
            for d50_model in MODEL_D50_LEVELS:
                for f_model in F_MODEL_LEVELS:
                    fit = selection_conditioned_smooth_hidden_host_log_likelihood(
                        z, gc, gm, h0_grid, distances,
                        sigma=SIGMA_MPC, f_missing=f_model,
                        d50_mpc=d50_model, width_mpc=WIDTH_MPC,
                        quadrature_order=QUADRATURE_ORDER,
                    )
                    post = posterior_summary(h0_grid, fit["log_likelihood"])
                    density = np.asarray(json.loads(post["posterior_density_grid_json"]), dtype=float)
                    cdf = np.r_[0.0, cumulative_trapezoid(density, h0_grid)]
                    cdf /= cdf[-1]
                    median = float(np.interp(0.5, cdf, h0_grid))
                    lo95, hi95 = float(np.interp(0.025, cdf, h0_grid)), float(np.interp(0.975, cdf, h0_grid))
                    fit_rows.append({
                        "replicate": rep, "data_seed": data_seed,
                        "d50_true_mpc": d50_true, "d50_model_mpc": d50_model,
                        "f_true_preselection": F_TRUE, "f_model_preselection": f_model,
                        "h0_true": H0_TRUE, "sigma_mpc": SIGMA_MPC,
                        "width_mpc": WIDTH_MPC, "n_events": n_events,
                        "posterior_mean": float(post["posterior_mean"]),
                        "posterior_median": median,
                        "ci68_low": float(post["ci68_low"]),
                        "ci68_high": float(post["ci68_high"]),
                        "ci95_low": lo95, "ci95_high": hi95,
                        "ci68_covered": bool(post["ci68_low"] <= H0_TRUE <= post["ci68_high"]),
                        "ci95_covered": bool(lo95 <= H0_TRUE <= hi95),
                        "ci68_width": float(post["ci68_high"] - post["ci68_low"]),
                        "posterior_integral": float(post["posterior_integral"]),
                        "h0_grid_json": json.dumps(h0_grid.tolist(), separators=(",", ":")),
                        "log_likelihood_grid_json": json.dumps(fit["log_likelihood"].tolist(), separators=(",", ":")),
                        "posterior_density_grid_json": post["posterior_density_grid_json"],
                        "beta_grid_json": json.dumps(fit["beta"].tolist(), separators=(",", ":")),
                    })

    events_df, fits_df = pd.DataFrame(event_rows), pd.DataFrame(fit_rows)
    aggregate = []
    for (dtrue, dmodel, fmodel), group in fits_df.groupby(
        ["d50_true_mpc", "d50_model_mpc", "f_model_preselection"], sort=True
    ):
        error = group["posterior_median"].to_numpy() - H0_TRUE
        aggregate.append({
            "d50_true_mpc": float(dtrue), "d50_model_mpc": float(dmodel),
            "f_model_preselection": float(fmodel), "n_replicates": int(len(group)),
            "bias_median": float(np.mean(error)),
            "rmse_median": float(np.sqrt(np.mean(error ** 2))),
            "coverage68": float(group["ci68_covered"].astype(bool).mean()),
            "coverage95": float(group["ci95_covered"].astype(bool).mean()),
            "mean_ci68_width": float(group["ci68_width"].mean()),
        })
    summary = {
        "status": "synthetic_smooth_selection_misspecification_pilot",
        "claim_boundary": "synthetic toy only; not real LVK data or an observational H0 measurement",
        "seed": int(seed), "replicates_per_true_curve": int(replicates),
        "n_events_per_dataset": int(n_events), "h0_true": H0_TRUE,
        "f_true_preselection": F_TRUE, "f_model_levels": list(F_MODEL_LEVELS),
        "d50_true_levels_mpc": list(TRUE_D50_LEVELS),
        "d50_model_levels_mpc": list(MODEL_D50_LEVELS),
        "width_mpc": WIDTH_MPC, "sigma_mpc": SIGMA_MPC,
        "redshift_support": [float(z[0]), float(z[-1]), int(z_points)],
        "h0_grid": [float(h0_grid[0]), float(h0_grid[-1]), int(h0_points)],
        "selection": "logistic probability expit(-(distance-d50)/width)",
        "component_densities": "normalized Gaussian toy densities on redshift grid: catalogue mean 0.12 SD 0.05; missing mean 0.25 SD 0.07",
        "aggregate_metrics": aggregate,
        "runtime_seconds": float(time.perf_counter() - started),
        "limitations": [
            "Synthetic Gaussian distance likelihood and logistic selection curve, not a detector response.",
            "No event-rate information, real galaxy catalogue, peculiar velocities, or redshift uncertainty.",
            "Fixed H0 truth; coverage is conditional on that truth and Monte Carlo uncertainty depends on replicate count.",
            "Selection-curve mismatch is a controlled sensitivity experiment, not evidence of real-world bias.",
        ],
    }
    return events_df, fits_df, summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--n-events", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--h0-points", type=int, default=61)
    parser.add_argument("--z-points", type=int, default=151)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    start = time.perf_counter()
    events, fits, summary = run_experiment(
        replicates=args.replicates, n_events=args.n_events, seed=args.seed,
        h0_points=args.h0_points, z_points=args.z_points,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    paths = {
        "events.csv": args.output_dir / "events.csv",
        "fits.csv": args.output_dir / "fits.csv",
        "summary.json": args.output_dir / "summary.json",
    }
    events.to_csv(paths["events.csv"], index=False)
    fits.to_csv(paths["fits.csv"], index=False)
    summary["runtime_seconds_wall_including_serialization"] = time.perf_counter() - start
    paths["summary.json"].write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    manifest = {
        "experiment": "smooth_selection_misspecification_v1",
        "git_commit": __import__("subprocess").run(
            ["git", "rev-parse", "HEAD"], capture_output=True, text=True
        ).stdout.strip() or "unavailable",
        "python": sys.version, "platform": platform.platform(),
        "command": " ".join(sys.argv),
        "configuration": {
            "replicates": args.replicates, "n_events": args.n_events,
            "seed": args.seed, "h0_points": args.h0_points,
            "z_points": args.z_points, "d50_true_levels_mpc": TRUE_D50_LEVELS,
            "d50_model_levels_mpc": MODEL_D50_LEVELS, "f_model_levels": F_MODEL_LEVELS,
        },
        "files": {
            name: {"sha256": sha256(path), "bytes": path.stat().st_size}
            for name, path in paths.items()
        },
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
