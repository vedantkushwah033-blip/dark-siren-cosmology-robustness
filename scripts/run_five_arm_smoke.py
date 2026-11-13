"""Run one paired five-arm smoke replicate; this is not confirmatory evidence."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from dark_siren_robustness.five_arm_comparison import (
    build_five_arm_weights,
    weight_concentration_diagnostics,
)
from dark_siren_robustness.five_arm_likelihood import matched_five_arm_selected_log_likelihood
from dark_siren_robustness.generative_selection import (
    generate_independent_validation_sample,
    generate_synthetic_population,
)


def posterior_summary(grid: np.ndarray, log_likelihood: np.ndarray) -> dict[str, float]:
    density = np.exp(np.asarray(log_likelihood, dtype=float) - np.max(log_likelihood))
    area = float(trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, np.cumsum(0.5 * (density[:-1] + density[1:]) * np.diff(grid))]
    cdf /= cdf[-1]
    mean = float(trapezoid(grid * density, grid))
    variance = float(trapezoid((grid - mean) ** 2 * density, grid))
    q025, median, q975 = np.interp([0.025, 0.5, 0.975], cdf, grid)
    return {
        "posterior_mean": mean,
        "posterior_sd": float(np.sqrt(max(variance, 0.0))),
        "posterior_median": float(median),
        "posterior_q025": float(q025),
        "posterior_q975": float(q975),
        "map_grid": float(grid[int(np.argmax(log_likelihood))]),
    }


def sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def run_smoke(output_dir: Path, seed: int = 20261113) -> dict[str, object]:
    output_dir.mkdir(parents=True, exist_ok=True)
    observable_dir = output_dir / "observable"
    truth_dir = output_dir / "truth"
    results_dir = output_dir / "results"
    for directory in (observable_dir, truth_dir, results_dir):
        directory.mkdir(parents=True, exist_ok=True)

    sizes = (100, 300, 1000)
    population = generate_synthetic_population(
        seed=seed, n_parent=4000, n_detected=30, h0_true=70.0,
        sigma_mpc=140.0, detection_threshold_mpc=250.0,
        z_min=0.01, z_max=0.40, catalogue_floor=0.10,
        catalogue_midpoint=0.22, catalogue_width=0.05,
    )
    validation = generate_independent_validation_sample(
        seed=seed + 1, n_validation_frame=4000, n_validation=1000,
        z_min=0.01, z_max=0.40, catalogue_floor=0.10,
        catalogue_midpoint=0.22, catalogue_width=0.05,
    )
    parent = population["parent_truth"]
    catalogue = population["catalogue_observed"]
    validation_frame = validation["truth"]
    parent_z = parent["redshift"].to_numpy(float)
    catalogue_z = catalogue["redshift"].to_numpy(float)
    base_parent = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    base_catalogue = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": 140.0}
        for row in population["gw_events_observed"].itertuples(index=False)
    ]
    h0_grid = np.linspace(55.0, 85.0, 301)
    summary_rows: list[dict[str, object]] = []
    ledger: list[dict[str, object]] = []
    curve_frames: list[pd.DataFrame] = []

    for n_validation in sizes:
        replicate = f"seed{seed}_validation{n_validation}"
        validation_observed = validation["validation_observed"].iloc[:n_validation].copy()
        try:
            arms = build_five_arm_weights(
                parent_z, catalogue_z,
                validation_observed["redshift"].to_numpy(float),
                validation_observed["retained"].to_numpy(np.int8),
                catalogue_floor=0.10, catalogue_midpoint=0.22,
                catalogue_width=0.05, degree=1,
                coordinate_center=0.205, coordinate_scale=0.195,
                quadrature_order=7,
                parent_base_weights=base_parent,
                catalogue_base_weights=base_catalogue,
            )
            likelihood = matched_five_arm_selected_log_likelihood(
                h0_grid, arms, events, detection_threshold_mpc=250.0
            )
            arm_likelihoods = {
                "A_catalogue_only": likelihood.arm_a_log_likelihood,
                "B_parent_oracle": likelihood.arm_b_log_likelihood,
                "C_known_probability_ipw": likelihood.arm_c_log_likelihood,
                "D_map_plugin_ipw": likelihood.arm_d_log_likelihood,
                "E_uncertainty_integrated_ipw": likelihood.arm_e_marginal_log_likelihood,
            }
            for arm, log_likelihood in arm_likelihoods.items():
                estimates = posterior_summary(h0_grid, log_likelihood)
                summary_rows.append({
                    "replicate_id": replicate, "validation_n": n_validation,
                    "arm": arm, "status": "SUCCESS", "h0_true": 70.0,
                    "mean_error": estimates["posterior_mean"] - 70.0,
                    "median_error": estimates["posterior_median"] - 70.0,
                    **estimates,
                })
                curve_frames.append(pd.DataFrame({
                    "replicate_id": replicate, "validation_n": n_validation,
                    "arm": arm, "H0": h0_grid, "log_likelihood": log_likelihood,
                }))
            for arm, weights in {
                "A_catalogue_only": arms.arm_a_catalogue_only,
                "B_parent_oracle": arms.arm_b_parent_oracle,
                "C_known_probability_ipw": arms.arm_c_known_probability,
                "D_map_plugin_ipw": arms.arm_d_map_plugin,
            }.items():
                ledger.append({
                    "replicate_id": replicate, "validation_n": n_validation,
                    "arm": arm, "status": "SUCCESS",
                    **weight_concentration_diagnostics(weights),
                })
            node_diagnostics = [
                weight_concentration_diagnostics(row)
                for row in arms.arm_e_node_weights
            ]
            ledger.append({
                "replicate_id": replicate, "validation_n": n_validation,
                "arm": "E_uncertainty_integrated_ipw", "status": "SUCCESS",
                "quadrature_ess": arms.inclusion_posterior.quadrature_ess,
                "minimum_node_kish_ess": min(x["kish_ess"] for x in node_diagnostics),
                "maximum_node_top_5_percent_mass": max(x["top_5_percent_mass"] for x in node_diagnostics),
                "quadrature_nodes": len(node_diagnostics),
            })
        except Exception as exc:
            ledger.append({
                "replicate_id": replicate, "validation_n": n_validation,
                "arm": "ALL", "status": "FAILED",
                "failure_type": type(exc).__name__, "failure_message": str(exc),
            })

    catalogue.to_csv(observable_dir / "catalogue_observed.csv", index=False)
    population["gw_events_observed"].to_csv(observable_dir / "gw_events_observed.csv", index=False)
    for n in sizes:
        validation["validation_observed"].iloc[:n].to_csv(
            observable_dir / f"validation_n{n}.csv", index=False
        )
    parent.to_csv(truth_dir / "parent_truth.csv", index=False)
    population["gw_event_truth"].to_csv(truth_dir / "gw_event_truth.csv", index=False)
    population["gw_generation_trace_truth"].to_csv(truth_dir / "gw_generation_trace_truth.csv", index=False)
    validation_frame.to_csv(truth_dir / "validation_frame_truth.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(results_dir / "posterior_summaries.csv", index=False)
    pd.DataFrame(ledger).to_csv(results_dir / "failure_and_diagnostic_ledger.csv", index=False)
    if curve_frames:
        pd.concat(curve_frames, ignore_index=True).to_csv(results_dir / "likelihood_curves.csv", index=False)

    config = {
        "experiment": "matched_five_arm_single_replicate_smoke",
        "interpretation": "engineering smoke only; not confirmatory, observational, or an estimator-superiority claim",
        "seed": seed, "validation_seed": seed + 1,
        "n_parent": 4000, "n_detected_events": 30,
        "n_catalogue_observed": int(len(catalogue)),
        "validation_frame_size": 4000, "validation_sizes": list(sizes),
        "h0_true": 70.0, "h0_grid": [55.0, 85.0, 301],
        "sigma_mpc": 140.0, "detection_threshold_mpc": 250.0,
        "quadrature_order": 7,
        "truth_used_for_non_oracle_inference": False,
        "summary_rows": len(summary_rows),
        "failed_validation_sizes": sum(row["status"] == "FAILED" for row in ledger),
    }
    (output_dir / "config.json").write_text(json.dumps(config, indent=2) + "\n", encoding="utf-8")
    files = sorted(p for p in output_dir.rglob("*") if p.is_file() and p.name != "manifest.json")
    manifest = {
        "status": "ARTIFACTS_WRITTEN",
        "files": {
            str(p.relative_to(output_dir)): {"bytes": p.stat().st_size, "sha256": sha256(p)}
            for p in files
        },
    }
    (output_dir / "manifest.json").write_text(json.dumps(manifest, indent=2) + "\n", encoding="utf-8")
    return {
        "status": "COMPLETED_WITH_FAILURES" if config["failed_validation_sizes"] else "COMPLETED",
        "summary_rows": len(summary_rows),
        "failed_validation_sizes": config["failed_validation_sizes"],
        "output_dir": str(output_dir),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261113)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run_smoke(args.output_dir, args.seed), indent=2))


if __name__ == "__main__":
    main()
