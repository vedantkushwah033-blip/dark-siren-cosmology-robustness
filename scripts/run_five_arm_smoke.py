"""Run a single matched five-arm synthetic dark-siren smoke replicate.

This is an orchestration/artifact smoke tool, not a confirmatory simulation.
It writes observable inputs separately from simulator truth and hashes every
artifact other than the manifest itself.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.five_arm_comparison import (
    build_five_arm_weights,
    weight_concentration_diagnostics,
)
from dark_siren_robustness.five_arm_likelihood import matched_five_arm_selected_log_likelihood
from dark_siren_robustness.generative_selection import (
    generate_independent_validation_sample,
    generate_synthetic_population,
)
from dark_siren_robustness.h0_posterior import normalize_h0_posterior


def _parse_sizes(raw: str) -> list[int]:
    try:
        values = [int(item.strip()) for item in raw.split(",") if item.strip()]
    except ValueError as exc:
        raise argparse.ArgumentTypeError("validation sizes must be comma-separated integers") from exc
    if not values or any(value < 5 for value in values) or sorted(set(values)) != values:
        raise argparse.ArgumentTypeError("validation sizes must be unique, increasing integers >= 5")
    return values


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _posterior_summary(grid: np.ndarray, log_likelihood: np.ndarray) -> dict[str, float]:
    posterior = normalize_h0_posterior(grid, log_likelihood, prior_min=50.0, prior_max=90.0)
    return {
        "median": posterior.quantile(0.50),
        "q025": posterior.quantile(0.025),
        "q05": posterior.quantile(0.05),
        "q16": posterior.quantile(0.16),
        "q84": posterior.quantile(0.84),
        "q95": posterior.quantile(0.95),
        "q975": posterior.quantile(0.975),
    }


def run(args: argparse.Namespace) -> dict:
    output = Path(args.output_dir).resolve()
    output.mkdir(parents=True, exist_ok=True)
    observable = output / "observable"
    truth_dir = output / "truth"
    results_dir = output / "results"
    for directory in (observable, truth_dir, results_dir):
        directory.mkdir(parents=True, exist_ok=True)

    sizes = args.validation_sizes
    if sizes[-1] > args.validation_frame:
        raise ValueError("largest validation size cannot exceed validation frame size")

    population = generate_synthetic_population(
        seed=args.seed,
        n_parent=args.n_parent,
        n_detected=args.n_detected,
        h0_true=args.h0_true,
        sigma_mpc=args.sigma_mpc,
        detection_threshold_mpc=args.detection_threshold_mpc,
        catalogue_floor=args.catalogue_floor,
    )
    validation_results = [
        generate_independent_validation_sample(
            seed=args.validation_seed,
            n_validation_frame=args.validation_frame,
            n_validation=size,
            catalogue_floor=args.catalogue_floor,
        )
        for size in sizes
    ]
    reference_truth = validation_results[-1]["truth"]
    truth_frame_columns = ["frame_id", "redshift", "inclusion_probability", "retained"]
    for smaller, larger in zip(validation_results, validation_results[1:]):
        if not smaller["observable"]["validation_id"].tolist() == larger["observable"]["validation_id"].iloc[:len(smaller["observable"])].tolist():
            raise AssertionError("validation samples are not nested prefixes")
        if not reference_truth[truth_frame_columns].equals(larger["truth"][truth_frame_columns]):
            raise AssertionError("validation parent frame changed across nested sample sizes")

    catalogue = population["catalogue_observed"]
    parent_truth = population["parent_truth"]
    events_frame = population["gw_events_observed"]
    parent_z = parent_truth["redshift"].to_numpy(float)
    catalogue_z = catalogue["redshift"].to_numpy(float)
    parent_base = 1.0 + 1.5 * (parent_z - 0.01) / 0.39
    catalogue_base = 1.0 + 1.5 * (catalogue_z - 0.01) / 0.39
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": float(args.sigma_mpc)}
        for row in events_frame.itertuples(index=False)
    ]
    h0_grid = np.linspace(50.0, 90.0, args.h0_points)
    summary_rows: list[dict] = []
    likelihood_arrays: dict[str, np.ndarray] = {}
    diagnostics: dict[str, object] = {}

    for validation_result, n_validation in zip(validation_results, sizes):
        validation = validation_result["observable"]
        arms = build_five_arm_weights(
            parent_z,
            catalogue_z,
            validation["redshift"].to_numpy(float),
            validation["retained"].to_numpy(np.int8),
            catalogue_floor=args.catalogue_floor,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
            quadrature_order=args.quadrature_order,
        )
        likelihood = matched_five_arm_selected_log_likelihood(
            h0_grid,
            arms,
            events,
            detection_threshold_mpc=args.detection_threshold_mpc,
        )
        arm_likelihoods = {
            "A_catalogue_only": likelihood.arm_a_log_likelihood,
            "B_parent_oracle": likelihood.arm_b_log_likelihood,
            "C_known_probability_ipw": likelihood.arm_c_log_likelihood,
            "D_map_plugin_ipw": likelihood.arm_d_log_likelihood,
            "E_uncertainty_integrated_ipw": likelihood.arm_e_marginal_log_likelihood,
        }
        for arm_name, log_likelihood in arm_likelihoods.items():
            likelihood_arrays[f"{arm_name}_n{n_validation}"] = np.asarray(log_likelihood)
            estimates = _posterior_summary(h0_grid, np.asarray(log_likelihood))
            summary_rows.append({
                "validation_size": n_validation,
                "arm": arm_name,
                "status": "SUCCESS",
                "h0_true": float(args.h0_true),
                "posterior_median": estimates["median"],
                "median_error": estimates["median"] - float(args.h0_true),
                **{key: value for key, value in estimates.items() if key != "median"},
            })
        likelihood_arrays[f"posterior_nodes_n{n_validation}"] = np.asarray(arms.inclusion_posterior.nodes)
        likelihood_arrays[f"posterior_mode_n{n_validation}"] = np.asarray(arms.inclusion_posterior.mode)
        likelihood_arrays[f"posterior_weights_n{n_validation}"] = np.asarray(arms.inclusion_posterior.weights)
        likelihood_arrays[f"E_node_log_likelihood_n{n_validation}"] = np.asarray(
            likelihood.arm_e_node_log_likelihood
        )
        likelihood_arrays[f"E_node_log_efficiency_n{n_validation}"] = np.asarray(
            likelihood.arm_e_node_log_detection_efficiency
        )

        arm_diag: dict[str, object] = {}
        for arm_name, weights in (
            ("A_catalogue_only", arms.arm_a_catalogue_only),
            ("B_parent_oracle", arms.arm_b_parent_oracle),
            ("C_known_probability_ipw", arms.arm_c_known_probability),
            ("D_map_plugin_ipw", arms.arm_d_map_plugin),
        ):
            arm_diag[arm_name] = weight_concentration_diagnostics(weights)
        node_diag = [
            weight_concentration_diagnostics(weights)
            for weights in arms.arm_e_node_weights
        ]
        arm_diag["E_uncertainty_integrated_ipw"] = {
            "quadrature_ess": arms.inclusion_posterior.quadrature_ess,
            "minimum_node_kish_ess": min(item["kish_ess"] for item in node_diag),
            "median_node_kish_ess": float(np.median([item["kish_ess"] for item in node_diag])),
            "maximum_node_top_5_percent_mass": max(item["top_5_percent_mass"] for item in node_diag),
            "number_of_nodes": len(node_diag),
        }
        diagnostics[str(n_validation)] = arm_diag
        validation.to_csv(observable / f"validation_n{n_validation}.csv", index=False)

    catalogue.to_csv(observable / "catalogue_observed.csv", index=False)
    events_frame.to_csv(observable / "gw_events_observed.csv", index=False)
    parent_truth.to_csv(truth_dir / "parent_truth.csv", index=False)
    population["gw_event_truth"].to_csv(truth_dir / "gw_event_truth.csv", index=False)
    population["gw_generation_trace_truth"].to_csv(truth_dir / "gw_generation_trace_truth.csv", index=False)
    reference_truth.to_csv(truth_dir / "validation_frame_truth.csv", index=False)
    pd.DataFrame(summary_rows).to_csv(results_dir / "posterior_summaries.csv", index=False)
    np.savez_compressed(
        results_dir / "log_likelihoods.npz",
        H0_grid=h0_grid,
        **likelihood_arrays,
    )
    (results_dir / "weight_diagnostics.json").write_text(
        json.dumps(diagnostics, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    config = {
        "experiment": "matched_five_arm_single_replicate_smoke",
        "status": "engineering_smoke_not_confirmatory_not_observational",
        "seed": args.seed,
        "validation_seed": args.validation_seed,
        "n_parent": args.n_parent,
        "n_catalogue_observed": int(len(catalogue)),
        "n_detected": args.n_detected,
        "validation_frame": args.validation_frame,
        "validation_sizes": sizes,
        "catalogue_floor": args.catalogue_floor,
        "h0_true": args.h0_true,
        "h0_grid_points": args.h0_points,
        "quadrature_order": args.quadrature_order,
        "sigma_mpc": args.sigma_mpc,
        "detection_threshold_mpc": args.detection_threshold_mpc,
        "arms": ["A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw", "D_map_plugin_ipw", "E_uncertainty_integrated_ipw"],
        "claim_limit": "One engineering replicate; no calibration, estimator superiority, or real-data claim.",
    }
    (output / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    files = sorted(path for path in output.rglob("*") if path.is_file() and path.name != "manifest.json")
    manifest = {
        "status": "PASS",
        "file_count": len(files),
        "files": {
            str(path.relative_to(output)): {
                "bytes": path.stat().st_size,
                "sha256": _sha256(path),
            }
            for path in files
        },
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"output_dir": str(output), "files": len(files), "summary_rows": len(summary_rows), "status": "PASS"}


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261113)
    parser.add_argument("--validation-seed", type=int, default=20261114)
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-detected", type=int, default=4)
    parser.add_argument("--validation-frame", type=int, default=1500)
    parser.add_argument("--validation-sizes", type=_parse_sizes, default=_parse_sizes("100,300,1000"))
    parser.add_argument("--catalogue-floor", type=float, default=0.10)
    parser.add_argument("--h0-true", type=float, default=70.0)
    parser.add_argument("--h0-points", type=int, default=41)
    parser.add_argument("--quadrature-order", type=int, default=5)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--detection-threshold-mpc", type=float, default=250.0)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
