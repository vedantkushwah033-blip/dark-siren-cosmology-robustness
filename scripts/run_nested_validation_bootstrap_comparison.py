"""Run the nested outer-SRSWOR / inner-bootstrap comparison.

Exploratory synthetic benchmark only. This runner does not establish validity
of any bootstrap interval for the nonlinear posterior-mean functional.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.bootstrap_calibration import summarize_bootstrap_calibration
from dark_siren_robustness.estimated_inclusion import estimate_host_density_from_validation
from dark_siren_robustness.finite_parent_sampling import nested_srswor_sample_indices
from dark_siren_robustness.generative_selection import generate_synthetic_population
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.posterior_validation_bootstrap import bootstrap_validation_h0_posterior
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood


def _seed(root_seed: int, *keys: int) -> int:
    state = np.random.SeedSequence([int(root_seed), *map(int, keys)])
    return int(state.generate_state(1, dtype=np.uint32)[0])


def _sha(values: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(values, dtype=np.float64))
    digest = hashlib.sha256()
    digest.update(np.asarray(arr.shape, dtype=np.int64).tobytes())
    digest.update(arr.tobytes())
    return digest.hexdigest()


def _posterior(z_grid, h0_grid, observed_z, validation_z, validation_y,
               distances, rate, sigma_mpc, threshold_mpc):
    fitted = estimate_host_density_from_validation(
        z_grid, observed_z, validation_z, validation_y, merger_weight=rate
    )
    likelihood = selection_aware_toy_likelihood(
        z_grid, fitted["density"], h0_grid, distances,
        sigma=sigma_mpc, detection_threshold=threshold_mpc,
    )
    return fitted, posterior_summary(h0_grid, likelihood["selected_log_likelihood"])


def run_comparison(*, scenarios: int, draws: int, n_parent: int,
                   n_detected: int, validation_sizes: tuple[int, ...],
                   bootstrap_replicates: int, seed: int,
                   h0_true: float = 70.0, sigma_mpc: float = 140.0,
                   threshold_mpc: float = 250.0) -> tuple[pd.DataFrame, pd.DataFrame]:
    integer_parameters = {
        "scenarios": scenarios, "draws": draws, "n_parent": n_parent,
        "n_detected": n_detected, "bootstrap_replicates": bootstrap_replicates,
        "seed": seed,
    }
    for name, value in integer_parameters.items():
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError(f"{name} must be an integer")
    if scenarios < 1 or draws < 2 or n_parent < 100 or n_detected < 1:
        raise ValueError("invalid scenario, draw, parent, or detected-event count")
    if bootstrap_replicates < 2 or seed < 0:
        raise ValueError("bootstrap_replicates must be >=2 and seed nonnegative")
    if not np.isfinite(h0_true) or h0_true <= 0:
        raise ValueError("h0_true must be finite and positive")
    if not np.isfinite(sigma_mpc) or sigma_mpc <= 0:
        raise ValueError("sigma_mpc must be finite and positive")
    if not np.isfinite(threshold_mpc) or threshold_mpc <= 0:
        raise ValueError("threshold_mpc must be finite and positive")
    raw_sizes = tuple(validation_sizes)
    for value in raw_sizes:
        if isinstance(value, (bool, np.bool_)) or not isinstance(value, (int, np.integer)):
            raise TypeError("validation sizes must be integers")
    sizes = tuple(sorted(set(int(n) for n in raw_sizes)))
    if not sizes or min(sizes) < 10 or max(sizes) > n_parent:
        raise ValueError("validation sizes must be unique integers in [10, n_parent]")
    if len(sizes) != len(raw_sizes):
        raise ValueError("validation sizes must not contain duplicates")

    z_grid = np.linspace(0.01, 0.40, 401)
    h0_grid = np.linspace(50.0, 90.0, 81)
    rate = 1.0 + 1.5 * (z_grid-z_grid[0])/(z_grid[-1]-z_grid[0])
    midpoint, half_range = 0.205, 0.195
    rows: list[dict] = []
    attempt_log: list[dict] = []

    for scenario in range(scenarios):
        root = np.random.SeedSequence([seed, scenario])
        sim_ss, thin_ss, outer_ss = root.spawn(3)
        sim_seed = int(sim_ss.generate_state(1, dtype=np.uint32)[0])
        thin_seed = int(thin_ss.generate_state(1, dtype=np.uint32)[0])
        outer_seed = int(outer_ss.generate_state(1, dtype=np.uint32)[0])
        sim = generate_synthetic_population(
            seed=sim_seed, n_parent=n_parent, n_detected=n_detected,
            h0_true=h0_true, sigma_mpc=sigma_mpc,
            detection_threshold_mpc=threshold_mpc,
            z_min=float(z_grid[0]), z_max=float(z_grid[-1]),
        )
        parent_z = sim["parent_truth"]["redshift"].to_numpy(dtype=float)
        distances = sim["gw_events_observed"]["distance_mpc"].to_numpy(dtype=float)
        x = (parent_z-midpoint)/half_range
        p_true = expit(0.25-1.2*x)
        retained = np.random.default_rng(thin_seed).random(n_parent) < p_true
        observed_z = parent_z[retained]
        parent_hash, catalogue_hash, gw_hash = _sha(parent_z), _sha(observed_z), _sha(distances)

        if observed_z.size < 5 or np.unique(observed_z).size < 3 or np.unique(retained).size != 2:
            raise RuntimeError(
                f"scenario {scenario} cannot define the prespecified finite-parent target: "
                f"catalogue_size={observed_z.size}, unique_catalogue_redshifts={np.unique(observed_z).size}, "
                f"retention_classes={np.unique(retained).size}"
            )
        _, full = _posterior(
            z_grid, h0_grid, observed_z, parent_z, retained.astype(int),
            distances, rate, sigma_mpc, threshold_mpc,
        )
        target = float(full["posterior_mean"])
        if not np.isfinite(target):
            raise ArithmeticError(f"scenario {scenario} produced a nonfinite finite-parent target")

        # One independent permutation per outer draw; its prefixes create
        # paired nested SRSWOR samples across the declared validation sizes.
        outer_rng = np.random.default_rng(outer_seed)
        for draw in range(draws):
            samples = nested_srswor_sample_indices(outer_rng, n_parent, sizes)
            for size in sizes:
                idx = samples[size]
                vz, vy = parent_z[idx], retained[idx].astype(int)
                outer_status, outer_mean = "failed", np.nan
                outer_failure = ""
                if not np.isfinite(target) or np.unique(vy).size != 2:
                    outer_failure = "finite_parent_target_unavailable_or_one_class_outer_sample"
                else:
                    try:
                        _, outer_post = _posterior(
                            z_grid, h0_grid, observed_z, vz, vy, distances,
                            rate, sigma_mpc, threshold_mpc,
                        )
                        outer_mean = float(outer_post["posterior_mean"])
                        outer_status = "success"
                    except Exception as exc:
                        outer_failure = f"{type(exc).__name__}: {exc}"

                method_results: dict[str, dict] = {}
                if outer_status == "success":
                    for method in ("paired_bootstrap", "design_aware_srswor"):
                        method_seed = _seed(seed, scenario, draw, size,
                                            1 if method == "paired_bootstrap" else 2)
                        try:
                            result = bootstrap_validation_h0_posterior(
                                z_grid, h0_grid, observed_z, vz, vy, distances,
                                sigma_mpc=sigma_mpc,
                                detection_threshold_mpc=threshold_mpc,
                                parent_population_size=n_parent,
                                method=method, n_bootstrap=bootstrap_replicates,
                                seed=method_seed, merger_weight=rate,
                            )
                            ok = int(result["n_successful"])
                            if ok < 2:
                                raise RuntimeError(
                                    f"fewer than two successful bootstrap replicates ({ok}/{bootstrap_replicates})"
                                )
                            lo, hi = map(float, result["posterior_mean_percentile_interval_95"])
                            method_results[method] = {
                                "bootstrap_status": "success",
                                "bootstrap_se": float(result["raw_bootstrap_standard_error"]),
                                "interval_low": lo, "interval_high": hi,
                                "n_bootstrap_successful": ok,
                                "n_bootstrap_attempted": bootstrap_replicates,
                                "bootstrap_seed": method_seed,
                                "failure_records": result["failure_records"],
                                "wrapper_failure_class": type(exc).__name__,
                                "wrapper_failure_message": str(exc),
                                "fpc_factor": float(result["fpc_factor"]),
                                "input_fingerprints": result["input_fingerprints"],
                                "success_mask": np.asarray(result["success_mask"], dtype=bool).tolist(),
                                "attempt_posterior_mean": np.asarray(result["attempt_posterior_mean"], dtype=float).tolist(),
                                "reference_posterior_mean_by_attempt": np.asarray(result["reference_posterior_mean_by_attempt"], dtype=float).tolist(),
                                "posterior_median_draws_by_attempt": np.asarray(result["posterior_median_draws_by_attempt"], dtype=float).tolist(),
                                "ci68_low_draws_by_attempt": np.asarray(result["ci68_low_draws_by_attempt"], dtype=float).tolist(),
                                "ci68_high_draws_by_attempt": np.asarray(result["ci68_high_draws_by_attempt"], dtype=float).tolist(),
                                "inverse_weight_ess_by_attempt": np.asarray(result["inverse_weight_ess_by_attempt"], dtype=float).tolist(),
                                "inverse_weight_maximum_by_attempt": np.asarray(result["inverse_weight_maximum_by_attempt"], dtype=float).tolist(),
                                "coefficient_boundary_limited_by_attempt": np.asarray(result["coefficient_boundary_limited_by_attempt"], dtype=bool).tolist(),
                                "posterior_density_draws": np.asarray(result["posterior_density_draws"], dtype=float).tolist(),
                                "fixed_full_validation_reference_mean": float(result["fixed_full_validation_reference_mean"]),
                                "sample_indices": [np.asarray(v, dtype=np.int64).tolist() for v in result["sample_indices"]],
                                "pseudo_population_indices": [np.asarray(v, dtype=np.int64).tolist() for v in result["pseudo_population_indices"]],
                                "selected_pseudo_population_positions": [np.asarray(v, dtype=np.int64).tolist() for v in result["selected_pseudo_population_positions"]],
                            }
                        except Exception as exc:
                            method_results[method] = {
                                "bootstrap_status": "failed", "bootstrap_se": np.nan,
                                "interval_low": np.nan, "interval_high": np.nan,
                                "n_bootstrap_successful": 0,
                                "n_bootstrap_attempted": bootstrap_replicates,
                                "bootstrap_seed": method_seed,
                                "failure_records": [{"reason": f"{type(exc).__name__}: {exc}"}],
                                "fpc_factor": float(np.sqrt(max(0.0, 1.0-size/n_parent))),
                            }
                    paired = method_results["paired_bootstrap"]
                    fpc_se = (paired["bootstrap_se"]*paired["fpc_factor"]
                              if paired["bootstrap_status"] == "success" else np.nan)
                    method_results["paired_bootstrap_fpc"] = {
                        "bootstrap_status": paired["bootstrap_status"],
                        "bootstrap_se": fpc_se,
                        "interval_low": outer_mean-1.96*fpc_se if np.isfinite(fpc_se) else np.nan,
                        "interval_high": outer_mean+1.96*fpc_se if np.isfinite(fpc_se) else np.nan,
                        "n_bootstrap_successful": paired["n_bootstrap_successful"],
                        "n_bootstrap_attempted": bootstrap_replicates,
                        "bootstrap_seed": paired["bootstrap_seed"],
                        "failure_records": paired["failure_records"],
                        "fpc_factor": paired["fpc_factor"],
                    }
                else:
                    for method in ("paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"):
                        method_results[method] = {
                            "bootstrap_status": "not_run", "bootstrap_se": np.nan,
                            "interval_low": np.nan, "interval_high": np.nan,
                            "n_bootstrap_successful": 0,
                            "n_bootstrap_attempted": 0, "bootstrap_seed": None,
                            "failure_records": [], "fpc_factor": float(np.sqrt(max(0.0, 1.0-size/n_parent))),
                        }

                for method, result in method_results.items():
                    rows.append({
                        "scenario": scenario, "validation_size": size, "outer_draw": draw,
                        "method": method, "outer_status": outer_status,
                        "bootstrap_status": result["bootstrap_status"],
                        "outer_posterior_mean": outer_mean,
                        "bootstrap_se": result["bootstrap_se"],
                        "interval_low": result["interval_low"],
                        "interval_high": result["interval_high"],
                        "finite_parent_target": target,
                    })
                    attempt_log.append({
                        "scenario": scenario, "validation_size": size, "outer_draw": draw,
                        "method": method, "outer_seed": outer_seed,
                        "simulation_seed": sim_seed, "thinning_seed": thin_seed,
                        "bootstrap_seed": result["bootstrap_seed"],
                        "outer_failure": outer_failure,
                        "bootstrap_failure_records": result["failure_records"],
                        "bootstrap_wrapper_failure_class": result.get("wrapper_failure_class"),
                        "bootstrap_wrapper_failure_message": result.get("wrapper_failure_message"),
                        "bootstrap_input_fingerprints": result.get("input_fingerprints"),
                        "bootstrap_success_mask": result.get("success_mask"),
                        "bootstrap_attempt_posterior_mean": result.get("attempt_posterior_mean"),
                        "bootstrap_reference_posterior_mean_by_attempt": result.get("reference_posterior_mean_by_attempt"),
                        "bootstrap_posterior_median_draws_by_attempt": result.get("posterior_median_draws_by_attempt"),
                        "bootstrap_ci68_low_draws_by_attempt": result.get("ci68_low_draws_by_attempt"),
                        "bootstrap_ci68_high_draws_by_attempt": result.get("ci68_high_draws_by_attempt"),
                        "bootstrap_inverse_weight_ess_by_attempt": result.get("inverse_weight_ess_by_attempt"),
                        "bootstrap_inverse_weight_maximum_by_attempt": result.get("inverse_weight_maximum_by_attempt"),
                        "bootstrap_coefficient_boundary_limited_by_attempt": result.get("coefficient_boundary_limited_by_attempt"),
                        "bootstrap_posterior_density_draws": result.get("posterior_density_draws"),
                        "bootstrap_fixed_full_validation_reference_mean": result.get("fixed_full_validation_reference_mean"),
                        "bootstrap_sample_indices": result.get("sample_indices"),
                        "bootstrap_pseudo_population_indices": result.get("pseudo_population_indices"),
                        "bootstrap_selected_pseudo_population_positions": result.get("selected_pseudo_population_positions"),
                        "n_bootstrap_attempted": result["n_bootstrap_attempted"],
                        "n_bootstrap_successful": result["n_bootstrap_successful"],
                        "fpc_factor": result["fpc_factor"],
                        "parent_frame_sha256": parent_hash,
                        "catalogue_sha256": catalogue_hash,
                        "gw_distance_sha256": gw_hash,
                        "validation_indices": idx.tolist(),
                    })
    return pd.DataFrame(rows), pd.DataFrame(attempt_log)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--scenarios", type=int, default=2)
    parser.add_argument("--draws", type=int, default=5)
    parser.add_argument("--n-parent", type=int, default=1000)
    parser.add_argument("--n-detected", type=int, default=10)
    parser.add_argument("--validation-sizes", default="100,300")
    parser.add_argument("--bootstrap-replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sizes = tuple(int(value) for value in args.validation_sizes.split(","))
    records, attempts = run_comparison(
        scenarios=args.scenarios, draws=args.draws, n_parent=args.n_parent,
        n_detected=args.n_detected, validation_sizes=sizes,
        bootstrap_replicates=args.bootstrap_replicates, seed=args.seed,
    )
    scenario_summary, aggregate = summarize_bootstrap_calibration(records)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records.to_csv(args.output_dir/"calibration_records.csv", index=False)
    attempts.to_json(args.output_dir/"attempt_audit.jsonl", orient="records", lines=True)
    scenario_summary.to_csv(args.output_dir/"scenario_summary.csv", index=False)
    aggregate.to_csv(args.output_dir/"scenario_balanced_summary.csv", index=False)
    manifest = {
        "schema": "nested_validation_bootstrap_comparison_smoke_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "seed": args.seed, "scenarios": args.scenarios, "outer_draws": args.draws,
        "parent_size": args.n_parent, "detected_event_count": args.n_detected,
        "validation_sizes": list(sizes),
        "model_parameters": {
            "h0_true_km_s_mpc": 70.0,
            "distance_noise_sigma_mpc": 140.0,
            "detection_threshold_mpc": 250.0,
            "redshift_grid": {"start": 0.01, "stop": 0.40, "points": 401},
            "h0_grid": {"start": 50.0, "stop": 90.0, "points": 81},
            "merger_weight": "linear from 1.0 at z=0.01 to 2.5 at z=0.40",
            "catalogue_retention_probability": "expit(0.25 - 1.2*x), x=(z-0.205)/0.195",
        },
        "bootstrap_replicates_per_method_per_outer_draw": args.bootstrap_replicates,
        "methods": ["paired_bootstrap", "paired_bootstrap_fpc", "design_aware_srswor"],
        "outer_design": "nested_simple_random_sampling_without_replacement",
        "target": "posterior_mean_from_same_pipeline_fitted_to_all_actual_parent_records",
        "limitations": [
            "Synthetic low-redshift Euclidean toy model; not real LVK observations or a new H0 measurement.",
            "Observed catalogue and GW distances are fixed within each parent scenario.",
            "The pseudo-population arm is exploratory and not studentized finite-population-mean inference.",
            "FPC-normal intervals are diagnostic comparators, not theoretically validated intervals for this nonlinear estimator.",
            "Outer draws are nested within parent scenarios; scenarios are the top-level independent units.",
            "Bootstrap intervals describe resampling variation of posterior means, not conditional H0 posterior credible intervals.",
        ],
    }
    (args.output_dir/"manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": manifest["status"], "records": len(records),
        "outer_attempts": args.scenarios*args.draws*len(sizes),
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
