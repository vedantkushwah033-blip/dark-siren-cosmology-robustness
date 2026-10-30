from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from dark_siren_robustness.finite_parent_sampling import run_finite_parent_benchmark


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a finite-parent SRSWOR validation benchmark.")
    parser.add_argument("--scenarios", type=int, default=10)
    parser.add_argument("--draws", type=int, default=20)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-sizes", default="100,300")
    parser.add_argument("--seed", type=int, default=20261032)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    sizes = tuple(int(x) for x in args.validation_sizes.split(",") if x.strip())
    frame = run_finite_parent_benchmark(
        n_scenarios=args.scenarios, n_validation_draws=args.draws,
        n_parent=args.n_parent, n_detected=args.n_detected,
        validation_sizes=sizes, seed=args.seed,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "records.csv", index=False)
    successful = frame.loc[frame["status"].eq("success")]
    summary = {
        "schema": "finite_parent_srswor_smoke_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "n_attempted": int(len(frame)),
        "n_successful": int(len(successful)),
        "n_failed": int((frame["status"] != "success").sum()),
        "failure_fraction": float((frame["status"] != "success").mean()),
        "validation_sizes": list(sizes),
        "scenarios": int(args.scenarios),
        "draws_per_scenario_and_size": int(args.draws),
        "seed": int(args.seed),
        "h0_points": 81,
        "z_points": 401,
        "h0_grid_bounds": [50.0, 90.0],
        "redshift_grid_bounds": [0.01, 0.40],
        "metrics_by_validation_size": {},
        "limitations": [
            "Synthetic conditional design only; not real LVK data or a new H0 measurement.",
            "Catalogue and GW distances are fixed within each parent scenario.",
            "Validation draws are nested within parent scenarios and are not independent top-level units; scenario-balanced summaries are reported alongside pooled descriptive summaries.",
            "The finite-parent target is the posterior mean obtained by fitting the declared estimated-inclusion pipeline to all N actual paired parent records; it is not the generating H0 truth.",
            "The finite-parent target is scenario-specific and fixed across validation draws; its error is a separate estimand from error relative to generating H0.",
            "This runner estimates SRSWOR repeated-sampling behavior; paired-bootstrap versus design-aware bootstrap calibration is a subsequent stage.",
        ],
    }
    for size in sizes:
        group = successful.loc[successful["validation_size"].eq(size)]
        if len(group):
            summary["metrics_by_validation_size"][str(size)] = {
                "posterior_mean_bias": float((group["posterior_mean"]-group["h0_true"]).mean()),
                "posterior_mean_rmse": float((group["absolute_error"]**2).mean()**0.5),
                "finite_parent_target_bias": float(group["error_vs_finite_parent_target"].mean()),
                "finite_parent_target_rmse": float(
                    (group["error_vs_finite_parent_target"]**2).mean()**0.5
                ),
                "coverage_68": float(group["covered_68"].mean()),
                "mean_ci68_width": float(group["ci68_width"].mean()),
                "mean_density_ise": float(group["density_ise"].mean()),
                "mean_inclusion_probability_mae": float(group["inclusion_probability_mae"].mean()),
                "mean_weight_ess": float(group["weight_ess"].mean()),
                "mean_max_inverse_weight": float(group["max_inverse_weight"].mean()),
                "boundary_limited_fraction": float(group["fit_boundary_limited"].mean()),
                # Repeated validation draws are nested within parent scenarios.
                # Report scenario-balanced summaries so scenarios, not draws,
                # are the top-level units for uncertainty interpretation.
                "n_scenarios_with_successes": int(group["scenario"].nunique()),
                "scenario_balanced_bias_mean": float(
                    group.assign(error=group["posterior_mean"]-group["h0_true"])
                    .groupby("scenario")["error"].mean().mean()
                ),
                "scenario_balanced_rmse_mean": float(
                    group.assign(squared_error=(group["posterior_mean"]-group["h0_true"])**2)
                    .groupby("scenario")["squared_error"].mean().pow(0.5).mean()
                ),
                "scenario_balanced_coverage_mean": float(
                    group.groupby("scenario")["covered_68"].mean().mean()
                ),
                "scenario_coverage_sd": float(
                    group.groupby("scenario")["covered_68"].mean().std(ddof=1)
                ) if group["scenario"].nunique() > 1 else None,
            }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": summary["status"], "n_attempted": summary["n_attempted"],
        "n_successful": summary["n_successful"], "n_failed": summary["n_failed"],
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
