"""Paired multi-replicate pilot for quadratic-curvature H0 calibration.

This orchestration reuses the single-cell smoke runner, preserving each cell's
observable/truth-separated raw artifact. It is a synthetic pilot, not an
observational H0 analysis or an SBC claim.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from scripts.run_quadratic_curvature_h0_smoke import run as run_single


ARMS = {
    "degree_1": "density_d1_o7",
    "degree_2": "density_d2_o7",
    "catalogue_only": "density_catalogue_only",
    "oracle_true_inclusion": "density_oracle_true_inclusion",
}


def _seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([master, replicate, stream]).generate_state(
        1, dtype=np.uint32
    )[0])


def _posterior_metrics(grid: np.ndarray, density: np.ndarray, truth: float) -> dict:
    if (grid.ndim != 1 or density.shape != grid.shape or grid.size < 2
            or not np.all(np.isfinite(grid)) or not np.all(np.isfinite(density))
            or np.any(np.diff(grid) <= 0) or np.any(density < 0)):
        raise ValueError("invalid H0 posterior curve")
    area = float(np.trapezoid(density, grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("H0 posterior has invalid integrated mass")
    density = density / area
    increments = 0.5 * (density[:-1] + density[1:]) * np.diff(grid)
    cdf = np.r_[0.0, np.cumsum(increments)]
    cdf /= cdf[-1]
    if not grid[0] < truth < grid[-1]:
        raise ValueError("generating H0 must lie in the interior of posterior support")
    quantiles = {
        str(q): float(np.interp(q, cdf, grid))
        for q in (0.025, 0.05, 0.16, 0.50, 0.84, 0.95, 0.975)
    }
    return {
        "area_before_normalization": area,
        "pit": float(np.interp(truth, grid, cdf)),
        "h0_q025": quantiles["0.025"],
        "h0_q05": quantiles["0.05"],
        "h0_q16": quantiles["0.16"],
        "h0_median": quantiles["0.5"],
        "h0_q84": quantiles["0.84"],
        "h0_q95": quantiles["0.95"],
        "h0_q975": quantiles["0.975"],
        "covered_68": int(quantiles["0.16"] <= truth <= quantiles["0.84"]),
        "covered_90": int(quantiles["0.05"] <= truth <= quantiles["0.95"]),
        "covered_95": int(quantiles["0.025"] <= truth <= quantiles["0.975"]),
        "width_68": quantiles["0.84"] - quantiles["0.16"],
        "width_90": quantiles["0.95"] - quantiles["0.05"],
        "width_95": quantiles["0.975"] - quantiles["0.025"],
    }


def _wilson(successes: int, trials: int, z: float = 1.959963984540054):
    if trials <= 0:
        return [None, None]
    p = successes / trials
    d = 1 + z * z / trials
    center = (p + z * z / (2 * trials)) / d
    half = z * math.sqrt(p * (1-p) / trials + z*z/(4*trials*trials)) / d
    return [max(0.0, center-half), min(1.0, center+half)]


def _aggregate(metrics: pd.DataFrame) -> list[dict]:
    output = []
    for (gamma, n, arm), group in metrics.groupby(
        ["gamma", "n_validation", "arm"], sort=True
    ):
        valid = group[group.status == "SUCCESS"]
        n_success = len(valid)
        item = {
            "gamma": float(gamma), "n_validation": int(n), "arm": str(arm),
            "n_planned": int(len(group)), "n_success": int(n_success),
            "n_failed": int((group.status == "FAILED").sum()),
            "success_fraction": float(n_success / len(group)),
        }
        if n_success:
            errors = valid.h0_median.to_numpy(float) - valid.h0_true.to_numpy(float)
            item.update({
                "mean_posterior_median_bias": float(np.mean(errors)),
                "rmse_posterior_median": float(np.sqrt(np.mean(errors**2))),
                "mean_pit": float(valid.pit.mean()),
                "pit_mean_mc_se": float(valid.pit.std(ddof=1) / np.sqrt(n_success)) if n_success > 1 else None,
            })
            for level, col in ((68, "covered_68"), (90, "covered_90"), (95, "covered_95")):
                successes = int(valid[col].sum())
                item[f"coverage_{level}"] = successes / n_success
                item[f"coverage_{level}_wilson95"] = _wilson(successes, n_success)
                item[f"coverage_{level}_mc_se"] = (
                    math.sqrt((successes/n_success)*(1-successes/n_success)/n_success)
                )
            for level, col in ((68, "width_68"), (90, "width_90"), (95, "width_95")):
                item[f"median_width_{level}"] = float(valid[col].median())
            for col in ("quadrature_ess", "maximum_node_mass", "minimum_log_inclusion_probability"):
                if col in valid and valid[col].notna().any():
                    item[f"median_{col}"] = float(valid[col].median())
        else:
            for name in ("mean_posterior_median_bias", "rmse_posterior_median",
                         "mean_pit", "pit_mean_mc_se"):
                item[name] = None
        output.append(item)
    return output


def run_pilot(args) -> dict:
    output = Path(args.output_dir)
    raw_root = output / "raw"
    raw_root.mkdir(parents=True, exist_ok=True)
    metric_rows = []
    curve_rows = []
    convergence_rows = []
    failure_rows = []
    planned = 0

    for replicate in range(args.replicates):
        # The same stream seeds across gamma cells create paired parent,
        # validation-frame, and event random numbers where the DGP permits it.
        population_seed = _seed(args.seed, replicate, 0)
        validation_seed = _seed(args.seed, replicate, 1)
        event_seed = _seed(args.seed, replicate, 2)
        for gamma in args.gamma_values:
            cell_base = raw_root / f"gamma_{gamma:+.2f}" / f"replicate_{replicate:04d}"
            for n_validation in args.validation_sizes:
                planned += len(ARMS)
                cell_dir = cell_base / f"n_validation_{n_validation:04d}"
                cell_args = SimpleNamespace(
                    seed=population_seed,
                    validation_seed=validation_seed,
                    event_seed=event_seed,
                    n_parent=args.n_parent,
                    n_validation_frame=args.n_validation_frame,
                    n_validation=n_validation,
                    n_events=args.n_events,
                    gamma=gamma,
                    beta=args.beta,
                    h0_true=args.h0_true,
                    sigma_mpc=args.sigma_mpc,
                    threshold=args.threshold,
                    z_min=args.z_min,
                    z_max=args.z_max,
                    h0_min=args.h0_min,
                    h0_max=args.h0_max,
                    h0_points=args.h0_points,
                    orders=(5, 7),
                    output_dir=str(cell_dir),
                )
                try:
                    summary = run_single(cell_args)
                    curves = pd.read_csv(cell_dir / "h0_posterior_curves.csv")
                    grid = curves.H0.to_numpy(float)
                    for arm, column in ARMS.items():
                        row = {
                            "replicate": replicate, "gamma": gamma,
                            "n_validation": n_validation, "arm": arm,
                            "status": "SUCCESS", "error": "",
                            "h0_true": args.h0_true,
                        }
                        if arm == "degree_1":
                            row.update(summary["results"]["degree_1"])
                            row["quadrature_ess"] = summary["results"]["degree_1"]["quadrature_ess"]
                            row["maximum_node_mass"] = summary["results"]["degree_1"]["maximum_node_mass"]
                            row["minimum_log_inclusion_probability"] = summary["results"]["degree_1"]["minimum_log_inclusion_probability"]
                        elif arm == "degree_2":
                            row.update(summary["results"]["degree_2"])
                            row["quadrature_ess"] = summary["results"]["degree_2"]["quadrature_ess"]
                            row["maximum_node_mass"] = summary["results"]["degree_2"]["maximum_node_mass"]
                            row["minimum_log_inclusion_probability"] = summary["results"]["degree_2"]["minimum_log_inclusion_probability"]
                        else:
                            row.update(summary["results"][arm])
                        posterior = _posterior_metrics(
                            grid, curves[column].to_numpy(float), args.h0_true
                        )
                        row.update(posterior)
                        metric_rows.append(row)
                        for h0, density in zip(grid, curves[column].to_numpy(float)):
                            curve_rows.append({
                                "replicate": replicate, "gamma": gamma,
                                "n_validation": n_validation, "arm": arm,
                                "H0": float(h0), "density": float(density),
                            })
                    for degree in (1, 2):
                        d = summary["results"][f"degree_{degree}"]
                        convergence_rows.append({
                            "replicate": replicate, "gamma": gamma,
                            "n_validation": n_validation, "degree": degree,
                            "from_order": 5, "to_order": 7,
                            "density_l1": d["successive_order_density_l1"],
                            "max_quantile_shift": d["successive_order_max_quantile_shift"],
                            "quadrature_ess": d["quadrature_ess"],
                            "maximum_node_mass": d["maximum_node_mass"],
                        })
                except Exception as exc:
                    message = f"{type(exc).__name__}: {exc}"
                    failure_rows.append({
                        "replicate": replicate, "gamma": gamma,
                        "n_validation": n_validation, "status": "FAILED",
                        "error": message,
                    })
                    for arm in ARMS:
                        metric_rows.append({
                            "replicate": replicate, "gamma": gamma,
                            "n_validation": n_validation, "arm": arm,
                            "status": "FAILED", "error": message,
                            "h0_true": args.h0_true,
                        })

    metrics = pd.DataFrame(metric_rows)
    curves = pd.DataFrame(curve_rows)
    metrics_path = output / "replicate_metrics.csv"
    curves_path = output / "posterior_curves.csv"
    convergence_path = output / "quadrature_convergence.csv"
    failure_path = output / "failure_ledger.csv"
    metrics.to_csv(metrics_path, index=False, float_format="%.17g")
    curves.to_csv(curves_path, index=False, float_format="%.17g")
    pd.DataFrame(convergence_rows).to_csv(convergence_path, index=False, float_format="%.17g")
    metrics[[
        "replicate", "gamma", "n_validation", "arm", "status", "error"
    ]].to_csv(failure_path, index=False)
    aggregate = _aggregate(metrics)
    (output / "aggregate_summary.json").write_text(
        json.dumps(aggregate, indent=2, sort_keys=True) + "\n"
    )
    config = {
        "experiment": "paired_quadratic_curvature_h0_calibration_pilot",
        "status": "synthetic_pilot_not_sbc_or_observational_measurement",
        "replicates": args.replicates, "master_seed": args.seed,
        "gamma_values": list(args.gamma_values),
        "validation_sizes": list(args.validation_sizes),
        "n_parent": args.n_parent,
        "n_validation_frame": args.n_validation_frame,
        "n_events": args.n_events,
        "h0_true": args.h0_true,
        "h0_grid": [args.h0_min, args.h0_max, args.h0_points],
        "quadrature_orders": [5, 7],
        "arms": ARMS,
        "paired_seed_streams": {
            "population": "SeedSequence([master_seed, replicate, 0])",
            "validation": "SeedSequence([master_seed, replicate, 1])",
            "events": "SeedSequence([master_seed, replicate, 2])",
        },
        "primary_order": 7,
        "limitations": [
            "fixed-truth repeated-sampling PIT is not prior-predictive SBC",
            "single pilot replicate count is not sufficient for strong calibration claims",
            "oracle inclusion is a simulation diagnostic only",
            "low-z Euclidean toy with common distance noise and one-sided threshold",
        ],
        "planned_metric_rows": planned,
        "actual_metric_rows": int(len(metrics)),
        "status_counts": metrics.status.value_counts().to_dict(),
    }
    config_path = output / "config.json"
    config_path.write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    manifest = {}
    for path in sorted(output.rglob("*")):
        if path.is_file() and path.name != "manifest.json":
            data = path.read_bytes()
            manifest[str(path.relative_to(output))] = {
                "bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()
            }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n"
    )
    return {
        "status": "PASS" if len(metrics) == planned and (metrics.status == "SUCCESS").all() else "INCOMPLETE",
        "planned_metric_rows": planned,
        "actual_metric_rows": int(len(metrics)),
        "success_rows": int((metrics.status == "SUCCESS").sum()),
        "failed_rows": int((metrics.status == "FAILED").sum()),
        "artifact_files": len(manifest),
        "output_dir": str(output),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=2)
    parser.add_argument("--seed", type=int, default=2026111900)
    parser.add_argument("--gamma-values", default="0,1")
    parser.add_argument("--validation-sizes", default="100,300,1000")
    parser.add_argument("--n-parent", type=int, default=1200)
    parser.add_argument("--n-validation-frame", type=int, default=1000)
    parser.add_argument("--n-events", type=int, default=4)
    parser.add_argument("--beta", type=float, default=-3.9)
    parser.add_argument("--h0-true", type=float, default=70.0)
    parser.add_argument("--sigma-mpc", type=float, default=140.0)
    parser.add_argument("--threshold", type=float, default=250.0)
    parser.add_argument("--z-min", type=float, default=0.01)
    parser.add_argument("--z-max", type=float, default=0.40)
    parser.add_argument("--h0-min", type=float, default=55.0)
    parser.add_argument("--h0-max", type=float, default=85.0)
    parser.add_argument("--h0-points", type=int, default=61)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/quadratic_h0_calibration_pilot"))
    args = parser.parse_args()
    try:
        args.gamma_values = tuple(float(x) for x in args.gamma_values.split(",") if x.strip())
        args.validation_sizes = tuple(int(x) for x in args.validation_sizes.split(",") if x.strip())
    except ValueError as exc:
        parser.error(str(exc))
    if args.replicates < 1 or args.seed < 0:
        parser.error("replicates must be positive and seed non-negative")
    if (not args.gamma_values or len(set(args.gamma_values)) != len(args.gamma_values)
            or any(not np.isfinite(x) for x in args.gamma_values)):
        parser.error("gamma values must be distinct and finite")
    if (not args.validation_sizes or tuple(sorted(set(args.validation_sizes))) != args.validation_sizes
            or args.validation_sizes[0] < 5 or args.validation_sizes[-1] > args.n_validation_frame):
        parser.error("validation sizes must be increasing, distinct, and within the validation frame")
    if args.n_parent < 100 or args.n_validation_frame < 5 or args.n_events < 1:
        parser.error("invalid population, validation-frame, or event count")
    if args.h0_points < 11 or not args.h0_min < args.h0_true < args.h0_max:
        parser.error("H0 grid must contain the generating truth in its interior")
    result = run_pilot(args)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
