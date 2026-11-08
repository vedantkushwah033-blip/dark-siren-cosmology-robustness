"""Run paired synthetic SBC/PIT and fixed-truth coverage experiments.

Each replicate holds its analysis catalogue and GW events fixed across nested
validation sizes. Every cell is executed through the audited single-replicate
runner; failed and non-evaluable attempts remain in the ledger.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd


def _wilson_interval(successes: int, trials: int, z: float = 1.959963984540054):
    if trials <= 0:
        return (float("nan"), float("nan"))
    p = successes / trials
    denom = 1.0 + z * z / trials
    center = (p + z * z / (2.0 * trials)) / denom
    half = z * math.sqrt(p * (1.0 - p) / trials + z * z / (4.0 * trials * trials)) / denom
    return (max(0.0, center - half), min(1.0, center + half))


def _posterior_pit(density_frame: pd.DataFrame, truth: float, density_column: str) -> float:
    grid = density_frame["H0"].to_numpy(dtype=float)
    density = density_frame[density_column].to_numpy(dtype=float)
    if (
        grid.ndim != 1 or density.shape != grid.shape or grid.size < 2
        or np.any(~np.isfinite(grid)) or np.any(~np.isfinite(density))
        or np.any(np.diff(grid) <= 0) or np.any(density < 0)
    ):
        raise ValueError("invalid posterior grid/density")
    if truth < grid[0] or truth > grid[-1]:
        raise ValueError("generating H0 lies outside posterior support")
    increments = 0.5 * (density[1:] + density[:-1]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ValueError("posterior density has invalid integrated mass")
    cdf /= cdf[-1]
    return float(np.interp(truth, grid, cdf))


def _seed(master: int, mode_id: int, floor_id: int, replicate: int, stream: int) -> int:
    state = np.random.SeedSequence([master, mode_id, floor_id, replicate, stream])
    return int(state.generate_state(1, dtype=np.uint32)[0])


def _summarize(replicates: pd.DataFrame, bins: int) -> dict:
    result = {}
    grouped = replicates.groupby(["experiment", "catalogue_floor", "n_validation"], dropna=False)
    for key, frame in grouped:
        valid = frame[frame["status"] == "SUCCESS"]
        item = {
            "experiment": key[0],
            "catalogue_floor": float(key[1]),
            "n_validation": int(key[2]),
            "attempted": int(len(frame)),
            "success": int(len(valid)),
            "non_evaluable": int((frame["status"] == "NON_EVALUABLE").sum()),
            "failed": int((frame["status"] == "FAILED").sum()),
            "not_run": int((frame["status"] == "NOT_RUN").sum()),
            "evaluable_fraction": float(len(valid) / len(frame)) if len(frame) else float("nan"),
        }
        if key[0] == "SBC" and len(valid):
            pit = valid["pit"].to_numpy(dtype=float)
            counts, _ = np.histogram(pit, bins=np.linspace(0.0, 1.0, bins + 1))
            item.update({
                "pit_mean": float(np.mean(pit)),
                "pit_variance": float(np.var(pit, ddof=1)) if len(pit) > 1 else None,
                "pit_bin_counts": counts.tolist(),
                "pit_bin_expected": float(len(pit) / bins),
                "pit_mean_mc_se": float(np.std(pit, ddof=1) / math.sqrt(len(pit))) if len(pit) > 1 else None,
            })
            for level, column in ((68, "sbc_included_68"), (90, "sbc_included_90"), (95, "sbc_included_95")):
                successes = int(valid[column].sum())
                lo, hi = _wilson_interval(successes, len(valid))
                item[f"prior_predictive_interval_inclusion_{level}"] = successes / len(valid)
                item[f"prior_predictive_interval_inclusion_{level}_wilson95"] = [lo, hi]
        elif key[0] == "FIXED_TRUTH" and len(valid):
            for level, column in ((68, "covered_68"), (90, "covered_90"), (95, "covered_95")):
                successes = int(valid[column].sum())
                lo, hi = _wilson_interval(successes, len(valid))
                item[f"coverage_{level}"] = successes / len(valid)
                item[f"coverage_{level}_wilson95"] = [lo, hi]
                item[f"coverage_{level}_mc_se"] = math.sqrt(
                    (successes / len(valid)) * (1.0 - successes / len(valid)) / len(valid)
                )
                width_col = f"width_{level}"
                item[f"median_interval_width_{level}"] = float(valid[width_col].median())
            item["mean_posterior_median_bias"] = float((valid["h0_median"] - valid["h0_true"]).mean())
            item["median_absolute_posterior_median_bias"] = float((valid["h0_median"] - valid["h0_true"]).abs().median())
            for diagnostic in ("coefficient_quadrature_ess", "minimum_catalogue_kish_ess", "maximum_catalogue_normalized_weight"):
                item[f"median_{diagnostic}"] = float(valid[diagnostic].median())
        result[f"{key[0]}|floor={key[1]:g}|n={key[2]}"] = item
    return result


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--experiments", default="SBC,FIXED_TRUTH")
    parser.add_argument("--floors", default="0.0,0.10,0.25")
    parser.add_argument("--validation-sizes", default="100,300,1000")
    parser.add_argument("--n-validation-frame", type=int, default=4000)
    parser.add_argument("--n-parent", type=int, default=2400)
    parser.add_argument("--n-detected", type=int, default=8)
    parser.add_argument("--h0-points", type=int, default=81)
    parser.add_argument("--orders", default="7,11,15,19")
    parser.add_argument("--degree", type=int, choices=(1, 2), default=1)
    parser.add_argument("--catalogue-midpoint", type=float, default=0.22)
    parser.add_argument("--catalogue-width", type=float, default=0.05)
    parser.add_argument("--prior-slope-multipliers", default="0.5,1.0,2.0")
    parser.add_argument("--density-tolerance", type=float, default=1e-3)
    parser.add_argument("--quantile-tolerance", type=float, default=0.05)
    parser.add_argument("--consecutive-passes", type=int, default=2)
    parser.add_argument("--prior-min", type=float, default=50.0)
    parser.add_argument("--prior-max", type=float, default=90.0)
    parser.add_argument("--fixed-h0", type=float, default=70.0)
    parser.add_argument("--pit-bins", type=int, default=10)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/inclusion_posterior_calibration"))
    parser.add_argument("--stop-on-error", action="store_true")
    args = parser.parse_args()

    experiments = tuple(x.strip().upper() for x in args.experiments.split(",") if x.strip())
    floors = tuple(float(x) for x in args.floors.split(",") if x.strip())
    sizes = tuple(int(x) for x in args.validation_sizes.split(",") if x.strip())
    if args.replicates < 1 or args.seed < 0 or not experiments:
        parser.error("replicates and seed must be valid and experiments non-empty")
    if set(experiments) - {"SBC", "FIXED_TRUTH"}:
        parser.error("experiments may contain only SBC and FIXED_TRUTH")
    if len(set(floors)) != len(floors) or any(not np.isfinite(x) or x < 0 or x >= 1 for x in floors):
        parser.error("floors must be distinct values in [0,1)")
    if len(set(sizes)) != len(sizes) or any(x < 5 or x > args.n_validation_frame for x in sizes):
        parser.error("validation sizes must be distinct and within [5, frame size]")
    if tuple(sorted(sizes)) != sizes:
        parser.error("validation sizes must be increasing so the paired design is explicit")
    if args.prior_min >= args.prior_max or not args.prior_min < args.fixed_h0 < args.prior_max:
        parser.error("H0 prior bounds must contain the fixed truth in their interior")
    if args.pit_bins < 2:
        parser.error("pit-bins must be at least 2")
    if not np.isfinite(args.catalogue_midpoint) or args.catalogue_midpoint <= 0:
        parser.error("catalogue midpoint must be finite and positive")
    if not np.isfinite(args.catalogue_width) or args.catalogue_width <= 0:
        parser.error("catalogue width must be finite and positive")
    if args.consecutive_passes < 1 or args.density_tolerance <= 0 or args.quantile_tolerance <= 0:
        parser.error("quadrature convergence tolerances and pass count must be positive")

    root = Path(__file__).resolve().parents[1]
    smoke_script = root / "scripts" / "run_independent_validation_inclusion_h0_smoke.py"
    audit_script = root / "scripts" / "audit_independent_validation_inclusion_h0_smoke.py"
    output = args.output_dir
    output.mkdir(parents=True, exist_ok=True)
    if any(output.iterdir()):
        parser.error("output directory must be empty; use a new path to preserve prior run artifacts")
    ledger = []
    planned_rows = []
    for mode_id, experiment in enumerate(experiments):
        for floor_id, floor in enumerate(floors):
            for replicate in range(args.replicates):
                analysis_seed = _seed(args.seed, mode_id, floor_id, replicate, 0)
                validation_seed = _seed(args.seed, mode_id, floor_id, replicate, 1)
                if analysis_seed == validation_seed:
                    validation_seed = (validation_seed + 1) % (2**32)
                if experiment == "SBC":
                    rng = np.random.default_rng(_seed(args.seed, mode_id, floor_id, replicate, 2))
                    h0_true = float(rng.uniform(args.prior_min, args.prior_max))
                else:
                    h0_true = float(args.fixed_h0)
                for n_validation in sizes:
                    cell = f"{experiment.lower()}_floor_{floor:g}_n_{n_validation}_rep_{replicate:04d}"
                    planned_rows.append({
                        "cell_id": cell, "experiment": experiment,
                        "catalogue_floor": floor, "replicate": replicate,
                        "n_validation": n_validation, "analysis_seed": analysis_seed,
                        "validation_seed": validation_seed, "h0_true": h0_true,
                        "status": "NOT_RUN",
                    })
    pd.DataFrame(planned_rows).to_csv(output / "planned_cells.csv", index=False)
    planned = len(planned_rows)
    attempt = 0

    for mode_id, experiment in enumerate(experiments):
        for floor_id, floor in enumerate(floors):
            for replicate in range(args.replicates):
                analysis_seed = _seed(args.seed, mode_id, floor_id, replicate, 0)
                validation_seed = _seed(args.seed, mode_id, floor_id, replicate, 1)
                if analysis_seed == validation_seed:
                    validation_seed = (validation_seed + 1) % (2**32)
                if experiment == "SBC":
                    rng = np.random.default_rng(_seed(args.seed, mode_id, floor_id, replicate, 2))
                    h0_true = float(rng.uniform(args.prior_min, args.prior_max))
                else:
                    h0_true = float(args.fixed_h0)

                for n_validation in sizes:
                    attempt += 1
                    cell = f"{experiment.lower()}_floor_{floor:g}_n_{n_validation}_rep_{replicate:04d}"
                    artifact_dir = output / "replicates" / cell
                    command = [
                        sys.executable, str(smoke_script),
                        "--seed", str(analysis_seed),
                        "--validation-seed", str(validation_seed),
                        "--n-parent", str(args.n_parent),
                        "--n-detected", str(args.n_detected),
                        "--h0-true", repr(h0_true),
                        "--n-validation-frame", str(args.n_validation_frame),
                        "--n-validation", str(n_validation),
                        "--catalogue-floor", repr(floor),
                        "--catalogue-midpoint", repr(args.catalogue_midpoint),
                        "--catalogue-width", repr(args.catalogue_width),
                        "--degree", str(args.degree),
                        "--prior-slope-multipliers", args.prior_slope_multipliers,
                        "--density-tolerance", repr(args.density_tolerance),
                        "--quantile-tolerance", repr(args.quantile_tolerance),
                        "--consecutive-passes", str(args.consecutive_passes),
                        "--orders", args.orders,
                        "--h0-points", str(args.h0_points),
                        "--prior-min", repr(args.prior_min),
                        "--prior-max", repr(args.prior_max),
                        "--output-dir", str(artifact_dir),
                    ]
                    row = {
                        "cell_id": cell,
                        "experiment": experiment,
                        "catalogue_floor": floor,
                        "replicate": replicate,
                        "n_validation": n_validation,
                        "analysis_seed": analysis_seed,
                        "validation_seed": validation_seed,
                        "h0_true": h0_true,
                        "artifact_dir": str(artifact_dir),
                        "status": "NOT_RUN",
                        "stage": "NOT_STARTED",
                        "error": "",
                        "pit": float("nan"),
                        "config_sha256": "",
                        "n_generation_attempts": None,
                        "n_analysis_catalogue": None,
                        "quadrature_orders_evaluated": "",
                        "converged_at_order": None,
                        "h0_grid_points": args.h0_points,
                    }
                    try:
                        artifact_dir.parent.mkdir(parents=True, exist_ok=True)
                        row["stage"] = "RUN_SMOKE"
                        run = subprocess.run(command, capture_output=True, text=True, check=False)
                        (artifact_dir.parent / f"{cell}.runner.stdout.txt").write_text(run.stdout, encoding="utf-8")
                        (artifact_dir.parent / f"{cell}.runner.stderr.txt").write_text(run.stderr, encoding="utf-8")
                        if run.returncode != 0:
                            row["status"] = "FAILED"
                            row["error"] = f"smoke runner exit={run.returncode}: {run.stderr[-2000:]}"
                            raise RuntimeError(row["error"])

                        config_path = artifact_dir / "config.json"
                        replicate_config = json.loads(config_path.read_text(encoding="utf-8"))
                        row["config_sha256"] = hashlib.sha256(config_path.read_bytes()).hexdigest()
                        row["n_generation_attempts"] = int(
                            replicate_config["analysis_generator_config"]["n_generation_attempts"]
                        )
                        row["n_analysis_catalogue"] = int(
                            replicate_config["analysis_generator_config"]["n_catalogue_observed"]
                        )
                        convergence_path = artifact_dir / "quadrature_convergence.json"
                        convergence = json.loads(convergence_path.read_text(encoding="utf-8"))
                        row["quadrature_status"] = convergence["status"]
                        row["quadrature_orders_evaluated"] = ",".join(
                            str(value) for value in convergence["orders_evaluated"]
                        )
                        row["converged_at_order"] = convergence["converged_at_order"]
                        row["stage"] = "INDEPENDENT_AUDIT"
                        audit_path = artifact_dir / "audit_report.json"
                        audit = subprocess.run(
                            [sys.executable, str(audit_script), str(artifact_dir), "--output", str(audit_path)],
                            capture_output=True, text=True, check=False,
                        )
                        (artifact_dir.parent / f"{cell}.auditor.stdout.txt").write_text(audit.stdout, encoding="utf-8")
                        (artifact_dir.parent / f"{cell}.auditor.stderr.txt").write_text(audit.stderr, encoding="utf-8")
                        if audit.returncode != 0:
                            row["status"] = "FAILED"
                            row["error"] = f"independent audit exit={audit.returncode}: {audit.stderr[-2000:]}"
                            raise RuntimeError(row["error"])
                        row["audit_status"] = "PASS"

                        if convergence["status"] != "CONVERGED":
                            row["status"] = "NON_EVALUABLE"
                            row["error"] = "predeclared quadrature-order convergence not achieved"
                        else:
                            density = pd.read_csv(artifact_dir / "h0_posterior_density_by_order.csv")
                            final_order = int(convergence["orders_evaluated"][-1])
                            column = f"density_order_{final_order}"
                            row["h0_q025"] = float(convergence["posterior_by_order"][-1]["h0_q025"])
                            row["h0_q05"] = float(convergence["posterior_by_order"][-1]["h0_q05"])
                            row["h0_q16"] = float(convergence["posterior_by_order"][-1]["h0_q16"])
                            row["h0_median"] = float(convergence["posterior_by_order"][-1]["h0_median"])
                            row["h0_q84"] = float(convergence["posterior_by_order"][-1]["h0_q84"])
                            row["h0_q95"] = float(convergence["posterior_by_order"][-1]["h0_q95"])
                            row["h0_q975"] = float(convergence["posterior_by_order"][-1]["h0_q975"])
                            sensitivity = json.loads((artifact_dir / "prior_scale_sensitivity.json").read_text(encoding="utf-8"))
                            reference = next(record for record in sensitivity["records"] if np.isclose(record["slope_multiplier"], 1.0))
                            row["coefficient_quadrature_ess"] = float(convergence["posterior_by_order"][-1]["quadrature_ess"])
                            row["minimum_catalogue_kish_ess"] = float(reference["minimum_catalogue_kish_ess"])
                            row["maximum_catalogue_normalized_weight"] = float(reference["maximum_catalogue_normalized_weight"])
                            if experiment == "SBC":
                                row["pit"] = _posterior_pit(density, h0_true, column)
                                for level, lower, upper in (
                                    (68, "h0_q16", "h0_q84"),
                                    (90, "h0_q05", "h0_q95"),
                                    (95, "h0_q025", "h0_q975"),
                                ):
                                    row[f"sbc_included_{level}"] = int(row[lower] <= h0_true <= row[upper])
                            else:
                                for level, lower, upper in (
                                    (68, "h0_q16", "h0_q84"),
                                    (90, "h0_q05", "h0_q95"),
                                    (95, "h0_q025", "h0_q975"),
                                ):
                                    row[f"covered_{level}"] = int(row[lower] <= h0_true <= row[upper])
                                    row[f"width_{level}"] = float(row[upper] - row[lower])
                            row["status"] = "SUCCESS"
                        row["stage"] = "COMPLETE"
                    except Exception as exc:
                        if row["status"] == "NOT_RUN":
                            row["status"] = "FAILED"
                        row["error"] = row["error"] or f"{type(exc).__name__}: {exc}"
                        if args.stop_on_error:
                            ledger.append(row)
                            pd.DataFrame(ledger).to_csv(output / "replicate_ledger.csv", index=False)
                            raise
                    ledger.append(row)
                    pd.DataFrame(ledger).to_csv(output / "replicate_ledger.csv", index=False)
                    print(f"[{attempt}/{planned}] {cell}: {row['status']} ({row['stage']})", flush=True)

    completed_cells = {row["cell_id"] for row in ledger}
    for planned_row in planned_rows:
        if planned_row["cell_id"] not in completed_cells:
            ledger.append({
                **planned_row,
                "stage": "NOT_STARTED",
                "error": "planned but not attempted before orchestration stopped",
                "pit": float("nan"),
            })
    frame = pd.DataFrame(ledger)
    frame.to_csv(output / "replicate_ledger.csv", index=False)
    summary = {
        "experiment": "synthetic_independent_validation_inclusion_posterior_calibration",
        "design": {
            "replicates_per_experiment_floor_cell": args.replicates,
            "experiments": list(experiments),
            "catalogue_floors": list(floors),
            "validation_sizes": list(sizes),
            "validation_frame_size": args.n_validation_frame,
            "n_parent": args.n_parent,
            "n_detected": args.n_detected,
            "h0_prior": [args.prior_min, args.prior_max],
            "fixed_h0_truth": args.fixed_h0,
            "quadrature_orders": args.orders,
            "quadrature_density_tolerance": args.density_tolerance,
            "quadrature_quantile_tolerance": args.quantile_tolerance,
            "quadrature_required_consecutive_passes": args.consecutive_passes,
            "inclusion_degree": args.degree,
            "catalogue_midpoint": args.catalogue_midpoint,
            "catalogue_width": args.catalogue_width,
            "prior_slope_multipliers": args.prior_slope_multipliers,
            "h0_grid_points": args.h0_points,
            "master_seed": args.seed,
            "paired_design": "same analysis and validation seeds within replicate/floor/experiment; validation samples are permutation prefixes",
        },
        "planned_cells": int(len(frame)),
        "attempted": int((frame["status"] != "NOT_RUN").sum()),
        "status_counts": frame["status"].value_counts().to_dict(),
        "summary_by_cell": _summarize(frame, args.pit_bins),
        "warning": "Pilot Monte Carlo summaries are descriptive. SBC/PIT and fixed-truth coverage are distinct; no real-data inference is implied.",
    }
    (output / "calibration_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True, allow_nan=True) + "\n",
        encoding="utf-8",
    )
    print(json.dumps({"attempted": len(frame), "status_counts": summary["status_counts"], "summary_path": str(output / "calibration_summary.json")}, indent=2))


if __name__ == "__main__":
    main()
