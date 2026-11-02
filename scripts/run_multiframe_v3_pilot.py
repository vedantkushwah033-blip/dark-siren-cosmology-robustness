"""Run a pilot v3 paired multi-frame finite-population benchmark."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import platform
from pathlib import Path

import numpy as np
from scipy.stats import norm

from dark_siren_robustness.finite_population_mean_bootstrap import (
    studentized_pseudopopulation_srswor_bootstrap,
)
from dark_siren_robustness.multiframe_design import validate_multiframe_manifest, generate_finite_frame
from dark_siren_robustness.multiframe_sampling import generate_outer_sample
from dark_siren_robustness.multiframe_seed_schedule import method_seed, seed_schedule_metadata

METHODS = {"wald_srswor", "ordinary_percentile_bootstrap", "studentized_pseudopop_srswor"}


def _method_coordinate(method_id: str) -> int:
    """Stable method coordinate independent of manifest order."""
    return int.from_bytes(hashlib.sha256(method_id.encode("utf-8")).digest()[:4], "little")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def run_multiframe_v3(manifest: dict, output_dir: Path) -> dict:
    design = validate_multiframe_manifest(manifest)
    method_ids = [item["method_id"] for item in design["methods"]]
    if set(method_ids) != METHODS or len(method_ids) != len(METHODS):
        raise ValueError(f"v3 pilot requires exactly these method_id values: {sorted(METHODS)}")
    output_dir.mkdir(parents=True, exist_ok=True)

    frame_path = output_dir / "frames.csv"
    replicate_path = output_dir / "replicates.csv"
    summary_path = output_dir / "frame_summaries.csv"
    aggregate_path = output_dir / "equal_frame_summaries.csv"

    frames = {}
    with frame_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=[
            "scenario_id", "frame_id", "unit_id", "value", "frame_seed", "frame_sha256"
        ])
        writer.writeheader()
        for scenario in design["scenarios"]:
            for fid in range(design["frames_per_scenario"]):
                frame = generate_finite_frame(design, scenario["scenario_id"], fid)
                values = np.asarray(frame["values"], dtype=np.float64)
                digest = hashlib.sha256(
                    np.ascontiguousarray(values, dtype="<f8").tobytes()
                ).hexdigest()
                frames[(scenario["scenario_id"], fid)] = (values, float(frame["population_mean"]), digest, int(frame["seed"]))
                for unit, value in enumerate(values):
                    writer.writerow({
                        "scenario_id": scenario["scenario_id"], "frame_id": fid,
                        "unit_id": unit, "value": format(float(value), ".17g"),
                        "frame_seed": frame["seed"], "frame_sha256": digest,
                    })

    replicate_fields = [
        "scenario_id", "frame_id", "sample_size", "replicate_id", "method_id",
        "method_seed", "sample_seed", "sample_indices", "sample_indices_sha256",
        "frame_sha256", "finite_population_mean", "sample_mean", "estimate",
        "interval_lower", "interval_upper", "covered", "interval_success",
        "failure_code", "inner_attempts", "inner_successes", "inner_failures",
        "inner_failure_code_counts",
    ]
    frame_records = []
    with replicate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=replicate_fields)
        writer.writeheader()
        for scenario in design["scenarios"]:
            sid = scenario["scenario_id"]
            for fid in range(design["frames_per_scenario"]):
                for n in design["sample_sizes"]:
                    for rid in range(design["outer_repetitions"]):
                        sample = generate_outer_sample(design, sid, fid, n, rid)
                        x = sample["values"]
                        target = float(sample["finite_population_mean"])
                        for method_id in method_ids:
                            mseed = method_seed(
                                design["root_seed"], sid, fid, n, rid,
                                _method_coordinate(method_id),
                            )
                            lower = upper = float("nan")
                            failure = ""
                            inner = {
                                "inner_attempts": "NA", "inner_successes": "NA",
                                "inner_failures": "NA", "inner_failure_code_counts": "NA",
                            }
                            try:
                                if method_id == "wald_srswor":
                                    se = np.sqrt((1.0 - n / design["population_size"]) * np.var(x, ddof=1) / n)
                                    z = float(norm.ppf(0.5 + design["nominal_level"] / 2.0))
                                    lower, upper = float(np.mean(x) - z * se), float(np.mean(x) + z * se)
                                elif method_id == "ordinary_percentile_bootstrap":
                                    rng = np.random.default_rng(mseed)
                                    draws = rng.integers(0, n, size=(design["bootstrap_replicates"], n))
                                    means = np.mean(x[draws], axis=1)
                                    alpha = (1.0 - design["nominal_level"]) / 2.0
                                    lower, upper = map(float, np.quantile(means, [alpha, 1.0-alpha]))
                                else:
                                    result = studentized_pseudopopulation_srswor_bootstrap(
                                        x, design["population_size"],
                                        n_bootstrap=design["bootstrap_replicates"],
                                        seed=mseed, confidence=design["nominal_level"],
                                    )
                                    inner = {
                                        "inner_attempts": result["n_attempted"],
                                        "inner_successes": result["n_successful"],
                                        "inner_failures": result["n_failed"],
                                        "inner_failure_code_counts": json.dumps({
                                            str(k): int(v) for k, v in __import__("collections").Counter(
                                                str(code) for code in result["failure_codes"] if code
                                            ).items()
                                        }, sort_keys=True, separators=(",", ":")),
                                    }
                                    if result["status"] != "ok":
                                        raise ArithmeticError(str(result["status"]))
                                    lower, upper = map(float, result["interval"])
                                if not np.isfinite(lower) or not np.isfinite(upper) or lower > upper:
                                    raise ArithmeticError("invalid_interval_endpoints")
                            except (ValueError, ArithmeticError, FloatingPointError) as exc:
                                failure = type(exc).__name__ + ":" + str(exc)
                                lower = upper = float("nan")

                            success = not failure
                            covered = int(lower <= target <= upper) if success else ""
                            row = {
                                "scenario_id": sid, "frame_id": fid, "sample_size": n,
                                "replicate_id": rid, "method_id": method_id,
                                "method_seed": mseed, "sample_seed": sample["seed"],
                                "sample_indices": json.dumps(sample["indices"].tolist(), separators=(",", ":")),
                                "sample_indices_sha256": sample["indices_sha256"],
                                "frame_sha256": sample["frame_sha256"],
                                "finite_population_mean": format(target, ".17g"),
                                "sample_mean": format(float(np.mean(x)), ".17g"),
                                "estimate": format(float(np.mean(x)), ".17g"),
                                "interval_lower": lower, "interval_upper": upper,
                                "covered": covered, "interval_success": int(success),
                                "failure_code": failure, **inner,
                            }
                            writer.writerow(row)
                            frame_records.append({
                                "scenario_id": sid, "frame_id": fid, "sample_size": n,
                                "method_id": method_id, "success": success,
                                "covered": int(covered) if success else 0,
                                "width": upper-lower if success else float("nan"),
                                "target": target,
                            })

    def summarize(records, fields):
        grouped = {}
        for rec in records:
            key = tuple(rec[k] for k in fields)
            grouped.setdefault(key, []).append(rec)
        result = []
        for key, values in grouped.items():
            ok = [v for v in values if v["success"]]
            successes = len(ok)
            unconditional = sum(v["covered"] for v in values) / len(values)
            result.append(dict(zip(fields, key)) | {
                "outer_attempts": len(values), "successful_intervals": successes,
                "failed_intervals": len(values)-successes,
                "interval_success_rate": successes/len(values),
                "conditional_coverage": sum(v["covered"] for v in ok)/successes if successes else float("nan"),
                "unconditional_success_and_cover": unconditional,
                "mean_interval_width": float(np.mean([v["width"] for v in ok])) if ok else float("nan"),
                "finite_population_mean": values[0]["target"],
            })
        return result

    frame_summaries = summarize(frame_records, ["scenario_id", "frame_id", "sample_size", "method_id"])
    with summary_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(frame_summaries[0]))
        writer.writeheader(); writer.writerows(frame_summaries)

    aggregate_records = []
    for sid in [s["scenario_id"] for s in design["scenarios"]]:
        for n in design["sample_sizes"]:
            for method_id in method_ids:
                selected = [r for r in frame_summaries if r["scenario_id"] == sid and r["sample_size"] == n and r["method_id"] == method_id]
                aggregate_records.append({
                    "scenario_id": sid, "sample_size": n, "method_id": method_id,
                    "frames": len(selected),
                    "equal_frame_weighted_conditional_coverage": float(np.mean([r["conditional_coverage"] for r in selected])),
                    "equal_frame_weighted_success_rate": float(np.mean([r["interval_success_rate"] for r in selected])),
                    "equal_frame_weighted_success_and_cover": float(np.mean([r["unconditional_success_and_cover"] for r in selected])),
                    "equal_frame_weighted_mean_width": float(np.mean([r["mean_interval_width"] for r in selected])),
                })
    with aggregate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(aggregate_records[0]))
        writer.writeheader(); writer.writerows(aggregate_records)

    artifact = {
        "schema": "linear_finite_population_multiframe_pilot_v3",
        "status": "PILOT_NOT_CONFIRMATORY",
        "manifest": design,
        "seed_schedule": seed_schedule_metadata(),
        "methods": method_ids,
        "python": platform.python_version(), "numpy": np.__version__,
        "frame_rows": design["population_size"] * design["frames_per_scenario"] * len(design["scenarios"]),
        "outer_sample_keys": len(design["scenarios"]) * design["frames_per_scenario"] * len(design["sample_sizes"]) * design["outer_repetitions"],
        "replicate_rows": len(frame_records),
        "frame_csv_sha256": _sha256(frame_path),
        "replicate_csv_sha256": _sha256(replicate_path),
        "frame_summary_csv_sha256": _sha256(summary_path),
        "equal_frame_summary_csv_sha256": _sha256(aggregate_path),
        "scientific_scope": "synthetic linear finite-population arithmetic mean only; no inference for nonlinear dark-siren posterior functionals",
    }
    (output_dir / "manifest.json").write_text(json.dumps(artifact, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    return artifact


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--manifest", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    result = run_multiframe_v3(json.loads(args.manifest.read_text(encoding="utf-8")), args.output_dir)
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
