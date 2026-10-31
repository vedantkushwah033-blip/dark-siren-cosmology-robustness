"""Independent integrity and numerical auditor for the v3 multi-frame pilot."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import Counter, defaultdict
from pathlib import Path

import numpy as np
from scipy.stats import norm

from dark_siren_robustness.multiframe_design import (
    generate_finite_frame, validate_multiframe_manifest,
)
from dark_siren_robustness.multiframe_seed_schedule import (
    frame_seed, method_seed, outer_sample_seed,
)

METHODS = {"wald_srswor", "ordinary_percentile_bootstrap", "studentized_pseudopop_srswor"}


def _method_coordinate(method_id: str) -> int:
    return int.from_bytes(hashlib.sha256(method_id.encode("utf-8")).digest()[:4], "little")


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _frame_hash(values: np.ndarray) -> str:
    return hashlib.sha256(np.ascontiguousarray(values, dtype="<f8").tobytes()).hexdigest()


def _indices_hash(indices: np.ndarray) -> str:
    canonical = np.ascontiguousarray(np.asarray(indices, dtype="<i8"))
    digest = hashlib.sha256()
    digest.update(str(canonical.shape).encode("ascii"))
    digest.update(canonical.tobytes())
    return digest.hexdigest()


def _studentized_replay(sample: np.ndarray, population_size: int, replicates: int,
                        seed: int, confidence: float) -> tuple[float, float, int, int, dict[str, int]]:
    """Reimplement the pseudo-population algorithm without importing the runner."""
    n = len(sample)
    rng = np.random.default_rng(seed)
    valid_t = []
    failures: Counter[str] = Counter()
    for _ in range(replicates):
        multiplicities = rng.multinomial(population_size, np.full(n, 1.0 / n))
        pseudo = np.repeat(sample, multiplicities)
        draw = pseudo[rng.choice(population_size, size=n, replace=False)]
        se = float(np.sqrt((1.0 - n / population_size) * np.var(draw, ddof=1) / n))
        if not np.isfinite(se) or se <= 0:
            failures["zero_or_nonfinite_bootstrap_variance"] += 1
        else:
            valid_t.append((float(np.mean(draw)) - float(np.mean(pseudo))) / se)
    if len(valid_t) < 2:
        raise ArithmeticError("insufficient_valid_bootstrap_replicates")
    alpha = (1.0 - confidence) / 2.0
    qlo, qhi = np.quantile(np.asarray(valid_t), [alpha, 1.0 - alpha])
    observed_se = float(np.sqrt((1.0 - n / population_size) * np.var(sample, ddof=1) / n))
    estimate = float(np.mean(sample))
    return (estimate - qhi * observed_se, estimate - qlo * observed_se,
            replicates, len(valid_t), dict(sorted(failures.items())))


def audit_artifact(output_dir: Path) -> dict[str, object]:
    output_dir = Path(output_dir)
    artifact = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    if artifact.get("schema") != "linear_finite_population_multiframe_pilot_v3":
        raise ValueError("unexpected v3 artifact schema")
    if artifact.get("status") != "PILOT_NOT_CONFIRMATORY":
        raise ValueError("v3 artifact must remain pilot-only")
    design = validate_multiframe_manifest(artifact.get("manifest", {}))
    methods = [m["method_id"] for m in design["methods"]]
    if set(methods) != METHODS or len(methods) != len(METHODS):
        raise ValueError("method set differs from the locked v3 pilot contract")
    files = {
        "frame_csv_sha256": "frames.csv",
        "replicate_csv_sha256": "replicates.csv",
        "frame_summary_csv_sha256": "frame_summaries.csv",
        "equal_frame_summary_csv_sha256": "equal_frame_summaries.csv",
    }
    for key, name in files.items():
        if artifact.get(key) != _sha256(output_dir / name):
            raise ValueError(f"checksum mismatch: {name}")

    frames: dict[tuple[int, int], np.ndarray] = {}
    frame_labels: dict[tuple[int, int], tuple[int, str]] = {}
    with (output_dir / "frames.csv").open(newline="", encoding="utf-8") as f:
        grouped: dict[tuple[int, int], list[tuple[int, float]]] = defaultdict(list)
        for row in csv.DictReader(f):
            key = (int(row["scenario_id"]), int(row["frame_id"]))
            grouped[key].append((int(row["unit_id"]), float(row["value"])))
            label = (int(row["frame_seed"]), row["frame_sha256"])
            if key in frame_labels and frame_labels[key] != label:
                raise ValueError("inconsistent frame provenance labels")
            frame_labels[key] = label

    expected_frames = {(s["scenario_id"], fid) for s in design["scenarios"]
                       for fid in range(design["frames_per_scenario"])}
    if set(grouped) != expected_frames:
        raise ValueError("finite-frame grid is incomplete or unexpected")
    for scenario in design["scenarios"]:
        sid = scenario["scenario_id"]
        for fid in range(design["frames_per_scenario"]):
            key = (sid, fid)
            rows = sorted(grouped[key])
            if [i for i, _ in rows] != list(range(design["population_size"])):
                raise ValueError("frame unit IDs are duplicated, missing, or out of order")
            stored = np.asarray([v for _, v in rows], dtype=np.float64)
            regenerated = np.asarray(generate_finite_frame(design, sid, fid)["values"])
            expected_seed = frame_seed(design["root_seed"], sid, fid)
            if not np.array_equal(stored, regenerated):
                raise ValueError("finite-frame values do not independently regenerate")
            if frame_labels[key] != (expected_seed, _frame_hash(regenerated)):
                raise ValueError("frame seed or hash mismatch")
            frames[key] = regenerated

    expected_outer = {
        (s["scenario_id"], fid, n, rid)
        for s in design["scenarios"]
        for fid in range(design["frames_per_scenario"])
        for n in design["sample_sizes"]
        for rid in range(design["outer_repetitions"])
    }
    by_outer: dict[tuple[int, int, int, int], list[dict[str, str]]] = defaultdict(list)
    seen: set[tuple[int, int, int, int, str]] = set()
    with (output_dir / "replicates.csv").open(newline="", encoding="utf-8") as f:
        for row in csv.DictReader(f):
            key = (int(row["scenario_id"]), int(row["frame_id"]),
                   int(row["sample_size"]), int(row["replicate_id"]))
            method = row["method_id"]
            if key not in expected_outer or method not in METHODS:
                raise ValueError("unexpected replicate key or method")
            if (*key, method) in seen:
                raise ValueError("duplicate method replicate row")
            seen.add((*key, method))
            by_outer[key].append(row)
    expected_rows = len(expected_outer) * len(METHODS)
    if len(seen) != expected_rows or int(artifact.get("replicate_rows", -1)) != expected_rows:
        raise ValueError("replicate grid is incomplete or row count is wrong")

    recomputed_records = []
    for key in sorted(expected_outer):
        sid, fid, n, rid = key
        frame = frames[(sid, fid)]
        sample_seed = outer_sample_seed(design["root_seed"], sid, fid, n, rid)
        indices = np.random.default_rng(sample_seed).choice(len(frame), size=n, replace=False)
        if len(np.unique(indices)) != n:
            raise ValueError("replayed SRSWOR indices are not unique")
        sample = frame[indices]
        paired_indices: set[str] = set()
        if len(by_outer[key]) != len(METHODS):
            raise ValueError("outer sample key does not have exactly one row per method")
        for row in by_outer[key]:
            method = row["method_id"]
            stored_indices = np.asarray(json.loads(row["sample_indices"]), dtype=np.int64)
            paired_indices.add(row["sample_indices"])
            if not np.array_equal(stored_indices, indices):
                raise ValueError("stored sample indices fail exact seed replay")
            if row["sample_seed"] != str(sample_seed) or row["sample_indices_sha256"] != _indices_hash(indices):
                raise ValueError("sample seed or index digest mismatch")
            frame_digest = _frame_hash(frame)
            if row["frame_sha256"] != frame_digest:
                raise ValueError("replicate frame digest mismatch")
            target = float(np.mean(frame))
            estimate = float(np.mean(sample))
            if not np.isclose(float(row["finite_population_mean"]), target, rtol=0, atol=1e-12):
                raise ValueError("finite-population target mismatch")
            if not np.isclose(float(row["sample_mean"]), estimate, rtol=0, atol=1e-12) or not np.isclose(float(row["estimate"]), estimate, rtol=0, atol=1e-12):
                raise ValueError("sample estimate mismatch")
            mseed = method_seed(design["root_seed"], sid, fid, n, rid, _method_coordinate(method))
            if row["method_seed"] != str(mseed):
                raise ValueError("method seed does not match domain-separated schedule")

            lower = upper = float("nan")
            failure = ""
            inner = ("NA", "NA", "NA", "NA")
            if method == "wald_srswor":
                se = np.sqrt((1.0 - n / len(frame)) * np.var(sample, ddof=1) / n)
                z = float(norm.ppf(0.5 + design["nominal_level"] / 2.0))
                lower, upper = estimate - z * se, estimate + z * se
            elif method == "ordinary_percentile_bootstrap":
                draws = np.random.default_rng(mseed).integers(0, n, size=(design["bootstrap_replicates"], n))
                means = np.mean(sample[draws], axis=1)
                alpha = (1.0 - design["nominal_level"]) / 2.0
                lower, upper = np.quantile(means, [alpha, 1.0 - alpha])
            else:
                if not np.isfinite(np.var(sample, ddof=1)) or np.var(sample, ddof=1) <= 0:
                    failure = "ValueError:sample must have positive finite variance"
                else:
                    try:
                        lower, upper, attempts, successes, failures = _studentized_replay(
                            sample, len(frame), design["bootstrap_replicates"], mseed,
                            design["nominal_level"])
                        inner = (str(attempts), str(successes), str(attempts-successes),
                                 json.dumps(failures, sort_keys=True, separators=(",", ":")))
                    except ArithmeticError as exc:
                        failure = "ArithmeticError:" + str(exc)
                        # Replay once more to independently derive the failure accounting.
                        rng = np.random.default_rng(mseed)
                        bad = Counter()
                        for _ in range(design["bootstrap_replicates"]):
                            counts = rng.multinomial(len(frame), np.full(n, 1.0/n))
                            pseudo = np.repeat(sample, counts)
                            draw = pseudo[rng.choice(len(frame), size=n, replace=False)]
                            se = np.sqrt((1.0-n/len(frame))*np.var(draw, ddof=1)/n)
                            if not np.isfinite(se) or se <= 0:
                                bad["zero_or_nonfinite_bootstrap_variance"] += 1
                        inner = (str(design["bootstrap_replicates"]),
                                 str(design["bootstrap_replicates"]-sum(bad.values())),
                                 str(sum(bad.values())),
                                 json.dumps(dict(sorted(bad.items())), sort_keys=True, separators=(",", ":")))
            if method != "studentized_pseudopop_srswor":
                inner = ("NA", "NA", "NA", "NA")
            expected_inner = dict(zip(("inner_attempts", "inner_successes", "inner_failures", "inner_failure_code_counts"), inner))
            for field, value in expected_inner.items():
                if row[field] != value:
                    raise ValueError(f"inner attempt accounting mismatch: {field}")
            if row["failure_code"] != failure:
                raise ValueError("interval failure classification mismatch")
            success = not failure
            covered = int(lower <= target <= upper) if success else 0
            if int(row["interval_success"]) != int(success):
                raise ValueError("interval success indicator mismatch")
            if success:
                if not np.isclose(float(row["interval_lower"]), lower, rtol=1e-10, atol=1e-12):
                    raise ValueError("lower endpoint fails independent reconstruction")
                if not np.isclose(float(row["interval_upper"]), upper, rtol=1e-10, atol=1e-12):
                    raise ValueError("upper endpoint fails independent reconstruction")
                if int(row["covered"]) != covered:
                    raise ValueError("coverage indicator mismatch")
            else:
                if row["covered"] != "" or not np.isnan(float(row["interval_lower"])) or not np.isnan(float(row["interval_upper"])):
                    raise ValueError("failed interval must have NaN endpoints and blank coverage")
            recomputed_records.append({
                "scenario_id": sid, "frame_id": fid, "sample_size": n,
                "method_id": method, "success": success, "covered": covered,
                "width": float(upper-lower) if success else float("nan"),
                "target": target,
            })
        if len(paired_indices) != 1:
            raise ValueError("methods do not share identical paired outer-sample indices")

    def summary(records: list[dict[str, object]], keys: list[str]) -> list[dict[str, object]]:
        buckets: dict[tuple, list[dict[str, object]]] = defaultdict(list)
        for record in records:
            buckets[tuple(record[k] for k in keys)].append(record)
        output = []
        for key, values in buckets.items():
            ok = [r for r in values if r["success"]]
            count = len(values)
            output.append(dict(zip(keys, key)) | {
                "outer_attempts": count, "successful_intervals": len(ok),
                "failed_intervals": count-len(ok),
                "interval_success_rate": len(ok)/count,
                "conditional_coverage": float(np.mean([r["covered"] for r in ok])) if ok else float("nan"),
                "unconditional_success_and_cover": sum(r["covered"] for r in values)/count,
                "mean_interval_width": float(np.mean([r["width"] for r in ok])) if ok else float("nan"),
                "finite_population_mean": values[0]["target"],
            })
        return output

    expected_frame = summary(recomputed_records, ["scenario_id", "frame_id", "sample_size", "method_id"])
    with (output_dir / "frame_summaries.csv").open(newline="", encoding="utf-8") as f:
        stored_frame = list(csv.DictReader(f))
    if len(stored_frame) != len(expected_frame):
        raise ValueError("frame summary grid has missing or duplicate rows")
    def numeric_equal(stored: str, expected: object) -> bool:
        try:
            return bool(np.isclose(float(stored), float(expected), rtol=1e-10, atol=1e-12, equal_nan=True))
        except (TypeError, ValueError):
            return False
    frame_fields = ["outer_attempts", "successful_intervals", "failed_intervals",
                    "interval_success_rate", "conditional_coverage",
                    "unconditional_success_and_cover", "mean_interval_width",
                    "finite_population_mean"]
    expected_by_key = {(r["scenario_id"], r["frame_id"], r["sample_size"], r["method_id"]): r for r in expected_frame}
    observed_keys = set()
    for row in stored_frame:
        k = (int(row["scenario_id"]), int(row["frame_id"]), int(row["sample_size"]), row["method_id"])
        if k in observed_keys or k not in expected_by_key:
            raise ValueError("frame summary contains duplicate or unexpected key")
        observed_keys.add(k)
        for field in frame_fields:
            if not numeric_equal(row[field], expected_by_key[k][field]):
                raise ValueError(f"frame summary value mismatch: {field}")
    if observed_keys != set(expected_by_key):
        raise ValueError("frame summary grid is incomplete")

    expected_aggregate = []
    for scenario in design["scenarios"]:
        sid = scenario["scenario_id"]
        for n in design["sample_sizes"]:
            for method in methods:
                cells = [r for r in expected_frame if r["scenario_id"] == sid and r["sample_size"] == n and r["method_id"] == method]
                expected_aggregate.append({
                    "scenario_id": sid, "sample_size": n, "method_id": method, "frames": len(cells),
                    "equal_frame_weighted_conditional_coverage": float(np.mean([r["conditional_coverage"] for r in cells])),
                    "equal_frame_weighted_success_rate": float(np.mean([r["interval_success_rate"] for r in cells])),
                    "equal_frame_weighted_success_and_cover": float(np.mean([r["unconditional_success_and_cover"] for r in cells])),
                    "equal_frame_weighted_mean_width": float(np.mean([r["mean_interval_width"] for r in cells])),
                })
    with (output_dir / "equal_frame_summaries.csv").open(newline="", encoding="utf-8") as f:
        stored_aggregate = list(csv.DictReader(f))
    if len(stored_aggregate) != len(expected_aggregate):
        raise ValueError("equal-frame summary grid has missing or duplicate rows")
    aggregate_fields = ["frames", "equal_frame_weighted_conditional_coverage",
                        "equal_frame_weighted_success_rate", "equal_frame_weighted_success_and_cover",
                        "equal_frame_weighted_mean_width"]
    aggregate_keys = set()
    for row in stored_aggregate:
        k = (int(row["scenario_id"]), int(row["sample_size"]), row["method_id"])
        if k in aggregate_keys:
            raise ValueError("duplicate equal-frame summary row")
        aggregate_keys.add(k)
        match = next((r for r in expected_aggregate if (r["scenario_id"], r["sample_size"], r["method_id"]) == k), None)
        if match is None:
            raise ValueError("unexpected equal-frame summary key")
        for field in aggregate_fields:
            if not numeric_equal(row[field], match[field]):
                raise ValueError(f"equal-frame summary mismatch: {field}")
    if len(aggregate_keys) != len(expected_aggregate):
        raise ValueError("equal-frame summary grid is incomplete")

    return {
        "status": "PASS", "artifact_schema": artifact["schema"],
        "frames_verified": len(frames), "outer_sample_keys_verified": len(expected_outer),
        "replicate_rows_independently_reconstructed": expected_rows,
        "frame_summary_cells_verified": len(expected_frame),
        "equal_frame_summary_cells_verified": len(expected_aggregate),
        "checksums": "PASS", "seed_and_sampling_replay": "PASS",
        "paired_method_inputs": "PASS", "interval_endpoints": "PASS",
        "inner_attempt_accounting": "PASS", "summary_reconstruction": "PASS",
        "scientific_scope": artifact.get("scientific_scope"),
        "interpretation": "pilot software-validation evidence only; not confirmatory scientific evidence",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = audit_artifact(args.output_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
