from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


def _hash_frame(values: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    digest = hashlib.sha256()
    digest.update(str(arr.shape).encode("ascii"))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def _studentized_interval_independent(
    sample: np.ndarray, population_size: int, replicates: int, seed: int
) -> tuple[float, float]:
    n = sample.size
    fpc = 1.0 - n / population_size
    observed_se = float(np.sqrt(fpc * np.var(sample, ddof=1) / n))
    rng = np.random.default_rng(seed)
    t_values = []
    probabilities = np.full(n, 1.0 / n)
    for _ in range(replicates):
        counts = rng.multinomial(population_size, probabilities)
        pseudo = np.repeat(sample, counts)
        selected = rng.choice(population_size, size=n, replace=False)
        draw = pseudo[selected]
        draw_se = float(np.sqrt(fpc * np.var(draw, ddof=1) / n))
        if np.isfinite(draw_se) and draw_se > 0:
            t_values.append((float(np.mean(draw)) - float(np.mean(pseudo))) / draw_se)
    if len(t_values) < 2 or observed_se <= 0 or not np.isfinite(observed_se):
        raise ArithmeticError("insufficient valid studentized replicates")
    qlo, qhi = np.quantile(np.asarray(t_values), [0.025, 0.975])
    estimate = float(np.mean(sample))
    return float(estimate - qhi * observed_se), float(estimate - qlo * observed_se)


def audit_artifact(output_dir: Path) -> dict[str, object]:
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "linear_finite_population_bootstrap_benchmark_v1":
        raise ValueError("unexpected artifact schema")
    if manifest.get("status") != "PILOT_NOT_CONFIRMATORY":
        raise ValueError("artifact status must remain pilot-only")
    frame_path = output_dir / "finite_population_frames.csv"
    replicate_path = output_dir / "linear_finite_population_bootstrap_replicates.csv"
    summary_path = output_dir / "linear_finite_population_bootstrap.csv"
    for key, path in (
        ("frame_csv_sha256", frame_path),
        ("replicate_csv_sha256", replicate_path),
        ("summary_csv_sha256", summary_path),
    ):
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != manifest[key]:
            raise ValueError(f"checksum mismatch: {key}")

    frames: dict[int, list[tuple[int, float]]] = defaultdict(list)
    frame_hashes: dict[int, str] = {}
    with frame_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            s = int(row["scenario"])
            frames[s].append((int(row["unit"]), float(row["value"])))
            frame_hashes[s] = row["frame_sha256"]

    populations: dict[int, np.ndarray] = {}
    targets: dict[int, float] = {}
    for s, items in frames.items():
        items.sort()
        if [i for i, _ in items] != list(range(len(items))):
            raise ValueError(f"frame indices are not complete for scenario {s}")
        pop = np.asarray([v for _, v in items], dtype=float)
        if _hash_frame(pop) != frame_hashes[s]:
            raise ValueError(f"finite-frame hash mismatch for scenario {s}")
        populations[s] = pop
        targets[s] = float(np.mean(pop))

    grouped: dict[tuple[int, int, str], list[dict[str, object]]] = defaultdict(list)
    checked = 0
    with replicate_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            s, n, repetition = int(row["scenario"]), int(row["sample_size"]), int(row["repetition"])
            method = row["method"]
            pop = populations[s]
            if int(row["population_size"]) != pop.size or row["frame_sha256"] != frame_hashes[s]:
                raise ValueError("replicate finite-frame provenance mismatch")
            target = targets[s]
            if not np.isclose(float(row["target"]), target, rtol=0, atol=1e-12):
                raise ValueError("replicate target differs from exact finite-frame mean")
            indices = np.asarray(json.loads(row["sample_indices"]), dtype=int)
            replay_indices = np.random.default_rng(int(row["sample_seed"])).choice(pop.size, size=n, replace=False)
            if not np.array_equal(indices, replay_indices):
                raise ValueError("stored SRSWOR indices do not replay from the recorded seed")
            if indices.shape != (n,) or np.any(indices < 0) or np.any(indices >= pop.size):
                raise ValueError("invalid sample index vector")
            if np.unique(indices).size != n:
                raise ValueError("SRSWOR sample contains duplicate indices")
            sample = pop[indices]
            estimate = float(np.mean(sample))
            if not np.isclose(float(row["estimate"]), estimate, rtol=0, atol=1e-12):
                raise ValueError("sample estimate mismatch")
            lower = upper = float("nan")
            failure = ""
            if method == "wald_srswor":
                se = float(np.sqrt((1.0 - n / pop.size) * np.var(sample, ddof=1) / n))
                z = 1.959963984540054
                lower, upper = estimate - z * se, estimate + z * se
            elif method == "ordinary_percentile_bootstrap":
                rng = np.random.default_rng(int(row["ordinary_bootstrap_seed"]))
                b = int(manifest["bootstrap_replicates"])
                indices_b = rng.integers(0, n, size=(b, n))
                means = np.mean(sample[indices_b], axis=1)
                lower, upper = map(float, np.quantile(means, [0.025, 0.975]))
            elif method == "studentized_pseudopop_srswor":
                try:
                    lower, upper = _studentized_interval_independent(
                        sample, pop.size, int(manifest["bootstrap_replicates"]),
                        int(row["studentized_bootstrap_seed"]),
                    )
                except (ValueError, ArithmeticError, FloatingPointError) as exc:
                    failure = type(exc).__name__ + ":" + str(exc)
            else:
                raise ValueError(f"unknown method {method}")
            stored_failure = row["failure_code"]
            if failure != stored_failure:
                raise ValueError("failure classification mismatch")
            if not failure:
                if not np.isclose(float(row["interval_lower"]), lower, rtol=1e-11, atol=1e-12):
                    raise ValueError("lower endpoint does not reconstruct")
                if not np.isclose(float(row["interval_upper"]), upper, rtol=1e-11, atol=1e-12):
                    raise ValueError("upper endpoint does not reconstruct")
                covered = int(lower <= target <= upper)
                if int(row["covered"]) != covered:
                    raise ValueError("coverage indicator mismatch")
            grouped[(s, n, method)].append({
                "covered": None if failure else int(lower <= target <= upper),
                "width": None if failure else upper - lower,
                "failed": bool(failure),
            })
            checked += 1

    with summary_path.open(newline="", encoding="utf-8") as handle:
        summaries = list(csv.DictReader(handle))
    if len(summaries) != len(grouped):
        raise ValueError("summary cell count mismatch")
    for row in summaries:
        key = (int(row["scenario"]), int(row["sample_size"]), row["method"])
        records = grouped[key]
        successes = [r for r in records if not r["failed"]]
        failures = len(records) - len(successes)
        coverage = float(np.mean([r["covered"] for r in successes])) if successes else float("nan")
        width = float(np.mean([r["width"] for r in successes])) if successes else float("nan")
        if int(row["successful_intervals"]) != len(successes) or int(row["failed_intervals"]) != failures:
            raise ValueError("summary attempt accounting mismatch")
        if not np.isclose(float(row["coverage"]), coverage, rtol=1e-11, atol=1e-12, equal_nan=True):
            raise ValueError("summary coverage mismatch")
        if not np.isclose(float(row["mean_interval_width"]), width, rtol=1e-11, atol=1e-12, equal_nan=True):
            raise ValueError("summary interval-width mismatch")

    if checked != int(manifest["replicate_rows"]):
        raise ValueError("replicate row count mismatch")
    return {
        "status": "PASS",
        "replicates_independently_reconstructed": checked,
        "finite_frames_verified": len(populations),
        "summary_cells_verified": len(summaries),
        "scientific_scope": "linear finite-population mean pilot only; no transfer to nonlinear dark-siren posterior",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Independently audit linear finite-population bootstrap artifacts")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--report", type=Path, required=True)
    args = parser.parse_args()
    report = audit_artifact(args.output_dir)
    args.report.parent.mkdir(parents=True, exist_ok=True)
    args.report.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
