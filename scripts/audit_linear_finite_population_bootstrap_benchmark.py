from __future__ import annotations

import argparse
import csv
import hashlib
import json
from collections import defaultdict
from pathlib import Path

import numpy as np


METHODS = (
    "studentized_pseudopop_srswor",
    "wald_srswor",
    "ordinary_percentile_bootstrap",
)
POPULATION_KINDS = ("normal", "gamma", "lognormal")


def _hash_frame(values: np.ndarray) -> str:
    arr = np.ascontiguousarray(np.asarray(values, dtype="<f8"))
    digest = hashlib.sha256()
    digest.update(str(arr.shape).encode("ascii"))
    digest.update(arr.tobytes())
    return digest.hexdigest()


def _population(rng: np.random.Generator, kind: str, size: int) -> np.ndarray:
    if kind == "normal":
        return rng.normal(0.0, 1.0, size)
    if kind == "gamma":
        return rng.gamma(shape=0.7, scale=1.0 / 0.7, size=size)
    if kind == "lognormal":
        return rng.lognormal(mean=0.0, sigma=1.0, size=size)
    raise ValueError(f"unknown population kind: {kind}")


def _studentized_interval_independent(
    sample: np.ndarray, population_size: int, replicates: int, seed: int
) -> tuple[float, float]:
    n = sample.size
    fpc = 1.0 - n / population_size
    observed_se = float(np.sqrt(fpc * np.var(sample, ddof=1) / n))
    if not np.isfinite(observed_se) or observed_se <= 0:
        raise ValueError("observed sample must have positive finite variance")
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
    if len(t_values) < 2:
        raise ArithmeticError("studentized bootstrap did not yield a finite interval")
    qlo, qhi = np.quantile(np.asarray(t_values), [0.025, 0.975])
    estimate = float(np.mean(sample))
    return float(estimate - qhi * observed_se), float(estimate - qlo * observed_se)



def _independent_inner_attempt_accounting(
    sample: np.ndarray, population_size: int, n_bootstrap: int, seed: int
) -> dict[str, object]:
    """Replay inner multinomial and SRSWOR attempts without trusting artifact fields."""
    n = int(sample.size)
    rng = np.random.default_rng(int(seed))
    success_count = 0
    failures: dict[str, int] = {}
    probabilities = np.full(n, 1.0 / n, dtype=float)
    fpc_variance = 1.0 - n / population_size
    for _ in range(n_bootstrap):
        counts = rng.multinomial(population_size, probabilities)
        pseudo_population = np.repeat(sample, counts)
        positions = rng.choice(population_size, size=n, replace=False)
        bootstrap_sample = pseudo_population[positions]
        bootstrap_variance = float(np.var(bootstrap_sample, ddof=1))
        bootstrap_se = float(np.sqrt(fpc_variance * bootstrap_variance / n))
        if not np.isfinite(bootstrap_se) or bootstrap_se <= 0.0:
            code = "zero_or_nonfinite_bootstrap_variance"
            failures[code] = failures.get(code, 0) + 1
        else:
            success_count += 1
    return {
        "inner_attempts": n_bootstrap,
        "inner_successes": success_count,
        "inner_failures": n_bootstrap - success_count,
        "inner_failure_code_counts": json.dumps(dict(sorted(failures.items())), sort_keys=True, separators=(",", ":")),
    }


def audit_artifact(output_dir: Path) -> dict[str, object]:
    manifest = json.loads((output_dir / "manifest.json").read_text(encoding="utf-8"))
    if manifest.get("schema") != "linear_finite_population_bootstrap_benchmark_v1":
        raise ValueError("unexpected artifact schema")
    if manifest.get("status") != "PILOT_NOT_CONFIRMATORY":
        raise ValueError("artifact status must remain pilot-only")

    scenarios = int(manifest["scenarios"])
    repetitions = int(manifest["repetitions_per_scenario_and_sample_size"])
    population_size = int(manifest["population_size"])
    sample_sizes = tuple(int(n) for n in manifest["sample_sizes"])
    bootstrap_replicates = int(manifest["bootstrap_replicates"])
    root_seed = int(manifest["seed"])
    if scenarios < 1 or repetitions < 2 or population_size < 20:
        raise ValueError("manifest contains invalid design dimensions")
    if not sample_sizes or len(set(sample_sizes)) != len(sample_sizes):
        raise ValueError("manifest sample sizes must be nonempty and unique")
    if any(n < 4 or n >= population_size for n in sample_sizes):
        raise ValueError("manifest contains invalid SRSWOR sample size")
    if bootstrap_replicates < 20:
        raise ValueError("manifest bootstrap replicate count is invalid")
    if tuple(manifest.get("methods", ())) != tuple(sorted(METHODS)):
        raise ValueError("manifest method set differs from the prespecified method set")
    if manifest.get("schema") != "linear_finite_population_bootstrap_benchmark_v2":
        raise ValueError("auditor requires the v2 inner-attempt provenance schema")

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
            if s in frame_hashes and frame_hashes[s] != row["frame_sha256"]:
                raise ValueError(f"inconsistent frame hash labels for scenario {s}")
            frame_hashes[s] = row["frame_sha256"]

    expected_scenarios = set(range(scenarios))
    if set(frames) != expected_scenarios:
        raise ValueError("finite-frame scenario set is incomplete or unexpected")
    if int(manifest.get("frame_rows", -1)) != scenarios * population_size:
        raise ValueError("manifest finite-frame row count is inconsistent")

    # Recreate every frame from the recorded root seed, independently of the CSV.
    scenario_seed_sequences = np.random.SeedSequence(root_seed).spawn(scenarios)
    populations: dict[int, np.ndarray] = {}
    targets: dict[int, float] = {}
    for s in range(scenarios):
        items = sorted(frames[s])
        if [i for i, _ in items] != list(range(population_size)):
            raise ValueError(f"frame indices are not complete for scenario {s}")
        pop = np.asarray([v for _, v in items], dtype=float)
        if pop.shape != (population_size,) or not np.all(np.isfinite(pop)):
            raise ValueError(f"frame values are invalid for scenario {s}")
        kind = POPULATION_KINDS[s % len(POPULATION_KINDS)]
        expected_pop = _population(np.random.default_rng(scenario_seed_sequences[s]), kind, population_size)
        if not np.array_equal(pop, expected_pop):
            raise ValueError(f"finite frame does not regenerate from protocol seed for scenario {s}")
        actual_hash = _hash_frame(pop)
        if actual_hash != frame_hashes[s]:
            raise ValueError(f"finite-frame hash mismatch for scenario {s}")
        populations[s] = pop
        targets[s] = float(np.mean(pop))

    expected_attempts = {
        (s, n, r, method)
        for s in range(scenarios)
        for n in sample_sizes
        for r in range(repetitions)
        for method in METHODS
    }
    expected_sample_seeds: dict[tuple[int, int, int], int] = {}
    for s, child in enumerate(scenario_seed_sequences):
        seeds = child.spawn(repetitions * len(sample_sizes))
        k = 0
        for n in sample_sizes:
            for r in range(repetitions):
                expected_sample_seeds[(s, n, r)] = int(seeds[k].generate_state(1)[0])
                k += 1

    grouped: dict[tuple[int, int, str], list[dict[str, object]]] = defaultdict(list)
    checked = 0
    seen_attempts: set[tuple[int, int, int, str]] = set()
    with replicate_path.open(newline="", encoding="utf-8") as handle:
        for row in csv.DictReader(handle):
            s, n, repetition = int(row["scenario"]), int(row["sample_size"]), int(row["repetition"])
            method = row["method"]
            attempt_key = (s, n, repetition, method)
            if attempt_key in seen_attempts:
                raise ValueError("duplicate replicate/method attempt")
            seen_attempts.add(attempt_key)
            if attempt_key not in expected_attempts:
                raise ValueError("unexpected replicate/method attempt")
            pop = populations[s]
            kind = POPULATION_KINDS[s % len(POPULATION_KINDS)]
            if row["population_kind"] != kind:
                raise ValueError("population regime label mismatch")
            if int(row["population_size"]) != pop.size or row["frame_sha256"] != frame_hashes[s]:
                raise ValueError("replicate finite-frame provenance mismatch")
            target = targets[s]
            if not np.isclose(float(row["target"]), target, rtol=0, atol=1e-12):
                raise ValueError("replicate target differs from exact finite-frame mean")

            sample_seed = int(row["sample_seed"])
            if sample_seed != expected_sample_seeds[(s, n, repetition)]:
                raise ValueError("sample seed does not match the prespecified seed sequence")
            indices = np.asarray(json.loads(row["sample_indices"]), dtype=int)
            replay_indices = np.random.default_rng(sample_seed).choice(pop.size, size=n, replace=False)
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
            expected_ordinary_seed = int(
                np.random.SeedSequence([root_seed, s, n, repetition, 41031]).generate_state(1)[0]
            )
            expected_studentized_seed = int(
                np.random.SeedSequence([root_seed, s, n, repetition, 41032]).generate_state(1)[0]
            )
            if method == "wald_srswor":
                if any(row.get(k) != "NA" for k in ("inner_attempts", "inner_successes", "inner_failures", "inner_failure_code_counts")):
                    raise ValueError("non-studentized method must mark inner attempt fields NA")
                if row["ordinary_bootstrap_seed"] or row["studentized_bootstrap_seed"]:
                    raise ValueError("Wald method contains an unexpected bootstrap seed")
                se = float(np.sqrt((1.0 - n / pop.size) * np.var(sample, ddof=1) / n))
                z = 1.959963984540054
                lower, upper = estimate - z * se, estimate + z * se
            elif method == "ordinary_percentile_bootstrap":
                if any(row.get(k) != "NA" for k in ("inner_attempts", "inner_successes", "inner_failures", "inner_failure_code_counts")):
                    raise ValueError("non-studentized method must mark inner attempt fields NA")
                if int(row["ordinary_bootstrap_seed"]) != expected_ordinary_seed:
                    raise ValueError("ordinary bootstrap seed mismatch")
                if row["studentized_bootstrap_seed"]:
                    raise ValueError("ordinary bootstrap row contains unexpected studentized seed")
                rng = np.random.default_rng(expected_ordinary_seed)
                indices_b = rng.integers(0, n, size=(bootstrap_replicates, n))
                means = np.mean(sample[indices_b], axis=1)
                lower, upper = map(float, np.quantile(means, [0.025, 0.975]))
            elif method == "studentized_pseudopop_srswor":
                if int(row["studentized_bootstrap_seed"]) != expected_studentized_seed:
                    raise ValueError("studentized bootstrap seed mismatch")
                if int(row["ordinary_bootstrap_seed"]) != expected_ordinary_seed:
                    raise ValueError("studentized row paired ordinary seed mismatch")
                # A degenerate observed sample fails before inner draws; the runner records NA.
                observed_variance = float(np.var(sample, ddof=1))
                if not np.isfinite(observed_variance) or observed_variance <= 0.0:
                    expected_inner = {"inner_attempts": "NA", "inner_successes": "NA", "inner_failures": "NA", "inner_failure_code_counts": "NA"}
                else:
                    expected_inner = _independent_inner_attempt_accounting(sample, pop.size, bootstrap_replicates, expected_studentized_seed)
                for field, expected_value in expected_inner.items():
                    if row.get(field) != str(expected_value):
                        raise ValueError(f"inner attempt provenance mismatch: {field}")
                try:
                    lower, upper = _studentized_interval_independent(
                        sample, pop.size, bootstrap_replicates, expected_studentized_seed
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
            else:
                if row["covered"] or not np.isnan(float(row["interval_lower"])) or not np.isnan(float(row["interval_upper"])):
                    raise ValueError("failed attempt must have blank coverage and NaN endpoints")

            grouped[(s, n, method)].append({
                "covered": None if failure else int(lower <= target <= upper),
                "width": None if failure else upper - lower,
                "failed": bool(failure),
            })
            checked += 1

    if seen_attempts != expected_attempts:
        raise ValueError("replicate attempt grid is incomplete")
    expected_replicate_rows = scenarios * len(sample_sizes) * repetitions * len(METHODS)
    if checked != expected_replicate_rows or int(manifest.get("replicate_rows", -1)) != expected_replicate_rows:
        raise ValueError("replicate row count differs from the complete prespecified design")

    with summary_path.open(newline="", encoding="utf-8") as handle:
        summaries = list(csv.DictReader(handle))
    expected_summary_keys = {
        (s, n, method) for s in range(scenarios) for n in sample_sizes for method in METHODS
    }
    observed_summary_keys = {
        (int(row["scenario"]), int(row["sample_size"]), row["method"]) for row in summaries
    }
    if observed_summary_keys != expected_summary_keys:
        raise ValueError("summary cell grid is incomplete or contains unexpected cells")
    if len(summaries) != len(expected_summary_keys) or int(manifest.get("rows", -1)) != len(expected_summary_keys):
        raise ValueError("summary row count mismatch")

    for row in summaries:
        key = (int(row["scenario"]), int(row["sample_size"]), row["method"])
        records = grouped[key]
        successes = [r for r in records if not r["failed"]]
        failures = len(records) - len(successes)
        coverage = float(np.mean([r["covered"] for r in successes])) if successes else float("nan")
        width = float(np.mean([r["width"] for r in successes])) if successes else float("nan")
        if int(row["repetitions"]) != repetitions or int(row["successful_intervals"]) != len(successes) or int(row["failed_intervals"]) != failures:
            raise ValueError("summary attempt accounting mismatch")
        if not np.isclose(float(row["sampling_fraction"]), key[1] / population_size, rtol=0, atol=1e-14):
            raise ValueError("summary sampling fraction mismatch")
        if not np.isclose(float(row["target_finite_population_mean"]), targets[key[0]], rtol=0, atol=1e-12):
            raise ValueError("summary finite-frame target mismatch")
        if not np.isclose(float(row["coverage"]), coverage, rtol=1e-11, atol=1e-12, equal_nan=True):
            raise ValueError("summary coverage mismatch")
        if not np.isclose(float(row["mean_interval_width"]), width, rtol=1e-11, atol=1e-12, equal_nan=True):
            raise ValueError("summary interval-width mismatch")
        mcse = float(np.sqrt(coverage * (1.0 - coverage) / len(successes))) if successes else float("nan")
        if not np.isclose(float(row["coverage_mcse_conditional_on_success"]), mcse, rtol=1e-11, atol=1e-12, equal_nan=True):
            raise ValueError("summary coverage Monte Carlo SE mismatch")

    return {
        "status": "PASS",
        "replicates_independently_reconstructed": checked,
        "finite_frames_verified": len(populations),
        "summary_cells_verified": len(summaries),
        "design_completeness": "PASS",
        "frame_seed_reconstruction": "PASS",
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
