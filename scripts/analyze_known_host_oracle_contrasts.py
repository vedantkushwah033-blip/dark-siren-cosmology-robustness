"""Paired comparison of full-population and known-true-host oracle posteriors."""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


SHARED_INPUTS = ("seed", "host_redshift", "observed_distance_mpc", "distance_sigma_mpc")
POSTERIOR_FIELDS = ("median", "mean", "posterior_sd", "lower_68", "upper_68", "lower_95", "upper_95")


def load_rows(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise ValueError(f"Missing or empty input: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No records in {path}")
    return rows


def paired_delta(full_rows: list[dict[str, str]], known_rows: list[dict[str, str]]) -> dict:
    """Validate event-level pairing and return per-run metrics (known-host minus full-population)."""
    if len(full_rows) != len(known_rows):
        raise ValueError("oracle input groups have different replicate counts")
    full = {int(row["run"]): row for row in full_rows}
    known = {int(row["run"]): row for row in known_rows}
    if len(full) != len(full_rows) or len(known) != len(known_rows):
        raise ValueError("run indices must be unique")
    if set(full) != set(known):
        raise ValueError("oracle inputs have different run-index sets")

    for run in sorted(full):
        left, right = full[run], known[run]
        for field in SHARED_INPUTS:
            if field not in left or field not in right:
                raise ValueError(f"missing shared event field {field} at run {run}")
            a, b = float(left[field]), float(right[field])
            if not math.isfinite(a) or not math.isfinite(b) or not math.isclose(a, b, rel_tol=0.0, abs_tol=1e-12):
                raise ValueError(f"shared event input {field} differs at run {run}: {a} vs {b}")

    def metrics(row: dict[str, str]) -> dict[str, float]:
        median = float(row["median"])
        lo68, hi68 = float(row["lower_68"]), float(row["upper_68"])
        lo95, hi95 = float(row["lower_95"]), float(row["upper_95"])
        if not all(math.isfinite(x) for x in (median, lo68, hi68, lo95, hi95)):
            raise ValueError("non-finite posterior summary in paired oracle inputs")
        if lo68 > hi68 or lo95 > hi95 or lo95 > lo68 or hi68 > hi95:
            raise ValueError("invalid or non-nested credible intervals in paired oracle inputs")
        error = median - 70.0
        return {
            "coverage_68": float(lo68 <= 70.0 <= hi68),
            "coverage_95": float(lo95 <= 70.0 <= hi95),
            "width_68": hi68 - lo68,
            "width_95": hi95 - lo95,
            "median_error": error,
            "absolute_median_error": abs(error),
            "squared_median_error": error * error,
        }

    result = {}
    for metric in metrics(known[sorted(known)[0]]):
        values = []
        for run in sorted(full):
            a = metrics(full[run])[metric]
            b = metrics(known[run])[metric]
            values.append(b - a)
        result[metric] = values
    return result


def paired_bootstrap_mean_difference(values: list[float], repetitions: int, seed: int) -> dict:
    """Bootstrap paired per-run differences using vectorized batches."""
    import numpy as np

    if not values:
        raise ValueError("paired differences must not be empty")
    if repetitions < 1:
        raise ValueError("bootstrap repetitions must be positive")
    array = np.asarray(values, dtype=float)
    if not np.isfinite(array).all():
        raise ValueError("paired differences must all be finite")
    rng = np.random.default_rng(seed)
    n = len(array)
    observed = float(array.mean())
    draws = np.empty(repetitions, dtype=float)
    batch_size = max(1, min(1000, repetitions))
    cursor = 0
    while cursor < repetitions:
        size = min(batch_size, repetitions - cursor)
        indices = rng.integers(0, n, size=(size, n))
        draws[cursor:cursor + size] = array[indices].mean(axis=1)
        cursor += size
    return {
        "difference_known_host_minus_full_population": observed,
        "paired_bootstrap_95_interval": [float(x) for x in np.quantile(draws, [0.025, 0.975])],
        "n_pairs": n,
        "bootstrap_replicates": repetitions,
        "bootstrap_seed": seed,
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--full-population-input", default="results/full-population-oracle/per_run_results.csv")
    parser.add_argument("--known-host-input", default="results/known-host-oracle/per_run_results.csv")
    parser.add_argument("--output", default="results/known-host-oracle/paired_contrasts.json")
    parser.add_argument("--bootstrap", type=int, default=50000)
    parser.add_argument("--seed", type=int, default=20261025)
    args = parser.parse_args()
    if args.bootstrap < 100:
        parser.error("--bootstrap must be at least 100")

    full_rows = load_rows(Path(args.full_population_input))
    # The full-population artifact has two masks; the posterior is mask-invariant.
    full_rows = [row for row in full_rows if row.get("scenario") == "random70_full_population_oracle"]
    known_rows = load_rows(Path(args.known_host_input))
    if any(float(row.get("H0_true", "nan")) != 70.0 for row in full_rows + known_rows):
        raise SystemExit("Paired oracle comparison requires true H0=70 in all records")
    deltas = paired_delta(full_rows, known_rows)
    output = {
        "true_h0_km_s_mpc": 70.0,
        "comparison": "known true host redshift minus full-population oracle",
        "warning": "Synthetic paired calibration diagnostic only. Both oracle methods use hidden simulation truth; neither is directly available for a real dark siren with an unknown host.",
        "bootstrap_method": "paired replicate bootstrap of per-run metric differences",
        "bootstrap_replicates": args.bootstrap,
        "bootstrap_seed": args.seed,
        "metrics": {},
    }
    for index, (metric, values) in enumerate(deltas.items()):
        result = paired_bootstrap_mean_difference(values, args.bootstrap, args.seed + index)
        result["full_population_mean"] = float(sum(
            (float(row["lower_68"]) <= 70.0 <= float(row["upper_68"])) if metric == "coverage_68" else
            (float(row["lower_95"]) <= 70.0 <= float(row["upper_95"])) if metric == "coverage_95" else
            float(row["upper_68"]) - float(row["lower_68"]) if metric == "width_68" else
            float(row["upper_95"]) - float(row["lower_95"]) if metric == "width_95" else
            float(row["median"]) - 70.0 if metric == "median_error" else
            abs(float(row["median"]) - 70.0) if metric == "absolute_median_error" else
            (float(row["median"]) - 70.0) ** 2
            for row in full_rows
        ) / len(full_rows))
        output["metrics"][metric] = result

    out = Path(args.output)
    out.parent.mkdir(parents=True, exist_ok=True)
    out.write_text(json.dumps(output, indent=2), encoding="utf-8")
    print(f"Wrote paired oracle calibration contrasts for {len(full_rows)} matched replicates to {out}")


if __name__ == "__main__":
    main()
