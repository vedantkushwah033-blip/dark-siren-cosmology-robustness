"""Reconcile fixed-truth coverage by separating prior-boundary and interior truths.

This is an independent summary utility for the auxiliary known-host distance
uncertainty experiment. It does not import the experiment runner or likelihood
implementation and does not turn the fixed-truth grid into prior-predictive SBC.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np

ARMS = (
    "A_fixed_realized_sigma",
    "B_fractional_sigma_H0",
    "C_fixed_reference_sigma",
)
EXPECTED_TRUTHS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
EXPECTED_FRACTIONS = (0.05, 0.10, 0.15, 0.20)
REQUIRED = {
    "true_h0", "fractional_error", "replicate", "arm",
    "lower_68", "upper_68", "lower_95", "upper_95",
}


def _wilson(successes: int, n: int) -> list[float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("Wilson interval requires n > 0 and 0 <= successes <= n")
    z = 1.959963984540054
    p = successes / n
    den = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / den
    half = z * np.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / den
    return [float(max(0.0, center - half)), float(min(1.0, center + half))]


def _rate(values: np.ndarray) -> dict:
    successes = int(np.sum(values))
    n = int(values.size)
    return {
        "successes": successes,
        "n": n,
        "coverage": successes / n if n else None,
        "wilson_95": _wilson(successes, n) if n else None,
    }


def _stratified_bootstrap(
    cells: list[np.ndarray], *, seed: int, resamples: int
) -> list[float]:
    if not cells or resamples < 100:
        raise ValueError("bootstrap requires cells and at least 100 resamples")
    sizes = {int(cell.size) for cell in cells}
    if len(sizes) != 1 or 0 in sizes:
        raise ValueError("all bootstrap strata must have the same positive size")
    n = next(iter(sizes))
    rng = np.random.default_rng(seed)
    means = np.zeros(resamples, dtype=float)
    for values in cells:
        indices = rng.integers(0, n, size=(resamples, n))
        means += values[indices].mean(axis=1) / len(cells)
    # Coverage rates are probabilities. Floating-point summation can put a
    # degenerate bootstrap endpoint a few ulps outside [0, 1] (for example,
    # 0.9999999999999997 when every resample has coverage exactly one).
    return [
        float(np.clip(np.quantile(means, 0.025), 0.0, 1.0)),
        float(np.clip(np.quantile(means, 0.975), 0.0, 1.0)),
    ]


def reconcile(
    raw_csv: str | Path,
    *,
    interior_min: float = 55.0,
    interior_max: float = 85.0,
    bootstrap_seed: int = 20261121,
    bootstrap_resamples: int = 10000,
) -> dict:
    path = Path(raw_csv)
    raw_bytes = path.read_bytes()
    with path.open(newline="", encoding="utf-8-sig") as stream:
        reader = csv.DictReader(stream)
        if not REQUIRED.issubset(reader.fieldnames or []):
            raise ValueError(f"raw CSV is missing columns: {sorted(REQUIRED - set(reader.fieldnames or []))}")
        rows = list(reader)
    if not rows:
        raise ValueError("raw CSV is empty")
    if not np.isfinite(interior_min) or not np.isfinite(interior_max) or interior_min >= interior_max:
        raise ValueError("interior bounds must be finite and increasing")
    if bootstrap_seed < 0 or bootstrap_resamples < 100:
        raise ValueError("bootstrap seed must be nonnegative and resamples >= 100")

    indexed = {}
    truths, fractions, arms, replicates = set(), set(), set(), set()
    for line, row in enumerate(rows, start=2):
        try:
            truth = float(row["true_h0"])
            fraction = float(row["fractional_error"])
            replicate = int(row["replicate"])
            arm = row["arm"]
            bounds = {
                68: (float(row["lower_68"]), float(row["upper_68"])),
                95: (float(row["lower_95"]), float(row["upper_95"])),
            }
        except (TypeError, ValueError) as exc:
            raise ValueError(f"invalid row at CSV line {line}") from exc
        numeric = [truth, fraction, *bounds[68], *bounds[95]]
        if not all(np.isfinite(value) for value in numeric) or replicate < 0:
            raise ValueError(f"non-finite or invalid value at CSV line {line}")
        if arm not in ARMS or fraction <= 0:
            raise ValueError(f"unexpected arm or fractional error at CSV line {line}")
        for level, (lower, upper) in bounds.items():
            if lower > upper:
                raise ValueError(f"reversed {level}% interval at CSV line {line}")
        if bounds[95][0] > bounds[68][0] or bounds[95][1] < bounds[68][1]:
            raise ValueError(f"95% interval does not contain 68% interval at CSV line {line}")
        key = (truth, fraction, replicate, arm)
        if key in indexed:
            raise ValueError(f"duplicate raw key at CSV line {line}: {key}")
        indexed[key] = {
            "covered_68": int(bounds[68][0] <= truth <= bounds[68][1]),
            "covered_95": int(bounds[95][0] <= truth <= bounds[95][1]),
        }
        truths.add(truth)
        fractions.add(fraction)
        arms.add(arm)
        replicates.add(replicate)

    if tuple(sorted(truths)) != EXPECTED_TRUTHS:
        raise ValueError("truth grid differs from the frozen nine-point design")
    if tuple(sorted(fractions)) != EXPECTED_FRACTIONS:
        raise ValueError("fractional-error grid differs from the frozen four-point design")
    if arms != set(ARMS):
        raise ValueError("raw CSV does not contain exactly the three expected arms")
    expected_replicates = set(range(max(replicates) + 1))
    if replicates != expected_replicates:
        raise ValueError("replicate IDs must be contiguous from zero")

    expected_keys = {
        (truth, fraction, replicate, arm)
        for truth in truths for fraction in fractions
        for replicate in replicates for arm in ARMS
    }
    if set(indexed) != expected_keys:
        raise ValueError(
            f"raw design incomplete: missing={len(expected_keys-set(indexed))}, "
            f"extra={len(set(indexed)-expected_keys)}"
        )

    interior_truths = tuple(
        truth for truth in sorted(truths) if interior_min <= truth <= interior_max
    )
    if not interior_truths:
        raise ValueError("interior truth selection is empty")
    cells_per_arm = len(interior_truths) * len(fractions)
    results = {}
    for arm in ARMS:
        all_values = {
            level: np.asarray([
                indexed[(truth, fraction, replicate, arm)][f"covered_{level}"]
                for truth in sorted(truths) for fraction in sorted(fractions)
                for replicate in sorted(replicates)
            ], dtype=int)
            for level in (68, 95)
        }
        interior_values = {
            level: np.asarray([
                indexed[(truth, fraction, replicate, arm)][f"covered_{level}"]
                for truth in interior_truths for fraction in sorted(fractions)
                for replicate in sorted(replicates)
            ], dtype=int)
            for level in (68, 95)
        }
        item = {
            "all_nine_truth_grid": {
                str(level): _rate(all_values[level]) for level in (68, 95)
            },
            "interior_truth_grid": {
                str(level): _rate(interior_values[level]) for level in (68, 95)
            },
            "by_truth": {},
            "by_fraction_interior": {},
            "stratified_bootstrap": {},
        }
        for truth in sorted(truths):
            item["by_truth"][str(int(truth))] = {
                str(level): _rate(np.asarray([
                    indexed[(truth, fraction, replicate, arm)][f"covered_{level}"]
                    for fraction in sorted(fractions) for replicate in sorted(replicates)
                ], dtype=int))
                for level in (68, 95)
            }
        for fraction in sorted(fractions):
            item["by_fraction_interior"][str(fraction)] = {
                str(level): _rate(np.asarray([
                    indexed[(truth, fraction, replicate, arm)][f"covered_{level}"]
                    for truth in interior_truths for replicate in sorted(replicates)
                ], dtype=int))
                for level in (68, 95)
            }
        for level in (68, 95):
            cells = [
                np.asarray([
                    indexed[(truth, fraction, replicate, arm)][f"covered_{level}"]
                    for replicate in sorted(replicates)
                ], dtype=float)
                for truth in interior_truths for fraction in sorted(fractions)
            ]
            item["stratified_bootstrap"][str(level)] = {
                "estimand": "equal-weight mean of fixed truth x fractional-error cell coverage rates",
                "truth_cells": list(interior_truths),
                "fractional_errors": list(sorted(fractions)),
                "strata": cells_per_arm,
                "replicates_per_stratum": len(replicates),
                "resamples": bootstrap_resamples,
                "seed": bootstrap_seed + level,
                "percentile_95_interval": _stratified_bootstrap(
                    cells, seed=bootstrap_seed + level, resamples=bootstrap_resamples
                ),
            }
        results[arm] = item

    return {
        "status": "PASS",
        "scope": "auxiliary synthetic known-host low-redshift fixed-truth coverage reconciliation only",
        "not_prior_predictive_calibration": True,
        "raw_record_count": len(rows),
        "raw_csv_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "truth_grid": list(sorted(truths)),
        "fractional_error_grid": list(sorted(fractions)),
        "replicates_per_truth_fraction_arm": len(replicates),
        "prior_boundary_truths": [min(truths), max(truths)],
        "interior_truths": list(interior_truths),
        "results": results,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("raw_csv", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    parser.add_argument("--interior-min", type=float, default=55.0)
    parser.add_argument("--interior-max", type=float, default=85.0)
    parser.add_argument("--bootstrap-seed", type=int, default=20261121)
    parser.add_argument("--bootstrap-resamples", type=int, default=10000)
    args = parser.parse_args()
    report = reconcile(
        args.raw_csv,
        interior_min=args.interior_min,
        interior_max=args.interior_max,
        bootstrap_seed=args.bootstrap_seed,
        bootstrap_resamples=args.bootstrap_resamples,
    )
    rendered = json.dumps(report, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "raw_record_count": report["raw_record_count"],
        "raw_csv_sha256": report["raw_csv_sha256"],
        "interior_truths": report["interior_truths"],
        "arms": {
            arm: {
                "coverage_68": result["interior_truth_grid"]["68"]["coverage"],
                "coverage_95": result["interior_truth_grid"]["95"]["coverage"],
                "coverage_68_bootstrap_95": result["stratified_bootstrap"]["68"]["percentile_95_interval"],
                "coverage_95_bootstrap_95": result["stratified_bootstrap"]["95"]["percentile_95_interval"],
            }
            for arm, result in report["results"].items()
        },
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
