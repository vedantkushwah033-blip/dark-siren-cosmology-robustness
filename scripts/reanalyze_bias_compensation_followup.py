"""Independent paired reanalysis for the synthetic bias-compensation follow-up.

The input CSV must contain replicate, arm, h0_true, and posterior_mean columns.
Each replicate must contain exactly one row for each requested arm. The primary
contrast is Delta_RMSE = RMSE(correct arm) - RMSE(misspecified arm); positive
values mean the misspecified arm has lower RMSE in that cell.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path
from typing import Iterable


def _finite(value: object, label: str) -> float:
    try:
        number = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(number):
        raise ValueError(f"{label} must be finite")
    return number


def _rmse(errors: list[float]) -> float:
    return math.sqrt(sum(x * x for x in errors) / len(errors))


def _variance(errors: list[float]) -> float | None:
    if len(errors) < 2:
        return None
    mean = sum(errors) / len(errors)
    return sum((x - mean) ** 2 for x in errors) / (len(errors) - 1)


def analyze_rows(
    rows: Iterable[dict[str, str]], *, true_h0: float | None = None,
    correct_arm: str = "correct_weighting", misspecified_arm: str = "shallow_slope",
    seed: int = 20261029, n_boot: int = 10_000,
) -> dict:
    import numpy as np

    if correct_arm == misspecified_arm:
        raise ValueError("correct_arm and misspecified_arm must differ")
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 1000:
        raise ValueError("n_boot must be an integer of at least 1000")
    if isinstance(seed, bool) or not isinstance(seed, int):
        raise ValueError("seed must be an integer")

    by_rep: dict[str, dict[str, tuple[float, float]]] = {}
    for line_number, row in enumerate(rows, start=2):
        required = {"replicate", "arm", "h0_true", "posterior_mean"}
        if not required.issubset(row):
            raise ValueError("CSV is missing one or more required columns: replicate, arm, h0_true, posterior_mean")
        rep = str(row["replicate"]).strip()
        arm = str(row["arm"]).strip()
        if not rep:
            raise ValueError(f"empty replicate key on CSV line {line_number}")
        if arm not in {correct_arm, misspecified_arm}:
            continue
        if arm in by_rep.setdefault(rep, {}):
            raise ValueError(f"duplicate row for replicate={rep}, arm={arm}")
        truth = _finite(row["h0_true"], "h0_true")
        mean = _finite(row["posterior_mean"], "posterior_mean")
        if true_h0 is not None and not math.isclose(truth, true_h0, rel_tol=0.0, abs_tol=1e-10):
            raise ValueError(f"h0_true mismatch for replicate={rep}: {truth} != {true_h0}")
        by_rep[rep][arm] = (truth, mean)

    if not by_rep:
        raise ValueError("no rows found for the requested pair of arms")
    required_arms = {correct_arm, misspecified_arm}
    for rep, arm_rows in by_rep.items():
        if set(arm_rows) != required_arms:
            raise ValueError(f"replicate={rep} is incomplete for requested arms")
        truth_values = {arm_rows[a][0] for a in required_arms}
        if len(truth_values) != 1:
            raise ValueError(f"paired rows disagree on true H0 for replicate={rep}")

    # Stable deterministic ordering, so input CSV row order does not change results.
    replicate_ids = sorted(by_rep, key=lambda x: (not x.isdigit(), int(x) if x.isdigit() else x))
    truths = [by_rep[r][correct_arm][0] for r in replicate_ids]
    inferred_truth = truths[0]
    if any(not math.isclose(t, inferred_truth, rel_tol=0.0, abs_tol=1e-10) for t in truths):
        raise ValueError("true H0 varies across replicates; analyze each truth cell separately")
    if true_h0 is not None and not math.isclose(inferred_truth, true_h0, rel_tol=0.0, abs_tol=1e-10):
        raise ValueError("inferred true H0 does not match requested true_h0")

    err_correct = [by_rep[r][correct_arm][1] - inferred_truth for r in replicate_ids]
    err_miss = [by_rep[r][misspecified_arm][1] - inferred_truth for r in replicate_ids]
    n = len(replicate_ids)
    if n < 2:
        raise ValueError("at least two paired replicates are required")

    def arm_summary(errors: list[float]) -> dict:
        bias = sum(errors) / n
        variance = _variance(errors)
        rmse = _rmse(errors)
        return {
            "n_replicates": n,
            "bias": bias,
            "bias_squared": bias * bias,
            "sample_variance": variance,
            "rmse": rmse,
            "rmse_squared": rmse * rmse,
        }

    point = _rmse(err_correct) - _rmse(err_miss)
    rng = np.random.default_rng(seed)
    ec = np.asarray(err_correct, dtype=float)
    em = np.asarray(err_miss, dtype=float)
    boot = np.empty(n_boot, dtype=float)
    # Resample paired replicate indices, preserving within-replicate dependence.
    for start in range(0, n_boot, 500):
        stop = min(start + 500, n_boot)
        idx = rng.integers(0, n, size=(stop - start, n))
        boot[start:stop] = np.sqrt(np.mean(ec[idx] ** 2, axis=1)) - np.sqrt(np.mean(em[idx] ** 2, axis=1))
    low, high = (float(x) for x in np.quantile(boot, [0.025, 0.975]))
    return {
        "status": "synthetic_followup_reanalysis_not_real_data_validation",
        "estimand": "RMSE(correct weighting) - RMSE(shallow/misspecified arm)",
        "positive_delta_means": f"{misspecified_arm} has lower RMSE in this cell",
        "true_h0": inferred_truth,
        "n_paired_replicates": n,
        "replicate_ids": replicate_ids,
        "correct_arm": correct_arm,
        "misspecified_arm": misspecified_arm,
        "correct_arm_summary": arm_summary(err_correct),
        "misspecified_arm_summary": arm_summary(err_miss),
        "delta_rmse": point,
        "paired_bootstrap_95_interval": [low, high],
        "bootstrap_resamples": n_boot,
        "bootstrap_seed": seed,
        "interpretation": (
            "An interval spanning zero is inconclusive. A positive contrast in one truth cell "
            "does not establish general superiority of a misspecified model."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("csv", type=Path, help="raw replicate-level CSV")
    parser.add_argument("--true-h0", type=float, default=None)
    parser.add_argument("--correct-arm", default="correct_weighting")
    parser.add_argument("--misspecified-arm", default="shallow_slope")
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--n-boot", type=int, default=10_000)
    parser.add_argument("--out", type=Path, default=None)
    args = parser.parse_args()
    with args.csv.open(newline="", encoding="utf-8") as stream:
        result = analyze_rows(csv.DictReader(stream), true_h0=args.true_h0,
                              correct_arm=args.correct_arm,
                              misspecified_arm=args.misspecified_arm,
                              seed=args.seed, n_boot=args.n_boot)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.out:
        args.out.parent.mkdir(parents=True, exist_ok=True)
        args.out.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")


if __name__ == "__main__":
    main()
