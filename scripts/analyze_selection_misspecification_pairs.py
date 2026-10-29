"""Independent paired-contrast analysis for saved selection-misspecification artifacts.

Run the primary artifact auditor before this supplementary post-run analysis.
All pairwise contrasts are exploratory secondary analyses, except the primary
contrast already locked and audited by the simulation runner.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ARMS = (
    "analytic_truth", "oracle", "catalogue_only", "correct_weighting",
    "shifted_transition", "shallow_slope",
)
PAIRS = (
    ("correct_weighting", "catalogue_only"),
    ("oracle", "analytic_truth"),
    ("correct_weighting", "analytic_truth"),
    ("shifted_transition", "correct_weighting"),
    ("shallow_slope", "correct_weighting"),
    ("shallow_slope", "analytic_truth"),
    ("shifted_transition", "analytic_truth"),
)
REQUIRED = {
    "replicate", "seed", "arm", "h0_true", "posterior_mean",
    "ci68_width", "covered_68", "density_ise", "n_parent", "n_catalogue_observed", "n_detected",
}


def _num(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _bool(value):
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError("covered_68 must be boolean")


def _wilson(successes, trials):
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid Wilson interval counts")
    z = 1.959963984540054
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / denominator
    return [max(0.0, center - half), min(1.0, center + half)]


def _bootstrap_stats(values, *, seed, n_boot):
    """Deterministic percentile bootstrap using bounded temporary arrays."""
    import numpy as np

    x = np.asarray(values, dtype=float)
    if x.ndim != 1 or x.size < 1 or not np.all(np.isfinite(x)):
        raise ValueError("paired differences must be a nonempty finite vector")
    if isinstance(n_boot, bool) or not isinstance(n_boot, int) or n_boot < 1000:
        raise ValueError("n_boot must be an integer of at least 1000")
    mean = float(np.mean(x))
    if x.size < 2:
        sd = se = normal_low = normal_high = None
    else:
        sd = float(np.std(x, ddof=1))
        se = sd / math.sqrt(x.size)
        normal_low, normal_high = mean - 1.96 * se, mean + 1.96 * se
    rng = np.random.default_rng(seed)
    boot_means = np.empty(n_boot, dtype=float)
    for start in range(0, n_boot, 500):
        stop = min(start + 500, n_boot)
        indices = rng.integers(0, x.size, size=(stop - start, x.size))
        boot_means[start:stop] = np.mean(x[indices], axis=1)
    low, high = (float(v) for v in np.quantile(boot_means, [0.025, 0.975]))
    return {
        "mean": mean,
        "standard_deviation": sd,
        "standard_error": se,
        "normal_95_interval": None if normal_low is None else [normal_low, normal_high],
        "normal_95_low": normal_low,
        "normal_95_high": normal_high,
        "bootstrap_95_interval": [low, high],
        "bootstrap_95_low": low,
        "bootstrap_95_high": high,
        "bootstrap_resamples": int(n_boot),
        "analysis_seed": int(seed),
    }


def analyze_records(rows, *, n_replicates, seed_start, h0_true, n_boot=10_000,
                    contrast_seed_offset=223607):
    """Analyze validated paired CSV rows without relying on summary aggregates."""
    import numpy as np

    if isinstance(n_replicates, bool) or not isinstance(n_replicates, int) or n_replicates < 2:
        raise ValueError("n_replicates must be an integer of at least two")
    if len(rows) != n_replicates * len(ARMS):
        raise ValueError("raw record count does not match replicates times six arms")
    by_key = {}
    for row in rows:
        if not REQUIRED.issubset(row):
            raise ValueError("raw rows are missing required columns")
        rep, seed, arm = int(row["replicate"]), int(row["seed"]), row["arm"]
        key = (rep, arm)
        if key in by_key:
            raise ValueError(f"duplicate replicate-arm record: {key}")
        if not 0 <= rep < n_replicates or arm not in ARMS:
            raise ValueError(f"unexpected replicate-arm key: {key}")
        if seed != seed_start + rep:
            raise ValueError(f"seed mismatch at {key}")
        if not math.isclose(_num(row["h0_true"], "h0_true"), h0_true, rel_tol=0, abs_tol=1e-10):
            raise ValueError(f"true H0 mismatch at {key}")
        _num(row["posterior_mean"], "posterior_mean")
        for column in ("ci68_width", "density_ise"):
            if _num(row[column], column) < 0:
                raise ValueError(f"{column} must be non-negative at {key}")
        if _num(row["ci68_width"], "ci68_width") <= 0:
            raise ValueError(f"ci68_width must be positive at {key}")
        by_key[key] = row
    for rep in range(n_replicates):
        if {(r, arm) for r, arm in by_key if r == rep} != {(rep, arm) for arm in ARMS}:
            raise ValueError(f"replicate {rep} is not complete across all six arms")
        for column in ("seed", "h0_true", "n_parent", "n_catalogue_observed", "n_detected"):
            vals = {str(by_key[(rep, arm)].get(column)) for arm in ARMS}
            if len(vals) != 1:
                raise ValueError(f"paired rows disagree on {column} for replicate {rep}")

    errors = {
        arm: np.asarray([
            _num(by_key[(rep, arm)]["posterior_mean"], "posterior_mean") - h0_true
            for rep in range(n_replicates)
        ], dtype=float)
        for arm in ARMS
    }
    coverage = {
        arm: np.asarray([_bool(by_key[(rep, arm)]["covered_68"]) for rep in range(n_replicates)])
        for arm in ARMS
    }
    arm_summary = {}
    for arm in ARMS:
        err = errors[arm]
        hits = int(coverage[arm].sum())
        arm_summary[arm] = {
            "n_replicates": n_replicates,
            "bias": float(err.mean()),
            "bias_standard_error": float(err.std(ddof=1) / math.sqrt(n_replicates)),
            "rmse": float(np.sqrt(np.mean(err ** 2))),
            "coverage_68": hits / n_replicates,
            "coverage_68_numerator": hits,
            "coverage_68_wilson_95_interval": _wilson(hits, n_replicates),
            "mean_68_interval_width": float(np.mean([
                float(by_key[(rep, arm)]["ci68_width"]) for rep in range(n_replicates)
            ])),
            "mean_density_ise": float(np.mean([
                float(by_key[(rep, arm)]["density_ise"]) for rep in range(n_replicates)
            ])),
        }

    # Recompute the locked primary contrast and use the same seed convention as
    # the runner/auditor; positive values indicate error reduction by weighting.
    primary_values = np.abs(errors["catalogue_only"]) - np.abs(errors["correct_weighting"])
    primary = _bootstrap_stats(
        primary_values, seed=int(seed_start) + 104729, n_boot=n_boot
    )
    primary.update({
        "definition": "abs(error_catalogue_only) - abs(error_correct_weighting)",
        "positive_favors": "correct_weighting",
        "n_replicates": n_replicates,
    })

    # These are supplemental post-run comparisons. A negative difference means
    # the first-listed arm has smaller absolute posterior-mean error.
    rng = np.random.default_rng(int(seed_start) + int(contrast_seed_offset))
    contrasts = []
    for first, second in PAIRS:
        difference = np.abs(errors[first]) - np.abs(errors[second])
        boot = np.empty(n_boot, dtype=float)
        for start in range(0, n_boot, 500):
            stop = min(start + 500, n_boot)
            indices = rng.integers(0, n_replicates, size=(stop - start, n_replicates))
            boot[start:stop] = np.mean(difference[indices], axis=1)
        lower, upper = (float(v) for v in np.quantile(boot, [0.025, 0.975]))
        wins = int(np.sum(np.abs(errors[first]) < np.abs(errors[second]))))
        contrasts.append({
            "first_arm": first,
            "second_arm": second,
            "estimand": "abs_error_first_minus_abs_error_second; negative favors first arm",
            "mean_difference": float(difference.mean()),
            "standard_error": float(difference.std(ddof=1) / math.sqrt(n_replicates)),
            "normal_95_interval": [
                float(difference.mean() - 1.96 * difference.std(ddof=1) / math.sqrt(n_replicates)),
                float(difference.mean() + 1.96 * difference.std(ddof=1) / math.sqrt(n_replicates)),
            ],
            "paired_bootstrap_95_interval": [lower, upper],
            "first_arm_lower_absolute_error_count": wins,
            "first_arm_lower_absolute_error_fraction": wins / n_replicates,
            "first_arm_win_fraction_wilson_95_interval": _wilson(wins, n_replicates),
            "bootstrap_resamples": int(n_boot),
            "analysis_seed_start": int(seed_start) + int(contrast_seed_offset),
        })
    return {
        "scope": "supplementary analysis of synthetic fixed-H0 selection-misspecification experiment",
        "n_replicates": n_replicates,
        "seed_start": int(seed_start),
        "seed_end": int(seed_start) + n_replicates - 1,
        "h0_true": float(h0_true),
        "arm_summary": arm_summary,
        "primary_estimand": primary,
        "secondary_paired_absolute_error_contrasts": contrasts,
        "interpretation_limit": (
            "These are repeated-sampling results for the declared toy model and one fixed H0. "
            "They do not validate real detector selection or a real H0 measurement."
        ),
    }


def analyze_artifact(output_dir, *, n_boot=10_000, contrast_seed_offset=223607):
    """Verify manifest hashes, analyze raw records, and cross-check the locked primary."""
    folder = Path(output_dir)
    raw_path, summary_path, manifest_path = (
        folder / "replicate_metrics.csv", folder / "summary.json", folder / "manifest.json"
    )
    if not all(path.is_file() for path in (raw_path, summary_path, manifest_path)):
        raise ValueError("replicate CSV, summary JSON, and manifest JSON are required")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for path in (raw_path, summary_path):
        expected = manifest.get("files", {}).get(path.name, {}).get("sha256")
        observed = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected != observed:
            raise ValueError(f"SHA-256 mismatch for {path.name}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "synthetic_selection_misspecification_pilot_not_real_H0_measurement":
        raise ValueError("unexpected scientific scope label")
    if manifest.get("git_commit") != summary.get("git_commit", manifest.get("git_commit")):
        raise ValueError("manifest/summary commit mismatch")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not REQUIRED.issubset(reader.fieldnames):
            raise ValueError("raw CSV is missing required columns")
        rows = list(reader)
    result = analyze_records(
        rows, n_replicates=summary["n_replicates"],
        seed_start=summary["seed_start"], h0_true=summary["h0_true_km_s_mpc"],
        n_boot=n_boot, contrast_seed_offset=contrast_seed_offset,
    )
    saved = summary.get("primary_estimand")
    if not isinstance(saved, dict):
        raise ValueError("locked primary estimand missing from saved summary")
    for key in ("mean", "standard_error", "normal_95_low", "normal_95_high",
                "bootstrap_95_low", "bootstrap_95_high", "bootstrap_resamples", "analysis_seed"):
        actual = result["primary_estimand"].get(key)
        expected = saved.get(key)
        if actual is None or expected is None:
            if actual != expected:
                raise ValueError(f"saved primary estimand mismatch for {key}")
        elif not math.isclose(float(actual), float(expected), rel_tol=0, abs_tol=1e-7):
            raise ValueError(f"saved primary estimand mismatch for {key}")
    result["artifact_sha256"] = {
        "replicate_metrics.csv": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
        "summary.json": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
    }
    result["locked_primary_verified_against_saved_summary"] = True
    result["source_commit"] = manifest.get("git_commit")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path, help="Unpacked workflow-artifact directory")
    parser.add_argument("--n-bootstrap", type=int, default=10_000)
    parser.add_argument("--contrast-seed-offset", type=int, default=223607)
    parser.add_argument("--output", type=Path, default=None, help="Optional JSON output path")
    args = parser.parse_args()
    result = analyze_artifact(
        args.output_dir, n_boot=args.n_bootstrap, contrast_seed_offset=args.contrast_seed_offset
    )
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output is None:
        print(rendered, end="")
    else:
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_text(rendered, encoding="utf-8")
        print(json.dumps({"status": "PASS", "output": str(args.output)}, indent=2))


if __name__ == "__main__":
    main()
