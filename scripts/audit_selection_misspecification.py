"""Independent consistency audit for selection-misspecification pilot artifacts."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ARMS = ("analytic_truth", "oracle", "catalogue_only", "correct_weighting", "shifted_transition", "shallow_slope")
REQUIRED = {
    "replicate", "seed", "arm", "h0_true", "posterior_mean", "posterior_median",
    "ci68_low", "ci68_high", "ci68_width", "posterior_integral",
    "posterior_density_grid_json", "covered_68", "density_ise", "weight_ess",
    "weight_ess_fraction", "max_normalized_weight", "q95_normalized_weight",
    "min_assumed_inclusion_probability", "n_parent", "n_catalogue_observed",
}


def _num(value, label):
    try:
        x = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(x):
        raise ValueError(f"{label} must be finite")
    return x


def _close(actual, expected, label, atol=1e-6):
    if not math.isclose(_num(actual, label), expected, rel_tol=1e-7, abs_tol=atol):
        raise ValueError(f"{label} mismatch: {actual!r} != {expected!r}")


def _trap(values, grid):
    return sum(0.5 * (values[i] + values[i + 1]) * (grid[i + 1] - grid[i])
               for i in range(len(grid) - 1))


def _quantile(grid, cdf, q):
    for i in range(1, len(cdf)):
        if cdf[i] >= q:
            delta = cdf[i] - cdf[i - 1]
            return grid[i] if delta <= 0 else grid[i - 1] + (q - cdf[i - 1]) / delta * (grid[i] - grid[i - 1])
    return grid[-1]


def _posterior(row, grid):
    try:
        density = json.loads(row["posterior_density_grid_json"])
    except (TypeError, json.JSONDecodeError) as exc:
        raise ValueError("posterior density grid is invalid JSON") from exc
    if not isinstance(density, list) or len(density) != len(grid):
        raise ValueError("posterior density grid length mismatch")
    density = [_num(x, "posterior density") for x in density]
    if any(x < 0 for x in density):
        raise ValueError("posterior density cannot be negative")
    area = _trap(density, grid)
    if area <= 0:
        raise ValueError("posterior density has non-positive integral")
    density = [x / area for x in density]
    mean = _trap([x * p for x, p in zip(grid, density)], grid)
    cdf = [0.0]
    for i in range(len(grid) - 1):
        cdf.append(cdf[-1] + 0.5 * (density[i] + density[i + 1]) * (grid[i + 1] - grid[i]))
    cdf = [x / cdf[-1] for x in cdf]
    low, median, high = (_quantile(grid, cdf, q) for q in (0.16, 0.5, 0.84))
    return {"posterior_integral": area, "posterior_mean": mean, "posterior_median": median,
            "ci68_low": low, "ci68_high": high, "ci68_width": high - low}


def _wilson(successes, trials):
    z = 1.959963984540054
    p = successes / trials
    z2 = z * z
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half = z * math.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def _paired_summary(values):
    mean = sum(values) / len(values)
    if len(values) < 2:
        return {"mean": mean, "standard_error": None, "normal_95_low": None, "normal_95_high": None}
    variance = sum((x - mean) ** 2 for x in values) / (len(values) - 1)
    se = math.sqrt(variance / len(values))
    return {"mean": mean, "standard_error": se,
            "normal_95_low": mean - 1.96 * se, "normal_95_high": mean + 1.96 * se}


def audit_artifacts(output_dir: str | Path) -> dict:
    folder = Path(output_dir)
    raw_path, summary_path, manifest_path = (
        folder / "replicate_metrics.csv", folder / "summary.json", folder / "manifest.json"
    )
    if not all(p.is_file() for p in (raw_path, summary_path, manifest_path)):
        raise ValueError("raw CSV, summary JSON, and manifest JSON are all required")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    for path in (raw_path, summary_path):
        expected_hash = manifest.get("files", {}).get(path.name, {}).get("sha256")
        actual_hash = hashlib.sha256(path.read_bytes()).hexdigest()
        if expected_hash != actual_hash:
            raise ValueError(f"SHA-256 mismatch for {path.name}")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    if summary.get("status") != "synthetic_selection_misspecification_pilot_not_real_H0_measurement":
        raise ValueError("unexpected or missing scientific scope label")
    n_replicates = summary.get("n_replicates")
    if isinstance(n_replicates, bool) or not isinstance(n_replicates, int) or n_replicates < 1:
        raise ValueError("invalid n_replicates")
    h0_min, h0_max, h0_points = summary["h0_grid"]
    if isinstance(h0_points, bool) or int(h0_points) < 31 or not 0 < float(h0_min) < float(h0_max):
        raise ValueError("invalid H0 grid")
    grid = [float(h0_min) + i * (float(h0_max) - float(h0_min)) / (int(h0_points) - 1)
            for i in range(int(h0_points))]
    h0_true = _num(summary["h0_true_km_s_mpc"], "true H0")
    if not float(h0_min) <= h0_true <= float(h0_max):
        raise ValueError("true H0 outside grid")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not REQUIRED.issubset(reader.fieldnames):
            raise ValueError("raw CSV missing required columns")
        rows = list(reader)
    if len(rows) != n_replicates * len(ARMS):
        raise ValueError("raw record count mismatch")
    by_key = {}
    for row in rows:
        try:
            rep, seed = int(row["replicate"]), int(row["seed"])
        except (TypeError, ValueError) as exc:
            raise ValueError("replicate and seed must be integers") from exc
        arm = row["arm"]
        key = (rep, arm)
        if key in by_key:
            raise ValueError(f"duplicate replicate-arm record: {key}")
        if not 0 <= rep < n_replicates or arm not in ARMS:
            raise ValueError(f"unexpected replicate-arm key: {key}")
        if seed != int(summary["seed_start"]) + rep:
            raise ValueError(f"seed mismatch at {key}")
        _close(row["h0_true"], h0_true, "row true H0")
        n_parent = int(row["n_parent"])
        n_catalogue = int(row["n_catalogue_observed"])
        if n_parent < 100 or not 1 <= n_catalogue <= n_parent:
            raise ValueError(f"invalid parent/catalogue counts at {key}")
        sample_size = n_parent if arm in {"analytic_truth", "oracle"} else n_catalogue
        ess = _num(row["weight_ess"], "weight_ess")
        ess_fraction = _num(row["weight_ess_fraction"], "weight_ess_fraction")
        max_weight = _num(row["max_normalized_weight"], "max_normalized_weight")
        q95_weight = _num(row["q95_normalized_weight"], "q95_normalized_weight")
        min_probability = _num(
            row["min_assumed_inclusion_probability"], "min_assumed_inclusion_probability"
        )
        if not 0 < ess <= sample_size + 1e-8:
            raise ValueError(f"weight ESS outside (0, sample size] at {key}")
        _close(ess_fraction, ess / sample_size, "weight ESS fraction", atol=1e-8)
        if not 0 < max_weight <= 1 or not 0 < q95_weight <= max_weight + 1e-9:
            raise ValueError(f"invalid normalized-weight diagnostics at {key}")
        if not 0 < min_probability <= 1:
            raise ValueError(f"minimum assumed inclusion probability outside (0, 1] at {key}")
        reconstructed = _posterior(row, grid)
        for metric, value in reconstructed.items():
            _close(row[metric], value, f"reconstructed {metric}")
        if not float(h0_min) <= reconstructed["ci68_low"] <= reconstructed["posterior_median"] <= reconstructed["ci68_high"] <= float(h0_max):
            raise ValueError(f"invalid interval ordering at {key}")
        if _bool(row["covered_68"]) != (reconstructed["ci68_low"] <= h0_true <= reconstructed["ci68_high"]):
            raise ValueError(f"coverage flag mismatch at {key}")
        for name in ("density_ise", "weight_ess", "weight_ess_fraction", "max_normalized_weight",
                     "q95_normalized_weight", "min_assumed_inclusion_probability"):
            _num(row[name], name)
        if _num(row["density_ise"], "density_ise") < 0:
            raise ValueError("density ISE cannot be negative")
        if not 0 < _num(row["weight_ess_fraction"], "weight_ess_fraction") <= 1 + 1e-9:
            raise ValueError("weight ESS fraction must lie in (0, 1]")
        by_key[key] = row
    for rep in range(n_replicates):
        if {arm for r, arm in by_key if r == rep} != set(ARMS):
            raise ValueError(f"replicate {rep} is not paired across all arms")
        paired = [by_key[(rep, arm)] for arm in ARMS]
        for column in ("seed", "h0_true", "n_parent", "n_catalogue_observed"):
            if len({str(row[column]) for row in paired}) != 1:
                raise ValueError(f"paired arms disagree on {column} at replicate {rep}")
    aggregates = summary.get("aggregate_metrics")
    aggregate_by_arm = {item.get("arm"): item for item in aggregates if isinstance(item, dict)}
    if set(aggregate_by_arm) != set(ARMS):
        raise ValueError("aggregate arm set mismatch")
    for arm in ARMS:
        group = [by_key[(rep, arm)] for rep in range(n_replicates)]
        errors = [_num(row["posterior_mean"], "posterior_mean") - h0_true for row in group]
        covered = [_bool(row["covered_68"]) for row in group]
        wilson_low, wilson_high = _wilson(sum(covered), n_replicates)
        expected = {
            "n_replicates": n_replicates,
            "mean_bias": sum(errors) / n_replicates,
            "rmse": math.sqrt(sum(x*x for x in errors) / n_replicates),
            "empirical_68_coverage": sum(covered) / n_replicates,
            "coverage_numerator": sum(covered),
            "coverage_denominator": n_replicates,
            "coverage_wilson_95_low": wilson_low,
            "coverage_wilson_95_high": wilson_high,
            "mean_68_interval_width": sum(float(row["ci68_width"]) for row in group) / n_replicates,
            "mean_density_ise": sum(float(row["density_ise"]) for row in group) / n_replicates,
            "mean_weight_ess_fraction": sum(float(row["weight_ess_fraction"]) for row in group) / n_replicates,
            "mean_max_normalized_weight": sum(float(row["max_normalized_weight"]) for row in group) / n_replicates,
        }
        for metric, value in expected.items():
            if metric == "n_replicates" or metric in {"coverage_numerator", "coverage_denominator"}:
                if aggregate_by_arm[arm].get(metric) != value:
                    raise ValueError(f"{arm} aggregate {metric} mismatch")
            else:
                _close(aggregate_by_arm[arm].get(metric), value, f"{arm} aggregate {metric}")
    baseline = [by_key[(rep, "catalogue_only")] for rep in range(n_replicates)]
    reconstructed_contrasts = {}
    for arm in ARMS:
        if arm == "catalogue_only":
            continue
        comparison = [by_key[(rep, arm)] for rep in range(n_replicates)]
        squared_differences = [
            (_num(row["posterior_mean"], "posterior_mean") - h0_true) ** 2
            - (_num(base["posterior_mean"], "posterior_mean") - h0_true) ** 2
            for row, base in zip(comparison, baseline)
        ]
        density_differences = [
            _num(row["density_ise"], "density_ise") - _num(base["density_ise"], "density_ise")
            for row, base in zip(comparison, baseline)
        ]
        reconstructed_contrasts[arm] = {
            "squared_error_difference": _paired_summary(squared_differences),
            "density_ise_difference": _paired_summary(density_differences),
        }
    saved_contrasts = summary.get("paired_contrasts_vs_catalogue_only")
    if not isinstance(saved_contrasts, list) or len(saved_contrasts) != len(ARMS) - 1:
        raise ValueError("paired contrast summaries missing or incomplete")
    saved_by_arm = {item.get("arm"): item for item in saved_contrasts if isinstance(item, dict)}
    if set(saved_by_arm) != set(reconstructed_contrasts):
        raise ValueError("paired contrast arm labels mismatch")
    for arm, contrasts in reconstructed_contrasts.items():
        for metric, fields in contrasts.items():
            saved = saved_by_arm[arm].get(metric, {})
            for field, value in fields.items():
                if value is None:
                    if saved.get(field) is not None:
                        raise ValueError(f"{arm} paired {metric} {field} should be null for one replicate")
                else:
                    _close(saved.get(field), value, f"{arm} paired {metric} {field}", atol=1e-7)
    return {"status": "PASS", "raw_records": len(rows), "replicates": n_replicates,
            "posterior_grids_reconstructed": len(rows), "aggregate_arms_reconstructed": len(ARMS),
            "replicate_metrics_sha256": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            "summary_sha256": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
            "scope": "synthetic toy artifact consistency; not scientific model validation"}


def _bool(value):
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError("coverage flag must be boolean")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_artifacts(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
