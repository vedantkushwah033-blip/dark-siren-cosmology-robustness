"""Consistency audit for host-density-to-H0 benchmark summary artifacts.

Important: the current schema stores posterior summaries and a claimed integral,
not the posterior grid. The auditor can check that claimed scalar but cannot
independently reconstruct posterior normalization or quantiles from the grid.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

ARMS = ("oracle", "catalogue_only", "known_completeness")
REQUIRED_COLUMNS = {
    "replicate", "seed", "arm", "h0_true", "n_parent", "n_catalogue_observed",
    "n_detected", "min_detection_efficiency", "posterior_mean", "posterior_median",
    "ci68_low", "ci68_high", "ci68_width", "posterior_integral", "covered_68",
}


def _number(value, label):
    try:
        result = float(value)
    except (TypeError, ValueError) as exc:
        raise ValueError(f"{label} must be numeric") from exc
    if not math.isfinite(result):
        raise ValueError(f"{label} must be finite")
    return result


def _close(actual, expected, label, atol=1e-8):
    if not math.isclose(_number(actual, label), _number(expected, label), rel_tol=1e-8, abs_tol=atol):
        raise ValueError(f"{label} mismatch: {actual!r} != {expected!r}")


def _bool(value, label):
    normalized = str(value).strip().lower()
    if normalized in {"true", "1"}:
        return True
    if normalized in {"false", "0"}:
        return False
    raise ValueError(f"{label} must be a boolean")


def audit_artifacts(output_dir: str | Path) -> dict:
    """Recompute aggregate metrics and validate paired replicate records."""
    folder = Path(output_dir)
    raw_path, summary_path = folder / "replicate_metrics.csv", folder / "summary.json"
    if not raw_path.is_file() or not summary_path.is_file():
        raise ValueError("replicate_metrics.csv and summary.json must both exist")
    raw_bytes, summary_bytes = raw_path.read_bytes(), summary_path.read_bytes()
    metadata = json.loads(summary_bytes)
    if metadata.get("status") != "synthetic_exploratory_benchmark_not_real_H0_measurement":
        raise ValueError("missing or unexpected scientific-scope status")
    n_replicates = metadata.get("n_replicates")
    seed_start = metadata.get("seed_start")
    h0_true = _number(metadata.get("h0_true_km_s_mpc"), "h0_true_km_s_mpc")
    h0_grid = metadata.get("h0_grid")
    if isinstance(n_replicates, bool) or not isinstance(n_replicates, int) or n_replicates < 1:
        raise ValueError("invalid n_replicates metadata")
    if isinstance(seed_start, bool) or not isinstance(seed_start, int) or seed_start < 0:
        raise ValueError("invalid seed_start metadata")
    if not isinstance(h0_grid, list) or len(h0_grid) != 3:
        raise ValueError("invalid H0 grid metadata")
    h0_min, h0_max, h0_points = h0_grid
    if not (0 < float(h0_min) < float(h0_max)) or isinstance(h0_points, bool) or int(h0_points) < 31:
        raise ValueError("invalid H0 grid bounds or resolution")
    if not float(h0_min) <= h0_true <= float(h0_max):
        raise ValueError("true H0 lies outside declared H0 grid")
    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not REQUIRED_COLUMNS.issubset(reader.fieldnames):
            raise ValueError("raw replicate CSV is missing required columns")
        rows = list(reader)
    if len(rows) != n_replicates * len(ARMS):
        raise ValueError("raw record count does not equal n_replicates × three arms")
    by_key = {}
    for row in rows:
        rep_raw = row["replicate"]
        seed_raw = row["seed"]
        try:
            rep, seed = int(rep_raw), int(seed_raw)
        except (TypeError, ValueError) as exc:
            raise ValueError("replicate and seed must be integers") from exc
        arm = row["arm"]
        key = (rep, arm)
        if key in by_key:
            raise ValueError(f"duplicate replicate-arm row: {key}")
        if not 0 <= rep < n_replicates or arm not in ARMS:
            raise ValueError(f"unexpected replicate/arm key: {key}")
        if seed != seed_start + rep:
            raise ValueError(f"seed mismatch for replicate {rep}: {seed}")
        _close(row["h0_true"], h0_true, "row H0 truth")
        for col in ("n_parent", "n_catalogue_observed", "n_detected"):
            value = int(row[col])
            if value <= 0:
                raise ValueError(f"{col} must be positive")
        efficiency = _number(row["min_detection_efficiency"], "min_detection_efficiency")
        if not 0 < efficiency <= 1:
            raise ValueError("minimum detection efficiency must lie in (0, 1]")
        mean = _number(row["posterior_mean"], "posterior_mean")
        median = _number(row["posterior_median"], "posterior_median")
        low = _number(row["ci68_low"], "ci68_low")
        high = _number(row["ci68_high"], "ci68_high")
        width = _number(row["ci68_width"], "ci68_width")
        integral = _number(row["posterior_integral"], "posterior_integral")
        if not float(h0_min) <= low <= median <= high <= float(h0_max):
            raise ValueError(f"posterior interval/median outside declared grid or misordered: {key}")
        if not float(h0_min) <= mean <= float(h0_max):
            raise ValueError(f"posterior mean outside declared grid: {key}")
        if width <= 0:
            raise ValueError(f"non-positive interval width: {key}")
        _close(width, high - low, "68% interval width")
        _close(integral, 1.0, "posterior integral", atol=1e-6)
        covered = _bool(row["covered_68"], "covered_68")
        expected_covered = low <= h0_true <= high
        if covered != expected_covered:
            raise ValueError(f"coverage indicator inconsistent with interval: {key}")
        by_key[key] = row
    for rep in range(n_replicates):
        present = {arm for (r, arm) in by_key if r == rep}
        if present != set(ARMS):
            raise ValueError(f"replicate {rep} is not paired across all arms")
        records = [by_key[(rep, arm)] for arm in ARMS]
        for col in ("n_parent", "n_detected", "h0_true", "seed"):
            if len({str(record[col]) for record in records}) != 1:
                raise ValueError(f"paired arms disagree on {col} for replicate {rep}")
    aggregates = metadata.get("aggregate_metrics")
    if not isinstance(aggregates, list) or len(aggregates) != len(ARMS):
        raise ValueError("aggregate_metrics must contain exactly three arm summaries")
    aggregate_by_arm = {item.get("arm"): item for item in aggregates if isinstance(item, dict)}
    if set(aggregate_by_arm) != set(ARMS):
        raise ValueError("aggregate arm labels do not match the raw records")
    for arm in ARMS:
        group = [by_key[(rep, arm)] for rep in range(n_replicates)]
        errors = [_number(row["posterior_mean"], "posterior_mean") - h0_true for row in group]
        expected = {
            "n_replicates": n_replicates,
            "mean_bias": sum(errors) / n_replicates,
            "rmse": math.sqrt(sum(error * error for error in errors) / n_replicates),
            "empirical_68_coverage": sum(_bool(row["covered_68"], "covered_68") for row in group) / n_replicates,
            "mean_68_interval_width": sum(float(row["ci68_width"]) for row in group) / n_replicates,
            "mean_posterior_mean": sum(float(row["posterior_mean"]) for row in group) / n_replicates,
        }
        observed = aggregate_by_arm[arm]
        for key, value in expected.items():
            if key == "n_replicates":
                if observed.get(key) != value:
                    raise ValueError(f"{arm} {key} mismatch")
            else:
                _close(observed.get(key), value, f"{arm} aggregate {key}")
    return {
        "status": "PASS",
        "raw_record_count": len(rows),
        "replicates_verified": n_replicates,
        "paired_arm_groups_verified": n_replicates,
        "recorded_posterior_integral_claims_checked": len(rows),
        "posterior_grid_independently_reconstructed": False,
        "aggregate_metrics_reconstructed": len(ARMS),
        "replicate_metrics_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "summary_json_sha256": hashlib.sha256(summary_bytes).hexdigest(),
        "scope": "synthetic toy benchmark only; not a real H0 measurement",
        "audit_boundary": "posterior_integral is a recorded scalar; posterior grid, mean, and quantiles are not independently reconstructed from grid samples",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    args = parser.parse_args()
    print(json.dumps(audit_artifacts(args.output_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
