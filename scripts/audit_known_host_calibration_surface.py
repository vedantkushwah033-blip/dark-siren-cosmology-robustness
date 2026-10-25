"""Independent raw-record audit for the known-host calibration surface.

This deliberately recomputes cell metrics from CSV records rather than calling
the experiment's summarize_cell implementation.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
import hashlib
import os
import importlib.metadata
import platform
from collections import defaultdict
from pathlib import Path
from statistics import median


def wilson(successes: int, n: int, z: float = 1.959963984540054) -> tuple[float, float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("invalid binomial counts")
    p = successes / n
    z2 = z * z
    denominator = 1 + z2 / n
    center = (p + z2 / (2 * n)) / denominator
    half = z * math.sqrt(p * (1 - p) / n + z2 / (4 * n * n)) / denominator
    return max(0.0, center - half), min(1.0, center + half)


def close(actual: float, expected: float, label: str, tol: float = 1e-8) -> None:
    if not math.isfinite(actual) or not math.isclose(actual, expected, rel_tol=tol, abs_tol=tol):
        raise ValueError(f"{label}: stored={actual!r}, recomputed={expected!r}")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_base_seed(value: object) -> int:
    """Require a JSON integer seed; reject booleans and coercible non-integers."""
    if isinstance(value, bool) or not isinstance(value, int) or value < 0:
        raise ValueError("base seed must be a non-negative integer")
    return value


def parse_nonnegative_integer(value: object, label: str) -> int:
    """Parse an integer-valued CSV field without silently truncating/coercing it."""
    if not isinstance(value, str) or not value.isascii() or not value.isdigit():
        raise ValueError(f"{label} must be a non-negative integer")
    return int(value)


def validate_metadata_seed_and_design(metadata: dict, replicates: int) -> int:
    """Validate protocol metadata types and fixed design before auditing records."""
    if not isinstance(metadata, dict):
        raise ValueError("metadata must be a JSON object")
    seed = validate_base_seed(metadata.get("seed"))
    if isinstance(metadata.get("raw_row_count"), bool) or not isinstance(metadata.get("raw_row_count"), int):
        raise ValueError("metadata raw_row_count must be an integer")
    if isinstance(metadata.get("cell_count"), bool) or not isinstance(metadata.get("cell_count"), int):
        raise ValueError("metadata cell_count must be an integer")
    if isinstance(metadata.get("replicates_per_cell"), bool) or not isinstance(metadata.get("replicates_per_cell"), int):
        raise ValueError("metadata replicates_per_cell must be an integer")
    if metadata["raw_row_count"] != 9 * 4 * replicates or metadata["cell_count"] != 36:
        raise ValueError("metadata row/cell counts do not match the design")
    if metadata["replicates_per_cell"] != replicates:
        raise ValueError("metadata replicate count does not match the requested design")
    if metadata.get("prior") != [50.0, 90.0] or metadata.get("posterior_grid_points") != 1601:
        raise ValueError("unexpected prior or posterior grid in metadata")
    if metadata.get("original_repository") != "vedantkushwah033-blip/bayesian-dark-siren-cosmology":
        raise ValueError("unexpected source simulator repository")
    if metadata.get("original_commit") != "b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f":
        raise ValueError("source simulator commit does not match the pinned protocol")
    if metadata.get("true_h0_values") != [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]:
        raise ValueError("metadata truth grid does not match the prespecified protocol")
    if metadata.get("fractional_distance_errors") != [0.05, 0.10, 0.15, 0.20]:
        raise ValueError("metadata noise grid does not match the prespecified protocol")
    cells = metadata.get("cells")
    if not isinstance(cells, list) or len(cells) != 36:
        raise ValueError("metadata must contain exactly 36 cell summaries")
    if any(not isinstance(cell, dict) for cell in cells):
        raise ValueError("each JSON metadata cell must be an object")
    return seed


def validate_metadata_cell_counts(metadata_cell: dict, expected: dict, key: tuple[float, float]) -> None:
    """Reject malformed JSON count/flag types before comparing their values."""
    for field in ("n", "coverage_68_successes", "coverage_95_successes"):
        value = metadata_cell[field]
        # JSON booleans are Python ints by subclassing, so reject them explicitly.
        if isinstance(value, bool) or not isinstance(value, int):
            raise ValueError(f"JSON metadata {key} count {field} must be an integer")
        if value != int(expected[field]):
            raise ValueError(f"JSON metadata {key} count {field} mismatch")
    boundary_flag = metadata_cell["truth_on_prior_boundary"]
    if not isinstance(boundary_flag, bool):
        raise ValueError(f"JSON metadata boundary-truth flag must be boolean for cell {key}")
    if boundary_flag != expected["truth_on_prior_boundary"]:
        raise ValueError(f"JSON metadata boundary-truth flag mismatch for cell {key}")


def recompute_posterior_from_record(row: dict[str, str]) -> dict[str, float]:
    """Independently reconstruct one-host posterior summaries from raw data."""
    import numpy as np

    grid = np.linspace(50.0, 90.0, 1601)
    z = float(row["host_redshift"])
    observed = float(row["observed_distance_mpc"])
    sigma = float(row["distance_sigma_mpc"])
    if (
        not all(math.isfinite(v) for v in (z, observed, sigma))
        or z < 0.005 or z > 0.08
        or observed <= 0
        or sigma <= 0
    ):
        raise ValueError("invalid raw observation for posterior reconstruction")
    model_distance = 299792.458 * z / grid
    log_likelihood = -0.5 * ((observed - model_distance) / sigma) ** 2
    log_likelihood -= float(np.max(log_likelihood))
    posterior = np.exp(log_likelihood)
    area = float(np.trapezoid(posterior, grid))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("independent posterior normalization failed")
    posterior = posterior / area
    mean = float(np.trapezoid(grid * posterior, grid))
    variance = float(np.trapezoid((grid - mean) ** 2 * posterior, grid))
    sd = math.sqrt(max(variance, 0.0))
    cdf = np.concatenate([[0.0], np.cumsum((posterior[1:] + posterior[:-1]) * np.diff(grid) / 2.0)])
    cdf = cdf / cdf[-1]
    return {
        "mean": mean,
        "median": float(np.interp(0.5, cdf, grid)),
        "posterior_sd": sd,
        "lower_68": float(np.interp(0.16, cdf, grid)),
        "upper_68": float(np.interp(0.84, cdf, grid)),
        "lower_95": float(np.interp(0.025, cdf, grid)),
        "upper_95": float(np.interp(0.975, cdf, grid)),
    }


def recompute(rows: list[dict[str, str]]) -> dict:
    n = len(rows)
    if n == 0:
        raise ValueError("empty cell")
    truth = float(rows[0]["true_h0"])
    frac = float(rows[0]["fractional_distance_error"])
    medians = [float(r["median"]) for r in rows]
    means = [float(r["mean"]) for r in rows]
    sd = [float(r["posterior_sd"]) for r in rows]
    l68 = [float(r["lower_68"]) for r in rows]
    u68 = [float(r["upper_68"]) for r in rows]
    l95 = [float(r["lower_95"]) for r in rows]
    u95 = [float(r["upper_95"]) for r in rows]
    numeric_fields = (
        "true_h0", "fractional_distance_error", "host_redshift",
        "true_distance_mpc", "observed_distance_mpc", "distance_sigma_mpc",
        "median", "mean", "posterior_sd", "lower_68", "upper_68",
        "lower_95", "upper_95", "width_68", "width_95", "pull_median",
    )
    for i, row in enumerate(rows):
        for field in numeric_fields:
            x = float(row[field])
            if not math.isfinite(x):
                raise ValueError(f"non-finite {field} in cell record {i}")
        if float(row["true_h0"]) != truth:
            raise ValueError(f"truth changes within a summary cell at record {i}")
        if float(row["fractional_distance_error"]) != frac:
            raise ValueError(f"fractional distance error changes within a summary cell at record {i}")
        if truth not in (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0):
            raise ValueError(f"unexpected true H0 in cell: {truth}")
        if frac not in (0.05, 0.10, 0.15, 0.20):
            raise ValueError(f"unexpected fractional distance error in cell: {frac}")
        if not 50 <= float(row["median"]) <= 90 or not 50 <= float(row["mean"]) <= 90:
            raise ValueError(f"posterior point summary outside prior in cell record {i}")
        if not (50 <= float(row["lower_95"]) <= float(row["lower_68"]) <= float(row["median"]) <= float(row["upper_68"]) <= float(row["upper_95"]) <= 90):
            raise ValueError(f"credible intervals invalid, fail to contain the median, or lie outside prior in cell record {i}")
        if float(row["posterior_sd"]) <= 0:
            raise ValueError(f"non-positive posterior SD in cell record {i}")
        if float(row["width_68"]) <= 0 or float(row["width_95"]) <= 0:
            raise ValueError(f"non-positive credible interval width in cell record {i}")
        if float(row["host_redshift"]) < 0.005 or float(row["host_redshift"]) > 0.08:
            raise ValueError(f"host redshift outside the prespecified catalogue range in cell record {i}")
        if float(row["observed_distance_mpc"]) <= 0:
            raise ValueError(f"non-positive observed distance in cell record {i}")
        if float(row["true_distance_mpc"]) <= 0 or float(row["distance_sigma_mpc"]) <= 0:
            raise ValueError(f"non-positive true distance or distance uncertainty in cell record {i}")
        close(float(row["width_68"]), float(row["upper_68"]) - float(row["lower_68"]), f"width_68 row {i}")
        close(float(row["width_95"]), float(row["upper_95"]) - float(row["lower_95"]), f"width_95 row {i}")
        close(float(row["pull_median"]), (float(row["median"]) - truth) / float(row["posterior_sd"]), f"pull row {i}")
        expected_sigma = frac * float(row["true_distance_mpc"])
        close(float(row["distance_sigma_mpc"]), expected_sigma, f"distance sigma row {i}")
        close(float(row["true_distance_mpc"]), 299792.458 * float(row["host_redshift"]) / truth, f"low-z truth distance row {i}", tol=1e-7)
    err = [x - truth for x in medians]
    mean_err = [x - truth for x in means]
    c68 = sum(a <= truth <= b for a, b in zip(l68, u68))
    c95 = sum(a <= truth <= b for a, b in zip(l95, u95))
    ci68 = wilson(c68, n)
    ci95 = wilson(c95, n)
    pulls = [e / s for e, s in zip(err, sd)]
    return {
        "true_h0": truth,
        "fractional_distance_error": frac,
        "n": n,
        "coverage_68": c68 / n,
        "coverage_68_successes": c68,
        "coverage_68_wilson_low": ci68[0],
        "coverage_68_wilson_high": ci68[1],
        "coverage_68_mc_se": math.sqrt((c68 / n) * (1 - c68 / n) / n),
        "coverage_95": c95 / n,
        "coverage_95_successes": c95,
        "coverage_95_wilson_low": ci95[0],
        "coverage_95_wilson_high": ci95[1],
        "coverage_95_mc_se": math.sqrt((c95 / n) * (1 - c95 / n) / n),
        "median_bias": median(err),
        "mean_median_error": math.fsum(err) / n,
        "mean_posterior_mean_error": math.fsum(mean_err) / n,
        "rmse_median": math.sqrt(math.fsum(e * e for e in err) / n),
        "rmse_mean": math.sqrt(math.fsum(e * e for e in mean_err) / n),
        "mean_width_68": math.fsum(b - a for a, b in zip(l68, u68)) / n,
        "mean_width_95": math.fsum(b - a for a, b in zip(l95, u95)) / n,
        "mean_posterior_sd": math.fsum(sd) / n,
        "mean_median_pull": math.fsum(pulls) / n,
        "sd_median_pull": math.sqrt(math.fsum((x - math.fsum(pulls) / n) ** 2 for x in pulls) / (n - 1)) if n > 1 else 0.0,
        "median_near_prior_edge_fraction": sum((x - 50 <= 0.5) or (90 - x <= 0.5) for x in medians) / n,
        "truth_on_prior_boundary": truth in (50.0, 90.0),
    }


def validate_paired_group(key: tuple[float, int], group: list[dict[str, str]]) -> None:
    """Validate common host and standardized-noise draws across uncertainty levels."""
    if len(group) != 4:
        raise ValueError(f"paired truth/replicate {key} has {len(group)} noise levels")
    hosts = [float(r["host_redshift"]) for r in group]
    distances = [float(r["true_distance_mpc"]) for r in group]
    residuals = [
        (float(r["observed_distance_mpc"]) - float(r["true_distance_mpc"]))
        / float(r["distance_sigma_mpc"])
        for r in group
    ]
    if max(hosts) - min(hosts) > 1e-12 or max(distances) - min(distances) > 1e-9:
        raise ValueError(f"host/population not paired across noise levels for {key}")
    if max(residuals) - min(residuals) > 1e-9:
        raise ValueError(f"standardized noise is not paired across noise levels for {key}")


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="results/known-host-calibration-surface")
    parser.add_argument("--replicates", type=int, required=True)
    args = parser.parse_args()
    if not 100 <= args.replicates <= 2000:
        parser.error("--replicates must be between 100 and 2000")
    out = Path(args.output_dir)
    raw = read_csv(out / "per_run_results.csv")
    stored = read_csv(out / "cell_summary.csv")
    metadata = json.loads((out / "cell_summary.json").read_text(encoding="utf-8"))
    expected_count = 9 * 4 * args.replicates
    if len(raw) != expected_count:
        raise ValueError(f"expected {expected_count} raw rows, found {len(raw)}")
    if len(stored) != 36:
        raise ValueError("expected 36 cell summaries in CSV")
    base_seed = validate_metadata_seed_and_design(metadata, args.replicates)
    for i, row in enumerate(raw):
        independently_reconstructed = recompute_posterior_from_record(row)
        for field, value in independently_reconstructed.items():
            close(float(row[field]), value, f"independent posterior {field} raw row {i}", tol=2e-6)
    expected_truths = [50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0]
    expected_noise = [0.05, 0.10, 0.15, 0.20]

    grouped: dict[tuple[float, float], list[dict[str, str]]] = defaultdict(list)
    seen = set()
    paired: dict[tuple[float, int], list[dict[str, str]]] = defaultdict(list)
    for row in raw:
        run_index = parse_nonnegative_integer(row["run"], "raw replicate index")
        row_seed = parse_nonnegative_integer(row["seed"], "raw replicate seed")
        key = (float(row["true_h0"]), float(row["fractional_distance_error"]), run_index)
        if key in seen:
            raise ValueError(f"duplicate raw replicate key: {key}")
        if not 0 <= key[2] < args.replicates:
            raise ValueError(f"replicate index outside requested range: {key}")
        if row_seed != base_seed + key[2]:
            raise ValueError(f"seed does not match replicate index for {key}")
        seen.add(key)
        grouped[key[:2]].append(row)
        paired[(key[0], key[2])].append(row)
    expected_cells = {(truth, noise) for truth in expected_truths for noise in expected_noise}
    if set(grouped) != expected_cells:
        raise ValueError(f"raw cell grid mismatch; missing={sorted(expected_cells - set(grouped))}, extra={sorted(set(grouped) - expected_cells)}")
    for key, cell in grouped.items():
        if len(cell) != args.replicates:
            raise ValueError(f"cell {key} has {len(cell)} records, expected {args.replicates}")
    for key, group in paired.items():
        validate_paired_group(key, group)

    stored_by_key = {(float(r["true_h0"]), float(r["fractional_distance_error"])): r for r in stored}
    if len(stored_by_key) != 36:
        raise ValueError("duplicate summary-cell keys")
    fields = (
        "coverage_68", "coverage_68_wilson_low", "coverage_68_wilson_high",
        "coverage_68_mc_se", "coverage_95", "coverage_95_wilson_low",
        "coverage_95_wilson_high", "coverage_95_mc_se", "median_bias",
        "mean_median_error", "mean_posterior_mean_error", "rmse_median",
        "rmse_mean", "mean_width_68", "mean_width_95", "mean_posterior_sd",
        "mean_median_pull", "sd_median_pull", "median_near_prior_edge_fraction",
    )
    metadata_by_key = {
        (float(cell["true_h0"]), float(cell["fractional_distance_error"])): cell
        for cell in metadata["cells"]
    }
    if len(metadata_by_key) != 36:
        raise ValueError("duplicate summary-cell keys in JSON metadata")
    for key, cell_rows in grouped.items():
        expected = recompute(cell_rows)
        actual = stored_by_key.get(key)
        metadata_cell = metadata_by_key.get(key)
        if metadata_cell is None:
            raise ValueError(f"JSON metadata missing cell {key}")
        if actual is None:
            raise ValueError(f"summary missing cell {key}")
        if int(actual["n"]) != expected["n"] or int(actual["coverage_68_successes"]) != expected["coverage_68_successes"] or int(actual["coverage_95_successes"]) != expected["coverage_95_successes"]:
            raise ValueError(f"summary counts disagree for cell {key}")
        for field in fields:
            close(float(actual[field]), float(expected[field]), f"summary {key} {field}", tol=1e-7)
        if actual["truth_on_prior_boundary"].lower() != str(expected["truth_on_prior_boundary"]).lower():
            raise ValueError(f"boundary-truth flag mismatch for cell {key}")
        for field in fields:
            close(float(metadata_cell[field]), float(expected[field]), f"JSON metadata {key} {field}", tol=1e-7)
        validate_metadata_cell_counts(metadata_cell, expected, key)

    audit_report = {
        "audit": "PASS",
        "raw_records": len(raw),
        "summary_cells": len(stored),
        "replicates_per_cell": args.replicates,
        "paired_truth_replicate_groups": len(paired),
        "summary_metrics_independently_recomputed": len(fields),
        "posterior_summaries_independently_reconstructed_from_raw_observations": True,
        "pairing_checks": "host redshift, true distance, and standardized noise across four noise levels",
        "metadata_json_crosschecked": True,
        "source_simulator_commit": metadata.get("original_commit"),
        "workflow_run_id": os.environ.get("GITHUB_RUN_ID"),
        "workflow_commit_sha": os.environ.get("GITHUB_SHA"),
        "workflow_ref": os.environ.get("GITHUB_REF"),
        "auditor_script_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
        "base_seed": base_seed,
        "prior": metadata["prior"],
        "python_version": platform.python_version(),
        "input_sha256": {
            filename: hashlib.sha256((out / filename).read_bytes()).hexdigest()
            for filename in ("per_run_results.csv", "cell_summary.csv", "cell_summary.json")
        },
    }
    for package in ("numpy", "scipy", "pandas", "matplotlib", "pytest"):
        try:
            audit_report[f"{package}_version"] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            audit_report[f"{package}_version"] = "not-installed"
    (out / "audit_report.json").write_text(json.dumps(audit_report, indent=2), encoding="utf-8")
    print(json.dumps(audit_report, indent=2))


if __name__ == "__main__":
    main()
