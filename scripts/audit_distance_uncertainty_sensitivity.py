"""Independent raw-record audit for the distance-uncertainty sensitivity experiment.

This module deliberately recomputes likelihoods and summaries from raw records
without importing the experiment runner or its posterior-summary implementation.
"""
from __future__ import annotations

import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

C_KM_S = 299_792.458
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
ARMS = ("A_fixed_realized_sigma", "B_fractional_sigma_H0", "C_fixed_reference_sigma")
EXPECTED_EQUATIONS = {
    "A_fixed_realized_sigma": "logL=-0.5*((d_obs-c*z/H0)/sigma_realized)^2-log(sigma_realized)-0.5*log(2*pi)",
    "B_fractional_sigma_H0": "sigma(H0)=f*c*z/H0; logL=-0.5*((d_obs-c*z/H0)/sigma(H0))^2-log(sigma(H0))-0.5*log(2*pi)",
    "C_fixed_reference_sigma": "sigma_ref=f*c*z/70; logL=-0.5*((d_obs-c*z/H0)/sigma_ref)^2-log(sigma_ref)-0.5*log(2*pi)",
}
RAW_COLUMNS = {
    "true_h0", "fractional_error", "replicate", "arm", "host_redshift",
    "standard_normal", "true_distance_mpc", "observed_distance_mpc",
    "generated_sigma_mpc", "pull_median", "median", "mean", "posterior_sd", "lower_68",
    "upper_68", "lower_95", "upper_95",
}
SUMMARY_COLUMNS = {
    "true_h0", "fractional_error", "arm", "n", "coverage_68",
    "coverage_68_successes", "coverage_68_wilson_low", "coverage_68_wilson_high",
    "coverage_95", "coverage_95_successes", "coverage_95_wilson_low", "coverage_95_wilson_high",
    "mean_pull_median", "sd_pull_median",
    "mean_median_error", "rmse_median", "mean_posterior_mean_error",
    "rmse_mean", "mean_width_68", "mean_width_95", "truth_on_prior_boundary",
}


def _posterior_from_formula(observed, z, frac, sigma_fixed, arm, grid):
    mean_distance = C_KM_S * z / grid
    if arm == ARMS[1]:
        sigma = frac * mean_distance
    else:
        sigma = np.full_like(grid, sigma_fixed, dtype=float)
    log_likelihood = (
        -0.5 * ((observed - mean_distance) / sigma) ** 2
        - np.log(sigma)
        - 0.5 * math.log(2.0 * math.pi)
    )
    weights = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(np.trapezoid(weights, grid))
    if not math.isfinite(area) or area <= 0:
        raise ValueError("posterior normalization failed")
    density = weights / area
    cdf = np.r_[0.0, np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))]
    cdf /= cdf[-1]
    mean = float(np.trapezoid(grid * density, grid))
    variance = float(np.trapezoid((grid - mean) ** 2 * density, grid))
    if not math.isfinite(variance) or variance <= 0:
        raise ValueError("posterior variance invalid")
    quantile = lambda p: float(np.interp(p, cdf, grid))
    return {
        "median": quantile(.5), "mean": mean, "posterior_sd": math.sqrt(variance),
        "lower_68": quantile(.16), "upper_68": quantile(.84),
        "lower_95": quantile(.025), "upper_95": quantile(.975),
    }


def _wilson_interval(successes, n, z=1.959963984540054):
    p = successes / n
    denominator = 1.0 + z*z/n
    center = (p + z*z/(2*n)) / denominator
    half = z * math.sqrt(p*(1-p)/n + z*z/(4*n*n)) / denominator
    return max(0.0, center-half), min(1.0, center+half)


def _close(actual, expected, label, atol=2e-8):
    if not math.isfinite(float(actual)) or not math.isclose(
        float(actual), float(expected), rel_tol=2e-9, abs_tol=atol
    ):
        raise ValueError(f"{label} mismatch: observed={actual!r}, reconstructed={expected!r}")


def audit_artifacts(output_dir: str | Path) -> dict:
    """Audit raw records, model assignments, posterior summaries, and cell summaries."""
    folder = Path(output_dir)
    raw_path, summary_path = folder / "replicates.csv", folder / "summary.json"
    if not raw_path.is_file() or not summary_path.is_file():
        raise ValueError("replicates.csv and summary.json must both exist")
    raw_bytes, summary_bytes = raw_path.read_bytes(), summary_path.read_bytes()
    metadata = json.loads(summary_bytes)
    if metadata.get("status") != "synthetic sensitivity diagnostic; not an observational GW analysis":
        raise ValueError("unexpected or missing scientific-scope metadata")
    if metadata.get("prior") != [50.0, 90.0]:
        raise ValueError("unexpected prior metadata")
    if metadata.get("arms") != list(ARMS):
        raise ValueError("unexpected arm metadata")
    if metadata.get("boundary_truths_must_be_interpreted_separately") is not True:
        raise ValueError("boundary interpretation warning missing")
    if not isinstance(metadata.get("code_revision"), str) or not metadata["code_revision"]:
        raise ValueError("code revision provenance missing")
    for field in ("python_version", "numpy_version", "scipy_version"):
        if not isinstance(metadata.get(field), str) or not metadata[field]:
            raise ValueError(f"{field} provenance missing")
    if metadata.get("failed_run_count") != 0:
        raise ValueError("completed-run metadata must record zero failed replicates")
    equations = metadata.get("likelihood_equations")
    if not isinstance(equations, dict) or set(equations) != set(ARMS):
        raise ValueError("likelihood equations metadata missing or incomplete")
    if equations != EXPECTED_EQUATIONS:
        raise ValueError("likelihood equation metadata does not match the audited model definitions")
    truths = metadata.get("truths")
    fractions = metadata.get("fractional_errors")
    reps = metadata.get("replicates_per_cell")
    grid_points = metadata.get("grid_points")
    seed = metadata.get("seed")
    bootstrap_resamples = metadata.get("bootstrap_resamples")
    if truths != [50, 55, 60, 65, 70, 75, 80, 85, 90]:
        raise ValueError("truth grid metadata mismatch")
    if fractions != [0.05, 0.1, 0.15, 0.2]:
        raise ValueError("fractional-error grid metadata mismatch")
    if isinstance(reps, bool) or not isinstance(reps, int) or reps < 2:
        raise ValueError("invalid replicate count metadata")
    if isinstance(grid_points, bool) or not isinstance(grid_points, int) or grid_points < 101:
        raise ValueError("invalid grid resolution metadata")
    if isinstance(seed, bool) or not isinstance(seed, int) or seed < 0:
        raise ValueError("invalid seed metadata")
    if isinstance(bootstrap_resamples, bool) or not isinstance(bootstrap_resamples, int) or bootstrap_resamples < 100:
        raise ValueError("invalid bootstrap resample metadata")

    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        if not reader.fieldnames or not RAW_COLUMNS.issubset(reader.fieldnames):
            raise ValueError("raw CSV missing required columns")
        rows = list(reader)
    expected_rows = len(truths) * len(fractions) * reps * len(ARMS)
    if len(rows) != expected_rows or metadata.get("raw_record_count") != expected_rows:
        raise ValueError(f"raw record count mismatch: {len(rows)} vs {expected_rows}")
    grid = np.linspace(PRIOR_MIN, PRIOR_MAX, grid_points)
    seen, grouped = set(), {}
    records_by_key = {}
    for row in rows:
        truth, frac, rep, arm = float(row["true_h0"]), float(row["fractional_error"]), int(row["replicate"]), row["arm"]
        key = (truth, frac, rep, arm)
        if key in seen:
            raise ValueError(f"duplicate raw key: {key}")
        seen.add(key)
        if truth not in truths or frac not in fractions or not 0 <= rep < reps or arm not in ARMS:
            raise ValueError(f"unexpected design key: {key}")
        numeric = {name: float(row[name]) for name in RAW_COLUMNS - {"arm"}}
        if not all(math.isfinite(value) for value in numeric.values()):
            raise ValueError(f"non-finite raw value for {key}")
        z, eps = numeric["host_redshift"], numeric["standard_normal"]
        distance = C_KM_S * z / truth
        _close(numeric["true_distance_mpc"], distance, f"{key}: true distance")
        sigma_a = frac * distance
        sigma_c = frac * C_KM_S * z / 70.0
        if arm == ARMS[2]:
            expected_sigma, expected_obs = sigma_c, distance + sigma_c * eps
        else:
            expected_sigma, expected_obs = sigma_a, distance + sigma_a * eps
        _close(numeric["generated_sigma_mpc"], expected_sigma, f"{key}: generated sigma")
        _close(numeric["observed_distance_mpc"], expected_obs, f"{key}: observed distance")
        summary = _posterior_from_formula(
            numeric["observed_distance_mpc"], z, frac, expected_sigma, arm, grid
        )
        for name, value in summary.items():
            _close(numeric[name], value, f"{key}: posterior {name}")
        _close(numeric["pull_median"], (summary["median"] - truth) / summary["posterior_sd"],
               f"{key}: posterior-median pull")
        records_by_key[key] = row
        grouped.setdefault((truth, frac, arm), []).append(numeric)

    expected_keys = {
        (float(t), float(f), r, arm)
        for t in truths for f in fractions for r in range(reps) for arm in ARMS
    }
    if seen != expected_keys:
        raise ValueError("raw replicate key grid is incomplete or contains unexpected keys")
    for truth in truths:
        for frac in fractions:
            for rep in range(reps):
                a = records_by_key[(float(truth), float(frac), rep, ARMS[0])]
                b = records_by_key[(float(truth), float(frac), rep, ARMS[1])]
                c = records_by_key[(float(truth), float(frac), rep, ARMS[2])]
                # All arms are intended to use common redshift and standard-normal
                # draws; only Arm C's generated absolute noise scale differs.
                for col in ("host_redshift", "standard_normal", "true_distance_mpc"):
                    _close(a[col], b[col], f"paired A/B {truth}/{frac}/{rep}: {col}")
                    _close(a[col], c[col], f"paired A/C {truth}/{frac}/{rep}: {col}")
                for col in ("observed_distance_mpc", "generated_sigma_mpc"):
                    _close(a[col], b[col], f"paired A/B {truth}/{frac}/{rep}: {col}")

    summary_rows = metadata.get("summary")
    if not isinstance(summary_rows, list) or len(summary_rows) != len(truths) * len(fractions) * len(ARMS):
        raise ValueError("summary cell count mismatch")
    if metadata.get("summary_cell_count") != len(summary_rows):
        raise ValueError("summary cell-count metadata mismatch")
    summary_by_key = {}
    for item in summary_rows:
        if not SUMMARY_COLUMNS.issubset(item):
            raise ValueError("summary row missing required fields")
        key = (float(item["true_h0"]), float(item["fractional_error"]), item["arm"])
        if key in summary_by_key:
            raise ValueError(f"duplicate summary key: {key}")
        summary_by_key[key] = item
    if set(summary_by_key) != set(grouped):
        raise ValueError("summary cell grid does not match raw records")

    for key, cell in grouped.items():
        truth, frac, arm = key
        n = len(cell)
        median_errors = np.array([r["median"] - truth for r in cell])
        mean_errors = np.array([r["mean"] - truth for r in cell])
        c68 = sum(r["lower_68"] <= truth <= r["upper_68"] for r in cell)
        c95 = sum(r["lower_95"] <= truth <= r["upper_95"] for r in cell)
        c68_low, c68_high = _wilson_interval(c68, n)
        c95_low, c95_high = _wilson_interval(c95, n)
        pulls = np.array([r["pull_median"] for r in cell])
        expected = {
            "n": n, "coverage_68": c68/n, "coverage_68_successes": c68,
            "coverage_68_wilson_low": c68_low, "coverage_68_wilson_high": c68_high,
            "coverage_95": c95/n, "coverage_95_successes": c95,
            "coverage_95_wilson_low": c95_low, "coverage_95_wilson_high": c95_high,
            "mean_pull_median": float(pulls.mean()),
            "sd_pull_median": float(pulls.std(ddof=1)) if n > 1 else 0.0,
            "mean_median_error": float(median_errors.mean()),
            "rmse_median": float(np.sqrt(np.mean(median_errors**2))),
            "mean_posterior_mean_error": float(mean_errors.mean()),
            "rmse_mean": float(np.sqrt(np.mean(mean_errors**2))),
            "mean_width_68": float(np.mean([r["upper_68"] - r["lower_68"] for r in cell])),
            "mean_width_95": float(np.mean([r["upper_95"] - r["lower_95"] for r in cell])),
            "truth_on_prior_boundary": truth in (PRIOR_MIN, PRIOR_MAX),
        }
        observed = summary_by_key[key]
        for name, value in expected.items():
            if isinstance(value, bool) or isinstance(value, int):
                if observed.get(name) != value:
                    raise ValueError(f"summary {key}: {name} mismatch")
            else:
                _close(observed.get(name), value, f"summary {key}: {name}")

    contrasts = metadata.get("paired_ab_contrasts")
    metric_names = ("coverage_68_indicator", "absolute_median_error", "squared_median_error",
                    "width_68", "width_95")
    expected_contrast_count = len(truths) * len(fractions) * len(metric_names)
    if not isinstance(contrasts, list) or len(contrasts) != expected_contrast_count:
        raise ValueError("paired A/B contrast count mismatch")
    if metadata.get("paired_ab_contrast_count") != expected_contrast_count:
        raise ValueError("paired A/B contrast-count metadata mismatch")
    contrasts_by_key = {}
    for item in contrasts:
        key = (float(item["true_h0"]), float(item["fractional_error"]), item["contrast"])
        if key in contrasts_by_key:
            raise ValueError(f"duplicate paired contrast: {key}")
        contrasts_by_key[key] = item
    for cell_index, (truth, frac) in enumerate((float(t), float(f)) for t in truths for f in fractions):
        pairs = []
        for rep in range(reps):
            a = records_by_key[(truth, frac, rep, ARMS[0])]
            b = records_by_key[(truth, frac, rep, ARMS[1])]
            a_med, b_med = float(a["median"]), float(b["median"])
            pairs.append({
                "coverage_68_indicator": float(float(a["lower_68"]) <= truth <= float(a["upper_68"])) -
                                          float(float(b["lower_68"]) <= truth <= float(b["upper_68"])),
                "absolute_median_error": abs(a_med-truth) - abs(b_med-truth),
                "squared_median_error": (a_med-truth)**2 - (b_med-truth)**2,
                "width_68": (float(a["upper_68"])-float(a["lower_68"])) -
                            (float(b["upper_68"])-float(b["lower_68"])),
                "width_95": (float(a["upper_95"])-float(a["lower_95"])) -
                            (float(b["upper_95"])-float(b["lower_95"])),
            })
        for metric_index, metric in enumerate(metric_names):
            diffs = np.asarray([p[metric] for p in pairs], dtype=float)
            rng_seed = (seed + 1_000_003*cell_index + 97_409*metric_index) % (2**32)
            rng = np.random.default_rng(rng_seed)
            indices = rng.integers(0, len(diffs), size=(bootstrap_resamples, len(diffs)))
            boot_means = diffs[indices].mean(axis=1)
            key = (truth, frac, metric)
            if key not in contrasts_by_key:
                raise ValueError(f"missing paired contrast {key}")
            item = contrasts_by_key[key]
            if item.get("direction") != "A_minus_B" or item.get("n_paired_replicates") != reps:
                raise ValueError(f"paired contrast metadata mismatch for {key}")
            if item.get("bootstrap_resamples") != bootstrap_resamples or item.get("bootstrap_seed") != rng_seed:
                raise ValueError(f"paired bootstrap metadata mismatch for {key}")
            _close(item.get("mean_difference"), float(diffs.mean()), f"{key}: mean difference")
            _close(item.get("bootstrap_95_low"), float(np.quantile(boot_means, .025)), f"{key}: bootstrap low")
            _close(item.get("bootstrap_95_high"), float(np.quantile(boot_means, .975)), f"{key}: bootstrap high")
    return {
        "status": "PASS",
        "scope": "synthetic known-host low-redshift sensitivity diagnostic only",
        "raw_record_count": len(rows),
        "summary_cell_count": len(summary_rows),
        "replicates_per_cell": reps,
        "grid_points": grid_points,
        "code_revision": metadata["code_revision"],
        "python_version": metadata["python_version"],
        "numpy_version": metadata["numpy_version"],
        "scipy_version": metadata["scipy_version"],
        "unique_raw_keys": len(seen),
        "a_b_pairs_verified": len(truths) * len(fractions) * reps,
        "posterior_rows_reconstructed": len(rows),
        "summary_cells_reconstructed": len(summary_rows),
        "paired_ab_contrast_count": len(contrasts),
        "paired_ab_bootstrap_intervals_reconstructed": len(contrasts),
        "boundary_truths_separated": True,
        "replicates_sha256": hashlib.sha256(raw_bytes).hexdigest(),
        "summary_sha256": hashlib.sha256(summary_bytes).hexdigest(),
        "audit_method": "independent formula reconstruction; does not import the experiment runner or likelihood module",
    }


def main() -> None:
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", nargs="?", default="artifacts/distance_uncertainty_sensitivity")
    args = parser.parse_args()
    report = audit_artifacts(args.output_dir)
    report_path = Path(args.output_dir) / "audit_report.json"
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
