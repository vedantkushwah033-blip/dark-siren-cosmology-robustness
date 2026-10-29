"""Reproducible exact-density oracle diagnostics for a simplified H0 toy model.

This is a controlled synthetic check, not a realistic gravitational-wave inference
pipeline. CLI runs save per-replicate data, replayable mode seeds, and a SHA-256
manifest so another implementation can audit both rows and aggregate summaries.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import os
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import chi2

from dark_siren_robustness.generative_selection import latent_galaxy_intensity, C_KM_S
from dark_siren_robustness.host_density_benchmark import posterior_summary
from dark_siren_robustness.selection_toy import selection_aware_toy_likelihood

SIGMA, THRESHOLD, N_DETECTED = 140.0, 250.0, 40
Z_MIN, Z_MAX, SAMPLER_GRID_POINTS = 0.01, 0.40, 40001
MODE_OFFSETS = {"prior": 0, "fixed_60": 100_000, "fixed_70": 200_000, "fixed_80": 300_000}
MODEL_SCOPE = "synthetic_low_redshift_exact_density_oracle_not_real_GW_inference"


def mode_label(h0_mode):
    if isinstance(h0_mode, str) and h0_mode == "prior":
        return "prior"
    value = float(h0_mode)
    if not np.isfinite(value) or value <= 0:
        raise ValueError("fixed H0 truth must be finite and positive")
    return f"fixed_{value:g}"


def exact_density(z):
    """Normalized declared target density, including the toy merger weight."""
    z = np.asarray(z, dtype=float)
    density = latent_galaxy_intensity(z) * (1.0 + 1.5 * (z - z[0]) / (z[-1] - z[0]))
    return density / trapezoid(density, z)


def make_sampler(zmax=Z_MAX, zmin=Z_MIN, n=SAMPLER_GRID_POINTS):
    grid = np.linspace(zmin, zmax, n)
    density = exact_density(grid)
    cdf = np.r_[0.0, np.cumsum(0.5 * (density[1:] + density[:-1]) * np.diff(grid))]
    cdf /= cdf[-1]
    return lambda rng, size: np.interp(rng.random(size), cdf, grid)


def simulate(rng, h0, sample_z):
    """Draw exactly N_DETECTED distances passing the fixed threshold."""
    distances = []
    attempts = 0
    while len(distances) < N_DETECTED:
        attempts += 1
        if attempts > 10000:
            raise RuntimeError("failed to generate enough detected events after 10,000 batches")
        z = sample_z(rng, 4 * N_DETECTED)
        d = rng.normal(C_KM_S * z / h0, SIGMA)
        distances.extend(d[d >= THRESHOLD].tolist())
    return np.asarray(distances[:N_DETECTED], dtype=float)


def run(
    h0_mode,
    reps,
    seed,
    z_points,
    h0_points,
    h0_lo=50.0,
    h0_hi=90.0,
    key="selected_log_likelihood",
    *,
    return_records=False,
):
    """Run one truth regime using a reproducible RNG stream.

    Existing programmatic callers receive the same summary mapping by default.
    Setting return_records=True additionally returns one audit row per replicate.
    """
    if isinstance(reps, bool) or int(reps) != reps or reps < 2:
        raise ValueError("reps must be an integer of at least two")
    if isinstance(seed, bool) or int(seed) != seed or seed < 0:
        raise ValueError("seed must be a non-negative integer")
    if int(z_points) != z_points or z_points < 101:
        raise ValueError("z_points must be an integer >= 101")
    if int(h0_points) != h0_points or h0_points < 31:
        raise ValueError("h0_points must be an integer >= 31")
    if not np.isfinite(h0_lo) or not np.isfinite(h0_hi) or not 0 < h0_lo < h0_hi:
        raise ValueError("H0 bounds must be finite, positive, and ordered")
    if key not in {"selected_log_likelihood", "unconditional_log_likelihood"}:
        raise ValueError("invalid likelihood key")

    mode = mode_label(h0_mode)
    rng = np.random.default_rng(int(seed))
    sample_z = make_sampler()
    z = np.linspace(Z_MIN, Z_MAX, int(z_points))
    density = exact_density(z)
    grid = np.linspace(float(h0_lo), float(h0_hi), int(h0_points))
    errors, median_errors, cover68, cover90, ranks = [], [], [], [], []
    records = []

    for replicate in range(int(reps)):
        h0 = float(rng.uniform(h0_lo, h0_hi)) if mode == "prior" else float(mode.removeprefix("fixed_"))
        distances = simulate(rng, h0, sample_z)
        likelihood = selection_aware_toy_likelihood(
            z, density, grid, distances, sigma=SIGMA, detection_threshold=THRESHOLD
        )
        loglike = likelihood[key]
        s = posterior_summary(grid, loglike)
        unnormalized = np.exp(loglike - np.max(loglike))
        cdf = np.r_[0.0, np.cumsum(0.5 * (unnormalized[1:] + unnormalized[:-1]) * np.diff(grid))]
        cdf /= cdf[-1]
        rank = float(np.interp(h0, grid, cdf))
        q05, q95 = (float(x) for x in np.interp([0.05, 0.95], cdf, grid))
        hit68 = bool(s["ci68_low"] <= h0 <= s["ci68_high"])
        hit90 = bool(q05 <= h0 <= q95)

        errors.append(s["posterior_mean"] - h0)
        median_errors.append(s["posterior_median"] - h0)
        cover68.append(hit68)
        cover90.append(hit90)
        ranks.append(rank)

        records.append({
            "mode": mode,
            "replicate": int(replicate),
            "mode_seed": int(seed),
            "h0_true": float(h0),
            "distances_json": json.dumps(distances.tolist(), separators=(",", ":")),
            "posterior_mean": float(s["posterior_mean"]),
            "posterior_median": float(s["posterior_median"]),
            "ci68_low": float(s["ci68_low"]),
            "ci68_high": float(s["ci68_high"]),
            "ci90_low": q05,
            "ci90_high": q95,
            "posterior_integral": float(s["posterior_integral"]),
            "rank": rank,
            "covered_68": hit68,
            "covered_90": hit90,
            "likelihood_component": key,
        })

    e = np.asarray(errors, dtype=float)
    med_error = np.asarray(median_errors, dtype=float)
    r = np.asarray(ranks, dtype=float)
    counts, _ = np.histogram(r, bins=10, range=(0.0, 1.0))
    x2 = float(np.sum((counts - reps / 10.0) ** 2 / (reps / 10.0)))
    result = {
        "h0": h0_mode,
        "mode": mode,
        "mode_seed": int(seed),
        "reps": int(reps),
        "z_points": int(z_points),
        "h0_points": int(h0_points),
        "mean_bias": float(e.mean()),
        "mean_bias_se": float(e.std(ddof=1) / np.sqrt(reps)),
        "median_bias": float(med_error.mean()),
        "median_bias_se": float(med_error.std(ddof=1) / np.sqrt(reps)),
        "rmse": float(np.sqrt(np.mean(e ** 2))),
        "cov68": float(np.mean(cover68)),
        "cov90": float(np.mean(cover90)),
        "pit_chi2_9df": x2,
        "pit_p": float(chi2.sf(x2, 9)),
        "pit_counts": counts.tolist(),
        "likelihood_component": key,
    }
    return (result, records) if return_records else result


def run_campaign(modes, reps, seed, z_points, h0_points, key="selected_log_likelihood"):
    """Run the four preregistered regimes with disjoint deterministic seed ranges."""
    normalized_modes = [mode_label(m) for m in modes]
    if len(set(normalized_modes)) != len(normalized_modes):
        raise ValueError("modes must be unique")
    if set(normalized_modes) != set(MODE_OFFSETS):
        raise ValueError(f"this locked campaign requires exactly these modes: {sorted(MODE_OFFSETS)}")
    mode_seeds = {mode: int(seed + offset) for mode, offset in MODE_OFFSETS.items()}
    values_by_label = {"prior": "prior", "fixed_60": 60.0, "fixed_70": 70.0, "fixed_80": 80.0}
    summaries, all_records = [], []
    for label in normalized_modes:
        summary, rows = run(
            values_by_label[label], reps, mode_seeds[label], z_points, h0_points,
            key=key, return_records=True
        )
        summaries.append(summary)
        all_records.extend(rows)
        print({k: summary[k] for k in ("mode", "mode_seed", "reps", "mean_bias", "mean_bias_se", "cov68", "cov90")}, flush=True)

    configuration = {
        "model_scope": MODEL_SCOPE,
        "seed_start": int(seed),
        "mode_seeds": mode_seeds,
        "replicates_per_mode": int(reps),
        "h0_prior": [50.0, 90.0],
        "h0_points": int(h0_points),
        "z_bounds": [Z_MIN, Z_MAX],
        "z_points": int(z_points),
        "sigma_mpc": SIGMA,
        "detection_threshold_mpc": THRESHOLD,
        "n_detected_events": N_DETECTED,
        "sampler_grid_points": SAMPLER_GRID_POINTS,
        "likelihood_component": key,
        "truth_modes": {"prior": None, "fixed_60": 60.0, "fixed_70": 70.0, "fixed_80": 80.0},
    }
    artifact_summary = {
        "status": "synthetic_exact_density_oracle_calibration_not_real_H0_measurement",
        "git_commit": os.environ.get("GITHUB_SHA"),
        "configuration": configuration,
        "results": summaries,
        "interpretation": (
            "Prior-drawn ranks and coverage assess internal calibration under the stipulated toy model. "
            "Fixed-truth regimes assess conditional bias/coverage. Neither validates realistic detector response."
        ),
    }
    return artifact_summary, all_records


def write_artifacts(summary, records, output_json):
    """Write JSON summary, raw per-replicate CSV, and SHA-256/byte-count manifest."""
    output_json = Path(output_json)
    output_json.parent.mkdir(parents=True, exist_ok=True)
    csv_path = output_json.with_name(output_json.stem + "_replicates.csv")
    manifest_path = output_json.with_name(output_json.stem + "_manifest.json")
    summary_path = output_json
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    if not records:
        raise ValueError("cannot write an empty replicate artifact")
    fieldnames = list(records[0])
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(records)

    files = {}
    for path in (summary_path, csv_path):
        data = path.read_bytes()
        files[path.name] = {"sha256": hashlib.sha256(data).hexdigest(), "bytes": len(data)}
    manifest = {
        "schema_version": 1,
        "git_commit": summary.get("git_commit"),
        "record_count": len(records),
        "configuration": summary["configuration"],
        "files": files,
    }
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    return {"summary": summary_path, "replicates": csv_path, "manifest": manifest_path}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--reps", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--out", type=Path, required=True)
    parser.add_argument("--modes", nargs="+", default=["70", "60", "80", "prior"])
    parser.add_argument("--z-points", type=int, default=401)
    parser.add_argument("--h0-points", type=int, default=81)
    args = parser.parse_args()
    parsed_modes = [m if m == "prior" else float(m) for m in args.modes]
    summary, records = run_campaign(parsed_modes, args.reps, args.seed, args.z_points, args.h0_points)
    paths = write_artifacts(summary, records, args.out)
    print(json.dumps({
        "status": summary["status"], "git_commit": summary["git_commit"],
        "replicate_rows": len(records), "files": {k: str(v) for k, v in paths.items()}
    }, indent=2), flush=True)


if __name__ == "__main__":
    main()
