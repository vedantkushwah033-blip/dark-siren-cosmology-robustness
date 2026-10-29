"""Independent replay and posterior/aggregate audit for the exact-density oracle study.

This module deliberately re-implements the density, sampler, likelihood, posterior
integration and summary calculations. It does not import the model runner or its
likelihood/posterior helpers.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
from pathlib import Path

import numpy as np
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import logsumexp
from scipy.stats import chi2, norm

C_KM_S = 299_792.458
SIGMA = 140.0
THRESHOLD = 250.0
N_DETECTED = 40
Z_MIN, Z_MAX, SAMPLER_GRID_POINTS = 0.01, 0.40, 40001
MODE_OFFSETS = {"prior": 0, "fixed_60": 100_000, "fixed_70": 200_000, "fixed_80": 300_000}
MODEL_SCOPE = "synthetic_low_redshift_exact_density_oracle_not_real_GW_inference"


def _label_truth(mode):
    if mode == "prior":
        return None
    if mode not in {"fixed_60", "fixed_70", "fixed_80"}:
        raise ValueError(f"unsupported mode {mode}")
    return float(mode.removeprefix("fixed_"))


def _target_unnormalized(z):
    # Independent copy of the declared base intensity times merger-weight shape.
    base = 1.0 + 2.2 * np.exp(-0.5 * ((z - 0.105) / 0.025) ** 2)
    base += 1.6 * np.exp(-0.5 * ((z - 0.275) / 0.040) ** 2)
    merger = 1.0 + 1.5 * (z - Z_MIN) / (Z_MAX - Z_MIN)
    return base * merger


def _normalized_target(z):
    unnormalized = _target_unnormalized(z)
    area = float(trapezoid(unnormalized, z))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("independent target-density normalization failed")
    return unnormalized / area


def _sample_cdf():
    grid = np.linspace(Z_MIN, Z_MAX, SAMPLER_GRID_POINTS)
    density = _normalized_target(grid)
    increments = 0.5 * (density[:-1] + density[1:]) * np.diff(grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf = cdf / cdf[-1]
    return grid, cdf


def _replay_distances(rng, h0_true, grid, cdf):
    accepted = []
    attempts = 0
    while len(accepted) < N_DETECTED:
        attempts += 1
        if attempts > 10000:
            raise ValueError("independent replay did not reach the detection-count target")
        zdraw = np.interp(rng.random(4 * N_DETECTED), cdf, grid)
        distance = rng.normal(C_KM_S * zdraw / h0_true, SIGMA)
        accepted.extend(distance[distance >= THRESHOLD].tolist())
    return np.asarray(accepted[:N_DETECTED], dtype=float)


def _quadrature_weights(grid):
    delta = np.diff(grid)
    weights = np.empty_like(grid)
    weights[0] = delta[0] / 2.0
    weights[-1] = delta[-1] / 2.0
    weights[1:-1] = (delta[:-1] + delta[1:]) / 2.0
    return weights


def _posterior_from_distances(distances, z, host_density, h0_grid, likelihood_component):
    """Independent likelihood implementation and posterior quadrature."""
    quad = _quadrature_weights(z)
    supported = (host_density > 0.0) & (quad > 0.0)
    log_mass = np.log(host_density[supported]) + np.log(quad[supported])
    z_supported = z[supported]
    means = C_KM_S * z_supported[None, :] / h0_grid[:, None]
    log_joint = np.zeros(h0_grid.size, dtype=float)
    for observed_distance in distances:
        terms = norm.logpdf(observed_distance, loc=means, scale=SIGMA) + log_mass[None, :]
        log_joint += logsumexp(terms, axis=1)
    beta_integrand = norm.cdf((C_KM_S * z[None, :] / h0_grid[:, None] - THRESHOLD) / SIGMA)
    beta = np.sum(host_density[None, :] * quad[None, :] * beta_integrand, axis=1)
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0.0) or np.any(beta > 1.0 + 1e-12):
        raise ValueError("independent selection efficiency is invalid")
    if likelihood_component == "selected_log_likelihood":
        log_likelihood = log_joint - distances.size * np.log(beta)
    elif likelihood_component == "unconditional_log_likelihood":
        log_likelihood = log_joint
    else:
        raise ValueError("unknown likelihood component in configuration")
    scaled = np.exp(log_likelihood - np.max(log_likelihood))
    area = float(trapezoid(scaled, h0_grid))
    if not np.isfinite(area) or area <= 0:
        raise ValueError("independent posterior normalization failed")
    posterior = scaled / area
    cdf = np.r_[0.0, cumulative_trapezoid(posterior, h0_grid)]
    cdf /= cdf[-1]
    q05, q16, q50, q84, q95 = np.interp([0.05, 0.16, 0.50, 0.84, 0.95], cdf, h0_grid)
    return {
        "posterior_mean": float(trapezoid(h0_grid * posterior, h0_grid)),
        "posterior_median": float(q50),
        "ci68_low": float(q16),
        "ci68_high": float(q84),
        "ci90_low": float(q05),
        "ci90_high": float(q95),
        "posterior_integral": float(trapezoid(posterior, h0_grid)),
        "cdf": cdf,
    }


def _as_bool(value, label):
    rendered = str(value).strip().lower()
    if rendered in {"true", "1"}:
        return True
    if rendered in {"false", "0"}:
        return False
    raise ValueError(f"{label} must be Boolean")


def _wilson(successes, trials):
    z = 1.959963984540054
    p = successes / trials
    z2 = z * z
    d = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / d
    half = z * np.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / d
    return [float(max(0.0, center - half)), float(min(1.0, center + half))]


def _recompute_summary(mode, rows, expected_n):
    h0 = np.asarray([row["h0_true"] for row in rows], dtype=float)
    mean_error = np.asarray([row["posterior_mean"] - row["h0_true"] for row in rows], dtype=float)
    median_error = np.asarray([row["posterior_median"] - row["h0_true"] for row in rows], dtype=float)
    covered68 = np.asarray([row["ci68_low"] <= row["h0_true"] <= row["ci68_high"] for row in rows])
    covered90 = np.asarray([row["ci90_low"] <= row["h0_true"] <= row["ci90_high"] for row in rows])
    ranks = np.asarray([row["rank"] for row in rows], dtype=float)
    counts, _ = np.histogram(ranks, bins=10, range=(0.0, 1.0))
    x2 = float(np.sum((counts - expected_n / 10.0) ** 2 / (expected_n / 10.0)))
    return {
        "mode": mode,
        "reps": expected_n,
        "mean_bias": float(mean_error.mean()),
        "mean_bias_se": float(mean_error.std(ddof=1) / np.sqrt(expected_n)),
        "median_bias": float(median_error.mean()),
        "median_bias_se": float(median_error.std(ddof=1) / np.sqrt(expected_n)),
        "rmse": float(np.sqrt(np.mean(mean_error ** 2))),
        "cov68": float(np.mean(covered68)),
        "cov90": float(np.mean(covered90)),
        "pit_chi2_9df": x2,
        "pit_p": float(chi2.sf(x2, 9)),
        "pit_counts": counts.tolist(),
    }


def audit(summary_path, raw_path=None, manifest_path=None):
    summary_path = Path(summary_path)
    raw_path = Path(raw_path) if raw_path else summary_path.with_name(summary_path.stem + "_replicates.csv")
    manifest_path = Path(manifest_path) if manifest_path else summary_path.with_name(summary_path.stem + "_manifest.json")
    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    config = summary.get("configuration")
    if not isinstance(config, dict) or config.get("model_scope") != MODEL_SCOPE:
        raise ValueError("unexpected oracle diagnostic model scope")
    if manifest.get("configuration") != config or manifest.get("git_commit") != summary.get("git_commit"):
        raise ValueError("manifest and summary provenance/configuration disagree")
    expected_files = {summary_path.name, raw_path.name}
    if set(manifest.get("files", {})) != expected_files:
        raise ValueError("manifest must hash exactly the summary and raw replicate CSV")
    for name, item in manifest["files"].items():
        path = summary_path if name == summary_path.name else raw_path
        data = path.read_bytes()
        if hashlib.sha256(data).hexdigest() != item.get("sha256") or len(data) != item.get("bytes"):
            raise ValueError(f"SHA-256/byte-count mismatch: {name}")

    reps = config.get("replicates_per_mode")
    base_seed = config.get("seed_start")
    if isinstance(reps, bool) or not isinstance(reps, int) or reps < 2:
        raise ValueError("replicates_per_mode must be an integer >= 2")
    if isinstance(base_seed, bool) or not isinstance(base_seed, int) or base_seed < 0:
        raise ValueError("seed_start must be a non-negative integer")
    if config.get("h0_prior") != [50.0, 90.0] or config.get("z_bounds") != [Z_MIN, Z_MAX]:
        raise ValueError("H0/redshift domain differs from the locked oracle protocol")
    if config.get("sigma_mpc") != SIGMA or config.get("detection_threshold_mpc") != THRESHOLD:
        raise ValueError("noise/threshold differs from the locked oracle protocol")
    if config.get("n_detected_events") != N_DETECTED or config.get("sampler_grid_points") != SAMPLER_GRID_POINTS:
        raise ValueError("event count or sampler resolution differs from the locked oracle protocol")
    if config.get("likelihood_component") not in {"selected_log_likelihood", "unconditional_log_likelihood"}:
        raise ValueError("invalid likelihood component")
    z_points, h0_points = config.get("z_points"), config.get("h0_points")
    if isinstance(z_points, bool) or not isinstance(z_points, int) or z_points < 101:
        raise ValueError("z_points must be an integer >= 101")
    if isinstance(h0_points, bool) or not isinstance(h0_points, int) or h0_points < 31:
        raise ValueError("h0_points must be an integer >= 31")
    expected_seeds = {mode: base_seed + offset for mode, offset in MODE_OFFSETS.items()}
    if config.get("mode_seeds") != expected_seeds:
        raise ValueError("mode-specific seeds do not match the disjoint preregistered offsets")
    if manifest.get("record_count") != reps * len(MODE_OFFSETS):
        raise ValueError("manifest record count does not match the protocol")

    with raw_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        required = {
            "mode", "replicate", "mode_seed", "h0_true", "distances_json",
            "posterior_mean", "posterior_median", "ci68_low", "ci68_high",
            "ci90_low", "ci90_high", "posterior_integral", "rank",
            "covered_68", "covered_90", "likelihood_component",
        }
        if reader.fieldnames is None or not required.issubset(reader.fieldnames):
            raise ValueError("raw replicate CSV is missing required fields")
        raw_rows = list(reader)
    if len(raw_rows) != reps * len(MODE_OFFSETS):
        raise ValueError("raw replicate count is incomplete")
    expected_keys = {(mode, rep) for mode in MODE_OFFSETS for rep in range(reps)}
    seen = set()
    by_mode = {mode: [] for mode in MODE_OFFSETS}
    max_error = {key: 0.0 for key in (
        "h0_true", "distances", "posterior_mean", "posterior_median",
        "ci68_low", "ci68_high", "ci90_low", "ci90_high", "rank",
    )}
    sampler_grid, sampler_cdf = _sample_cdf()
    z = np.linspace(Z_MIN, Z_MAX, z_points)
    density = _normalized_target(z)
    h0_grid = np.linspace(50.0, 90.0, h0_points)
    summary_by_mode = {item.get("mode"): item for item in summary.get("results", [])}
    if set(summary_by_mode) != set(MODE_OFFSETS):
        raise ValueError("summary does not contain exactly the four locked truth regimes")

    # Sort before sequential replay. Each mode has a unique independent RNG stream.
    parsed = []
    for row in raw_rows:
        mode = row.get("mode")
        try:
            rep = int(row["replicate"])
            seed = int(row["mode_seed"])
        except (KeyError, ValueError) as exc:
            raise ValueError("invalid replicate index or seed") from exc
        key = (mode, rep)
        if key not in expected_keys or key in seen:
            raise ValueError(f"invalid/duplicate replicate key: {key}")
        seen.add(key)
        if seed != expected_seeds[mode]:
            raise ValueError(f"mode seed mismatch at {key}")
        if row.get("likelihood_component") != config["likelihood_component"]:
            raise ValueError(f"likelihood component mismatch at {key}")
        try:
            true = float(row["h0_true"])
            distances = np.asarray(json.loads(row["distances_json"]), dtype=float)
            saved = {key_name: float(row[key_name]) for key_name in max_error}
        except (TypeError, ValueError, KeyError, json.JSONDecodeError) as exc:
            raise ValueError(f"invalid numeric/data field at {key}") from exc
        if not np.isfinite(true) or not np.all(np.isfinite(distances)) or distances.shape != (N_DETECTED,):
            raise ValueError(f"non-finite truth/distances or wrong event count at {key}")
        if np.any(distances < THRESHOLD):
            raise ValueError(f"saved distance does not satisfy the detection rule at {key}")
        if not all(np.isfinite(v) for v in saved.values()):
            raise ValueError(f"non-finite saved posterior metric at {key}")
        if abs(float(row["posterior_integral"]) - 1.0) > 2e-10:
            raise ValueError(f"posterior integral failure at {key}")
        parsed.append((mode, rep, row, true, distances, saved))

    if seen != expected_keys:
        raise ValueError("raw replicate table is not complete across mode/replicate keys")
    parsed.sort(key=lambda item: (list(MODE_OFFSETS).index(item[0]), item[1]))
    # Replay each regime independently and independently reconstruct the posterior.
    rng_by_mode = {mode: np.random.default_rng(expected_seeds[mode]) for mode in MODE_OFFSETS}
    for mode, rep, row, saved_true, saved_distances, saved in parsed:
        rng = rng_by_mode[mode]
        expected_true = float(rng.uniform(50.0, 90.0)) if mode == "prior" else float(_label_truth(mode))
        replayed_distances = _replay_distances(rng, expected_true, sampler_grid, sampler_cdf)
        delta_true = abs(saved_true - expected_true)
        delta_distance = float(np.max(np.abs(saved_distances - replayed_distances)))
        max_error["h0_true"] = max(max_error["h0_true"], delta_true)
        max_error["distances"] = max(max_error["distances"], delta_distance)
        if delta_true > 1e-12 or delta_distance > 1e-10:
            raise ValueError(f"data replay mismatch at {(mode, rep)}")
        reconstructed = _posterior_from_distances(
            replayed_distances, z, density, h0_grid, config["likelihood_component"]
        )
        rank = float(np.interp(expected_true, h0_grid, reconstructed["cdf"]))
        reconstructed["rank"] = rank
        reconstructed["h0_true"] = expected_true
        reconstructed["mode"] = mode
        reconstructed["covered_68"] = reconstructed["ci68_low"] <= expected_true <= reconstructed["ci68_high"]
        reconstructed["covered_90"] = reconstructed["ci90_low"] <= expected_true <= reconstructed["ci90_high"]
        for field in ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci90_low", "ci90_high", "rank"):
            delta = abs(float(saved[field]) - reconstructed[field])
            max_error[field] = max(max_error[field], delta)
            if not np.isfinite(delta) or delta > 1e-7:
                raise ValueError(f"{field} reconstruction mismatch at {(mode, rep)}")
        if _as_bool(row["covered_68"], "covered_68") != reconstructed["covered_68"]:
            raise ValueError(f"covered_68 flag mismatch at {(mode, rep)}")
        if _as_bool(row["covered_90"], "covered_90") != reconstructed["covered_90"]:
            raise ValueError(f"covered_90 flag mismatch at {(mode, rep)}")
        by_mode[mode].append(reconstructed)

    summary_numeric = (
        "mean_bias", "mean_bias_se", "median_bias", "median_bias_se",
        "rmse", "cov68", "cov90", "pit_chi2_9df", "pit_p",
    )
    for mode, rows in by_mode.items():
        rows.sort(key=lambda item: item["h0_true"] if False else item["posterior_mean"])  # sorting is not used for the statistics below
        expected = _recompute_summary(mode, rows, reps)
        saved_summary = summary_by_mode[mode]
        if saved_summary.get("mode_seed") != expected_seeds[mode] or saved_summary.get("reps") != reps:
            raise ValueError(f"saved summary configuration mismatch for {mode}")
        for field in summary_numeric:
            value = float(saved_summary.get(field, float("nan")))
            if not np.isfinite(value) or not np.isclose(value, expected[field], rtol=1e-9, atol=1e-9):
                raise ValueError(f"summary aggregate {field} mismatch for {mode}")
        if saved_summary.get("pit_counts") != expected["pit_counts"]:
            raise ValueError(f"summary aggregate pit_counts mismatch for {mode}")

    return {
        "status": "PASS",
        "raw_records": len(raw_rows),
        "replicates_per_mode": reps,
        "rows_replayed": len(raw_rows),
        "posteriors_independently_reconstructed": len(raw_rows),
        "summary_aggregates_reconstructed": True,
        "files_hash_checked": len(manifest["files"]),
        "maximum_absolute_reconstruction_error": max_error,
        "claim_scope": "synthetic_model_internal_validation_only",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("summary", type=Path)
    parser.add_argument("--raw", type=Path, default=None)
    parser.add_argument("--manifest", type=Path, default=None)
    parser.add_argument("--report", type=Path, default=None)
    args = parser.parse_args()
    report = audit(args.summary, args.raw, args.manifest)
    report_path = args.report or args.summary.with_name(args.summary.stem + "_audit.json")
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True), flush=True)


if __name__ == "__main__":
    main()
