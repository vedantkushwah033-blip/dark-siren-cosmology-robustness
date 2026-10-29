"""Independent generator replay and likelihood reconstruction for H0-noise SBC artifacts."""
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
PRIOR_MIN, PRIOR_MAX = 50.0, 90.0
FRACTIONAL_SIGMA, SIGMA_FLOOR_MPC = 0.12, 8.0
THRESHOLD_MPC, N_EVENTS, H0_POINTS = 300.0, 3, 161
MODES = {"prior": None, "fixed_60": 60.0, "fixed_70": 70.0, "fixed_80": 80.0}


def normalized_weights(raw):
    values = np.asarray(raw, dtype=float)
    scaled = values / np.max(values)
    return scaled / np.sum(scaled)


def replay_detected_distances(rng, h0_true, redshifts, weights, n_events=N_EVENTS):
    """Separate implementation of the discrete-host rejection sampler."""
    means = C_KM_S * redshifts / h0_true
    sigmas = np.hypot(SIGMA_FLOOR_MPC, FRACTIONAL_SIGMA * means)
    accepted = []
    proposals = 0
    maximum = max(100_000, n_events * 20_000)
    while len(accepted) < n_events and proposals < maximum:
        remaining = n_events - len(accepted)
        batch = min(max(256, 2 * remaining), maximum - proposals)
        host_idx = rng.choice(len(redshifts), size=batch, replace=True, p=weights)
        values = rng.normal(means[host_idx], sigmas[host_idx])
        hits = values >= THRESHOLD_MPC
        if np.any(hits):
            accepted.extend(float(x) for x in values[hits][:remaining])
        proposals += batch
    if len(accepted) != n_events:
        raise RuntimeError("independent simulator did not generate the requested detected events")
    return accepted


def independently_reconstruct_posterior(distances, redshifts, raw_weights, h0_true):
    """Independent Normal-mixture and beta calculation, not the production likelihood."""
    weights = normalized_weights(raw_weights)
    # The event candidate weights are the normalized population weights returned by the simulator.
    event_weights = normalized_weights(weights)
    log_event_weights = np.log(event_weights)
    h0_grid = np.linspace(PRIOR_MIN, PRIOR_MAX, H0_POINTS)
    host_means = C_KM_S * redshifts[None, :] / h0_grid[:, None]
    host_sigma = np.hypot(SIGMA_FLOOR_MPC, FRACTIONAL_SIGMA * host_means)
    total_loglike = np.zeros(h0_grid.size)
    for distance in distances:
        components = norm.logpdf(distance, loc=host_means, scale=host_sigma) + log_event_weights[None, :]
        total_loglike += logsumexp(components, axis=1)
    beta = np.empty(h0_grid.size)
    for i in range(h0_grid.size):
        beta[i] = np.sum(weights * norm.sf((THRESHOLD_MPC - host_means[i]) / host_sigma[i]))
    if not np.all(np.isfinite(beta)) or np.any(beta <= 0):
        raise ArithmeticError("independent detection efficiency invalid")
    selected = total_loglike - len(distances) * np.log(beta)
    density = np.exp(selected - np.max(selected))
    density /= trapezoid(density, h0_grid)
    cdf = np.r_[0.0, cumulative_trapezoid(density, h0_grid)]
    cdf /= cdf[-1]
    q05, q16, med, q84, q95 = np.interp([.05, .16, .50, .84, .95], cdf, h0_grid)
    return {
        "posterior_mean": float(trapezoid(h0_grid*density, h0_grid)),
        "posterior_median": float(med), "ci68_low": float(q16), "ci68_high": float(q84),
        "ci90_low": float(q05), "ci90_high": float(q95),
        "rank": float(np.interp(h0_true, h0_grid, cdf)),
        "posterior_integral": float(trapezoid(density, h0_grid)),
    }



def _as_bool(value, name):
    """Parse booleans as stored by csv.DictWriter and fail closed on unknown values."""
    if isinstance(value, (bool, np.bool_)):
        return bool(value)
    rendered = str(value).strip().lower()
    if rendered in {"true", "1"}:
        return True
    if rendered in {"false", "0"}:
        return False
    raise ValueError(f"{name} must be a Boolean value.")


def _wilson_interval(successes, trials):
    """Independent Wilson score interval used for coverage summaries."""
    if trials < 1 or not 0 <= successes <= trials:
        raise ValueError("invalid coverage counts in reconstructed summary")
    z_score = 1.959963984540054
    p = successes / trials
    z2 = z_score * z_score
    denominator = 1.0 + z2 / trials
    center = (p + z2 / (2.0 * trials)) / denominator
    half = z_score * np.sqrt(p * (1.0 - p) / trials + z2 / (4.0 * trials * trials)) / denominator
    return [float(max(0.0, center - half)), float(min(1.0, center + half))]


def _validate_locked_configuration(config, manifest):
    """Ensure the artifact really represents the protocol this auditor implements."""
    expected_z = np.linspace(0.025, 0.16, 12)
    expected_w = np.linspace(1.0, 2.5, 12)
    if config.get("prior") != [PRIOR_MIN, PRIOR_MAX]:
        raise ValueError("configuration prior does not match the locked audit contract")
    if int(config.get("h0_points", -1)) != H0_POINTS:
        raise ValueError("configuration H0 grid does not match the locked audit contract")
    if int(config.get("n_detected_events", -1)) != N_EVENTS:
        raise ValueError("configuration event count does not match the locked audit contract")
    if not np.isclose(float(config.get("fractional_sigma", np.nan)), FRACTIONAL_SIGMA, rtol=0, atol=0):
        raise ValueError("configuration fractional noise does not match the locked audit contract")
    if not np.isclose(float(config.get("sigma_floor_mpc", np.nan)), SIGMA_FLOOR_MPC, rtol=0, atol=0):
        raise ValueError("configuration noise floor does not match the locked audit contract")
    if not np.isclose(float(config.get("detection_distance_min_mpc", np.nan)), THRESHOLD_MPC, rtol=0, atol=0):
        raise ValueError("configuration detection threshold does not match the locked audit contract")
    z = np.asarray(config.get("population_redshifts", []), dtype=float)
    w = np.asarray(config.get("population_weights_raw", []), dtype=float)
    if z.shape != expected_z.shape or not np.allclose(z, expected_z, rtol=0, atol=1e-14):
        raise ValueError("configuration population redshifts do not match the locked audit contract")
    if w.shape != expected_w.shape or not np.allclose(w, expected_w, rtol=0, atol=1e-14):
        raise ValueError("configuration population weights do not match the locked audit contract")
    expected_modes = {name: value for name, value in MODES.items()}
    if config.get("modes") != expected_modes:
        raise ValueError("configuration truth modes do not match the locked audit contract")
    replicates = config.get("replicates_per_mode")
    seed_start = config.get("seed_start")
    if isinstance(replicates, bool) or not isinstance(replicates, int) or replicates < 30:
        raise ValueError("configuration replicates_per_mode must be an integer >= 30")
    if isinstance(seed_start, bool) or not isinstance(seed_start, int) or seed_start < 0:
        raise ValueError("configuration seed_start must be a non-negative integer")
    expected_ends = {
        mode: seed_start + index * 100_000 + replicates - 1
        for index, mode in enumerate(MODES)
    }
    if config.get("seed_end_by_mode") != expected_ends:
        raise ValueError("configuration seed ranges do not match the locked offset scheme")
    if config.get("model_scope") != "synthetic_discrete_host_distance_only_model_not_real_GW_inference":
        raise ValueError("configuration model scope does not match the audit contract")
    if manifest.get("schema_version") != 1:
        raise ValueError("unsupported artifact manifest schema")
    if manifest.get("replicate_count") != replicates * len(MODES):
        raise ValueError("manifest replicate count does not match the locked configuration")
    if manifest.get("git_commit") != manifest.get("_summary_git_commit"):
        raise ValueError("manifest and summary commit identifiers disagree")


def _validate_summary_aggregates(summary, reconstructed_rows, expected_n):
    """Recompute every headline aggregate from independently reconstructed rows."""
    summary_rows = summary.get("summaries")
    if not isinstance(summary_rows, list):
        raise ValueError("summary.summaries must be a list")
    by_mode = {}
    for item in summary_rows:
        if not isinstance(item, dict) or item.get("mode") not in MODES or item["mode"] in by_mode:
            raise ValueError("summary contains an invalid or duplicate mode")
        by_mode[item["mode"]] = item
    if set(by_mode) != set(MODES):
        raise ValueError("summary is missing a truth mode")

    numeric_fields = (
        "bias", "bias_standard_error", "rmse", "coverage68", "coverage90",
        "mean_rank", "rank_standard_error", "rank_chi2_9df",
        "rank_chi2_p_descriptive", "mean_interval68_width", "mean_interval90_width",
    )
    for mode in MODES:
        records = [row for row in reconstructed_rows if row["mode"] == mode]
        if len(records) != expected_n:
            raise ValueError(f"independent reconstruction has wrong count for mode {mode}")
        truth = np.asarray([row["h0_true"] for row in records], dtype=float)
        means = np.asarray([row["posterior_mean"] for row in records], dtype=float)
        errors = means - truth
        hit68 = np.asarray([row["ci68_low"] <= row["h0_true"] <= row["ci68_high"] for row in records])
        hit90 = np.asarray([row["ci90_low"] <= row["h0_true"] <= row["ci90_high"] for row in records])
        ranks = np.asarray([row["rank"] for row in records], dtype=float)
        counts, _ = np.histogram(ranks, bins=10, range=(0.0, 1.0))
        expected = {
            "mode": mode,
            "n_replicates": expected_n,
            "bias": float(np.mean(errors)),
            "bias_standard_error": float(np.std(errors, ddof=1) / np.sqrt(expected_n)),
            "rmse": float(np.sqrt(np.mean(errors ** 2))),
            "coverage68": float(np.mean(hit68)),
            "coverage68_count": int(np.sum(hit68)),
            "coverage68_wilson95": _wilson_interval(int(np.sum(hit68)), expected_n),
            "coverage90": float(np.mean(hit90)),
            "coverage90_count": int(np.sum(hit90)),
            "coverage90_wilson95": _wilson_interval(int(np.sum(hit90)), expected_n),
            "mean_rank": float(np.mean(ranks)),
            "rank_standard_error": float(np.std(ranks, ddof=1) / np.sqrt(expected_n)),
            "rank_decile_counts": counts.tolist(),
            "rank_chi2_9df": float(np.sum((counts - expected_n / 10.0) ** 2 / (expected_n / 10.0))),
            "rank_chi2_p_descriptive": float(chi2.sf(
                np.sum((counts - expected_n / 10.0) ** 2 / (expected_n / 10.0)), 9
            )),
            "mean_interval68_width": float(np.mean([row["ci68_high"] - row["ci68_low"] for row in records])),
            "mean_interval90_width": float(np.mean([row["ci90_high"] - row["ci90_low"] for row in records])),
        }
        saved = by_mode[mode]
        for field in numeric_fields:
            try:
                actual = float(saved[field])
            except (KeyError, TypeError, ValueError) as exc:
                raise ValueError(f"summary field {field} is absent/invalid for mode {mode}") from exc
            if not np.isfinite(actual) or not np.isclose(actual, expected[field], rtol=1e-9, atol=1e-9):
                raise ValueError(f"summary aggregate {field} mismatch for mode {mode}")
        for field in ("n_replicates", "coverage68_count", "coverage90_count"):
            if saved.get(field) != expected[field]:
                raise ValueError(f"summary aggregate {field} mismatch for mode {mode}")
        if saved.get("rank_decile_counts") != expected["rank_decile_counts"]:
            raise ValueError(f"summary aggregate rank_decile_counts mismatch for mode {mode}")
        for field in ("coverage68_wilson95", "coverage90_wilson95"):
            if not np.allclose(np.asarray(saved.get(field), dtype=float), expected[field], rtol=1e-9, atol=1e-9):
                raise ValueError(f"summary aggregate {field} mismatch for mode {mode}")



def audit(output_dir):
    folder = Path(output_dir)
    manifest = json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    expected_files = {"replicates.csv", "summary.json"}
    if set(manifest.get("files", {})) != expected_files:
        raise ValueError("manifest must cover exactly replicates.csv and summary.json")
    for name, item in manifest["files"].items():
        path = folder/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"SHA-256 mismatch/missing file: {name}")
        if int(item.get("bytes", -1)) != path.stat().st_size:
            raise ValueError(f"manifest byte count mismatch: {name}")
    summary = json.loads((folder/"summary.json").read_text(encoding="utf-8"))
    config = summary["configuration"]
    if summary.get("status") != "synthetic_h0_dependent_noise_selection_sbc_not_real_h0_measurement":
        raise ValueError("unexpected scientific scope label")
    if manifest.get("configuration") != config:
        raise ValueError("manifest and summary configurations disagree")
    if manifest.get("git_commit") != summary.get("git_commit"):
        raise ValueError("manifest and summary commit identifiers disagree")
    manifest_for_lock = dict(manifest)
    manifest_for_lock["_summary_git_commit"] = summary.get("git_commit")
    _validate_locked_configuration(config, manifest_for_lock)
    z = np.asarray(config["population_redshifts"], dtype=float)
    raw_w = np.asarray(config["population_weights_raw"], dtype=float)
    rows = list(csv.DictReader((folder/"replicates.csv").open(newline="", encoding="utf-8")))
    expected_n = config["replicates_per_mode"]
    if len(rows) != expected_n * len(MODES):
        raise ValueError("unexpected raw replicate count")
    seen = set()
    reconstructed_rows = []
    maximum_error = {key: 0.0 for key in ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci90_low", "ci90_high", "rank")}
    for row in rows:
        mode, rep = row["mode"], int(row["replicate"])
        if mode not in MODES or rep < 0 or rep >= expected_n or (mode, rep) in seen:
            raise ValueError("invalid/duplicate replicate key")
        seen.add((mode, rep))
        mode_index = list(MODES).index(mode)
        expected_seed = int(config["seed_start"]) + mode_index * 100_000 + rep
        seed = int(row["seed"])
        if seed != expected_seed:
            raise ValueError(f"seed does not match the locked range for {(mode, rep)}")
        rng = np.random.default_rng(seed)
        fixed = MODES[mode]
        true = float(rng.uniform(PRIOR_MIN, PRIOR_MAX)) if fixed is None else float(fixed)
        if not np.isclose(float(row["h0_true"]), true, rtol=0, atol=1e-12):
            raise ValueError(f"truth replay failed at {(mode, rep)}")
        # Reproduce the population normalization used by the sampler, independent of production helper.
        q = normalized_weights(raw_w)
        distances = replay_detected_distances(rng, true, z, q)
        saved_distances = json.loads(row["distances_mpc_json"])
        if len(saved_distances) != N_EVENTS or not np.array_equal(
            np.asarray(saved_distances, dtype=float), np.asarray(distances, dtype=float)
        ):
            raise ValueError(f"distance replay failed at {(mode, rep)}")
        computed = independently_reconstruct_posterior(saved_distances, z, raw_w, true)
        for key in maximum_error:
            delta = abs(float(row[key]) - computed[key])
            maximum_error[key] = max(maximum_error[key], delta)
            if delta > 1e-7:
                raise ValueError(f"{key} mismatch at {(mode, rep)}: {delta}")
        if abs(float(row["posterior_integral"]) - 1.0) > 2e-10:
            raise ValueError(f"posterior integral failure at {(mode, rep)}")
        expected_covered_68 = computed["ci68_low"] <= true <= computed["ci68_high"]
        expected_covered_90 = computed["ci90_low"] <= true <= computed["ci90_high"]
        if _as_bool(row["covered_68"], "covered_68") != expected_covered_68:
            raise ValueError(f"covered_68 flag mismatch at {(mode, rep)}")
        if _as_bool(row["covered_90"], "covered_90") != expected_covered_90:
            raise ValueError(f"covered_90 flag mismatch at {(mode, rep)}")
        reconstructed_rows.append({"mode": mode, "h0_true": true, **computed})
    expected_keys = {(mode, rep) for mode in MODES for rep in range(expected_n)}
    if seen != expected_keys:
        raise ValueError("replicate table is not complete across modes")
    if int(manifest.get("replicate_count", -1)) != len(rows):
        raise ValueError("manifest replicate count does not match the raw CSV")
    _validate_summary_aggregates(summary, reconstructed_rows, expected_n)
    return {"status": "PASS", "raw_records": len(rows), "replicates_per_mode": expected_n,
            "files_hash_checked": len(manifest["files"]), "rows_reconstructed": len(rows),
            "summary_aggregates_reconstructed": True,
            "maximum_absolute_reconstruction_error": maximum_error,
            "claim_scope": "synthetic_model_internal_validation_only"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output_dir", type=Path)
    parser.add_argument("--output", type=Path, default=None)
    args = parser.parse_args()
    result = audit(args.output_dir)
    rendered = json.dumps(result, indent=2, sort_keys=True) + "\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    print(rendered, end="")
    raise SystemExit(0 if result["status"] == "PASS" else 1)


if __name__ == "__main__":
    main()
