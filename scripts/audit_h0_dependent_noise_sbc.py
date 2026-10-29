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
from scipy.stats import norm

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


def audit(output_dir):
    folder = Path(output_dir)
    manifest = json.loads((folder/"manifest.json").read_text(encoding="utf-8"))
    for name, item in manifest["files"].items():
        path = folder/name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != item["sha256"]:
            raise ValueError(f"SHA-256 mismatch/missing file: {name}")
    summary = json.loads((folder/"summary.json").read_text(encoding="utf-8"))
    config = summary["configuration"]
    if summary.get("status") != "synthetic_h0_dependent_noise_selection_sbc_not_real_h0_measurement":
        raise ValueError("unexpected scientific scope label")
    if manifest.get("configuration") != config:
        raise ValueError("manifest and summary configurations disagree")
    z = np.asarray(config["population_redshifts"], dtype=float)
    raw_w = np.asarray(config["population_weights_raw"], dtype=float)
    rows = list(csv.DictReader((folder/"replicates.csv").open(newline="", encoding="utf-8")))
    expected_n = config["replicates_per_mode"]
    if len(rows) != expected_n * len(MODES):
        raise ValueError("unexpected raw replicate count")
    seen = set()
    maximum_error = {key: 0.0 for key in ("posterior_mean", "posterior_median", "ci68_low", "ci68_high", "ci90_low", "ci90_high", "rank")}
    for row in rows:
        mode, rep = row["mode"], int(row["replicate"])
        if mode not in MODES or rep < 0 or rep >= expected_n or (mode, rep) in seen:
            raise ValueError("invalid/duplicate replicate key")
        seen.add((mode, rep))
        seed = int(row["seed"])
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
    expected_keys = {(mode, rep) for mode in MODES for rep in range(expected_n)}
    if seen != expected_keys:
        raise ValueError("replicate table is not complete across modes")
    return {"status": "PASS", "raw_records": len(rows), "replicates_per_mode": expected_n,
            "files_hash_checked": len(manifest["files"]), "rows_reconstructed": len(rows),
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
