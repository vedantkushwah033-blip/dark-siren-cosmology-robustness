"""Targeted grid-resolution audit for difficult replicates in the 1,000-replicate SBC source.

The fixed replicate list was selected after inspecting source posterior widths and
low-H0 truths. This is a post-hoc numerical stress subset, not a random sample
or a full 1,000-replicate convergence certificate.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299_792.458
TARGET_REPLICATES = (208, 622, 538, 17, 158, 883, 89, 847, 705, 224)
QUANTILE_PROBS = np.array([.025, .05, .10, .16, .25, .50, .75, .84, .90, .95, .975])


def seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate), int(stream)]).generate_state(1, dtype=np.uint32)[0])


def posterior(z, host_weights, distances, grid, sigma, threshold, chunk_size=128):
    z = np.asarray(z, float)
    w = np.asarray(host_weights, float)
    d = np.asarray(distances, float)
    h = np.asarray(grid, float)
    logw = np.log(w) - logsumexp(np.log(w))
    logpost = np.empty(len(h), dtype=float)
    log_norm = np.log(sigma * np.sqrt(2.0 * np.pi))
    for start in range(0, len(h), chunk_size):
        hh = h[start:start + chunk_size]
        mu = C_KM_S * z[None, :] / hh[:, None]
        log_beta = logsumexp(log_ndtr((mu - threshold) / sigma) + logw[None, :], axis=1)
        log_num = np.zeros(len(hh), dtype=float)
        for distance in d:
            logpdf = -0.5 * ((distance - mu) / sigma) ** 2 - log_norm
            log_num += logsumexp(logpdf + logw[None, :], axis=1)
        logpost[start:start + len(hh)] = log_num - len(d) * log_beta
    density = np.exp(logpost - np.max(logpost))
    density /= trapezoid(density, h)
    cdf = np.r_[0.0, cumulative_trapezoid(density, h)]
    cdf /= cdf[-1]
    quantiles = np.interp(QUANTILE_PROBS, cdf, h)
    return density, cdf, quantiles


def density_l1(grid_a, density_a, grid_b, density_b):
    interp = np.interp(grid_b, grid_a, density_a)
    interp /= trapezoid(interp, grid_b)
    return float(trapezoid(np.abs(density_b - interp), grid_b))


def audit(source_dir: Path, output_dir: Path, intermediate_points=2241,
          fine_points=4481, extra_points=8961) -> dict:
    source_dir = Path(source_dir).resolve()
    output_dir = Path(output_dir).resolve()
    output_dir.mkdir(parents=True, exist_ok=True)
    manifest = json.loads((source_dir / "manifest.json").read_text())
    for name, digest in manifest.items():
        path = source_dir / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError(f"source manifest mismatch: {name}")
    design = json.loads((source_dir / "design.json").read_text())
    independent_audit_path = source_dir / "independent_audit.json"
    if independent_audit_path.is_file():
        independent_audit = json.loads(independent_audit_path.read_text())
        if independent_audit.get("status") != "PASS" or int(independent_audit.get("replicates_regenerated", -1)) != 1000:
            raise AssertionError("source independent audit is not a complete 1000-replicate PASS")
    else:
        raise AssertionError("source artifact omits the independent 1000-replicate audit")
    raw = pd.read_csv(source_dir / "sbc_replicates.csv")
    if design.get("status") != "COMPLETE" or int(design.get("replicates", -1)) != 1000 or len(raw) != 1000:
        raise AssertionError("expected the complete immutable 1,000-replicate source")
    source_points = int(design["h0_grid_points"])
    if not source_points < intermediate_points < fine_points < extra_points:
        raise ValueError("grid points must be strictly increasing")
    lo, hi = map(float, design["h0_prior"])
    grids = {
        source_points: np.linspace(lo, hi, source_points),
        intermediate_points: np.linspace(lo, hi, intermediate_points),
        fine_points: np.linspace(lo, hi, fine_points),
        extra_points: np.linspace(lo, hi, extra_points),
    }
    rows = []
    for rep in TARGET_REPLICATES:
        row = raw.loc[raw.replicate.astype(int) == rep].iloc[0]
        truth_seed = seed(int(design["master_seed"]), rep, 0)
        sim_seed = seed(int(design["master_seed"]), rep, 1)
        truth = float(np.random.default_rng(truth_seed).uniform(lo, hi))
        if int(row.truth_seed) != truth_seed or int(row.simulation_seed) != sim_seed:
            raise AssertionError(f"seed provenance mismatch for replicate {rep}")
        if not np.isclose(float(row.h0_true), truth, rtol=0, atol=1e-12):
            raise AssertionError(f"prior draw mismatch for replicate {rep}")
        sim = generate_synthetic_population(
            seed=sim_seed, n_parent=int(design["n_parent"]),
            n_detected=int(design["n_detected"]), h0_true=truth,
        )
        parent = sim["parent_truth"]
        events = sim["gw_events_observed"]
        cfg = sim["config"]
        z = parent.redshift.to_numpy(float)
        zmin, zmax = map(float, cfg["redshift_domain"])
        host_weights = 1.0 + 1.5 * (z - zmin) / (zmax - zmin)
        distances = events.distance_mpc.to_numpy(float)
        results = {}
        for points, grid in grids.items():
            density, cdf, quantiles = posterior(
                z, host_weights, distances, grid,
                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]),
            )
            results[points] = (grid, density, cdf, quantiles)
        source_grid, source_density, source_cdf, source_q = results[source_points]
        stored_q = np.array([
            row.q025, row.q05, row.q10, row.q16, row.q25, row.h0_median,
            row.q75, row.q84, row.q90, row.q95, row.q975,
        ], dtype=float)
        source_rank = float(np.interp(truth, source_grid, source_cdf))
        if not np.allclose(source_q, stored_q, rtol=0, atol=2e-8):
            raise AssertionError(f"source quantile reproduction failed for replicate {rep}")
        if not np.isclose(source_rank, float(row.posterior_rank), rtol=0, atol=2e-8):
            raise AssertionError(f"source rank reproduction failed for replicate {rep}")
        record = {
            "replicate": rep, "h0_true": truth,
            "source_quantile_max_abs_error": float(np.max(np.abs(source_q - stored_q))),
            "source_rank_abs_error": abs(source_rank - float(row.posterior_rank)),
            "source_q95_width": float(row.q975 - row.q025),
        }
        ordered = sorted(results)
        for a, b in zip(ordered[:-1], ordered[1:]):
            ga, da, ca, qa = results[a]
            gb, db, cb, qb = results[b]
            record[f"{a}_to_{b}_density_l1"] = density_l1(ga, da, gb, db)
            record[f"{a}_to_{b}_max_quantile_shift"] = float(np.max(np.abs(qb - qa)))
            record[f"{a}_to_{b}_rank_shift"] = abs(float(np.interp(truth, gb, cb)) - float(np.interp(truth, ga, ca)))
        rows.append(record)
    frame = pd.DataFrame(rows)
    summary = {
        "status": "DIAGNOSTIC_COMPLETE",
        "scope": "Post-hoc targeted numerical stress subset from the 1,000-replicate synthetic full-parent SBC source; not a full-source convergence audit.",
        "source_manifest_entries_verified": len(manifest),
        "source_replicates": 1000,
        "targeted_replicates": list(TARGET_REPLICATES),
        "source_grid_points": source_points,
        "grid_points": [source_points, intermediate_points, fine_points, extra_points],
        "source_reconstruction_max_quantile_error": float(frame.source_quantile_max_abs_error.max()),
        "source_reconstruction_max_rank_error": float(frame.source_rank_abs_error.max()),
        "comparisons": {},
    }
    ordered = sorted(grids)
    for a, b in zip(ordered[:-1], ordered[1:]):
        key = f"{a}_to_{b}"
        summary["comparisons"][key] = {
            "max_density_l1": float(frame[f"{key}_density_l1"].max()),
            "max_quantile_shift_km_s_Mpc": float(frame[f"{key}_max_quantile_shift"].max()),
            "max_rank_shift": float(frame[f"{key}_rank_shift"].max()),
            "worst_density_replicate": int(frame.loc[frame[f"{key}_density_l1"].idxmax(), "replicate"]),
        }
    frame.to_csv(output_dir / "sbc_1000_boundary_grid_replicates.csv", index=False)
    (output_dir / "sbc_1000_boundary_grid_summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    return summary


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--intermediate-points", type=int, default=2241)
    parser.add_argument("--fine-points", type=int, default=4481)
    parser.add_argument("--extra-points", type=int, default=8961)
    args = parser.parse_args()
    audit(args.source_dir, args.output_dir, args.intermediate_points, args.fine_points, args.extra_points)


if __name__ == "__main__":
    main()
