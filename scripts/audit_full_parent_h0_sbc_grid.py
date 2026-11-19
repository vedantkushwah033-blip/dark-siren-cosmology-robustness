"""Independent H0-grid sensitivity audit for prior-predictive full-parent SBC."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.special import logsumexp
from scipy.stats import norm

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299792.458


def _seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate), int(stream)]).generate_state(1, dtype=np.uint32)[0])


def _posterior(z, weights, distances, grid, sigma, threshold):
    z = np.asarray(z, float)
    w = np.asarray(weights, float)
    d = np.asarray(distances, float)
    h = np.asarray(grid, float)
    logw = np.log(w) - logsumexp(np.log(w))
    mu = C_KM_S * z[None, :] / h[:, None]
    log_beta = logsumexp(norm.logcdf((mu-threshold)/sigma) + logw[None, :], axis=1)
    log_num = np.zeros(len(h), dtype=float)
    for distance in d:
        log_num += logsumexp(norm.logpdf(distance, loc=mu, scale=sigma) + logw[None, :], axis=1)
    log_post = log_num - len(d)*log_beta
    density = np.exp(log_post - np.max(log_post))
    area = float(trapezoid(density, h))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h)]
    cdf /= cdf[-1]
    probs = (0.025, 0.05, 0.10, 0.16, 0.25, 0.50, 0.75, 0.84, 0.90, 0.95, 0.975)
    quantiles = np.interp(probs, cdf, h)
    return density, cdf, quantiles


def audit(root: Path, fine_points: int = 1121, max_quantile_shift: float = 0.05,
          max_rank_shift: float = 0.002, max_density_l1: float = 0.001,
          output: Path | None = None) -> dict:
    root = Path(root).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    if not {"sbc_replicates.csv", "sbc_summary.json", "design.json"}.issubset(manifest):
        raise AssertionError("source manifest omits required SBC files")
    for name, digest in manifest.items():
        path = root / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError(f"source manifest mismatch: {name}")
    design = json.loads((root / "design.json").read_text())
    raw = pd.read_csv(root / "sbc_replicates.csv")
    source_summary = json.loads((root / "sbc_summary.json").read_text())
    independent_audit_path = root / "independent_audit.json"
    if independent_audit_path.is_file():
        independent_audit = json.loads(independent_audit_path.read_text())
        if independent_audit.get("status") != "PASS":
            raise AssertionError("source SBC independent reconstruction audit did not pass")
        if int(independent_audit.get("replicates_regenerated", -1)) != 200:
            raise AssertionError("source SBC independent audit has an incomplete denominator")
    n = int(design["replicates"])
    if design.get("status") != "COMPLETE" or n != 200 or len(raw) != n:
        raise AssertionError("expected the complete 200-replicate SBC source artifact")
    required_columns = {
        "replicate", "truth_seed", "simulation_seed", "h0_true", "n_parent",
        "n_detected", "posterior_rank", "h0_median", "median_bias",
        "q025", "q05", "q10", "q16", "q25", "q75", "q84", "q90", "q95", "q975",
        "covered50", "covered68", "covered80", "covered90", "covered95",
    }
    if required_columns - set(raw.columns):
        raise AssertionError("source SBC replicate table omits required rank/quantile/coverage columns")
    if source_summary.get("status") != "COMPLETE" or int(source_summary.get("replicates", -1)) != n:
        raise AssertionError("source SBC summary is incomplete or inconsistent")
    if int(source_summary.get("h0_grid_points", -1)) != int(design["h0_grid_points"]):
        raise AssertionError("source SBC summary grid metadata mismatch")
    for level in (50, 68, 80, 90, 95):
        empirical = float(raw[f"covered{level}"].mean())
        if not np.isclose(float(source_summary[f"coverage{level}"]), empirical, rtol=0, atol=1e-12):
            raise AssertionError(f"source SBC coverage{level} summary mismatch")
    if set(raw.replicate.astype(int)) != set(range(n)) or raw.replicate.duplicated().any():
        raise AssertionError("SBC replicate inventory is incomplete or duplicated")
    coarse_points = int(design["h0_grid_points"])
    if fine_points <= coarse_points:
        raise ValueError("fine grid must have more points than the source grid")
    lo, hi = map(float, design["h0_prior"])
    coarse = np.linspace(lo, hi, coarse_points)
    fine = np.linspace(lo, hi, fine_points)
    master = int(design["master_seed"])
    n_parent = int(design["n_parent"])
    n_detected = int(design["n_detected"])
    rows = []
    for rep in range(n):
        row = raw.loc[raw.replicate == rep].iloc[0]
        truth_seed = _seed(master, rep, 0)
        sim_seed = _seed(master, rep, 1)
        if int(row.truth_seed) != truth_seed or int(row.simulation_seed) != sim_seed:
            raise AssertionError(f"seed provenance mismatch at replicate {rep}")
        truth = float(np.random.default_rng(truth_seed).uniform(lo, hi))
        if not np.isclose(float(row.h0_true), truth, rtol=0, atol=1e-12):
            raise AssertionError(f"prior draw mismatch at replicate {rep}")
        sim = generate_synthetic_population(seed=sim_seed, n_parent=n_parent,
                                             n_detected=n_detected, h0_true=truth)
        parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
        if len(parent) != n_parent or len(events) != n_detected:
            raise AssertionError(f"regenerated dimensions mismatch at replicate {rep}")
        z = parent.redshift.to_numpy(float)
        zmin, zmax = map(float, cfg["redshift_domain"])
        weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
        distances = events.distance_mpc.to_numpy(float)
        dc, cc, qc = _posterior(z, weights, distances, coarse,
                                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]))
        df, cf, qf = _posterior(z, weights, distances, fine,
                                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]))
        rank_c = float(np.interp(truth, coarse, cc))
        rank_f = float(np.interp(truth, fine, cf))
        stored_q = np.array([row.q025, row.q05, row.q10, row.q16, row.q25,
                             row.h0_median, row.q75, row.q84, row.q90, row.q95, row.q975], float)
        if not np.allclose(qc, stored_q, rtol=0, atol=2e-8):
            raise AssertionError(f"independent source-grid quantile mismatch at replicate {rep}")
        if not np.isclose(rank_c, float(row.posterior_rank), rtol=0, atol=2e-8):
            raise AssertionError(f"independent source-grid rank mismatch at replicate {rep}")
        if not np.isclose(float(row.median_bias), float(qc[5]-truth), rtol=0, atol=2e-8):
            raise AssertionError(f"independent source-grid median-bias mismatch at replicate {rep}")
        for level, low_idx, high_idx in ((50, 4, 6), (68, 3, 7), (80, 2, 8),
                                         (90, 1, 9), (95, 0, 10)):
            expected_coverage = int(qc[low_idx] <= truth <= qc[high_idx])
            if int(row[f"covered{level}"]) != expected_coverage:
                raise AssertionError(f"source-grid coverage{level} mismatch at replicate {rep}")
        coarse_on_fine = np.interp(fine, coarse, dc)
        coarse_on_fine /= trapezoid(coarse_on_fine, fine)
        density_l1 = float(trapezoid(np.abs(df-coarse_on_fine), fine))
        rows.append({
            "replicate": rep,
            "source_grid_rank": rank_c,
            "fine_grid_rank": rank_f,
            "absolute_rank_shift": abs(rank_f-rank_c),
            "max_quantile_shift": float(np.max(np.abs(qf-qc))),
            "density_l1": density_l1,
            "source_grid_median": float(qc[5]),
            "fine_grid_median": float(qf[5]),
        })
    frame = pd.DataFrame(rows)
    summary = {
        "status": "PASS" if (
            frame.max_quantile_shift.max() <= max_quantile_shift
            and frame.absolute_rank_shift.max() <= max_rank_shift
            and frame.density_l1.max() <= max_density_l1
        ) else "FAIL",
        "source_manifest_entries_verified": len(manifest),
        "replicates_reconstructed": len(frame),
        "source_grid_points": coarse_points,
        "fine_grid_points": fine_points,
        "max_quantile_shift_tolerance_km_s_Mpc": max_quantile_shift,
        "max_rank_shift_tolerance": max_rank_shift,
        "max_density_L1_tolerance": max_density_l1,
        "max_quantile_shift_observed": float(frame.max_quantile_shift.max()),
        "median_quantile_shift_observed": float(frame.max_quantile_shift.median()),
        "max_absolute_rank_shift_observed": float(frame.absolute_rank_shift.max()),
        "median_absolute_rank_shift_observed": float(frame.absolute_rank_shift.median()),
        "max_density_L1_observed": float(frame.density_l1.max()),
        "median_density_L1_observed": float(frame.density_l1.median()),
        "scope": "Synthetic full-parent prior-predictive numerical grid audit; not catalogue-inclusion SBC or observational inference.",
    }
    worst_cols = ["replicate", "max_quantile_shift", "absolute_rank_shift", "density_l1"]
    summary["worst_quantile_shift_replicates"] = (
        frame.nlargest(10, "max_quantile_shift")[worst_cols].to_dict("records")
    )
    summary["worst_rank_shift_replicates"] = (
        frame.nlargest(10, "absolute_rank_shift")[worst_cols].to_dict("records")
    )
    summary["worst_density_l1_replicates"] = (
        frame.nlargest(10, "density_l1")[worst_cols].to_dict("records")
    )
    if output:
        output = Path(output)
        output.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(output.with_name("sbc_grid_audit_replicates.csv"), index=False)
        output.write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n")
    print(json.dumps(summary, indent=2, sort_keys=True))
    if summary["status"] != "PASS":
        raise AssertionError(json.dumps(summary, indent=2))
    return summary


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--fine-points", type=int, default=1121)
    p.add_argument("--max-quantile-shift", type=float, default=0.05)
    p.add_argument("--max-rank-shift", type=float, default=0.002)
    p.add_argument("--max-density-l1", type=float, default=0.001)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.artifact_dir, a.fine_points, a.max_quantile_shift,
                           a.max_rank_shift, a.max_density_l1, a.output), indent=2))


if __name__ == "__main__":
    main()
