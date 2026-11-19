"""Independent 561-to-1121 H0-grid stability audit for the event-count study."""
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


def _seed(master: int, replicate: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate)]).generate_state(1, dtype=np.uint32)[0])


def _posterior(z, host_weights, distances, h0, sigma, threshold):
    z = np.asarray(z, float)
    w = np.asarray(host_weights, float)
    d = np.asarray(distances, float)
    h = np.asarray(h0, float)
    logw = np.log(w)
    logw -= logsumexp(logw)
    mu = 299792.458 * z[None, :] / h[:, None]
    log_det = logsumexp(norm.logcdf((mu-threshold)/sigma) + logw[None, :], axis=1)
    log_num = np.zeros(len(h), dtype=float)
    for distance in d:
        log_num += logsumexp(norm.logpdf(distance, loc=mu, scale=sigma) + logw[None, :], axis=1)
    logp = log_num - len(d)*log_det
    density = np.exp(logp - np.max(logp))
    area = float(trapezoid(density, h))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior normalization failed")
    density /= area
    cdf = np.r_[0.0, cumulative_trapezoid(density, h)]
    cdf /= cdf[-1]
    quantiles = {name: float(np.interp(q, cdf, h)) for name, q in (
        ("q025", .025), ("q050", .05), ("q160", .16), ("q500", .50),
        ("q840", .84), ("q950", .95), ("q975", .975),
    )}
    return density, quantiles


def audit(artifact_dir: Path, sample_size: int = 50, fine_points: int = 1121,
          max_quantile_shift: float = 0.05, max_density_l1: float = 0.001,
          output: Path | None = None) -> dict:
    root = Path(artifact_dir)
    manifest = json.loads((root/"manifest.json").read_text())
    expected_files = {"replicates.csv", "summary.csv", "paired_contrast.csv", "design.json"}
    if set(manifest) != expected_files:
        raise AssertionError("source manifest does not contain the exact four expected files")
    for name, digest in manifest.items():
        if hashlib.sha256((root/name).read_bytes()).hexdigest() != digest:
            raise AssertionError(f"source manifest mismatch: {name}")
    design = json.loads((root/"design.json").read_text())
    raw = pd.read_csv(root/"replicates.csv")
    nrep = int(design["replicates"])
    if design.get("status") != "COMPLETE" or nrep != 500 or design.get("event_counts") != [4, 30]:
        raise AssertionError("unexpected source design; expected complete 500-replicate 4/30-event study")
    if len(raw) != 2*nrep or raw.duplicated(["replicate", "event_count"]).any():
        raise AssertionError("source replicate table is incomplete or duplicated")
    if raw.groupby("event_count").size().to_dict() != {4: nrep, 30: nrep}:
        raise AssertionError("source event-count arms have incomplete denominators")
    if set(raw.replicate.astype(int)) != set(range(nrep)):
        raise AssertionError("source replicate IDs do not match the planned range")
    if set(raw.event_count.astype(int)) != {4, 30}:
        raise AssertionError("source artifact contains unexpected event-count arms")
    if sample_size < 2 or sample_size > nrep or fine_points <= int(design["h0_grid_points"]):
        raise ValueError("invalid sample size or finer grid resolution")
    prior_min, prior_max = map(float, design["uniform_h0_prior"])
    coarse = np.linspace(prior_min, prior_max, int(design["h0_grid_points"]))
    fine = np.linspace(prior_min, prior_max, fine_points)
    truth = float(design["h0_true"])
    n_parent = int(design["n_parent"])
    master = int(design["master_seed"])
    reps = np.unique(np.linspace(0, nrep-1, sample_size, dtype=int))
    rows = []
    for rep in reps:
        block = raw[raw.replicate == rep]
        seed = _seed(master, int(rep))
        if block.simulation_seed.nunique() != 1 or int(block.simulation_seed.iloc[0]) != seed:
            raise AssertionError(f"seed mismatch at replicate {rep}")
        sim = generate_synthetic_population(seed=seed, n_parent=n_parent, n_detected=30, h0_true=truth)
        parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
        if len(parent) != n_parent or len(events) != 30:
            raise AssertionError(f"regenerated dimensions mismatch at replicate {rep}")
        z = parent.redshift.to_numpy(float)
        zmin, zmax = map(float, cfg["redshift_domain"])
        weights = 1.0 + 1.5*(z-zmin)/(zmax-zmin)
        distances = events.distance_mpc.to_numpy(float)
        for count in (4, 30):
            archived = block[block.event_count == count]
            if len(archived) != 1:
                raise AssertionError(f"missing event-count arm at replicate {rep}, count {count}")
            dc, qc = _posterior(z, weights, distances[:count], coarse,
                                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]))
            df, qf = _posterior(z, weights, distances[:count], fine,
                                float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"]))
            source = archived.iloc[0]
            source_quantiles = {
                "q025": "q025", "q050": "q050", "q160": "q160",
                "q500": "h0_median", "q840": "q840", "q950": "q950", "q975": "q975",
            }
            for quantile, column in source_quantiles.items():
                if not np.isclose(float(source[column]), qc[quantile], rtol=0, atol=2e-8):
                    raise AssertionError(
                        f"independent 561-grid {quantile} reconstruction mismatch "
                        f"at rep {rep}, count {count}"
                    )
            if not np.isclose(
                float(source.event_prefix_distance_sum), float(distances[:count].sum()),
                rtol=0, atol=1e-8
            ):
                raise AssertionError(f"event-distance prefix mismatch at rep {rep}, count {count}")
            shifts = {k: abs(qf[k]-qc[k]) for k in qc}
            coarse_on_fine = np.interp(fine, coarse, dc)
            coarse_on_fine /= trapezoid(coarse_on_fine, fine)
            l1 = float(trapezoid(np.abs(df-coarse_on_fine), fine))
            rows.append({"replicate": int(rep), "event_count": count,
                         "max_quantile_shift": max(shifts.values()), "density_l1": l1,
                         **{f"shift_{k}": v for k, v in shifts.items()}})
    table = pd.DataFrame(rows)
    summary = {}
    for count in (4, 30):
        g = table[table.event_count == count]
        summary[str(count)] = {
            "posterior_curves": int(len(g)),
            "max_quantile_shift_km_s_Mpc": float(g.max_quantile_shift.max()),
            "median_quantile_shift_km_s_Mpc": float(g.max_quantile_shift.median()),
            "max_density_L1": float(g.density_l1.max()),
            "median_density_L1": float(g.density_l1.median()),
        }
    result = {
        "status": "PASS" if all(
            v["max_quantile_shift_km_s_Mpc"] <= max_quantile_shift and v["max_density_L1"] <= max_density_l1
            for v in summary.values()
        ) else "FAIL",
        "source_manifest_entries_verified": len(manifest),
        "source_replicates": nrep,
        "replicates_audited": len(reps),
        "posterior_curves_reconstructed": len(table),
        "coarse_grid_points": len(coarse),
        "fine_grid_points": len(fine),
        "max_quantile_shift_tolerance_km_s_Mpc": max_quantile_shift,
        "max_density_L1_tolerance": max_density_l1,
        "by_event_count": summary,
        "scope": "Synthetic full-parent fixed-truth numerical grid-stability diagnostic; not SBC or observational H0 inference.",
    }
    if output:
        Path(output).parent.mkdir(parents=True, exist_ok=True)
        Path(output).write_text(json.dumps(result, indent=2, sort_keys=True)+"\n")
    if result["status"] != "PASS":
        raise AssertionError(json.dumps(result, indent=2))
    return result


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir", type=Path)
    p.add_argument("--sample-size", type=int, default=50)
    p.add_argument("--fine-points", type=int, default=1121)
    p.add_argument("--max-quantile-shift", type=float, default=0.05)
    p.add_argument("--max-density-l1", type=float, default=0.001)
    p.add_argument("--output", type=Path, required=True)
    a = p.parse_args()
    print(json.dumps(audit(a.artifact_dir, a.sample_size, a.fine_points,
                           a.max_quantile_shift, a.max_density_l1, a.output), indent=2))


if __name__ == "__main__":
    main()
