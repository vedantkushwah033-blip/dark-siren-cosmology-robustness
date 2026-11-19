"""Adaptive-CDF reference audit for boundary-concentrated full-parent H0 SBC draws.

This targeted diagnostic independently regenerates selected prior-predictive
replicates and compares the uniform-grid posterior ranks/quantiles with
adaptive quadrature. It does not change the production likelihood or tolerances.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad, trapezoid
from scipy.optimize import brentq, minimize_scalar
from scipy.special import log_ndtr, logsumexp
from scipy.stats import norm

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299792.458
PROBS = np.array([0.025, 0.05, 0.10, 0.16, 0.25, 0.50, 0.75, 0.84, 0.90, 0.95, 0.975])


def _seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([int(master), int(replicate), int(stream)]).generate_state(1, dtype=np.uint32)[0])


def _log_posterior(h0: float, z: np.ndarray, logw: np.ndarray,
                   distances: np.ndarray, sigma: float, threshold: float) -> float:
    if h0 <= 0:
        return -np.inf
    mu = C_KM_S * z / h0
    log_beta = logsumexp(log_ndtr((mu - threshold) / sigma) + logw)
    log_num = 0.0
    for distance in distances:
        log_num += logsumexp(norm.logpdf(distance, loc=mu, scale=sigma) + logw)
    return float(log_num - len(distances) * log_beta)


def _adaptive_reference(log_density, lower: float, upper: float, truth: float,
                        coarse_grid: np.ndarray) -> dict:
    """Integrate a positive unnormalized density and invert its CDF adaptively."""
    values = np.asarray([log_density(float(x)) for x in coarse_grid])
    if not np.isfinite(values).any():
        raise ArithmeticError("posterior log density is non-finite on the reference grid")
    imax = int(np.nanargmax(values))
    candidates = [(float(values[0]), lower), (float(values[-1]), upper),
                  (float(values[imax]), float(coarse_grid[imax]))]
    left = float(coarse_grid[max(0, imax - 1)])
    right = float(coarse_grid[min(len(coarse_grid) - 1, imax + 1)])
    if right > left:
        opt = minimize_scalar(lambda x: -log_density(float(x)), bounds=(left, right),
                              method="bounded", options={"xatol": 1e-10})
        candidates.append((-float(opt.fun), float(opt.x)))
    log_scale, mode = max(candidates)
    if not np.isfinite(log_scale):
        raise ArithmeticError("could not identify a finite posterior scale")

    # Fine knots near both prior boundaries, plus local refinement around every
    # coarse-grid local maximum. The broad interior is integrated adaptively.
    near_lo = np.linspace(lower, min(lower + 5.0, upper), 201)
    near_hi = np.linspace(max(upper - 5.0, lower), upper, 201)
    broad = np.linspace(lower, upper, 141)
    local_peaks = np.flatnonzero((values[1:-1] >= values[:-2]) &
                                 (values[1:-1] >= values[2:])) + 1
    local_knots = []
    for idx in local_peaks:
        center = float(coarse_grid[idx])
        local_knots.extend(np.linspace(max(lower, center - 0.5),
                                       min(upper, center + 0.5), 41))
    knots = np.unique(np.concatenate([near_lo, near_hi, broad,
                                      np.asarray(local_knots, dtype=float)]))
    knots = knots[(knots >= lower) & (knots <= upper)]
    knots[0], knots[-1] = lower, upper

    def scaled(x: float) -> float:
        value = log_density(float(x)) - log_scale
        return float(np.exp(value)) if value > -745 else 0.0

    masses = np.empty(len(knots) - 1)
    errors = np.empty(len(knots) - 1)
    for i, (a, b) in enumerate(zip(knots[:-1], knots[1:])):
        masses[i], errors[i] = quad(scaled, float(a), float(b),
                                    epsabs=2e-13, epsrel=2e-10, limit=100)
    cumulative = np.concatenate([[0.0], np.cumsum(masses)])
    total = float(cumulative[-1])
    total_error = float(errors.sum())
    if not np.isfinite(total) or total <= 0:
        raise ArithmeticError("adaptive posterior normalization failed")

    def integral_to(x: float) -> float:
        x = float(np.clip(x, lower, upper))
        if x <= lower:
            return 0.0
        if x >= upper:
            return total
        idx = int(np.searchsorted(knots, x, side="right") - 1)
        partial, _ = quad(scaled, float(knots[idx]), x,
                          epsabs=2e-13, epsrel=2e-10, limit=100)
        return float(cumulative[idx] + partial)

    quantiles = []
    for probability in PROBS:
        target = probability * total
        idx = min(int(np.searchsorted(cumulative, target, side="right") - 1),
                  len(masses) - 1)
        a, b = float(knots[idx]), float(knots[idx + 1])
        root = brentq(lambda x: cumulative[idx] + quad(
            scaled, a, float(x), epsabs=2e-13, epsrel=2e-10, limit=100
        )[0] - target, a, b, xtol=2e-9, rtol=1e-12)
        quantiles.append(float(root))

    boundary_masses = {}
    for fraction in (0.01, 0.025, 0.05):
        width = (upper - lower) * fraction
        boundary_masses[f"mass_lower_{fraction:.3f}"] = integral_to(lower + width) / total
        boundary_masses[f"mass_upper_{fraction:.3f}"] = 1.0 - integral_to(upper - width) / total

    return {
        "adaptive_rank": integral_to(truth) / total,
        "adaptive_quantiles": quantiles,
        "adaptive_mode_h0": mode,
        "adaptive_log_scale": log_scale,
        "adaptive_normalization_scaled": total,
        "adaptive_relative_quad_error_bound": total_error / total,
        **boundary_masses,
    }


def audit(artifact_dir: Path, replicates=(100, 174, 156, 62), grid_points=2241) -> dict:
    root = Path(artifact_dir).resolve()
    manifest = json.loads((root / "manifest.json").read_text())
    required = {"sbc_replicates.csv", "sbc_summary.json", "design.json"}
    if not required.issubset(manifest):
        raise AssertionError("source manifest omits required SBC files")
    for name, digest in manifest.items():
        path = root / name
        if not path.is_file() or hashlib.sha256(path.read_bytes()).hexdigest() != digest:
            raise AssertionError(f"source manifest mismatch: {name}")
    design = json.loads((root / "design.json").read_text())
    raw = pd.read_csv(root / "sbc_replicates.csv")
    n = int(design["replicates"])
    if design.get("status") != "COMPLETE" or not set(replicates).issubset(set(raw.replicate.astype(int))):
        raise AssertionError("source artifact is incomplete or omits a requested replicate")
    lower, upper = map(float, design["h0_prior"])
    master = int(design["master_seed"])
    n_parent, n_detected = int(design["n_parent"]), int(design["n_detected"])
    grid = np.linspace(lower, upper, grid_points)
    rows = []
    for rep in replicates:
        archived = raw.loc[raw.replicate == rep].iloc[0]
        truth_seed, sim_seed = _seed(master, rep, 0), _seed(master, rep, 1)
        if int(archived.truth_seed) != truth_seed or int(archived.simulation_seed) != sim_seed:
            raise AssertionError(f"seed provenance mismatch for replicate {rep}")
        truth = float(np.random.default_rng(truth_seed).uniform(lower, upper))
        if not np.isclose(float(archived.h0_true), truth, rtol=0, atol=1e-12):
            raise AssertionError(f"truth draw mismatch for replicate {rep}")
        sim = generate_synthetic_population(seed=sim_seed, n_parent=n_parent,
                                             n_detected=n_detected, h0_true=truth)
        parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
        z = parent.redshift.to_numpy(float)
        zmin, zmax = map(float, cfg["redshift_domain"])
        weights = 1.0 + 1.5 * (z - zmin) / (zmax - zmin)
        logw = np.log(weights) - logsumexp(np.log(weights))
        distances = events.distance_mpc.to_numpy(float)
        sigma, threshold = float(cfg["sigma_mpc"]), float(cfg["detection_threshold_mpc"])
        logpost = lambda h: _log_posterior(h, z, logw, distances, sigma, threshold)
        reference = _adaptive_reference(logpost, lower, upper, truth, grid)

        # Independently reproduce the uniform-grid calculation at the requested
        # resolution, so the comparison is tied to the actual archived seed.
        log_grid = np.asarray([logpost(float(h)) for h in grid])
        density = np.exp(log_grid - np.max(log_grid))
        area = float(trapezoid(density, grid))
        density /= area
        cdf = np.concatenate([[0.0], np.cumsum((density[:-1] + density[1:]) *
                                                np.diff(grid) * 0.5)])
        cdf /= cdf[-1]
        grid_rank = float(np.interp(truth, grid, cdf))
        grid_quantiles = np.interp(PROBS, cdf, grid)
        row = {
            "replicate": int(rep), "truth_h0": truth,
            "grid_rank": grid_rank, "adaptive_rank": reference["adaptive_rank"],
            "absolute_rank_difference": abs(grid_rank - reference["adaptive_rank"]),
            "max_quantile_difference": float(np.max(np.abs(grid_quantiles - reference["adaptive_quantiles"]))),
            "grid_mode_h0": float(grid[int(np.argmax(log_grid))]),
            "adaptive_mode_h0": reference["adaptive_mode_h0"],
            "grid_median": float(grid_quantiles[5]),
            "adaptive_median": float(reference["adaptive_quantiles"][5]),
            "adaptive_relative_quad_error_bound": reference["adaptive_relative_quad_error_bound"],
            **{f"adaptive_q{int(round(p*1000)):04d}": q for p, q in zip(PROBS, reference["adaptive_quantiles"])},
            **{k: v for k, v in reference.items() if k.startswith("mass_")},
        }
        rows.append(row)
    frame = pd.DataFrame(rows)
    summary = {
        "status": "DIAGNOSTIC_COMPLETE",
        "source_replicates": n,
        "selected_replicates": list(map(int, replicates)),
        "source_manifest_entries_verified": len(manifest),
        "grid_points": int(grid_points),
        "max_absolute_rank_difference": float(frame.absolute_rank_difference.max()),
        "max_quantile_difference_km_s_Mpc": float(frame.max_quantile_difference.max()),
        "maximum_adaptive_relative_quad_error_bound": float(frame.adaptive_relative_quad_error_bound.max()),
        "scope": "Targeted adaptive-CDF diagnostic for selected synthetic full-parent SBC boundary cases; not a complete 200/1000-replicate grid-convergence pass.",
    }
    return {"summary": summary, "replicates": frame.to_dict("records")}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--replicates", type=int, nargs="+", default=[100, 174, 156, 62])
    parser.add_argument("--grid-points", type=int, default=2241)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.artifact_dir, tuple(args.replicates), args.grid_points)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\\n")
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
