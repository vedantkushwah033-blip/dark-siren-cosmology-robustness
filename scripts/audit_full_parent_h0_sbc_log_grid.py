"""Audit a geometric H0 grid against an independently reconstructed uniform-grid audit.

The prior remains Uniform(H0_min, H0_max) in H0. The geometric grid is only
a nonuniform set of H0 quadrature nodes; integration is still performed with
respect to dH0 using trapezoidal quadrature. No log-coordinate density is
formed, so no log-Jacobian is omitted or required.
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

from dark_siren_robustness.generative_selection import generate_synthetic_population

C_KM_S = 299_792.458
QUANTILE_PROBABILITIES = np.array(
    [0.025, 0.05, 0.10, 0.16, 0.25, 0.50, 0.75, 0.84, 0.90, 0.95, 0.975],
    dtype=float,
)


def verify_source_manifest(source_dir: Path) -> int:
    """Verify every file listed in the immutable source artifact manifest."""
    root = Path(source_dir).resolve()
    manifest_path = root / "manifest.json"
    if not manifest_path.is_file():
        raise FileNotFoundError(f"missing source manifest: {manifest_path}")
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if not manifest:
        raise AssertionError("source manifest is empty")
    for relative, expected in manifest.items():
        path = root / relative
        if not path.is_file():
            raise AssertionError(f"source manifest file is missing: {relative}")
        actual = hashlib.sha256(path.read_bytes()).hexdigest()
        if actual != expected:
            raise AssertionError(f"source manifest digest mismatch: {relative}")
    return len(manifest)


def normalize_posterior_on_h0_grid(
    h0_grid: np.ndarray, log_likelihood: np.ndarray
) -> tuple[np.ndarray, np.ndarray]:
    """Normalize likelihood under a uniform prior, integrating directly in H0."""
    grid = np.asarray(h0_grid, dtype=float)
    log_like = np.asarray(log_likelihood, dtype=float)
    if grid.ndim != 1 or len(grid) < 2 or not np.isfinite(grid).all():
        raise ValueError("H0 grid must be finite, one-dimensional and contain >=2 points")
    if np.any(grid <= 0) or np.any(np.diff(grid) <= 0):
        raise ValueError("H0 grid must be positive and strictly increasing")
    if log_like.shape != grid.shape or not np.isfinite(log_like).all():
        raise ValueError("log likelihood must be finite and match the H0 grid")
    scaled = np.exp(log_like - float(np.max(log_like)))
    area = float(trapezoid(scaled, grid))
    if not np.isfinite(area) or area <= 0:
        raise ArithmeticError("posterior normalization failed")
    density = scaled / area
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, grid)))
    if not np.isfinite(cdf[-1]) or cdf[-1] <= 0:
        raise ArithmeticError("posterior CDF normalization failed")
    cdf /= cdf[-1]
    return density, cdf


def log_posterior_on_grid(
    h0_grid: np.ndarray,
    redshift: np.ndarray,
    distances: np.ndarray,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    redshift_domain: tuple[float, float] = (0.01, 0.40),
) -> np.ndarray:
    """Vectorized full-parent selected-event log posterior on arbitrary H0 nodes."""
    grid = np.asarray(h0_grid, dtype=float)
    z = np.asarray(redshift, dtype=float)
    observed = np.asarray(distances, dtype=float)
    if z.ndim != 1 or len(z) < 1 or not np.isfinite(z).all() or np.any(z <= 0):
        raise ValueError("redshift must be a finite positive one-dimensional array")
    if observed.ndim != 1 or len(observed) < 1 or not np.isfinite(observed).all():
        raise ValueError("distances must be a finite one-dimensional array")
    if np.any(observed < detection_threshold_mpc):
        raise ValueError("all supplied events must satisfy the detection threshold")
    if not np.isfinite(sigma_mpc) or sigma_mpc <= 0:
        raise ValueError("sigma_mpc must be positive and finite")
    z_min, z_max = map(float, redshift_domain)
    if not 0 < z_min < z_max:
        raise ValueError("invalid redshift domain")
    host_weight = 1.0 + 1.5 * (z - z_min) / (z_max - z_min)
    log_weight = np.log(host_weight) - logsumexp(np.log(host_weight))
    mu = C_KM_S * z[None, :] / grid[:, None]
    log_beta = logsumexp(
        log_ndtr((mu - detection_threshold_mpc) / sigma_mpc) + log_weight[None, :],
        axis=1,
    )
    log_numerator = np.zeros(len(grid), dtype=float)
    log_norm = np.log(sigma_mpc * np.sqrt(2.0 * np.pi))
    for distance in observed:
        log_pdf = -0.5 * ((distance - mu) / sigma_mpc) ** 2 - log_norm
        log_numerator += logsumexp(log_pdf + log_weight[None, :], axis=1)
    return log_numerator - len(observed) * log_beta


def summarize_grid(
    h0_grid: np.ndarray,
    redshift: np.ndarray,
    distances: np.ndarray,
    truth: float,
    sigma_mpc: float,
    detection_threshold_mpc: float,
    redshift_domain: tuple[float, float] = (0.01, 0.40),
) -> dict[str, object]:
    log_post = log_posterior_on_grid(
        h0_grid, redshift, distances, sigma_mpc, detection_threshold_mpc,
        redshift_domain=redshift_domain,
    )
    density, cdf = normalize_posterior_on_h0_grid(h0_grid, log_post)
    quantiles = np.interp(QUANTILE_PROBABILITIES, cdf, h0_grid)
    return {
        "rank": float(np.interp(truth, h0_grid, cdf)),
        "quantiles": quantiles,
        "median": float(quantiles[5]),
        "density": density,
    }


def audit(
    source_dir: Path,
    baseline_csv: Path,
    output_dir: Path,
    *,
    grid_points: int = 2241,
) -> dict[str, object]:
    source = Path(source_dir).resolve()
    output = Path(output_dir).resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be empty")
    output.mkdir(parents=True, exist_ok=True)
    manifest_count = verify_source_manifest(source)
    design = json.loads((source / "design.json").read_text(encoding="utf-8"))
    raw = pd.read_csv(source / "sbc_replicates.csv")
    baseline = pd.read_csv(baseline_csv)

    required_design = {
        "status", "replicates", "master_seed", "n_parent", "n_detected",
        "h0_prior", "h0_grid_points",
    }
    if not required_design.issubset(design):
        raise AssertionError("source design omits required fields")
    if design["status"] != "COMPLETE":
        raise AssertionError("source SBC artifact is not complete")
    if len(raw) != int(design["replicates"]):
        raise AssertionError("source replicate count differs from design")
    if raw.replicate.astype(int).tolist() != list(range(len(raw))):
        raise AssertionError("source replicates must be complete and ordered from zero")
    baseline_required = {"replicate", "fine_grid_rank", "fine_grid_median"}
    if not baseline_required.issubset(baseline.columns):
        raise AssertionError("baseline audit CSV omits fine-grid rank or median")
    if baseline.replicate.astype(int).tolist() != list(range(len(raw))):
        raise AssertionError("baseline audit does not contain the same ordered replicate IDs")
    if len(baseline) != len(raw):
        raise AssertionError("baseline and source replicate counts differ")

    lower, upper = map(float, design["h0_prior"])
    if grid_points < 101:
        raise ValueError("grid_points must be at least 101")
    grid = np.geomspace(lower, upper, int(grid_points))
    rows = []
    for row in raw.itertuples(index=False):
        replicate = int(row.replicate)
        truth_seed = int(np.random.SeedSequence(
            [int(design["master_seed"]), replicate, 0]
        ).generate_state(1, dtype=np.uint32)[0])
        simulation_seed = int(np.random.SeedSequence(
            [int(design["master_seed"]), replicate, 1]
        ).generate_state(1, dtype=np.uint32)[0])
        if int(row.truth_seed) != truth_seed or int(row.simulation_seed) != simulation_seed:
            raise AssertionError(f"seed provenance mismatch in replicate {replicate}")
        truth = float(np.random.default_rng(truth_seed).uniform(lower, upper))
        if not np.isclose(float(row.h0_true), truth, rtol=0, atol=1e-12):
            raise AssertionError(f"truth draw mismatch in replicate {replicate}")
        simulation = generate_synthetic_population(
            seed=simulation_seed,
            n_parent=int(design["n_parent"]),
            n_detected=int(design["n_detected"]),
            h0_true=truth,
        )
        parent = simulation["parent_truth"]
        events = simulation["gw_events_observed"]
        config = simulation["config"]
        if len(parent) != int(design["n_parent"]) or len(events) != int(design["n_detected"]):
            raise AssertionError(f"regenerated dimensions mismatch in replicate {replicate}")
        if int(row.n_parent) != len(parent) or int(row.n_detected) != len(events):
            raise AssertionError(f"archived dimensions mismatch in replicate {replicate}")
        posterior = summarize_grid(
            grid,
            parent.redshift.to_numpy(float),
            events.distance_mpc.to_numpy(float),
            truth,
            float(config["sigma_mpc"]),
            float(config["detection_threshold_mpc"]),
            tuple(config["redshift_domain"]),
        )
        base = baseline.iloc[replicate]
        rows.append({
            "replicate": replicate,
            "h0_true": truth,
            "log_grid_rank": posterior["rank"],
            "uniform_4481_rank": float(base.fine_grid_rank),
            "abs_rank_difference": abs(posterior["rank"] - float(base.fine_grid_rank)),
            "log_grid_median": posterior["median"],
            "uniform_4481_median": float(base.fine_grid_median),
            "abs_median_difference": abs(posterior["median"] - float(base.fine_grid_median)),
            **{
                f"log_grid_q{int(round(p * 1000)):04d}": float(q)
                for p, q in zip(QUANTILE_PROBABILITIES, posterior["quantiles"])
            },
        })
    frame = pd.DataFrame(rows)
    frame.to_csv(output / "log_grid_replicates.csv", index=False)
    summary = {
        "status": "DIAGNOSTIC_COMPLETE_NOT_FULL_CONVERGENCE_PASS",
        "replicates": int(len(frame)),
        "grid_points": int(grid_points),
        "grid_type": "geometric H0 nodes; posterior integrated with respect to dH0",
        "prior": [lower, upper],
        "source_manifest_entries_verified": int(manifest_count),
        "max_abs_rank_difference_vs_uniform_4481": float(frame.abs_rank_difference.max()),
        "median_abs_rank_difference_vs_uniform_4481": float(frame.abs_rank_difference.median()),
        "max_abs_median_difference_vs_uniform_4481": float(frame.abs_median_difference.max()),
        "median_abs_median_difference_vs_uniform_4481": float(frame.abs_median_difference.median()),
        "max_rank_difference_tolerance": 0.002,
        "max_quantile_difference_tolerance_km_s_Mpc": 0.05,
        "density_l1_status": "NOT_EVALUATED_IN_THIS_FULL_SAMPLE_COMPARISON",
        "interpretation": (
            "A geometric H0 grid is a promising rank/median alternative on this "
            "200-replicate synthetic baseline. This comparison is against a "
            "4481-point uniform grid, not a complete independent adaptive "
            "reference for every replicate. It does not certify density-L1 "
            "convergence or catalogue-inclusion SBC."
        ),
    }
    (output / "log_grid_summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--source-dir", type=Path, required=True)
    parser.add_argument("--baseline-csv", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--grid-points", type=int, default=2241)
    args = parser.parse_args()
    print(json.dumps(audit(
        args.source_dir, args.baseline_csv, args.output_dir,
        grid_points=args.grid_points,
    ), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
