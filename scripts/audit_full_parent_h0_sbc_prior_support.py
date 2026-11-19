"""Prior-support sensitivity audit for low-H0 draws in the synthetic full-parent SBC run.

The archived truths were generated from Uniform(20, 300). This diagnostic keeps
those same data fixed and recomputes likelihood posteriors under wider uniform
supports [15, 300] and [10, 300]. These are posterior-sensitivity comparisons,
NOT SBC ranks under the original generating prior.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import quad

from audit_full_parent_h0_sbc_adaptive import (
    _adaptive_reference,
    _log_posterior,
    _seed,
)
from dark_siren_robustness.generative_selection import generate_synthetic_population


def _posterior_mass_below(log_density, reference: dict, lower: float, cutoff: float) -> float:
    """Integrate normalized posterior mass from lower to cutoff."""
    if cutoff <= lower:
        return 0.0
    if cutoff >= np.inf:
        return 1.0
    scale = float(reference["adaptive_log_scale"])
    total = float(reference["adaptive_normalization_scaled"])

    def scaled(x: float) -> float:
        value = float(log_density(float(x))) - scale
        return float(np.exp(value)) if value > -745 else 0.0

    numerator, _ = quad(scaled, lower, cutoff, epsabs=2e-13, epsrel=2e-10, limit=200)
    return float(np.clip(numerator / total, 0.0, 1.0))


def audit(artifact_dir: Path, output: Path, grid_points: int = 4481,
          selected_count: int = 10, alternative_lowers=(10.0, 15.0)) -> dict:
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
    if design.get("status") != "COMPLETE" or int(design["replicates"]) != 1000:
        raise AssertionError("expected the immutable complete 1000-replicate source")
    source_lower, upper = map(float, design["h0_prior"])
    if not np.isclose(source_lower, 20.0) or not np.isclose(upper, 300.0):
        raise AssertionError("unexpected source prior support; this audit is preregistered for [20, 300]")
    if "q975" not in raw:
        raise AssertionError("source replicate table lacks q975 for deterministic tail selection")

    selected = raw.nsmallest(selected_count, "q975").replicate.astype(int).tolist()
    master = int(design["master_seed"])
    n_parent, n_detected = int(design["n_parent"]), int(design["n_detected"])
    rows = []

    for rep in selected:
        archived = raw.loc[raw.replicate == rep].iloc[0]
        truth_seed = _seed(master, rep, 0)
        simulation_seed = _seed(master, rep, 1)
        if int(archived.truth_seed) != truth_seed or int(archived.simulation_seed) != simulation_seed:
            raise AssertionError(f"seed provenance mismatch for replicate {rep}")
        truth = float(np.random.default_rng(truth_seed).uniform(source_lower, upper))
        if not np.isclose(float(archived.h0_true), truth, rtol=0, atol=1e-12):
            raise AssertionError(f"truth draw mismatch for replicate {rep}")

        sim = generate_synthetic_population(
            seed=simulation_seed, n_parent=n_parent, n_detected=n_detected, h0_true=truth
        )
        parent, events, cfg = sim["parent_truth"], sim["gw_events_observed"], sim["config"]
        if len(parent) != n_parent or len(events) != n_detected:
            raise AssertionError(f"regenerated dimensions mismatch for replicate {rep}")
        z = parent.redshift.to_numpy(float)
        zmin, zmax = map(float, cfg["redshift_domain"])
        weights = 1.0 + 1.5 * (z - zmin) / (zmax - zmin)
        logw = np.log(weights)
        from scipy.special import logsumexp
        logw -= logsumexp(logw)
        distances = events.distance_mpc.to_numpy(float)
        sigma = float(cfg["sigma_mpc"])
        threshold = float(cfg["detection_threshold_mpc"])
        logpost = lambda h: _log_posterior(h, z, logw, distances, sigma, threshold)

        baseline_grid = np.linspace(source_lower, upper, grid_points)
        baseline = _adaptive_reference(logpost, source_lower, upper, truth, baseline_grid)
        baseline_q = np.asarray(baseline["adaptive_quantiles"], dtype=float)
        support_rows = []
        for lower in sorted(set(map(float, alternative_lowers)) | {source_lower}):
            if lower >= source_lower or lower <= 0 or lower >= upper:
                if np.isclose(lower, source_lower):
                    ref = baseline
                else:
                    raise ValueError("alternative lower bounds must be positive and below source lower bound")
            else:
                grid = np.linspace(lower, upper, grid_points)
                ref = _adaptive_reference(logpost, lower, upper, truth, grid)

            q = np.asarray(ref["adaptive_quantiles"], dtype=float)
            below20 = (
                _posterior_mass_below(logpost, ref, lower, source_lower)
                if lower < source_lower else 0.0
            )
            support_rows.append({
                "prior_support": [lower, upper],
                "posterior_median": float(q[5]),
                "posterior_q025": float(q[0]),
                "posterior_q975": float(q[-1]),
                "median_shift_vs_original": float(q[5] - baseline_q[5]),
                "q025_shift_vs_original": float(q[0] - baseline_q[0]),
                "q975_shift_vs_original": float(q[-1] - baseline_q[-1]),
                "posterior_mass_below_original_lower_20": below20,
                "adaptive_relative_quad_error_estimate": float(ref["adaptive_relative_quad_error_bound"]),
            })
        rows.append({
            "replicate": int(rep),
            "truth_h0": truth,
            "archived_q975": float(archived.q975),
            "original_support_posterior_median": float(baseline_q[5]),
            "original_support_posterior_q025": float(baseline_q[0]),
            "original_support_posterior_q975": float(baseline_q[-1]),
            "posterior_support_comparisons": support_rows,
        })

    comparison = [
        item for row in rows for item in row["posterior_support_comparisons"]
        if not np.isclose(item["prior_support"][0], source_lower)
    ]
    summary = {
        "status": "SUPPORT_SENSITIVITY_COMPLETE",
        "source_replicates": 1000,
        "selected_replicates": selected,
        "selected_count": len(selected),
        "source_prior": [source_lower, upper],
        "alternative_lower_bounds": sorted(set(map(float, alternative_lowers))),
        "grid_points_per_support": int(grid_points),
        "manifest_entries_verified": len(manifest),
        "max_absolute_median_shift_km_s_Mpc": max(abs(x["median_shift_vs_original"]) for x in comparison),
        "max_absolute_q025_shift_km_s_Mpc": max(abs(x["q025_shift_vs_original"]) for x in comparison),
        "max_absolute_q975_shift_km_s_Mpc": max(abs(x["q975_shift_vs_original"]) for x in comparison),
        "max_posterior_mass_below_20": max(x["posterior_mass_below_original_lower_20"] for x in comparison),
        "max_adaptive_relative_quad_error_estimate": max(x["adaptive_relative_quad_error_estimate"] for x in comparison),
        "interpretation_guardrail": (
            "The source truths were drawn from Uniform(20,300). Wider-support posteriors are "
            "conditional sensitivity analyses on the same data, not calibrated SBC posteriors; "
            "their posterior ranks must not be pooled with the original SBC ranks."
        ),
    }
    result = {"summary": summary, "replicates": rows}
    output = Path(output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text(json.dumps(result, indent=2, sort_keys=True) + "\n")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_dir", type=Path)
    parser.add_argument("--grid-points", type=int, default=4481)
    parser.add_argument("--selected-count", type=int, default=10)
    parser.add_argument("--alternative-lower-bounds", type=float, nargs="+", default=[10.0, 15.0])
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    result = audit(args.artifact_dir, args.output, args.grid_points,
                   args.selected_count, tuple(args.alternative_lower_bounds))
    print(json.dumps(result["summary"], indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
