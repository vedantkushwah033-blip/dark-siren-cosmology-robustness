"""Audit Arm-E posterior-mean KDE density convergence over coefficient quadrature order.

This is a density-only numerical audit on matched synthetic curvature-pilot inputs.
It does not evaluate H0 likelihoods, selection calibration, or coverage.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import expit

from dark_siren_robustness.catalogue_density import estimate_host_density
from dark_siren_robustness.generative_selection import latent_galaxy_intensity
from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature

Z_MIN, Z_MAX = 0.01, 0.40


def _find_cells(root: Path) -> dict[float, Path]:
    found = {}
    for path in root.iterdir():
        cfg_path = path / "config.json"
        if not cfg_path.is_file():
            continue
        cfg = json.loads(cfg_path.read_text())
        if cfg.get("catalogue_selection_model") != "quadratic_logit":
            continue
        gamma = float(cfg["gamma"])
        if gamma in found:
            raise AssertionError(f"duplicate gamma artifact: {gamma}")
        found[gamma] = path
    return found


def _density_for_order(validation: pd.DataFrame, catalogue_z: np.ndarray,
                       grid: np.ndarray, merger: np.ndarray, order: int,
                       width: float, center: float, scale: float) -> tuple[np.ndarray, dict]:
    fit = fit_inclusion_posterior_quadrature(
        validation.redshift.to_numpy(float),
        validation.retained.to_numpy(int),
        degree=1,
        coordinate_center=center,
        coordinate_scale=scale,
        order=order,
    )
    x = (catalogue_z - fit.center) / fit.scale
    design = np.column_stack((np.ones_like(x), x))
    probability = expit(design @ fit.nodes.T).T
    if probability.shape != (fit.nodes.shape[0], catalogue_z.size):
        raise AssertionError("node-specific inclusion probability shape mismatch")
    if not np.all(np.isfinite(probability)) or np.any(probability <= 0.0):
        raise ArithmeticError("non-finite or non-positive inclusion probability")
    node_densities = np.stack([
        estimate_host_density(
            grid, catalogue_z, sample_weights=1.0 / probability[j],
            merger_weight=merger, absolute_bandwidth_sd=width,
        )
        for j in range(probability.shape[0])
    ])
    mean_density = fit.weights @ node_densities
    area = float(trapezoid(mean_density, grid))
    if not np.isfinite(area) or not np.isclose(area, 1.0, atol=1e-9, rtol=1e-9):
        raise ArithmeticError(f"posterior-mean KDE does not integrate to one: {area}")
    metadata = {
        "order": order,
        "node_count": int(fit.nodes.shape[0]),
        "quadrature_ess": float(fit.quadrature_ess),
        "max_quadrature_weight": float(fit.maximum_normalized_weight),
        "log_evidence": float(fit.log_evidence_unnormalized),
        "density_integral": area,
    }
    return mean_density, metadata


def run(args) -> dict:
    root = args.artifact_root.resolve()
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise ValueError("output directory must be empty")
    out.mkdir(parents=True, exist_ok=True)
    cells = _find_cells(root)
    gammas = tuple(float(v) for v in args.gammas.split(","))
    if not set(gammas).issubset(cells):
        raise AssertionError(f"missing gamma artifacts: {sorted(set(gammas)-set(cells))}")
    orders = tuple(int(v) for v in args.orders.split(","))
    if len(orders) < 3 or tuple(sorted(set(orders))) != orders or orders[-1] < 11:
        raise ValueError("orders must be strictly increasing, unique, and include order >= 11")
    sizes = tuple(int(v) for v in args.validation_sizes.split(","))
    widths = tuple(float(v) for v in args.widths.split(","))
    if not sizes or not widths or any(w <= 0 or not np.isfinite(w) for w in widths):
        raise ValueError("validation sizes and positive finite widths are required")

    grid = np.linspace(Z_MIN, Z_MAX, args.grid_points)
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    target_raw = latent_galaxy_intensity(grid) * merger
    target = target_raw / float(trapezoid(target_raw, grid))
    rows = []
    curves = {}
    expected_reps = int(args.expected_replicates)

    for gamma in gammas:
        cell = cells[gamma]
        cfg = json.loads((cell / "config.json").read_text())
        if int(cfg["replicates"]) != expected_reps:
            raise AssertionError(f"gamma={gamma}: expected {expected_reps} replicates")
        center, scale = 0.205, 0.195
        for rid in range(args.replicates_to_audit):
            child = cell / "replicates" / f"replicate_{rid:04d}"
            catalogue = pd.read_csv(child / "observable" / "catalogue_observed.csv")
            catalogue_z = catalogue.redshift.to_numpy(float)
            if len(catalogue_z) < 5:
                raise AssertionError(f"gamma={gamma}, replicate={rid}: catalogue too small")
            for n in sizes:
                validation_path = child / "observable" / f"validation_n{n}.csv"
                if not validation_path.is_file():
                    raise FileNotFoundError(validation_path)
                validation = pd.read_csv(validation_path)
                if len(validation) != n:
                    raise AssertionError(f"validation size mismatch: {validation_path}")
                for width in widths:
                    density_by_order = {}
                    metadata_by_order = {}
                    for order in orders:
                        density, metadata = _density_for_order(
                            validation, catalogue_z, grid, merger, order, width, center, scale
                        )
                        density_by_order[order] = density
                        metadata_by_order[order] = metadata
                        curves[f"g{gamma:+.1f}_r{rid:02d}_n{n}_w{width:.4f}_q{order}"] = density
                    reference = density_by_order[orders[-1]]
                    for order in orders[:-1]:
                        candidate = density_by_order[order]
                        rows.append({
                            "gamma": gamma, "replicate": rid, "validation_size": n,
                            "absolute_bandwidth_sd": width, "order": order,
                            "reference_order": orders[-1],
                            "max_abs_density_difference": float(np.max(np.abs(candidate-reference))),
                            "density_tv_distance": float(0.5 * trapezoid(np.abs(candidate-reference), grid)),
                            "iae_to_truth": float(trapezoid(np.abs(candidate-target), grid)),
                            "quadrature_ess": metadata_by_order[order]["quadrature_ess"],
                            "max_quadrature_weight": metadata_by_order[order]["max_quadrature_weight"],
                            "node_count": metadata_by_order[order]["node_count"],
                        })
                    for lower, upper in zip(orders[:-1], orders[1:]):
                        delta = density_by_order[upper] - density_by_order[lower]
                        rows.append({
                            "gamma": gamma, "replicate": rid, "validation_size": n,
                            "absolute_bandwidth_sd": width, "order": upper,
                            "reference_order": lower,
                            "max_abs_density_difference": float(np.max(np.abs(delta))),
                            "density_tv_distance": float(0.5 * trapezoid(np.abs(delta), grid)),
                            "iae_to_truth": float(trapezoid(np.abs(density_by_order[upper]-target), grid)),
                            "quadrature_ess": metadata_by_order[upper]["quadrature_ess"],
                            "max_quadrature_weight": metadata_by_order[upper]["max_quadrature_weight"],
                            "node_count": metadata_by_order[upper]["node_count"],
                        })

    frame = pd.DataFrame(rows)
    frame.to_csv(out / "density_quadrature_convergence.csv", index=False)
    np.savez_compressed(out / "posterior_mean_density_curves.npz", redshift_grid=grid, **curves)
    summary = (frame.groupby(["gamma", "validation_size", "absolute_bandwidth_sd", "order", "reference_order"], as_index=False)
        .agg(count=("max_abs_density_difference", "count"),
             max_abs_density_difference_mean=("max_abs_density_difference", "mean"),
             max_abs_density_difference_p95=("max_abs_density_difference", lambda x: float(np.quantile(x, .95))),
             max_abs_density_difference_max=("max_abs_density_difference", "max"),
             tv_distance_mean=("density_tv_distance", "mean"),
             tv_distance_p95=("density_tv_distance", lambda x: float(np.quantile(x, .95))),
             tv_distance_max=("density_tv_distance", "max")))
    summary.to_csv(out / "density_quadrature_convergence_summary.csv", index=False)
    report = {
        "status": "PASS",
        "gammas": list(gammas), "replicates_per_gamma": expected_reps,
        "validation_sizes": list(sizes), "absolute_bandwidths": list(widths),
        "orders": list(orders), "grid_points": args.grid_points,
        "replicates_audited_per_gamma": args.replicates_to_audit,\n        "case_count": int(args.replicates_to_audit * len(gammas) * len(sizes) * len(widths)),
        "comparison_rows": int(len(frame)),
        "scope": "Arm-E posterior-mean KDE density convergence only; synthetic low-z density stage; no H0 likelihood, posterior calibration, PIT, or coverage.",
        "reference_order": orders[-1],
    }
    (out / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2, sort_keys=True))
    return report


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_root", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--gammas", default="-2,0,2")
    p.add_argument("--orders", default="5,7,9,11,13")
    p.add_argument("--validation-sizes", default="100,300,1000")
    p.add_argument("--widths", default="0.016,0.024,0.032")
    p.add_argument("--expected-replicates", type=int, default=20)\n    p.add_argument("--replicates-to-audit", type=int, default=5)
    p.add_argument("--grid-points", type=int, default=401)
    args = p.parse_args()
    run(args)


if __name__ == "__main__":
    main()
