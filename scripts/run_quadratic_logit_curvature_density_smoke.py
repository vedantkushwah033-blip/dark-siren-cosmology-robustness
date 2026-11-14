"""Density-only smoke runner for the parent-calibrated quadratic-logit inclusion DGP.

This is a synthetic curvature stress experiment. Arm C receives the true
quadratic inclusion probabilities as an oracle diagnostic; Arms D/E fit the
declared degree-one logistic model from independent validation outcomes.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from dark_siren_robustness.generative_selection import (
    generate_quadratic_logit_inclusion_samples,
    latent_galaxy_intensity,
)
from dark_siren_robustness.five_arm_comparison import (
    build_five_arm_weights,
    weight_concentration_diagnostics,
)
from dark_siren_robustness.catalogue_density import estimate_host_density, reflected_weighted_kde_bandwidth

Z_MIN, Z_MAX = 0.01, 0.40
DEFAULT_WIDTHS = (0.008, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080)
DEFAULT_GAMMAS = (-2.0, -1.0, -0.5, 0.0, 0.5, 1.0, 2.0)
INCLUSION_COORDINATE_CENTER = 0.205
INCLUSION_COORDINATE_SCALE = 0.195
INCLUSION_PRIOR_MEAN = (0.0, 0.0)
INCLUSION_PRIOR_SD = (5.0, 2.5)
BOUNDARY_LOWER, BOUNDARY_UPPER = 0.049, 0.361


def _hash(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1 << 20), b""):
            h.update(block)
    return h.hexdigest()


def _metric(density: np.ndarray, truth: np.ndarray, grid: np.ndarray) -> dict[str, float]:
    delta = np.asarray(density, dtype=float) - truth
    left = grid <= BOUNDARY_LOWER
    right = grid >= BOUNDARY_UPPER
    return {
        "iae": float(trapezoid(np.abs(delta), grid)),
        "ise": float(trapezoid(delta**2, grid)),
        "boundary_iae": float(
            trapezoid(np.abs(delta[left]), grid[left])
            + trapezoid(np.abs(delta[right]), grid[right])
        ),
    }


def run(args: argparse.Namespace) -> dict[str, object]:
    out = Path(args.output_dir).resolve()
    out.mkdir(parents=True, exist_ok=True)
    if any(out.iterdir()):
        raise ValueError("output directory must be empty")

    widths = tuple(sorted(set(float(x) for x in args.absolute_widths)))
    sizes = list(args.validation_sizes)
    if len(widths) < 2 or any(not np.isfinite(x) or x <= 0 for x in widths):
        raise ValueError("at least two distinct positive finite widths are required")
    if not sizes or sorted(set(sizes)) != sizes or sizes[0] < 5:
        raise ValueError("validation sizes must be strictly increasing and >= 5")
    if sizes[-1] > args.validation_frame:
        raise ValueError("validation size exceeds validation frame")
    if args.n_parent < 100 or args.quadrature_order < 3 or args.grid_points < 51:
        raise ValueError("invalid parent size, quadrature order, or grid resolution")
    if args.seed < 0 or args.validation_seed < 0:
        raise ValueError("seeds must be non-negative")
    if not np.isfinite(args.gamma):
        raise ValueError("gamma must be finite")

    for folder in ("observable", "truth", "results"):
        (out / folder).mkdir()

    sample = generate_quadratic_logit_inclusion_samples(
        population_seed=args.seed,
        validation_seed=args.validation_seed,
        n_parent=args.n_parent,
        n_validation_frame=args.validation_frame,
        n_validation=sizes[-1],
        gamma=args.gamma,
        beta=args.beta,
        z_min=Z_MIN,
        z_max=Z_MAX,
    )
    parent = sample["parent_truth"]
    catalogue = sample["catalogue_observed"]
    validation_max = sample["validation_observed"]
    validation_truth = sample["validation_frame_truth"]

    parent.to_csv(out / "truth" / "parent_truth.csv", index=False)
    validation_truth.to_csv(out / "truth" / "validation_frame_truth.csv", index=False)
    catalogue.to_csv(out / "observable" / "catalogue_observed.csv", index=False)
    validation_max.to_csv(out / "observable" / "validation_max.csv", index=False)

    pz = parent["redshift"].to_numpy(float)
    cz = catalogue["redshift"].to_numpy(float)
    # Catalogue observables intentionally omit inclusion probabilities. Recover
    # the oracle values only from the truth-only parent table, preserving the
    # retained-row order used to construct catalogue_observed.
    p_catalogue = parent.loc[parent["retained"].astype(bool), "inclusion_probability"].to_numpy(float)
    if len(p_catalogue) != len(catalogue):
        raise AssertionError("truth-only oracle probabilities do not align with catalogue rows")
    parent_base = 1.0 + 1.5 * (pz - Z_MIN) / (Z_MAX - Z_MIN)
    catalogue_base = 1.0 + 1.5 * (cz - Z_MIN) / (Z_MAX - Z_MIN)

    grid = np.linspace(Z_MIN, Z_MAX, args.grid_points)
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    truth_raw = latent_galaxy_intensity(grid) * merger
    truth_density = truth_raw / float(trapezoid(truth_raw, grid))

    archive: dict[str, np.ndarray] = {
        "redshift_grid": grid,
        "merger_weight_grid": merger,
        "truth_density": truth_density,
        "parent_redshift": pz,
        "catalogue_redshift": cz,
        "catalogue_true_inclusion_probability": p_catalogue,
        "parent_base_weights": parent_base,
        "catalogue_base_weights": catalogue_base,
        "absolute_widths": np.asarray(widths),
    }
    metric_rows: list[dict[str, object]] = []
    diagnostic_rows: list[dict[str, object]] = []
    scott_rows: list[dict[str, object]] = []

    for n in sizes:
        validation = validation_max.iloc[:n].reset_index(drop=True)
        validation.to_csv(out / "observable" / f"validation_n{n}.csv", index=False)
        arms = build_five_arm_weights(
            pz,
            cz,
            validation["redshift"].to_numpy(float),
            validation["retained"].to_numpy(int),
            catalogue_floor=None,
            degree=1,
            coordinate_center=INCLUSION_COORDINATE_CENTER,
            coordinate_scale=INCLUSION_COORDINATE_SCALE,
            prior_mean=INCLUSION_PRIOR_MEAN,
            prior_sd=INCLUSION_PRIOR_SD,
            quadrature_order=args.quadrature_order,
            parent_base_weights=parent_base,
            catalogue_base_weights=catalogue_base,
            known_inclusion_probability=p_catalogue,
        )
        archive[f"posterior_nodes_n{n}"] = arms.inclusion_posterior.nodes
        archive[f"posterior_weights_n{n}"] = arms.inclusion_posterior.weights
        archive[f"posterior_mode_n{n}"] = arms.inclusion_posterior.mode

        raw_inputs = {
            "A_catalogue_only": (cz, arms.arm_a_catalogue_only),
            "B_parent_oracle": (pz, arms.arm_b_parent_oracle),
            "C_known_probability_ipw": (cz, arms.arm_c_known_probability),
            "D_map_plugin_ipw": (cz, arms.arm_d_map_plugin),
        }
        for name, (sample_z, weights) in raw_inputs.items():
            scott_rows.append({
                "validation_size": n,
                "arm": name,
                "node": -1,
                "realized_bandwidth_redshift": reflected_weighted_kde_bandwidth(
                    sample_z, lower_bound=Z_MIN, upper_bound=Z_MAX,
                    sample_weights=weights, scott_multiplier=args.scott_multiplier,
                ),
            })
            diagnostic_rows.append({
                "validation_size": n, "arm": name, "node": -1,
                **weight_concentration_diagnostics(weights),
            })

        node_diagnostics = [
            weight_concentration_diagnostics(weights)
            for weights in arms.arm_e_node_weights
        ]
        for node, diag in enumerate(node_diagnostics):
            diagnostic_rows.append({
                "validation_size": n, "arm": "E_uncertainty_integrated_ipw",
                "node": node, **diag,
            })
            scott_rows.append({
                "validation_size": n, "arm": "E_uncertainty_integrated_ipw",
                "node": node,
                "realized_bandwidth_redshift": reflected_weighted_kde_bandwidth(
                    cz, lower_bound=Z_MIN, upper_bound=Z_MAX,
                    sample_weights=arms.arm_e_node_weights[node],
                    scott_multiplier=args.scott_multiplier,
                ),
            })
        diagnostic_rows.append({
            "validation_size": n, "arm": "E_uncertainty_integrated_ipw",
            "node": -1,
            "kish_ess": min(x["kish_ess"] for x in node_diagnostics),
            "kish_ess_fraction": min(x["kish_ess_fraction"] for x in node_diagnostics),
            "entropy_ess": min(x["entropy_ess"] for x in node_diagnostics),
            "maximum_normalized_weight": max(x["maximum_normalized_weight"] for x in node_diagnostics),
            "top_1_percent_mass": max(x["top_1_percent_mass"] for x in node_diagnostics),
            "top_5_percent_mass": max(x["top_5_percent_mass"] for x in node_diagnostics),
        })

        archive[f"weights_E_nodes_n{n}"] = arms.arm_e_node_weights
        archive[f"weights_E_quadrature_n{n}"] = arms.arm_e_quadrature_weights
        for name, (_, weights) in raw_inputs.items():
            archive[f"weights_{name}_n{n}"] = weights

        for width in widths:
            density_by_arm: dict[str, np.ndarray] = {}
            for name, (sample_z, weights) in raw_inputs.items():
                density = estimate_host_density(
                    grid, sample_z, sample_weights=weights,
                    merger_weight=merger, absolute_bandwidth_sd=width,
                )
                density_by_arm[name] = density
                archive[f"density_{name}_n{n}_w{width:.6f}"] = density
                metric_rows.append({
                    "validation_size": n, "absolute_bandwidth_sd": width,
                    "arm": name, "node": -1,
                    **_metric(density, truth_density, grid), "status": "SUCCESS",
                    "kish_ess": weight_concentration_diagnostics(weights)["kish_ess"],
                })

            node_densities = np.asarray([
                estimate_host_density(
                    grid, cz, sample_weights=weights, merger_weight=merger,
                    absolute_bandwidth_sd=width,
                )
                for weights in arms.arm_e_node_weights
            ])
            archive[f"density_E_nodes_n{n}_w{width:.6f}"] = node_densities
            q = arms.arm_e_quadrature_weights
            posterior_mean_density = np.tensordot(q, node_densities, axes=(0, 0))
            archive[f"density_E_posterior_mean_n{n}_w{width:.6f}"] = posterior_mean_density
            metric_rows.append({
                "validation_size": n, "absolute_bandwidth_sd": width,
                "arm": "E_posterior_mean_density", "node": -1,
                **_metric(posterior_mean_density, truth_density, grid),
                "status": "SUCCESS", "kish_ess": float("nan"),
            })
            node_metrics = [_metric(d, truth_density, grid) for d in node_densities]
            metric_rows.append({
                "validation_size": n, "absolute_bandwidth_sd": width,
                "arm": "E_posterior_weighted_node_metric", "node": -1,
                **{
                    key: float(np.dot(q, [m[key] for m in node_metrics]))
                    for key in ("iae", "ise", "boundary_iae")
                },
                "status": "SUCCESS", "kish_ess": float("nan"),
            })
            for node, weights in enumerate(arms.arm_e_node_weights):
                metric_rows.append({
                    "validation_size": n, "absolute_bandwidth_sd": width,
                    "arm": "E_uncertainty_integrated_ipw", "node": node,
                    **_metric(node_densities[node], truth_density, grid),
                    "status": "SUCCESS",
                    "kish_ess": node_diagnostics[node]["kish_ess"],
                })

    np.savez_compressed(out / "results" / "density_reconstruction_inputs.npz", **archive)
    pd.DataFrame(metric_rows).to_csv(out / "results" / "density_metrics.csv", index=False)
    pd.DataFrame(diagnostic_rows).to_csv(out / "results" / "weight_diagnostics.csv", index=False)
    pd.DataFrame(scott_rows).to_csv(out / "results" / "scott_widths.csv", index=False)

    config = {
        "schema_version": 1,
        "experiment": "quadratic_logit_curvature_five_arm_density_smoke",
        "status": "ENGINEERING_SMOKE_NOT_SELECTION_NOT_CONFIRMATORY",
        "population_seed": int(args.seed),
        "validation_seed": int(args.validation_seed),
        "n_parent": int(args.n_parent),
        "n_catalogue": int(len(catalogue)),
        "validation_frame": int(args.validation_frame),
        "validation_sizes": sizes,
        "gamma": float(args.gamma),
        "beta": float(args.beta),
        "quadratic_logit_calibration": sample["config"],
        "analysis_model": "degree_one_logistic",
        "arm_C": "true quadratic inclusion probabilities; oracle diagnostic only",
        "arm_D_E": "degree-one logistic fitted to independent validation outcomes",
        "quadrature_order": int(args.quadrature_order),
        "coordinate_center": INCLUSION_COORDINATE_CENTER,
        "coordinate_scale": INCLUSION_COORDINATE_SCALE,
        "coefficient_prior": {"family": "independent_gaussian", "mean": list(INCLUSION_PRIOR_MEAN), "sd": list(INCLUSION_PRIOR_SD)},
        "grid_points": int(args.grid_points),
        "redshift_bounds": [Z_MIN, Z_MAX],
        "absolute_bandwidth_candidates": list(widths),
        "scott_multiplier": float(args.scott_multiplier),
        "analytic_target": "normalized latent_galaxy_intensity(z) * (1+1.5*(z-z_min)/(z_max-z_min))",
        "arm_E_metric": "posterior-mean density loss is primary; posterior expected node loss is secondary",
        "limitations": [
            "one engineering replicate only",
            "no bandwidth selection or holdout validation",
            "no H0 likelihood, posterior, PIT, or coverage",
            "Arm C uses truth-only probabilities and is an oracle diagnostic",
            "synthetic low-redshift Euclidean methodology only",
        ],
    }
    (out / "config.json").write_text(json.dumps(config, indent=2, sort_keys=True) + "\n")
    files = sorted(p for p in out.rglob("*") if p.is_file() and p.name != "manifest.json")
    manifest = {
        "status": "PASS",
        "file_count": len(files),
        "files": {
            str(p.relative_to(out)): {"bytes": p.stat().st_size, "sha256": _hash(p)}
            for p in files
        },
    }
    (out / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return {
        "status": "PASS",
        "replicate_count": 1,
        "gamma": float(args.gamma),
        "validation_sizes": sizes,
        "width_count": len(widths),
        "metric_rows": len(metric_rows),
        "weight_diagnostic_rows": len(diagnostic_rows),
        "scott_width_rows": len(scott_rows),
        "manifest_files": len(files),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261115)
    parser.add_argument("--validation-seed", type=int, default=20261116)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--validation-frame", type=int, default=1500)
    parser.add_argument("--validation-sizes", type=lambda x: tuple(int(v) for v in x.split(",")), default=(100, 300, 1000))
    parser.add_argument("--gamma", type=float, default=0.0)
    parser.add_argument("--beta", type=float, default=-3.9)
    parser.add_argument("--quadrature-order", type=int, default=5)
    parser.add_argument("--grid-points", type=int, default=401)
    parser.add_argument("--scott-multiplier", type=float, default=0.7)
    parser.add_argument("--absolute-widths", type=lambda x: tuple(float(v) for v in x.split(",")), default=DEFAULT_WIDTHS)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
