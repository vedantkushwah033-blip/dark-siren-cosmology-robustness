"""Run a modular inclusion-posterior to selected-H0 smoke experiment.

This is a low-redshift Euclidean synthetic benchmark. Analysis population and
validation frame are independently generated with separate seeds. The fitter
sees only the validation SRS redshifts and retention labels.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.generative_selection import (
    generate_independent_validation_sample,
    generate_synthetic_population,
)
from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)
from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature
from dark_siren_robustness.inclusion_prior_sensitivity import (
    run_inclusion_prior_scale_sensitivity,
)
from dark_siren_robustness.shared_inclusion_likelihood import (
    shared_inclusion_selected_log_likelihood,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1 << 20), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--validation-seed", type=int, default=20261109)
    parser.add_argument("--n-parent", type=int, default=2400)
    parser.add_argument("--n-detected", type=int, default=8)
    parser.add_argument("--n-validation-frame", type=int, default=4000)
    parser.add_argument("--n-validation", type=int, default=300)
    parser.add_argument("--degree", type=int, choices=(1, 2), default=1)
    parser.add_argument("--orders", default="7,11,15,19")
    parser.add_argument("--h0-points", type=int, default=81)
    parser.add_argument("--prior-min", type=float, default=50.0)
    parser.add_argument("--prior-max", type=float, default=90.0)
    parser.add_argument("--density-tolerance", type=float, default=1e-3)
    parser.add_argument("--quantile-tolerance", type=float, default=0.05)
    parser.add_argument("--consecutive-passes", type=int, default=2)
    parser.add_argument("--prior-slope-multipliers", default="0.5,1.0,2.0")
    parser.add_argument("--catalogue-floor", type=float, default=0.10)
    parser.add_argument("--catalogue-midpoint", type=float, default=0.22)
    parser.add_argument("--catalogue-width", type=float, default=0.05)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/inclusion_h0_smoke"))
    args = parser.parse_args()

    orders = tuple(int(item) for item in args.orders.split(",") if item.strip())
    prior_slope_multipliers = tuple(
        float(item) for item in args.prior_slope_multipliers.split(",") if item.strip()
    )
    if (
        len(prior_slope_multipliers) < 2
        or len(set(prior_slope_multipliers)) != len(prior_slope_multipliers)
        or any(not np.isfinite(value) or value <= 0 for value in prior_slope_multipliers)
        or 1.0 not in prior_slope_multipliers
    ):
        parser.error("--prior-slope-multipliers must be distinct positive values including 1.0")
    if len(orders) < 2 or any(b <= a for a, b in zip(orders, orders[1:])):
        parser.error("--orders must contain at least two strictly increasing integers")
    if args.seed == args.validation_seed:
        parser.error("analysis and validation seeds must differ")
    if args.n_validation > args.n_validation_frame or args.n_validation < 5:
        parser.error("validation sizes must satisfy 5 <= n_validation <= n_validation_frame")
    if args.h0_points < 11 or args.prior_min >= args.prior_max:
        parser.error("invalid H0 grid or prior bounds")

    output = args.output_dir
    observable_dir = output / "observable"
    truth_dir = output / "truth"
    observable_dir.mkdir(parents=True, exist_ok=True)
    truth_dir.mkdir(parents=True, exist_ok=True)

    analysis = generate_synthetic_population(
        seed=args.seed, n_parent=args.n_parent, n_detected=args.n_detected,
        catalogue_floor=args.catalogue_floor,
        catalogue_midpoint=args.catalogue_midpoint,
        catalogue_width=args.catalogue_width,
    )
    validation = generate_independent_validation_sample(
        seed=args.validation_seed,
        n_validation_frame=args.n_validation_frame,
        n_validation=args.n_validation,
        z_min=analysis["config"]["redshift_domain"][0],
        z_max=analysis["config"]["redshift_domain"][1],
        catalogue_floor=args.catalogue_floor,
        catalogue_midpoint=args.catalogue_midpoint,
        catalogue_width=args.catalogue_width,
    )
    analysis_catalogue = analysis["catalogue_observed"]
    events_frame = analysis["gw_events_observed"]
    validation_frame = validation["validation_observed"]

    # Inference input tables are deliberately stripped of all truth-only columns.
    analysis_catalogue.to_csv(observable_dir / "analysis_catalogue.csv", index=False)
    events_frame.to_csv(observable_dir / "gw_events_observed.csv", index=False)
    validation_frame.to_csv(observable_dir / "validation_observed.csv", index=False)
    analysis["parent_truth"].to_csv(truth_dir / "analysis_parent_truth.csv", index=False)
    analysis["gw_event_truth"].to_csv(truth_dir / "gw_event_truth.csv", index=False)
    analysis["gw_generation_trace_truth"].to_csv(
        truth_dir / "gw_generation_trace_truth.csv", index=False
    )
    validation["validation_frame_truth"].to_csv(
        truth_dir / "validation_frame_truth.csv", index=False
    )

    z_validation = validation_frame["redshift"].to_numpy(dtype=float)
    retained_validation = validation_frame["retained"].to_numpy(dtype=int)
    z_catalogue = analysis_catalogue["redshift"].to_numpy(dtype=float)
    z_min, z_max = analysis["config"]["redshift_domain"]
    base_weights = 1.0 + 1.5 * (z_catalogue - z_min) / (z_max - z_min)
    events = [
        {
            "gw_distance": float(row.distance_mpc),
            "gw_sigma": float(analysis["config"]["sigma_mpc"]),
        }
        for row in events_frame.itertuples(index=False)
    ]
    h0_grid = np.linspace(args.prior_min, args.prior_max, args.h0_points)
    # Record the coefficient prior explicitly so the posterior target is fully
    # specified in the artifact rather than relying on library defaults.
    coefficient_prior_mean = np.zeros(args.degree + 1, dtype=float)
    coefficient_prior_sd = np.asarray(
        [5.0] + [2.5] * args.degree, dtype=float
    )

    posterior_records = []
    density_columns = {"H0": h0_grid.tolist()}
    convergence_records = []
    previous_posterior = None
    consecutive = 0
    converged_at = None
    last_fit = None

    for order in orders:
        fit = fit_inclusion_posterior_quadrature(
            z_validation, retained_validation,
            degree=args.degree, prior_mean=coefficient_prior_mean,
            prior_sd=coefficient_prior_sd, order=order,
        )
        likelihood = shared_inclusion_selected_log_likelihood(
            h0_grid, z_catalogue, events, fit,
            detection_threshold_mpc=float(analysis["config"]["detection_threshold_mpc"]),
            population_base_weights=base_weights,
        )
        coefficient_columns = {
            f"theta_{index}": fit.nodes[:, index]
            for index in range(fit.nodes.shape[1])
        }
        coefficient_columns["posterior_weight"] = fit.weights
        pd.DataFrame(coefficient_columns).to_csv(
            output / f"inclusion_quadrature_order_{order}.csv", index=False
        )
        node_count = fit.nodes.shape[0]
        h0_count = h0_grid.size
        likelihood_record = pd.DataFrame({
            "node_index": np.repeat(np.arange(node_count), h0_count),
            "H0": np.tile(h0_grid, node_count),
            "node_log_selected_likelihood": likelihood.node_log_selected_likelihood.reshape(-1),
            "node_log_detection_efficiency": likelihood.node_log_detection_efficiency.reshape(-1),
            "marginal_log_likelihood": np.tile(likelihood.marginal_log_likelihood, node_count),
        })
        likelihood_record.to_csv(
            output / f"selected_likelihood_order_{order}.csv", index=False
        )
        h0_posterior = normalize_h0_posterior(
            h0_grid, likelihood.marginal_log_likelihood,
            prior_min=args.prior_min, prior_max=args.prior_max,
        )
        density_columns[f"density_order_{order}"] = np.interp(
            h0_grid, h0_posterior.grid, h0_posterior.density
        ).tolist()
        posterior_records.append({
            "order": order,
            "validation_center": fit.center,
            "validation_scale": fit.scale,
            "posterior_mean_coefficients": fit.posterior_mean.tolist(),
            "posterior_covariance": fit.posterior_covariance.tolist(),
            "quadrature_ess": fit.quadrature_ess,
            "maximum_node_mass": fit.maximum_normalized_weight,
            "minimum_predicted_inclusion_probability": fit.min_predicted_probability,
            "maximum_predicted_inclusion_probability": fit.max_predicted_probability,
            "h0_q05": h0_posterior.quantile(0.05),
            "h0_q16": h0_posterior.quantile(0.16),
            "h0_median": h0_posterior.quantile(0.50),
            "h0_q84": h0_posterior.quantile(0.84),
            "h0_q95": h0_posterior.quantile(0.95),
            "h0_log_likelihood_integral": h0_posterior.log_likelihood_integral,
        })
        if previous_posterior is not None:
            comparison = compare_h0_grid_posteriors(
                previous_posterior, h0_posterior,
                density_tolerance=args.density_tolerance,
                quantile_tolerance=args.quantile_tolerance,
            )
            consecutive = consecutive + 1 if comparison.passed else 0
            convergence_records.append({
                "from_order": last_fit.order,
                "to_order": order,
                "max_abs_density_difference": comparison.max_abs_density_difference,
                "max_abs_quantile_difference": comparison.max_abs_quantile_difference,
                "passed": comparison.passed,
                "consecutive_passes": consecutive,
            })
            if consecutive >= args.consecutive_passes:
                converged_at = order
        previous_posterior = h0_posterior
        last_fit = fit
        if converged_at is not None:
            break

    if last_fit is None:
        raise RuntimeError("no inclusion posterior fit was completed")
    prior_sensitivity = run_inclusion_prior_scale_sensitivity(
        z_validation,
        retained_validation,
        h0_grid,
        z_catalogue,
        events,
        detection_threshold_mpc=float(analysis["config"]["detection_threshold_mpc"]),
        population_base_weights=base_weights,
        degree=args.degree,
        prior_mean=coefficient_prior_mean,
        prior_sd=coefficient_prior_sd,
        slope_multipliers=prior_slope_multipliers,
        order=last_fit.order,
        prior_min=args.prior_min,
        prior_max=args.prior_max,
    )
    sensitivity_density = {"H0": h0_grid.tolist()}
    sensitivity_records = []
    for record, density_difference, quantile_shift in zip(
        prior_sensitivity.records,
        prior_sensitivity.max_abs_density_difference_from_reference,
        prior_sensitivity.quantile_shift_from_reference,
    ):
        sensitivity_density[
            f"density_slope_multiplier_{record.slope_multiplier:g}"
        ] = np.interp(
            h0_grid,
            record.h0_posterior.grid,
            record.h0_posterior.density,
        ).tolist()
        sensitivity_records.append({
            "slope_multiplier": record.slope_multiplier,
            "prior_mean": record.prior_mean.tolist(),
            "prior_sd": record.prior_sd.tolist(),
            "posterior_mean_coefficients": record.inclusion_posterior.posterior_mean.tolist(),
            "posterior_covariance": record.inclusion_posterior.posterior_covariance.tolist(),
            "quadrature_ess": record.inclusion_posterior.quadrature_ess,
            "minimum_catalogue_kish_ess": record.minimum_catalogue_kish_ess,
            "median_catalogue_kish_ess": record.median_catalogue_kish_ess,
            "maximum_catalogue_normalized_weight": record.maximum_catalogue_normalized_weight,
            "h0_quantiles": record.h0_quantiles,
            "max_abs_density_difference_from_reference": density_difference,
            "quantile_shift_from_reference": quantile_shift,
        })
    pd.DataFrame(sensitivity_density).to_csv(
        output / "prior_scale_sensitivity_density.csv", index=False
    )
    prior_sensitivity_summary = {
        "reference_multiplier": prior_sensitivity.reference_multiplier,
        "quadrature_order": last_fit.order,
        "slope_multipliers": list(prior_slope_multipliers),
        "records": sensitivity_records,
        "interpretation": (
            "Prior sensitivity only: all observable data, catalogue, events, "
            "selection rule, H0 grid and base population weights are held fixed."
        ),
        "warning": "Sensitivity is not evidence of calibration or coverage.",
    }
    (output / "prior_scale_sensitivity.json").write_text(
        json.dumps(prior_sensitivity_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    status = "CONVERGED" if converged_at is not None else "NON_CONVERGED"
    pd.DataFrame(density_columns).to_csv(output / "h0_posterior_density_by_order.csv", index=False)
    posterior_summary = {
        "status": status,
        "converged_at_order": converged_at,
        "orders_evaluated": [item["order"] for item in posterior_records],
        "required_consecutive_passes": args.consecutive_passes,
        "density_tolerance": args.density_tolerance,
        "quantile_tolerance": args.quantile_tolerance,
        "posterior_by_order": posterior_records,
        "h0_posterior_convergence": convergence_records,
        "warning": "Numerical order stability is not SBC, calibration, coverage, or scientific validation.",
    }
    (output / "quadrature_convergence.json").write_text(
        json.dumps(posterior_summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )

    config = {
        "experiment": "independent_validation_shared_inclusion_selected_h0_smoke",
        "status": "synthetic_design_not_real_data_or_cosmological_measurement",
        "analysis_generator_config": analysis["config"],
        "validation_generator_config": validation["config"],
        "inference": {
            "degree": args.degree,
            "coefficient_prior": {
                "family": "independent_normal",
                "mean": coefficient_prior_mean.tolist(),
                "sd": coefficient_prior_sd.tolist(),
            },
            "prior_slope_multipliers": list(prior_slope_multipliers),
            "quadrature_orders_requested": list(orders),
            "h0_grid_points": args.h0_points,
            "uniform_h0_prior": [args.prior_min, args.prior_max],
            "selection_rule": "d_obs >= detection_threshold_mpc",
            "shared_inclusion_coefficients_across_events": True,
            "population_base_weight": "1 + 1.5*(z-z_min)/(z_max-z_min)",
        },
        "quadrature_convergence_status": status,
        "truth_firewall": "Only validation redshift/retention, observable catalogue redshift, and observed event distances enter inference.",
        "scientific_limit": "Low-redshift Euclidean synthetic toy; no real survey, LVK data, H0 measurement, SBC, or coverage claim.",
    }
    (output / "config.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )

    # Save checksums for every artifact present before writing the manifest itself.
    manifest = {
        str(path.relative_to(output)): _sha256(path)
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "manifest.json"
    }
    (output / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": status,
        "converged_at_order": converged_at,
        "orders_evaluated": posterior_summary["orders_evaluated"],
        "n_analysis_catalogue": len(analysis_catalogue),
        "n_validation": len(validation_frame),
        "n_detected_events": len(events_frame),
        "output_dir": str(output),
        "manifest_files": len(manifest),
    }, indent=2))


if __name__ == "__main__":
    main()
