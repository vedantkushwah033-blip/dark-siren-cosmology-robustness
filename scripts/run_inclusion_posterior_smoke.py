"""Run a small end-to-end finite-validation shared-inclusion H0 smoke experiment.

All data are synthetic. Observable inputs and simulation truth are written to
separate directories. This is a plumbing/numerical smoke runner, not a
calibration or coverage experiment.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

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
from dark_siren_robustness.inclusion_posterior import (
    fit_inclusion_posterior_quadrature,
)
from dark_siren_robustness.shared_inclusion_likelihood import (
    shared_inclusion_selected_log_likelihood,
)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def _evaluate(
    analysis: dict[str, object],
    validation: dict[str, object],
    *,
    order: int,
    h0_points: int,
) -> tuple[dict[str, object], object]:
    validation_observed = validation["validation_observed"]
    posterior = fit_inclusion_posterior_quadrature(
        validation_observed["redshift"].to_numpy(),
        validation_observed["retained"].to_numpy(),
        degree=1,
        order=order,
    )
    catalogue = analysis["catalogue_observed"]
    events_table = analysis["gw_events_observed"]
    sigma = float(analysis["config"]["sigma_mpc"])
    events = [
        {"gw_distance": float(row.distance_mpc), "gw_sigma": sigma}
        for row in events_table.itertuples(index=False)
    ]
    z_catalogue = catalogue["redshift"].to_numpy(dtype=float)
    z_min, z_max = analysis["config"]["redshift_domain"]
    base_weights = 1.0 + 1.5 * (z_catalogue - z_min) / (z_max - z_min)
    h0_grid = np.linspace(50.0, 90.0, h0_points)
    likelihood = shared_inclusion_selected_log_likelihood(
        h0_grid,
        z_catalogue,
        events,
        posterior,
        detection_threshold_mpc=float(
            analysis["config"]["detection_threshold_mpc"]
        ),
        population_base_weights=base_weights,
    )
    h0_posterior = normalize_h0_posterior(
        h0_grid, likelihood.marginal_log_likelihood,
        prior_min=50.0, prior_max=90.0,
    )
    result = {
        "order": int(order),
        "h0_points": int(h0_points),
        "posterior": posterior,
        "likelihood": likelihood,
        "h0_posterior": h0_posterior,
        "summary": {
            "quadrature_order": int(order),
            "quadrature_nodes": int(posterior.nodes.shape[0]),
            "quadrature_ess": float(posterior.quadrature_ess),
            "quadrature_maximum_normalized_weight": float(
                posterior.maximum_normalized_weight
            ),
            "validation_degree": 1,
            "validation_n": int(len(validation_observed)),
            "catalogue_n": int(len(catalogue)),
            "n_events": int(len(events)),
            "catalogue_weight_kish_ess_min": float(
                np.min(likelihood.node_kish_ess)
            ),
            "catalogue_weight_kish_ess_median": float(
                np.median(likelihood.node_kish_ess)
            ),
            "catalogue_weight_maximum_normalized_weight_max": float(
                np.max(likelihood.node_maximum_normalized_weight)
            ),
            "h0_q05": h0_posterior.quantile(0.05),
            "h0_q16": h0_posterior.quantile(0.16),
            "h0_median": h0_posterior.quantile(0.50),
            "h0_q84": h0_posterior.quantile(0.84),
            "h0_q95": h0_posterior.quantile(0.95),
            "h0_prior": [50.0, 90.0],
        },
    }
    return result, posterior


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-events", type=int, default=30)
    parser.add_argument("--n-validation-frame", type=int, default=4000)
    parser.add_argument("--n-validation", type=int, default=300)
    parser.add_argument("--output-dir", type=Path, default=Path(
        "artifacts/inclusion_posterior_smoke"
    ))
    args = parser.parse_args()
    if args.seed < 0:
        parser.error("--seed must be non-negative")
    analysis_seed = args.seed
    validation_seed = args.seed + 1
    if validation_seed == analysis_seed:
        parser.error("analysis and validation seeds must differ")

    analysis = generate_synthetic_population(
        seed=analysis_seed,
        n_parent=args.n_parent,
        n_detected=args.n_events,
        h0_true=70.0,
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
    )
    validation = generate_independent_validation_sample(
        seed=validation_seed,
        n_validation_frame=args.n_validation_frame,
        n_validation=args.n_validation,
        z_min=analysis["config"]["redshift_domain"][0],
        z_max=analysis["config"]["redshift_domain"][1],
    )

    # Three low-order fits are a smoke-level quadrature ladder only. The final
    # H0 likelihood itself is compared across these orders below.
    fits = [
        _evaluate(analysis, validation, order=order, h0_points=81)[0]
        for order in (5, 7, 9)
    ]
    fine = fits[-1]
    coarse_grid_result, _ = _evaluate(
        analysis, validation, order=9, h0_points=41
    )
    h0_grid_check = compare_h0_grid_posteriors(
        coarse_grid_result["h0_posterior"],
        fine["h0_posterior"],
        density_tolerance=0.02,
        quantile_tolerance=0.25,
    )

    # Compare likelihood shape after removing the arbitrary node-independent
    # additive constant; the evidence scale is not used as an H0 posterior.
    quadrature_comparisons = []
    reference = fine["likelihood"].marginal_log_likelihood
    for current in fits[:-1]:
        values = current["likelihood"].marginal_log_likelihood
        centered_difference = (
            values - values[0] - (reference - reference[0])
        )
        quadrature_comparisons.append({
            "order": current["order"],
            "reference_order": fine["order"],
            "max_abs_centered_log_likelihood_difference": float(
                np.max(np.abs(centered_difference))
            ),
        })

    output = args.output_dir
    observable_dir = output / "observable"
    truth_dir = output / "truth"
    analysis_dir = output / "analysis"
    for directory in (observable_dir, truth_dir, analysis_dir):
        directory.mkdir(parents=True, exist_ok=True)

    analysis["catalogue_observed"].to_csv(
        observable_dir / "catalogue_observed.csv", index=False
    )
    analysis["gw_events_observed"].to_csv(
        observable_dir / "gw_events_observed.csv", index=False
    )
    validation["validation_observed"].to_csv(
        observable_dir / "validation_observed.csv", index=False
    )
    analysis["parent_truth"].to_csv(
        truth_dir / "analysis_parent_truth.csv", index=False
    )
    analysis["gw_event_truth"].to_csv(
        truth_dir / "gw_event_truth.csv", index=False
    )
    validation["validation_frame_truth"].to_csv(
        truth_dir / "validation_frame_truth.csv", index=False
    )

    # The highest-order fit supplies the archived posterior curve.
    pd.DataFrame({
        "H0_km_s_Mpc": fine["h0_posterior"].grid,
        "posterior_density": fine["h0_posterior"].density,
        "posterior_cdf": fine["h0_posterior"].cdf,
    }).to_csv(analysis_dir / "h0_posterior.csv", index=False)

    summary = {
        "status": "SMOKE_COMPLETE_NOT_CALIBRATION",
        "scientific_scope": (
            "synthetic low-redshift Euclidean toy; no real survey, detector, "
            "LVK/GWTC, or observational H0 inference"
        ),
        "analysis_seed": analysis_seed,
        "validation_seed": validation_seed,
        "seeds_are_distinct": analysis_seed != validation_seed,
        "analysis_config": analysis["config"],
        "validation_config": validation["config"],
        "quadrature_runs": [fit["summary"] for fit in fits],
        "quadrature_h0_shape_comparisons": quadrature_comparisons,
        "h0_grid_comparison": {
            "coarse_points": 41,
            "fine_points": 81,
            "max_abs_density_difference": (
                h0_grid_check.max_abs_density_difference
            ),
            "max_abs_quantile_difference": (
                h0_grid_check.max_abs_quantile_difference
            ),
            "passed_declared_smoke_tolerances": h0_grid_check.passed,
            "tolerances_are_smoke_conventions_not_scientific_standards": True,
        },
        "selection_model": "one-sided d_obs >= threshold; common sigma",
        "population_weight_model": (
            "1 + 1.5*(z-z_min)/(z_max-z_min), passed identically to inference"
        ),
        "limitations": [
            "one synthetic replicate only",
            "no independent auditor yet",
            "no SBC or fixed-truth coverage claim",
            "quadrature orders 5/7/9 are smoke-level only",
            "no prior-scale sensitivity in this runner",
            "truth tables must never enter the fitter or likelihood",
        ],
    }
    (analysis_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    files = sorted(
        path for path in output.rglob("*")
        if path.is_file() and path.name != "sha256_manifest.json"
    )
    manifest = {
        str(path.relative_to(output)): _sha256(path)
        for path in files
    }
    (output / "sha256_manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    subprocess.run(
        [
            sys.executable,
            str(Path(__file__).with_name("audit_inclusion_posterior_smoke.py")),
            str(output),
        ],
        check=True,
    )
    print(json.dumps({
        "status": summary["status"],
        "output_dir": str(output),
        "analysis_seed": analysis_seed,
        "validation_seed": validation_seed,
        "n_catalogue": len(analysis["catalogue_observed"]),
        "n_events": len(analysis["gw_events_observed"]),
        "n_validation": len(validation["validation_observed"]),
        "h0_median_order9": fine["summary"]["h0_median"],
        "h0_grid_smoke_check_passed": h0_grid_check.passed,
        "files_manifested": len(manifest),
        "truth_is_separate": True,
    }, indent=2))


if __name__ == "__main__":
    main()
