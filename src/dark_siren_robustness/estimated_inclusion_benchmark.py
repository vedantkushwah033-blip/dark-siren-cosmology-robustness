"""Paired first-stage benchmark for finite-validation inclusion probabilities.

This module assesses probability estimation and IPW weight concentration only.
It does not construct H0 posteriors or claim downstream cosmological benefit.
Parent truth is used only to generate the synthetic frame and score predictions
on held-out objects; it is never passed to the fitter.
"""
from __future__ import annotations
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from dark_siren_robustness.estimated_inclusion import fit_inclusion_logistic
from dark_siren_robustness.inclusion_diagnostics import (
    inverse_probability_weight_diagnostics,
    score_inclusion_predictions,
)
from dark_siren_robustness.validation_design import split_disjoint_validation


def _draw_parent_redshifts(n_parent: int, rng: np.random.Generator):
    z_grid = np.linspace(0.01, 0.40, 20001)
    intensity = (
        1.0
        + 2.2 * np.exp(-0.5 * ((z_grid - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((z_grid - 0.275) / 0.040) ** 2)
    )
    increments = 0.5 * (intensity[:-1] + intensity[1:]) * np.diff(z_grid)
    cdf = np.concatenate(([0.0], np.cumsum(increments)))
    cdf /= cdf[-1]
    return np.interp(rng.random(n_parent), cdf, z_grid), z_grid, intensity


def _true_inclusion_probability(z, gamma, z_grid, intensity):
    x_grid = (z_grid - 0.01) / (0.40 - 0.01)
    center = float(np.trapezoid(intensity * x_grid**2, z_grid) /
                   np.trapezoid(intensity, z_grid))
    x = (z - 0.01) / (0.40 - 0.01)
    return expit(0.6 - 1.0 * x + gamma * (x**2 - center))


def run_probability_validation_benchmark(
    *,
    n_parent: int = 4000,
    validation_sizes=(100, 300, 1000),
    curvature_values=(-2.0, -1.2, 0.0, 1.2, 2.0),
    degrees=(1, 2),
    replicates: int = 20,
    seed: int = 20261108,
    raw_output: str | Path | None = None,
) -> pd.DataFrame:
    """Run paired validation-size/model-specification probability experiments."""
    if not isinstance(n_parent, int) or n_parent < 100:
        raise ValueError("n_parent must be an integer >= 100")
    if not isinstance(replicates, int) or replicates < 1:
        raise ValueError("replicates must be a positive integer")
    sizes = tuple(int(n) for n in validation_sizes)
    if not sizes or any(n < 5 or n >= n_parent for n in sizes):
        raise ValueError("validation_sizes must be between 5 and n_parent - 1")
    gammas = tuple(float(g) for g in curvature_values)
    model_degrees = tuple(int(d) for d in degrees)
    if not model_degrees or any(d not in (1, 2) for d in model_degrees):
        raise ValueError("degrees must contain only 1 or 2")

    rows = []
    raw_rows = []
    for gamma_index, gamma in enumerate(gammas):
        for replicate in range(replicates):
            scenario_seed = int(seed + 10007 * gamma_index + replicate * 1009)
            rng = np.random.default_rng(scenario_seed)
            z, z_grid, intensity = _draw_parent_redshifts(n_parent, rng)
            p_true = _true_inclusion_probability(z, gamma, z_grid, intensity)
            retained = rng.binomial(1, p_true)
            parent = pd.DataFrame({
                "object_id": np.arange(n_parent, dtype=np.int64),
                "redshift": z,
                "retained": retained,
                "inclusion_probability": p_true,
            })

            for n_validation in sizes:
                split_seed = scenario_seed + 37 * n_validation
                split = split_disjoint_validation(
                    parent, n_validation=n_validation, seed=split_seed
                )
                holdout = parent.loc[
                    ~parent["object_id"].isin(split.validation_object_ids)
                ].reset_index(drop=True)
                # Oracle weights are scoring-only diagnostics; never supplied to the fitter.
                truth_by_id = parent.set_index("object_id")
                if len(split.analysis_catalogue):
                    oracle_probability = truth_by_id.loc[
                        split.analysis_object_ids, "inclusion_probability"
                    ].to_numpy()
                    oracle_weight = inverse_probability_weight_diagnostics(oracle_probability)
                else:
                    oracle_weight = {
                        "kish_ess": np.nan, "entropy_ess": np.nan,
                        "max_normalized_weight": np.nan,
                    }

                for degree in model_degrees:
                    fit = fit_inclusion_logistic(
                        split.validation["redshift"].to_numpy(),
                        split.validation["retained"].to_numpy(),
                        degree=degree,
                    )
                    model_label = (
                        "correctly_specified"
                        if (gamma == 0.0 and degree == 1) or degree == 2
                        else "deliberately_misspecified"
                    )
                    row = {
                        "scenario_seed": scenario_seed,
                        "replicate": replicate,
                        "gamma": gamma,
                        "model_degree": degree,
                        "model_specification": model_label,
                        "n_parent": n_parent,
                        "n_validation": n_validation,
                        "validation_retained": fit.n_retained,
                        "validation_not_retained": fit.n_not_retained,
                        "n_analysis_catalogue": len(split.analysis_catalogue),
                        "oracle_kish_ess": oracle_weight["kish_ess"],
                        "oracle_entropy_ess": oracle_weight["entropy_ess"],
                        "oracle_max_normalized_weight": oracle_weight["max_normalized_weight"],
                        "fit_status": fit.status,
                        "fit_message": fit.message,
                        "design_rank": fit.design_rank,
                        "design_condition_number": fit.design_condition_number,
                        "gradient_inf_norm": fit.gradient_inf_norm,
                        "hessian_condition_number": fit.hessian_condition_number,
                        "max_abs_coefficient": fit.max_abs_coefficient,
                        "iterations": fit.iterations,
                        "fit_center": fit.center,
                        "fit_scale": fit.scale,
                        "coefficients_json": (
                            None if fit.coefficients is None
                            else json.dumps(fit.coefficients.tolist())
                        ),
                        "covariance_json": (
                            None if fit.covariance is None
                            else json.dumps(fit.covariance.tolist())
                        ),
                    }
                    if fit.status == "SUCCESS":
                        p_holdout = fit.predict(holdout["redshift"].to_numpy())
                        p_validation = fit.predict(
                            split.validation["redshift"].to_numpy()
                        )
                    else:
                        p_holdout = np.full(len(holdout), np.nan)
                        p_validation = np.full(len(split.validation), np.nan)

                    # Complete unit-level partition: validation, retained analysis
                    # catalogue, and non-retained holdout. IDs and roles are auditable.
                    raw_key = {
                        "scenario_seed": scenario_seed,
                        "replicate": replicate,
                        "gamma": gamma,
                        "n_validation": n_validation,
                        "model_degree": degree,
                    }
                    for position, unit in enumerate(split.validation.itertuples(index=False)):
                        raw_rows.append({
                            **raw_key,
                            "split_role": "validation",
                            "object_id": int(unit.object_id),
                            "redshift": float(unit.redshift),
                            "retained": int(unit.retained),
                            "p_true": float(truth_by_id.loc[unit.object_id, "inclusion_probability"]),
                            "p_hat": float(p_validation[position]),
                        })
                    for position, unit in enumerate(holdout.itertuples(index=False)):
                        raw_rows.append({
                            **raw_key,
                            "split_role": (
                                "analysis_catalogue" if int(unit.retained) == 1
                                else "nonretained_holdout"
                            ),
                            "object_id": int(unit.object_id),
                            "redshift": float(unit.redshift),
                            "retained": int(unit.retained),
                            "p_true": float(truth_by_id.loc[unit.object_id, "inclusion_probability"]),
                            "p_hat": float(p_holdout[position]),
                        })
                    if fit.status == "SUCCESS":
                        score = score_inclusion_predictions(
                            holdout["retained"].to_numpy(), p_holdout
                        )
                        if len(split.analysis_catalogue):
                            p_analysis = fit.predict(
                                split.analysis_catalogue["redshift"].to_numpy()
                            )
                            weight = inverse_probability_weight_diagnostics(p_analysis)
                        else:
                            weight = {
                                "kish_ess": np.nan, "entropy_ess": np.nan,
                                "max_normalized_weight": np.nan,
                                "top_1pct_weight_mass": np.nan,
                                "top_5pct_weight_mass": np.nan,
                            }
                        row.update({
                            "heldout_brier_score": score["brier_score"],
                            "heldout_log_loss": score["log_loss"],
                            "heldout_probability_rmse": float(np.sqrt(np.mean(
                                (p_holdout - holdout["inclusion_probability"].to_numpy())**2
                            ))),
                            "kish_ess": weight["kish_ess"],
                            "entropy_ess": weight["entropy_ess"],
                            "max_normalized_weight": weight["max_normalized_weight"],
                            "top_1pct_weight_mass": weight["top_1pct_weight_mass"],
                            "top_5pct_weight_mass": weight["top_5pct_weight_mass"],
                        })
                    else:
                        row.update({
                            "heldout_brier_score": np.nan,
                            "heldout_log_loss": np.nan,
                            "heldout_probability_rmse": np.nan,
                            "kish_ess": np.nan,
                            "entropy_ess": np.nan,
                            "max_normalized_weight": np.nan,
                            "top_1pct_weight_mass": np.nan,
                            "top_5pct_weight_mass": np.nan,
                        })
                    rows.append(row)
    summary = pd.DataFrame(rows)
    if raw_output is not None:
        raw_path = Path(raw_output)
        raw_path.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame(raw_rows).to_csv(raw_path, index=False)
    return summary
