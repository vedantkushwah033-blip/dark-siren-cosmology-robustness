"""Full-support inclusion-posterior quadrature and host-weight stability audit.

This is a numerical gate for the density-only curvature pilot. It does not
validate an H0 likelihood or posterior. It re-fits every persisted validation
sample at a predeclared Gauss-Hermite order ladder and measures downstream
catalogue-weight stability without changing the raw pilot artifacts.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.special import expit

from dark_siren_robustness.inclusion_posterior import (
    build_node_specific_host_weights,
    fit_inclusion_posterior_quadrature_ladder,
    posterior_predictive_inclusion,
)

Z_MIN, Z_MAX = 0.01, 0.40


def _mixture_host_weights(fit, catalogue_redshift: np.ndarray) -> np.ndarray:
    host = build_node_specific_host_weights(fit, catalogue_redshift)
    mixture = fit.weights @ host.normalized_weight
    if (mixture.ndim != 1 or not np.all(np.isfinite(mixture))
            or np.any(mixture < 0) or not np.isclose(mixture.sum(), 1.0)):
        raise ArithmeticError("posterior-mixture host weights are invalid")
    return mixture


def _weight_summary(weights: np.ndarray) -> dict[str, float]:
    w = np.asarray(weights, dtype=float)
    if w.ndim != 1 or w.size == 0 or not np.all(np.isfinite(w)) or np.any(w < 0):
        raise ValueError("mixture weights must be a finite non-negative vector")
    total = float(w.sum())
    if total <= 0:
        raise ValueError("mixture weights must have positive total mass")
    p = w / total
    positive = p[p > 0]
    ordered = np.sort(p)[::-1]
    n1 = max(1, int(np.ceil(0.01 * p.size)))
    n5 = max(1, int(np.ceil(0.05 * p.size)))
    return {
        "kish_ess": float(1.0 / np.sum(p * p)),
        "entropy_ess": float(np.exp(-np.sum(positive * np.log(positive)))),
        "maximum_normalized_weight": float(ordered[0]),
        "top_1_percent_mass": float(ordered[:n1].sum()),
        "top_5_percent_mass": float(ordered[:n5].sum()),
    }


def _outside_validation_mass(
    fit, catalogue_redshift: np.ndarray, validation_redshift: np.ndarray
) -> float:
    host = build_node_specific_host_weights(fit, catalogue_redshift)
    outside = ((catalogue_redshift < np.min(validation_redshift))
               | (catalogue_redshift > np.max(validation_redshift)))
    node_mass = host.normalized_weight[:, outside].sum(axis=1)
    return float(fit.weights @ node_mass)


def audit(
    artifact_root: Path,
    *,
    orders: tuple[int, ...] = (5, 7, 9, 11, 13),
    full_support_points: int = 401,
    mean_tolerance: float = 1e-3,
    covariance_tolerance: float = 1e-3,
    predictive_tolerance: float = 1e-3,
    log_evidence_tolerance: float = 1e-3,
    mixture_max_abs_tolerance: float = 1e-3,
    mixture_total_variation_tolerance: float = 5e-3,
    consecutive_passes: int = 2,
) -> dict:
    root = artifact_root.resolve()
    cfg = json.loads((root / "config.json").read_text())
    if cfg.get("catalogue_selection_model") != "quadratic_logit":
        raise ValueError("quadrature gate expects a quadratic-logit curvature artifact")
    expected_replicates = int(cfg["replicates"])
    validation_sizes = [int(x) for x in cfg["validation_sizes"]]
    full_grid = np.linspace(Z_MIN, Z_MAX, full_support_points)
    rows: list[dict] = []
    failures: list[str] = []

    for replicate_id in range(expected_replicates):
        child = root / "replicates" / f"replicate_{replicate_id:04d}"
        if not child.is_dir():
            failures.append(f"replicate {replicate_id}: missing attempt directory")
            continue
        child_cfg_path = child / "config.json"
        if not child_cfg_path.is_file():
            failures.append(f"replicate {replicate_id}: missing child config")
            continue
        child_cfg = json.loads(child_cfg_path.read_text())
        catalogue_path = child / "observable" / "catalogue_observed.csv"
        if not catalogue_path.is_file():
            failures.append(f"replicate {replicate_id}: missing observed catalogue")
            continue
        catalogue = pd.read_csv(catalogue_path).redshift.to_numpy(float)
        if catalogue.size == 0 or not np.all(np.isfinite(catalogue)):
            failures.append(f"replicate {replicate_id}: invalid observed catalogue")
            continue

        for n in validation_sizes:
            validation_path = child / "observable" / f"validation_n{n}.csv"
            if not validation_path.is_file():
                failures.append(f"replicate {replicate_id}, n={n}: missing validation sample")
                continue
            validation = pd.read_csv(validation_path)
            vz = validation.redshift.to_numpy(float)
            retained = validation.retained.to_numpy(int)
            if vz.size != n:
                failures.append(f"replicate {replicate_id}, n={n}: validation row count mismatch")
                continue

            try:
                ladder = fit_inclusion_posterior_quadrature_ladder(
                    vz,
                    retained,
                    degree=1,
                    coordinate_center=0.205,
                    coordinate_scale=0.195,
                    orders=orders,
                    predictive_redshift=full_grid,
                    mean_tolerance=mean_tolerance,
                    covariance_tolerance=covariance_tolerance,
                    predictive_tolerance=predictive_tolerance,
                    log_evidence_tolerance=log_evidence_tolerance,
                    consecutive_passes=consecutive_passes,
                )
                mixtures = []
                for fit in ladder.fits:
                    mixture = _mixture_host_weights(fit, catalogue)
                    mixtures.append(mixture)
                    outside_mass = _outside_validation_mass(fit, catalogue, vz)
                    x_full = (full_grid - fit.center) / fit.scale
                    full_design = np.column_stack((np.ones(full_grid.size), x_full))
                    node_probability = expit(fit.nodes @ full_design.T)
                    posterior_mean_probability = posterior_predictive_inclusion(fit, full_grid)
                    rows.append({
                        "replicate": replicate_id,
                        "validation_size": n,
                        "order": fit.order,
                        "node_count": int(fit.weights.size),
                        "quadrature_ess": fit.quadrature_ess,
                        "maximum_quadrature_node_weight": fit.maximum_normalized_weight,
                        "posterior_mean_intercept": float(fit.posterior_mean[0]),
                        "posterior_mean_slope": float(fit.posterior_mean[1]),
                        "posterior_sd_intercept": float(np.sqrt(fit.posterior_covariance[0, 0])),
                        "posterior_sd_slope": float(np.sqrt(fit.posterior_covariance[1, 1])),
                        "log_evidence_unnormalized": fit.log_evidence_unnormalized,
                        "full_support_min_node_inclusion": float(np.min(node_probability)),
                        "full_support_max_node_inclusion": float(np.max(node_probability)),
                        "full_support_min_posterior_mean_inclusion": float(np.min(posterior_mean_probability)),
                        "full_support_max_posterior_mean_inclusion": float(np.max(posterior_mean_probability)),
                        "catalogue_fraction_outside_validation_range": float(
                            np.mean((catalogue < vz.min()) | (catalogue > vz.max()))
                        ),
                        "posterior_mean_outside_validation_ipw_mass": outside_mass,
                        **_weight_summary(mixture),
                        "ladder_status": ladder.status,
                    })

                for j in range(1, len(mixtures)):
                    delta = np.abs(mixtures[j] - mixtures[j - 1])
                    max_abs = float(np.max(delta))
                    total_variation = float(0.5 * np.sum(delta))
                    rows[j + len(rows) - len(mixtures)].update({
                        "mixture_max_abs_change_from_previous_order": max_abs,
                        "mixture_total_variation_from_previous_order": total_variation,
                    })
                    if (max_abs > mixture_max_abs_tolerance
                            or total_variation > mixture_total_variation_tolerance):
                        failures.append(
                            f"replicate {replicate_id}, n={n}: host-weight mixture "
                            f"unstable between orders {ladder.fits[j-1].order} and "
                            f"{ladder.fits[j].order} (max_abs={max_abs:.6g}, "
                            f"TV={total_variation:.6g})"
                        )
                if ladder.status != "CONVERGED":
                    failures.append(
                        f"replicate {replicate_id}, n={n}: quadrature ladder "
                        f"did not converge by order {orders[-1]}"
                    )
            except Exception as exc:
                failures.append(
                    f"replicate {replicate_id}, n={n}: {type(exc).__name__}: {exc}"
                )

    expected_cases = expected_replicates * len(validation_sizes)
    completed_cases = sum(
        1 for rid in range(expected_replicates) for n in validation_sizes
        if any(row["replicate"] == rid and row["validation_size"] == n for row in rows)
    )
    status = "PASS" if completed_cases == expected_cases and not failures else "NOT_PASSED"
    return {
        "status": status,
        "scope": "full-support inclusion quadrature and posterior-mixture host-weight stability; no H0 validation",
        "artifact_root": str(root),
        "planned_replicates": expected_replicates,
        "validation_sizes": validation_sizes,
        "planned_replicate_validation_cases": expected_cases,
        "completed_replicate_validation_cases": completed_cases,
        "orders": list(orders),
        "full_support": [Z_MIN, Z_MAX],
        "full_support_grid_points": full_support_points,
        "tolerances": {
            "posterior_mean": mean_tolerance,
            "posterior_covariance": covariance_tolerance,
            "full_support_predictive_inclusion": predictive_tolerance,
            "log_evidence": log_evidence_tolerance,
            "host_mixture_max_abs": mixture_max_abs_tolerance,
            "host_mixture_total_variation": mixture_total_variation_tolerance,
            "consecutive_passes": consecutive_passes,
        },
        "failures": failures,
        "rows": rows,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("artifact_root", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--orders", type=lambda x: tuple(int(v) for v in x.split(",")),
                        default=(5, 7, 9, 11, 13))
    parser.add_argument("--full-support-points", type=int, default=401)
    parser.add_argument("--mean-tolerance", type=float, default=1e-3)
    parser.add_argument("--covariance-tolerance", type=float, default=1e-3)
    parser.add_argument("--predictive-tolerance", type=float, default=1e-3)
    parser.add_argument("--log-evidence-tolerance", type=float, default=1e-3)
    parser.add_argument("--mixture-max-abs-tolerance", type=float, default=1e-3)
    parser.add_argument("--mixture-total-variation-tolerance", type=float, default=5e-3)
    parser.add_argument("--consecutive-passes", type=int, default=2)
    parser.add_argument("--require-convergence", action="store_true")
    args = parser.parse_args()
    report = audit(
        args.artifact_root,
        orders=args.orders,
        full_support_points=args.full_support_points,
        mean_tolerance=args.mean_tolerance,
        covariance_tolerance=args.covariance_tolerance,
        predictive_tolerance=args.predictive_tolerance,
        log_evidence_tolerance=args.log_evidence_tolerance,
        mixture_max_abs_tolerance=args.mixture_max_abs_tolerance,
        mixture_total_variation_tolerance=args.mixture_total_variation_tolerance,
        consecutive_passes=args.consecutive_passes,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    if report["rows"]:
        pd.DataFrame(report["rows"]).to_csv(args.output.with_suffix(".csv"), index=False)
    print(json.dumps({k: v for k, v in report.items() if k != "rows"}, indent=2, sort_keys=True))
    if args.require_convergence and report["status"] != "PASS":
        raise SystemExit(2)


if __name__ == "__main__":
    main()
