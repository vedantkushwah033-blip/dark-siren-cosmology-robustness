"""Density-level posterior quadrature stability audit for corrected curvature artifacts."""
from __future__ import annotations
import argparse
import json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.special import expit
from dark_siren_robustness.inclusion_posterior import fit_inclusion_posterior_quadrature
from audit_absolute_bandwidth_density_smoke import _independent_kde

DEFAULT_REPLICATES = (0, 5, 10, 15, 19)
DEFAULT_WIDTHS = (0.012, 0.024, 0.040)
DEFAULT_ORDERS = (5, 9, 11)

def density_total_variation(a, b, grid):
    a, b, grid = np.asarray(a, float), np.asarray(b, float), np.asarray(grid, float)
    if a.ndim != 1 or a.shape != b.shape or a.shape != grid.shape:
        raise ValueError("densities and grid must be matching one-dimensional arrays")
    if not np.all(np.isfinite(a)) or not np.all(np.isfinite(b)) or np.any(a < 0) or np.any(b < 0):
        raise ValueError("density curves must be finite and non-negative")
    if np.any(np.diff(grid) <= 0):
        raise ValueError("grid must be strictly increasing")
    return float(0.5 * trapezoid(np.abs(a - b), grid))

def _posterior_mean_density(catalogue, validation, retained, grid, rate, width, order):
    fit = fit_inclusion_posterior_quadrature(
        validation, retained, degree=1, coordinate_center=0.205,
        coordinate_scale=0.195, order=int(order),
    )
    design = np.column_stack((np.ones(catalogue.size), (catalogue - 0.205) / 0.195))
    probabilities = expit(fit.nodes @ design.T)
    if not np.all(np.isfinite(probabilities)) or np.any(probabilities <= 0) or np.any(probabilities > 1):
        raise ArithmeticError("quadrature node produced unusable inclusion probabilities")
    curves = np.asarray([
        _independent_kde(grid, catalogue, 1.0 / p, rate, float(width))
        for p in probabilities
    ])
    q = np.asarray(fit.weights, float)
    if q.shape != (curves.shape[0],) or np.any(q < 0) or not np.isclose(q.sum(), 1.0, atol=1e-10):
        raise ArithmeticError("posterior quadrature weights are invalid")
    mean = np.tensordot(q, curves, axes=(0, 0))
    area = float(trapezoid(mean, grid))
    if not np.isclose(area, 1.0, rtol=1e-9, atol=1e-9):
        raise ArithmeticError(f"posterior-mean density is not normalized: {area}")
    return mean, int(q.size), float(fit.quadrature_ess)

def audit_cell(artifact_root, output_dir, replicate_ids=DEFAULT_REPLICATES,
               widths=DEFAULT_WIDTHS, orders=DEFAULT_ORDERS, tv_tolerance=0.01):
    root = Path(artifact_root).resolve()
    cfg = json.loads((root / "config.json").read_text())
    gamma = float(cfg["gamma"])
    if cfg.get("catalogue_selection_model") != "quadratic_logit":
        raise ValueError("expected a quadratic-logit curvature artifact")
    if int(cfg.get("quadrature_order", -1)) != int(orders[0]):
        raise ValueError("first order must match production density-pilot order")
    if len(orders) < 3 or any(b <= a for a, b in zip(orders, orders[1:])):
        raise ValueError("orders must be strictly increasing and contain at least three values")
    if len(set(replicate_ids)) != len(replicate_ids) or not replicate_ids:
        raise ValueError("replicate IDs must be unique and non-empty")
    if any(r < 0 or r >= int(cfg["replicates"]) for r in replicate_ids):
        raise ValueError("replicate ID outside planned cell")
    if any(w <= 0 or not np.isfinite(w) for w in widths) or len(set(widths)) != len(widths):
        raise ValueError("bandwidths must be unique, finite and positive")
    rows, failures = [], []
    validation_sizes = tuple(int(n) for n in cfg["validation_sizes"])
    expected = len(replicate_ids) * len(validation_sizes) * len(widths)
    for rid in replicate_ids:
        child = root / "replicates" / f"replicate_{rid:04d}"
        cat_path = child / "observable" / "catalogue_observed.csv"
        archive_path = child / "results" / "density_reconstruction_inputs.npz"
        if not cat_path.is_file() or not archive_path.is_file():
            failures.append(f"replicate {rid}: missing catalogue or density archive")
            continue
        catalogue = pd.read_csv(cat_path).redshift.to_numpy(float)
        with np.load(archive_path, allow_pickle=False) as z:
            grid = np.asarray(z["redshift_grid"], float)
            rate = np.asarray(z["merger_weight_grid"], float)
            truth = np.asarray(z["truth_density"], float)
        if (catalogue.size < 5 or not np.all(np.isfinite(catalogue))
                or rate.shape != grid.shape or truth.shape != grid.shape
                or not np.isclose(trapezoid(truth, grid), 1.0, atol=1e-9)):
            failures.append(f"replicate {rid}: invalid catalogue or density archive")
            continue
        for n in validation_sizes:
            vp = child / "observable" / f"validation_n{n}.csv"
            if not vp.is_file():
                failures.append(f"replicate {rid}, n={n}: missing validation file")
                continue
            v = pd.read_csv(vp)
            validation = v.redshift.to_numpy(float)
            retained = v.retained.to_numpy(int)
            if validation.size != n or not np.isin(retained, (0, 1)).all():
                failures.append(f"replicate {rid}, n={n}: invalid validation sample")
                continue
            for width in widths:
                try:
                    curves, nodes, qess = {}, {}, {}
                    for order in orders:
                        curves[int(order)], nodes[int(order)], qess[int(order)] = _posterior_mean_density(
                            catalogue, validation, retained, grid, rate, width, order
                        )
                    first, previous, final = int(orders[0]), int(orders[-2]), int(orders[-1])
                    tv_first = density_total_variation(curves[first], curves[final], grid)
                    tv_previous = density_total_variation(curves[previous], curves[final], grid)
                    iae = {o: float(trapezoid(np.abs(curves[o] - truth), grid)) for o in curves}
                    status = "PASS" if tv_first <= tv_tolerance and tv_previous <= tv_tolerance else "NOT_PASSED"
                    rows.append({
                        "gamma": gamma, "replicate": rid, "validation_size": n,
                        "absolute_bandwidth_sd": float(width),
                        "production_order": first, "previous_refined_order": previous,
                        "refined_order": final, "production_node_count": nodes[first],
                        "previous_refined_node_count": nodes[previous], "refined_node_count": nodes[final],
                        "production_quadrature_ess": qess[first],
                        "previous_refined_quadrature_ess": qess[previous],
                        "refined_quadrature_ess": qess[final],
                        "tv_production_vs_refined": tv_first,
                        "tv_previous_refined_vs_refined": tv_previous,
                        "iae_production": iae[first], "iae_previous_refined": iae[previous],
                        "iae_refined": iae[final],
                        "iae_change_production_vs_refined": iae[first] - iae[final],
                        "status": status,
                    })
                except Exception as exc:
                    failures.append(f"replicate {rid}, n={n}, width={width}: {type(exc).__name__}: {exc}")
    frame = pd.DataFrame(rows)
    completed = len(frame)
    failed = int((frame.status != "PASS").sum()) if completed else 0
    status = "PASS" if completed == expected and not failures and failed == 0 else "NOT_PASSED"
    report = {
        "status": status, "gamma": gamma,
        "scope": "density-level Arm-E quadrature stability on a predeclared replicate subset and bandwidth grid; no H0 calibration",
        "replicate_ids_audited": list(replicate_ids),
        "validation_sizes": list(validation_sizes), "absolute_bandwidths": list(widths),
        "orders": list(orders), "planned_cases": expected, "completed_cases": completed,
        "failed_cases": failed, "tv_tolerance": tv_tolerance,
        "maximum_tv_production_vs_refined": float(frame.tv_production_vs_refined.max()) if completed else None,
        "maximum_tv_previous_refined_vs_refined": float(frame.tv_previous_refined_vs_refined.max()) if completed else None,
        "failures": failures,
    }
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    frame.to_csv(out / "density_quadrature.csv", index=False)
    (out / "audit.json").write_text(json.dumps(report, indent=2, sort_keys=True) + "\n")
    return report

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_root", type=Path)
    p.add_argument("--output-dir", type=Path, required=True)
    p.add_argument("--replicate-ids", type=lambda x: tuple(int(v) for v in x.split(",")), default=DEFAULT_REPLICATES)
    p.add_argument("--widths", type=lambda x: tuple(float(v) for v in x.split(",")), default=DEFAULT_WIDTHS)
    p.add_argument("--orders", type=lambda x: tuple(int(v) for v in x.split(",")), default=DEFAULT_ORDERS)
    p.add_argument("--tv-tolerance", type=float, default=0.01)
    a = p.parse_args()
    report = audit_cell(a.artifact_root, a.output_dir, a.replicate_ids, a.widths, a.orders, a.tv_tolerance)
    print(json.dumps(report, indent=2, sort_keys=True))
    if report["status"] != "PASS":
        raise SystemExit(1)

if __name__ == "__main__":
    main()
