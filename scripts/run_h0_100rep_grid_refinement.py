"""Paired H0-grid refinement using archived fixed-truth pilot seeds.

Re-runs a small, predeclared subset at finer H0 grids while holding all
simulation seeds and observable/truth data fixed. This diagnoses numerical
grid sensitivity; it does not add calibration replicates.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import os
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.h0_posterior import (
    compare_h0_grid_posteriors,
    normalize_h0_posterior,
)


DATA_FILES = (
    "observable/analysis_catalogue.csv",
    "observable/gw_events_observed.csv",
    "observable/validation_observed.csv",
    "truth/analysis_parent_truth.csv",
    "truth/gw_event_truth.csv",
    "truth/gw_generation_trace_truth.csv",
    "truth/validation_frame_truth.csv",
)


def verify_manifest(root: Path) -> int:
    manifest = json.loads((root / "manifest.json").read_text())
    bad = [
        name for name, digest in manifest.items()
        if not (root / name).is_file()
        or hashlib.sha256((root / name).read_bytes()).hexdigest() != digest
    ]
    if bad:
        raise AssertionError(f"manifest integrity failures in {root}: {bad[:10]}")
    return len(manifest)


def _posterior_from_child(directory: Path):
    convergence = json.loads((directory / "quadrature_convergence.json").read_text())
    if convergence.get("status") != "CONVERGED":
        raise AssertionError(f"coefficient quadrature not converged: {directory}")
    final = convergence["posterior_by_order"][-1]
    order = int(final["order"])
    likelihood = pd.read_csv(directory / f"selected_likelihood_order_{order}.csv")
    unique = likelihood[["H0", "marginal_log_likelihood"]].drop_duplicates("H0")
    unique = unique.sort_values("H0")
    cfg = json.loads((directory / "config.json").read_text())
    prior_min, prior_max = cfg["inference"]["uniform_h0_prior"]
    posterior = normalize_h0_posterior(
        unique["H0"].to_numpy(dtype=float),
        unique["marginal_log_likelihood"].to_numpy(dtype=float),
        prior_min=float(prior_min),
        prior_max=float(prior_max),
    )
    return posterior, final, order


def compare_grid_pair(coarse, fine, *, density_tolerance=1e-3, quantile_tolerance=0.05):
    result = compare_h0_grid_posteriors(
        coarse, fine,
        density_tolerance=density_tolerance,
        quantile_tolerance=quantile_tolerance,
    )
    return {
        "max_abs_density_difference": result.max_abs_density_difference,
        "max_abs_quantile_difference": result.max_abs_quantile_difference,
        "passed_predeclared_tolerances": result.passed,
        "density_tolerance": density_tolerance,
        "quantile_tolerance": quantile_tolerance,
    }


def _assert_same_generated_data(
    baseline: Path,
    refined: Path,
    *,
    float_rtol: float = 1e-12,
    float_atol: float = 2e-12,
) -> None:
    """Require identical records, allowing only harmless floating-point roundoff.

    IDs, flags, column order, row order and non-numeric values remain exact.
    Floating columns may differ only within tight tolerances because the same
    deterministic generator can vary at machine precision across NumPy/SciPy
    builds. Material data drift still fails the audit.
    """
    for relative in DATA_FILES:
        left = pd.read_csv(baseline / relative)
        right = pd.read_csv(refined / relative)
        if left.shape != right.shape or not left.columns.equals(right.columns):
            raise AssertionError(f"generated-data shape/schema mismatch: {relative}")
        for column in left.columns:
            a, b = left[column], right[column]
            if (
                pd.api.types.is_bool_dtype(a.dtype)
                or pd.api.types.is_bool_dtype(b.dtype)
                or pd.api.types.is_integer_dtype(a.dtype)
                or pd.api.types.is_integer_dtype(b.dtype)
            ):
                if not np.array_equal(a.to_numpy(), b.to_numpy()):
                    raise AssertionError(f"generated-data exact mismatch: {relative}:{column}")
            elif pd.api.types.is_numeric_dtype(a.dtype) and pd.api.types.is_numeric_dtype(b.dtype):
                av = a.to_numpy(dtype=float)
                bv = b.to_numpy(dtype=float)
                if not np.allclose(av, bv, rtol=float_rtol, atol=float_atol, equal_nan=True):
                    error = float(np.nanmax(np.abs(av - bv)))
                    raise AssertionError(
                        f"generated-data numeric mismatch: {relative}:{column}; "
                        f"max_abs_difference={error:.6g}"
                    )
            elif not a.equals(b):
                raise AssertionError(f"generated-data exact mismatch: {relative}:{column}")


def run(args) -> dict:
    baseline = args.baseline_dir.resolve()
    output = args.output_dir.resolve()
    if output.exists() and any(output.iterdir()):
        raise ValueError("output directory must be empty")
    output.mkdir(parents=True, exist_ok=True)

    manifest_entries = verify_manifest(baseline)
    ledger = pd.read_csv(baseline / "replicate_ledger.csv")
    allowed_statuses = {"SUCCESS", "NON_EVALUABLE"}
    if ledger.empty or not set(ledger.status).issubset(allowed_statuses):
        raise AssertionError(f"baseline pilot has unsupported statuses: {sorted(set(ledger.status))}")
    planned_path = baseline / "planned_cells.csv"
    if not planned_path.is_file():
        raise AssertionError("baseline planned_cells.csv is missing")
    planned = pd.read_csv(planned_path)
    keys = ["cell_id", "replicate", "gamma", "n_validation", "analysis_seed", "validation_seed", "h0_true"]
    if planned.duplicated("cell_id").any() or ledger.duplicated("cell_id").any():
        raise AssertionError("duplicate cell IDs in planned cells or result ledger")
    left = planned[keys].sort_values(keys[:4]).reset_index(drop=True)
    right = ledger[keys].sort_values(keys[:4]).reset_index(drop=True)
    pd.testing.assert_frame_equal(left, right, check_dtype=False, check_exact=True)
    planned_subset = planned[planned.replicate < args.replicates].copy()
    if planned_subset.empty or planned_subset.replicate.nunique() != args.replicates:
        raise AssertionError("requested replicate subset is absent from baseline plan")
    selected = ledger[(ledger.replicate < args.replicates) & (ledger.status == "SUCCESS")].copy()
    if selected.empty or selected.replicate.nunique() != args.replicates:
        raise AssertionError("requested replicate subset has no successful child in one or more replicates")
    planned_per_rep = planned_subset.groupby("replicate").size()
    if not planned_per_rep.eq(planned_per_rep.iloc[0]).all():
        raise AssertionError("planned subset has unequal cell counts across replicates")
    grids = tuple(int(v) for v in args.grid_points.split(","))
    if len(grids) < 2 or tuple(sorted(set(grids))) != grids:
        raise ValueError("grid points must be distinct, increasing integers")
    if any(v < 141 for v in grids):
        raise ValueError("refinement grids must be at least 141 points")

    repo_root = Path(__file__).resolve().parents[1]
    runner = repo_root / "scripts/run_quadratic_logit_inclusion_h0_smoke.py"
    auditor = repo_root / "scripts/audit_quadratic_logit_inclusion_h0_smoke.py"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(repo_root / "src") + os.pathsep + env.get("PYTHONPATH", "")

    rows = []
    for cell in selected.itertuples(index=False):
        baseline_dir = baseline / "replicates" / cell.cell_id
        verify_manifest(baseline_dir)
        cfg = json.loads((baseline_dir / "config.json").read_text())
        analysis_cfg = cfg["analysis_generator_config"]
        validation_cfg = cfg["validation_generator_config"]
        inference_cfg = cfg["inference"]
        prior_min, prior_max = inference_cfg["uniform_h0_prior"]
        orders = ",".join(map(str, inference_cfg["quadrature_orders_requested"]))
        multipliers = ",".join(map(str, inference_cfg["prior_slope_multipliers"]))

        for points in grids:
            dest = output / "replicates" / cell.cell_id / f"grid_{points}"
            dest.mkdir(parents=True, exist_ok=False)
            cmd = [
                sys.executable, str(runner),
                "--input-dir", str(baseline_dir),
                "--seed", str(int(cell.analysis_seed)),
                "--validation-seed", str(int(cell.validation_seed)),
                "--n-parent", str(int(analysis_cfg["n_parent"])),
                "--n-detected", str(int(analysis_cfg["n_detected_events"])),
                "--h0-true", str(float(cell.h0_true)),
                "--n-validation-frame", str(int(validation_cfg["n_validation_frame"])),
                "--n-validation", str(int(cell.n_validation)),
                "--degree", str(int(inference_cfg["degree"])),
                "--orders", orders,
                "--h0-points", str(points),
                "--prior-min", str(float(prior_min)),
                "--prior-max", str(float(prior_max)),
                "--density-tolerance", str(float(cfg.get("density_tolerance", 0.001))),
                "--quantile-tolerance", str(float(cfg.get("quantile_tolerance", 0.05))),
                "--consecutive-passes", str(int(cfg.get("consecutive_passes", 2))),
                "--prior-slope-multipliers", multipliers,
                "--catalogue-model", "quadratic_logit",
                "--curvature-gamma", str(float(cell.gamma)),
                "--curvature-beta", str(float(inference_cfg["curvature_beta"])),
                "--output-dir", str(dest),
            ]
            record = {
                "cell_id": cell.cell_id,
                "replicate": int(cell.replicate),
                "gamma": float(cell.gamma),
                "n_validation": int(cell.n_validation),
                "analysis_seed": int(cell.analysis_seed),
                "validation_seed": int(cell.validation_seed),
                "grid_points": points,
                "status": "FAILED",
                "error": "",
            }
            try:
                proc = subprocess.run(cmd, cwd=repo_root, env=env, capture_output=True, text=True)
                (dest / "runner.stdout.txt").write_text(proc.stdout)
                (dest / "runner.stderr.txt").write_text(proc.stderr)
                if proc.returncode:
                    raise RuntimeError(f"runner exit {proc.returncode}: {proc.stderr[-2000:]}")
                audit = subprocess.run(
                    [sys.executable, str(auditor), str(dest), "--output", str(dest / "audit_report.json")],
                    cwd=repo_root, env=env, capture_output=True, text=True,
                )
                (dest / "auditor.stdout.txt").write_text(audit.stdout)
                (dest / "auditor.stderr.txt").write_text(audit.stderr)
                if audit.returncode:
                    raise RuntimeError(f"child audit exit {audit.returncode}: {audit.stderr[-2000:]}")
                child_report = json.loads((dest / "audit_report.json").read_text())
                if child_report.get("status") != "PASS":
                    raise AssertionError("child raw-data audit did not PASS")
                _assert_same_generated_data(baseline_dir, dest)
                verify_manifest(dest)
                posterior, final, order = _posterior_from_child(dest)
                record.update({
                    "status": "SUCCESS",
                    "audit_status": child_report["status"],
                    "quadrature_order": order,
                    "h0_median": float(final["h0_median"]),
                    "h0_q05": float(final["h0_q05"]),
                    "h0_q16": float(final["h0_q16"]),
                    "h0_q84": float(final["h0_q84"]),
                    "h0_q95": float(final["h0_q95"]),
                    "max_abs_density_difference_from_141": np.nan,
                    "max_abs_quantile_difference_from_141": np.nan,
                })
            except Exception as exc:
                record["error"] = str(exc)
            rows.append(record)
            pd.DataFrame(rows).to_csv(output / "grid_refinement_ledger.csv", index=False)

    results = pd.DataFrame(rows)
    comparisons = []
    for cell in selected.itertuples(index=False):
        base_dir = baseline / "replicates" / cell.cell_id
        base_post, base_final, base_order = _posterior_from_child(base_dir)
        by_grid = {}
        for points in grids:
            row = results[(results.cell_id == cell.cell_id) & (results.grid_points == points)]
            if len(row) != 1 or row.iloc[0].status != "SUCCESS":
                continue
            directory = output / "replicates" / cell.cell_id / f"grid_{points}"
            post, final, order = _posterior_from_child(directory)
            by_grid[points] = (post, final, order)
            comp = compare_grid_pair(base_post, post)
            comparisons.append({
                "cell_id": cell.cell_id, "replicate": int(cell.replicate),
                "gamma": float(cell.gamma), "n_validation": int(cell.n_validation),
                "from_grid_points": int(base_post.grid.size),
                "to_grid_points": points, **comp,
                "median_shift": float(final["h0_median"] - base_final["h0_median"]),
                "q16_shift": float(final["h0_q16"] - base_final["h0_q16"]),
                "q84_shift": float(final["h0_q84"] - base_final["h0_q84"]),
            })
        for grid_a, grid_b in zip(grids, grids[1:]):
            if grid_a not in by_grid or grid_b not in by_grid:
                continue
            post_a, final_a, _ = by_grid[grid_a]
            post_b, final_b, _ = by_grid[grid_b]
            comp = compare_grid_pair(post_a, post_b)
            comparisons.append({
                "cell_id": cell.cell_id, "replicate": int(cell.replicate),
                "gamma": float(cell.gamma), "n_validation": int(cell.n_validation),
                "from_grid_points": grid_a, "to_grid_points": grid_b, **comp,
                "median_shift": float(final_b["h0_median"] - final_a["h0_median"]),
                "q16_shift": float(final_b["h0_q16"] - final_a["h0_q16"]),
                "q84_shift": float(final_b["h0_q84"] - final_a["h0_q84"]),
            })
    comp_df = pd.DataFrame(comparisons)
    comp_df.to_csv(output / "grid_refinement_comparisons.csv", index=False)
    results.to_csv(output / "grid_refinement_ledger.csv", index=False)
    comparison_pass_count = (
        int(comp_df.passed_predeclared_tolerances.sum()) if len(comp_df) else 0
    )
    comparison_fail_count = int(len(comp_df) - comparison_pass_count)
    report = {
        # PASS means execution, child audits and artifact production completed.
        # Numerical grid stability is reported separately below.
        "status": "PASS" if (results.status == "SUCCESS").all() else "INCOMPLETE",
        "baseline_manifest_entries_verified": manifest_entries,
        "baseline_planned_cells_in_subset": int(len(planned_subset)),
        "baseline_successful_cells_selected": int(len(selected)),
        "baseline_non_evaluable_cells_in_subset": int((planned_subset.cell_id.isin(
            ledger.loc[ledger.status == "NON_EVALUABLE", "cell_id"]
        )).sum()),
        "planned_refined_runs": int(len(selected) * len(grids)),
        "successful_refined_runs": int((results.status == "SUCCESS").sum()),
        "failed_refined_runs": int((results.status != "SUCCESS").sum()),
        "grid_points": list(grids),
        "comparison_rows": int(len(comp_df)),
        "comparison_pass_count": comparison_pass_count,
        "comparison_fail_count": comparison_fail_count,
        "all_comparisons_within_tolerance": bool(
            len(comp_df) > 0 and comparison_fail_count == 0
        ),
        "interpretation": (
            "PASS denotes completed paired reruns and child audits, not posterior-grid "
            "convergence. Per-comparison tolerance outcomes are reported separately. "
            "This is numerical H0-grid sensitivity on fixed-truth pilot seeds, not a new "
            "calibration study, not SBC, and not observational inference."
        ),
    }
    (output / "grid_refinement_report.json").write_text(json.dumps(report, indent=2) + "\n")
    manifest = {
        str(path.relative_to(output)): hashlib.sha256(path.read_bytes()).hexdigest()
        for path in sorted(output.rglob("*"))
        if path.is_file() and path.name != "manifest.json"
    }
    (output / "manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    print(json.dumps(report, indent=2))
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--baseline-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--replicates", type=int, default=3)
    parser.add_argument("--grid-points", default="281,561")
    args = parser.parse_args()
    if args.replicates < 1:
        parser.error("--replicates must be positive")
    run(args)


if __name__ == "__main__":
    main()
