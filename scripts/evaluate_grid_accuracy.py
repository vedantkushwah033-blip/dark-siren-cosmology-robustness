"""Evaluate grid-refinement CSVs against the prospective accuracy protocol.

This evaluator intentionally separates cell-level convergence failures from
individual cases that require a finer-grid follow-up. It never silently treats
an outlier as resolved. It does not assess posterior normalization because the
existing grid-audit CSV does not store the density arrays needed for an
independent integral check.
"""
from __future__ import annotations

import argparse
import csv
import json
from pathlib import Path

import numpy as np

IDENTITY = ("H0_true", "z2", "noise_fraction", "replicate", "inference_weighting")
CELL_KEYS = ("H0_true", "z2", "noise_fraction", "inference_weighting")
SUMMARY_COLUMNS = {
    "posterior_median": "posterior_median",
    "lower_68": "lower_68",
    "upper_68": "upper_68",
    "lower_95": "lower_95",
    "upper_95": "upper_95",
}
P95_TOLERANCE = 0.01
CASE_GUARDRAIL = 0.05
REFINEMENT_TOLERANCE = 0.001


def load_and_pair(path: Path):
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    if not rows:
        raise ValueError(f"No rows found in {path}")
    paired = {}
    for line, row in enumerate(rows, 2):
        try:
            key = tuple(row[k] for k in IDENTITY)
            size = int(row["grid_size"])
            numeric = {k: float(v) for k, v in row.items()
                       if k not in IDENTITY + ("grid_size", "observed_distance_mpc")}
        except (KeyError, ValueError, TypeError) as exc:
            raise ValueError(f"Malformed row {line}: {exc}") from exc
        if size not in (1001, 2001, 4001):
            raise ValueError(f"Unexpected grid_size={size} at row {line}")
        interval_keys = ("lower_95", "lower_68", "posterior_median", "upper_68", "upper_95")
        interval = [float(row[name]) for name in interval_keys]
        if not np.all(np.isfinite(interval)):
            raise ValueError(f"Non-finite posterior summary at row {line}")
        if interval != sorted(interval):
            raise ValueError(f"Unordered credible-interval summaries at row {line}")
        if key not in paired:
            paired[key] = {"metadata": row, "grids": {}}
        if size in paired[key]["grids"]:
            raise ValueError(f"Duplicate grid size {size} for case {key}")
        paired[key]["grids"][size] = row
    observations_by_simulation = {}
    for key, item in paired.items():
        if set(item["grids"]) != {1001, 2001, 4001}:
            raise ValueError(f"Incomplete resolution triplet for case {key}")
        obs = [float(item["grids"][s]["observed_distance_mpc"]) for s in (1001, 2001, 4001)]
        if not np.allclose(obs, obs[0], rtol=0.0, atol=1e-12):
            raise ValueError(f"Observation differs across grids for case {key}")
        meta = item["metadata"]
        sim_key = tuple(meta[k] for k in ("H0_true", "z2", "noise_fraction", "replicate"))
        observations_by_simulation.setdefault(sim_key, set()).add(obs[0])
    if any(len(values) != 1 for values in observations_by_simulation.values()):
        raise ValueError("Simulated observations differ across inference-weighting conditions")
    return paired


def evaluate(per_run_path: Path, targeted_path: Path | None = None):
    paired = load_and_pair(per_run_path)
    metrics = {}
    case_rows = []
    for key, item in paired.items():
        meta = item["metadata"]
        for label, column in SUMMARY_COLUMNS.items():
            vals = {size: float(item["grids"][size][column]) for size in (1001, 2001, 4001)}
            if not np.all(np.isfinite(list(vals.values()))):
                raise ValueError(f"Non-finite {label} for case {key}")
            if not vals[1001] <= vals[4001] and label.startswith("lower"):
                # A lower quantile can move either way across grid sizes; do not
                # assume monotonicity of the endpoint itself.
                pass
            diff = abs(vals[1001] - vals[4001])
            cell = tuple(meta[k] for k in CELL_KEYS)
            metrics.setdefault((cell, label), []).append(diff)
            case_rows.append({
                **{k: meta[k] for k in CELL_KEYS},
                "replicate": meta["replicate"],
                "summary": label,
                "value_1001": vals[1001],
                "value_2001": vals[2001],
                "value_4001": vals[4001],
                "abs_diff_1001_vs_4001": diff,
                "abs_diff_2001_vs_4001": abs(vals[2001] - vals[4001]),
                "flag_for_targeted_refinement": diff > CASE_GUARDRAIL,
                "targeted_refinement_status": "required" if diff > CASE_GUARDRAIL else "not_flagged",
            })
    cell_rows = []
    for (cell, summary), diffs in sorted(metrics.items(), key=lambda x: (x[0][0], x[0][1])):
        a = np.asarray(diffs, dtype=float)
        p95 = float(np.quantile(a, 0.95))
        maximum = float(np.max(a))
        cell_rows.append({
            **dict(zip(CELL_KEYS, cell)),
            "summary": summary,
            "n_cases": len(a),
            "median_abs_difference": float(np.median(a)),
            "p95_abs_difference": p95,
            "max_abs_difference": maximum,
            "p95_tolerance": P95_TOLERANCE,
            "p95_pass": p95 <= P95_TOLERANCE,
            "case_guardrail": CASE_GUARDRAIL,
            "cases_flagged_for_refinement": int(np.count_nonzero(a > CASE_GUARDRAIL)),
        })
    failed_cells = [r for r in cell_rows if not r["p95_pass"]]
    flagged_cases = [r for r in case_rows if r["flag_for_targeted_refinement"]]
    targeted_lookup = {}
    if targeted_path is not None:
        with targeted_path.open(newline="", encoding="utf-8") as handle:
            targeted_rows = list(csv.DictReader(handle))
        grouped_targeted = {}
        for row in targeted_rows:
            case_key = tuple(row[k] for k in IDENTITY)
            grouped_targeted.setdefault(case_key, {})[int(row["grid_size"])] = row
        for case_key, grids in grouped_targeted.items():
            if set(grids) != {8001, 16001, 32001, 64001}:
                raise ValueError(f"Incomplete targeted refinement ladder for {case_key}")
            for label, column in SUMMARY_COLUMNS.items():
                values = [float(grids[size][column]) for size in (8001, 16001, 32001, 64001)]
                if not np.all(np.isfinite(values)):
                    raise ValueError(f"Non-finite targeted summary {label} for {case_key}")
                targeted_lookup[(case_key, label)] = {
                    "value_8001": values[0], "value_16001": values[1],
                    "value_32001": values[2], "value_64001": values[3],
                    "abs_diff_32001_vs_64001": abs(values[2] - values[3]),
                    "refinement_pass": abs(values[2] - values[3]) <= REFINEMENT_TOLERANCE,
                }
        for case in flagged_cases:
            case_key = tuple(str(case[k]) for k in IDENTITY)
            target = targeted_lookup.get((case_key, case["summary"]))
            if target is None:
                case["targeted_refinement_status"] = "unresolved_missing_targeted_result"
            else:
                case.update(target)
                case["targeted_refinement_status"] = "resolved_within_tolerance" if target["refinement_pass"] else "unresolved_tolerance_failed"
    unresolved_flagged = [r for r in flagged_cases if r["targeted_refinement_status"] != "resolved_within_tolerance"]
    assessment = {
        "input_file": str(per_run_path),
        "matched_cases": len(paired),
        "design_cells": len({tuple(v["metadata"][k] for k in CELL_KEYS) for v in paired.values()}),
        "summary_metrics": list(SUMMARY_COLUMNS),
        "grid_comparison": "1001 versus 4001; 2001 versus 4001 also reported",
        "p95_tolerance": P95_TOLERANCE,
        "individual_case_guardrail": CASE_GUARDRAIL,
        "targeted_refinement_tolerance": REFINEMENT_TOLERANCE,
        "cell_summary_tests": len(cell_rows),
        "failed_cell_summary_tests": len(failed_cells),
        "flagged_case_summary_pairs": len(flagged_cases),
        "targeted_refinement_input": str(targeted_path) if targeted_path else None,
        "targeted_refinement_resolved_pairs": sum(r["targeted_refinement_status"] == "resolved_within_tolerance" for r in flagged_cases),
        "targeted_refinement_unresolved_pairs": len(unresolved_flagged),
        "status": "FAIL" if failed_cells or unresolved_flagged else "PASS_FOR_GRID_DIFFERENCES_ONLY",
        "limitations": [
            "The grid audit does not store posterior density arrays, so this evaluator cannot independently verify normalization.",
            "Targeted flagged cases are unresolved unless the 32001-to-64001 endpoint difference is within tolerance.",
            "The cell-level 95th-percentile criterion still applies at the original resolutions; targeted refinement does not erase a cell-level failure.",
            "Passing these numerical checks does not validate the scientific model or observational realism.",
        ],
    }
    return assessment, cell_rows, case_rows


def write_csv(path, rows):
    if not rows:
        raise ValueError(f"Refusing to write empty report: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--per-run", required=True, help="Path to per_run.csv")
    parser.add_argument("--output-dir", default="results/grid-accuracy-assessment")
    parser.add_argument("--targeted-refinement", help="Optional targeted_refinement.csv from the high-resolution follow-up")
    args = parser.parse_args()
    source = Path(args.per_run)
    output = Path(args.output_dir)
    assessment, cells, cases = evaluate(source, Path(args.targeted_refinement) if args.targeted_refinement else None)
    write_csv(output / "cell_accuracy_assessment.csv", cells)
    write_csv(output / "case_accuracy_assessment.csv", cases)
    output.mkdir(parents=True, exist_ok=True)
    (output / "assessment.json").write_text(json.dumps(assessment, indent=2) + "\n", encoding="utf-8")
    lines = [
        "# Posterior-grid accuracy assessment",
        "",
        f"**Status: {assessment['status']}**",
        "",
        f"- Matched observations: {assessment['matched_cases']}",
        f"- Design cells: {assessment['design_cells']}",
        f"- Cell-summary tests: {assessment['cell_summary_tests']}",
        f"- Failed cell-summary tests: {assessment['failed_cell_summary_tests']}",
        f"- Flagged case-summary pairs requiring targeted refinement: {assessment['flagged_case_summary_pairs']}",
        f"- Flagged pairs resolved by targeted refinement: {assessment['targeted_refinement_resolved_pairs']}",
        f"- Flagged pairs unresolved: {assessment['targeted_refinement_unresolved_pairs']}",
        f"- Prespecified 95th-percentile tolerance: {P95_TOLERANCE} km s^-1 Mpc^-1",
        f"- Individual-case guardrail: {CASE_GUARDRAIL} km s^-1 Mpc^-1",
        "",
        "Targeted refinement must compare 32,001- and 64,001-point summaries. A FAIL result is not repaired by dropping outliers or by ignoring a cell-level percentile failure.",
        "",
        "## Failed cell-summary checks",
        "",
        "| H0 true | z2 | Noise fraction | Weighting | Summary | p95 absolute difference | Maximum | Flagged cases |",
        "|---:|---:|---:|---|---|---:|---:|---:|",
    ]
    for r in cells:
        if not r["p95_pass"] or r["cases_flagged_for_refinement"]:
            lines.append(f"| {r['H0_true']} | {r['z2']} | {r['noise_fraction']} | {r['inference_weighting']} | {r['summary']} | {r['p95_abs_difference']:.8g} | {r['max_abs_difference']:.8g} | {r['cases_flagged_for_refinement']} |")
    lines += ["", "## Limitations", ""]
    lines += [f"- {x}" for x in assessment["limitations"]]
    (output / "assessment.md").write_text("\n".join(lines) + "\n", encoding="utf-8")
    print(json.dumps(assessment, indent=2))
    print(f"Reports written to {output}")


if __name__ == "__main__":
    main()
