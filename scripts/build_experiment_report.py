"""Create a cautious, source-linked Markdown report from validated experiment artifacts.

This script summarizes generated outputs; it does not run simulations or claim
that the model is observationally realistic. Run only after the experiment
workflow has completed and its artifact bundle has been downloaded.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path


def load_json(path: Path) -> dict:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty JSON artifact: {path}")
    return json.loads(path.read_text(encoding="utf-8"))


def load_csv(path: Path) -> list[dict[str, str]]:
    if not path.is_file() or path.stat().st_size == 0:
        raise SystemExit(f"Missing or empty CSV artifact: {path}")
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def fmt(value: object, digits: int = 3) -> str:
    try:
        number = float(value)
    except (TypeError, ValueError):
        return str(value)
    if not math.isfinite(number):
        return "INVALID"
    return f"{number:.{digits}f}"


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("--results-root", default="results")
    parser.add_argument("--output", default="experiment_report.md")
    args = parser.parse_args()
    root = Path(args.results_root)

    factorial_root = root / "factorial-sensitivity"
    factorial = load_json(factorial_root / "scenario_summary.json")
    scenario_rows = load_csv(factorial_root / "scenario_summary.csv")
    per_run = load_csv(factorial_root / "per_run_results.csv")
    contrasts_path = factorial_root / "contrast_analysis.json"
    contrasts = load_json(contrasts_path) if contrasts_path.is_file() else None

    if len(scenario_rows) != 9 or len(factorial.get("scenarios", [])) != 9:
        raise SystemExit("Expected exactly nine main factorial scenarios.")
    n = int(factorial["requested_replicates_per_scenario"])
    if len(per_run) != 9 * n:
        raise SystemExit(f"Expected {9*n} main per-run rows; found {len(per_run)}.")
    names = {row["scenario"] for row in scenario_rows}
    if {row["scenario"] for row in per_run} != names:
        raise SystemExit("Scenario labels disagree between summary and per-run files.")

    lines = [
        "# Dark-siren synthetic experiment: results report",
        "",
        "> **Scope:** This report summarizes a synthetic simulation experiment. It is not an analysis of real gravitational-wave detections, a measurement of the real Hubble constant, or evidence that a method will work on a real galaxy survey.",
        "",
        "## Reproducibility metadata",
        "",
        f"- Original simulator commit: `{factorial.get('original_commit', 'not recorded')}`",
        f"- Base random seed: `{factorial.get('base_seed', 'not recorded')}`",
        f"- Requested replicates per scenario: **{n}**",
        f"- Main scenarios: **{len(scenario_rows)}**",
        f"- Main per-run records: **{len(per_run)}**",
        "",
        "## Main factorial results",
        "",
        "| Scenario | Completeness | Missingness | Host generation | Inference weighting | n | 68% coverage | 68% MC SE | 95% coverage | 95% MC SE | Median bias | Median RMSE | Mean 68% width | Mean 95% width |",
        "|---|---:|---|---|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for row in scenario_rows:
        lines.append(
            "| {scenario} | {completeness} | {missingness} | {host_generation} | {inference_weighting} | {n} | {c68} | {se68} | {c95} | {se95} | {bias} | {rmse} | {w68} | {w95} |".format(
                scenario=row["scenario"],
                completeness=fmt(row.get("target_completeness", "not recorded")),
                missingness=row.get("missingness", "not recorded"),
                host_generation=row.get("host_generation", "not recorded"),
                inference_weighting=row.get("inference_weighting", "not recorded"),
                n=row.get("n", row.get("requested_replicates", n)),
                c68=fmt(row.get("coverage_68")),
                se68=fmt(row.get("coverage_68_mc_se")),
                c95=fmt(row.get("coverage_95")),
                se95=fmt(row.get("coverage_95_mc_se")),
                bias=fmt(row.get("median_bias", row.get("median_error", row.get("bias_median")))),
                rmse=fmt(row.get("rmse_median")),
                w68=fmt(row.get("mean_width_68", row.get("width_68"))),
                w95=fmt(row.get("mean_width_95", row.get("width_95"))),
            )
        )

    lines.extend([
        "",
        "## Contrast analysis",
        "",
    ])
    if contrasts is None:
        lines.extend([
            "No `contrast_analysis.json` was found; contrast conclusions are omitted.",
            "",
        ])
    else:
        contrast_rows = contrasts.get("contrasts", [])
        contrast_status = "OK" if len(contrast_rows) == 56 else "CHECK"
        lines.extend([
            f"- Bootstrap repetitions: {contrasts.get('bootstrap_replicates', 'not recorded')}",
            f"- Contrasts in file: {len(contrast_rows)}",
            f"- Expected contrasts for the current design: 56; status: {contrast_status}",
            "- Each reported contrast is scenario B minus scenario A.",
            "- The intervals quantify Monte Carlo variability conditional on this synthetic design, not uncertainty in the real-world value of H0.",
            "- The contrast set is exploratory; do not treat unadjusted intervals as confirmatory tests.",
            "",
            "| Contrast | Metric | B − A | Bootstrap 95% interval |",
            "|---|---|---:|---:|",
        ])
        for row in contrast_rows:
            interval = row.get("bootstrap_95_interval", [])
            interval_text = (
                f"[{fmt(interval[0])}, {fmt(interval[1])}]"
                if isinstance(interval, list) and len(interval) == 2 else "not available"
            )
            label = str(row.get("label", "")).replace("|", "\\|")
            metric = str(row.get("metric", "")).replace("|", "\\|")
            lines.append(
                f"| {label} | {metric} | {fmt(row.get('difference_B_minus_A'))} | {interval_text} |"
            )

    lines.extend([
        "",
        "## Interpretation cautions",
        "",
        "- At 500 replicates and coverage near 0.68, the approximate binomial Monte Carlo standard error is about 0.021 (2.1 percentage points). Small differences should not be overstated.",
        "- The 70% catalogue scenarios use fixed-size masks retaining exactly 70% of the generated galaxies. This conditions on the retained count.",
        "- Redshift-dependent missingness is a synthetic stress test, not a calibrated galaxy-survey selection function.",
        "- The full-population oracle benchmark uses information unavailable in an incomplete real survey and must be described as a simulation-only reference.",
        "- The low-redshift distance approximation, generated host population, likelihood assumptions, and finite H0 grid limit generalization.",
        "- Report inconclusive or null contrasts as such. Do not select only favorable metrics or scenarios.",
        "",
        "## Required checks before external claims",
        "",
        "- Archive the full original artifact bundle with the GitHub Actions run ID and both repository commit SHAs.",
        "- Confirm the workflow succeeded and inspect failed/invalid replicate counts.",
        "- Verify the missing-host control and full-population oracle artifacts separately.",
        "- Review the assumptions and code independently before presenting scientific conclusions.",
        "",
        "## Source artifacts",
        "",
        "- `factorial-sensitivity/scenario_summary.json`",
        "- `factorial-sensitivity/scenario_summary.csv`",
        "- `factorial-sensitivity/per_run_results.csv`",
        "- `factorial-sensitivity/contrast_analysis.json` (if present)",
        "- `oracle-missing-host/` and `full-population-oracle/` control outputs (review separately)",
        "",
    ])
    output = Path(args.output)
    output.parent.mkdir(parents=True, exist_ok=True)
    output.write_text("\n".join(lines), encoding="utf-8")
    print(f"Wrote cautious report to {output}")


if __name__ == "__main__":
    main()
