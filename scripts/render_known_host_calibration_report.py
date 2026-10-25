"""Render a transparent Markdown report from audited calibration cell summaries.

This script is intentionally descriptive: it does not claim a universal coverage
rate by pooling different fixed truths, and it labels prior-boundary summaries
separately from interior-truth summaries.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

TRUTHS = (50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0)
NOISES = (0.05, 0.10, 0.15, 0.20)
NUMERIC = (
    "coverage_68", "coverage_68_wilson_low", "coverage_68_wilson_high",
    "coverage_95", "coverage_95_wilson_low", "coverage_95_wilson_high",
    "mean_median_error", "rmse_median", "mean_width_68", "mean_width_95",
    "median_near_prior_edge_fraction",
)


def read_summary(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def validate_summary(rows: list[dict[str, str]], metadata: dict) -> None:
    expected = {(t, n) for t in TRUTHS for n in NOISES}
    keyed = {}
    for i, row in enumerate(rows):
        truth, noise = float(row["true_h0"]), float(row["fractional_distance_error"])
        key = (truth, noise)
        if key in keyed:
            raise ValueError(f"duplicate summary cell {key}")
        keyed[key] = row
        for field in NUMERIC:
            value = float(row[field])
            if not math.isfinite(value):
                raise ValueError(f"non-finite {field} in summary cell {key}")
        n = int(row["n"])
        if n <= 0:
            raise ValueError(f"non-positive replicate count in summary cell {key}")
        for level in ("68", "95"):
            coverage = float(row[f"coverage_{level}"])
            low, high = float(row[f"coverage_{level}_wilson_low"]), float(row[f"coverage_{level}_wilson_high"])
            if not 0 <= low <= coverage <= high <= 1:
                raise ValueError(f"invalid coverage/Wilson interval for {level}% in cell {key}")
        if float(row["rmse_median"]) < 0 or float(row["mean_width_68"]) <= 0 or float(row["mean_width_95"]) <= 0:
            raise ValueError(f"invalid RMSE or interval widths in summary cell {key}")
    if set(keyed) != expected:
        raise ValueError(f"expected 36 design cells; missing={sorted(expected-set(keyed))}; extra={sorted(set(keyed)-expected)}")
    if int(metadata.get("cell_count", -1)) != 36 or int(metadata.get("raw_row_count", -1)) != sum(int(r["n"]) for r in rows):
        raise ValueError("metadata counts do not agree with summary table")
    if metadata.get("prior") != [50.0, 90.0]:
        raise ValueError("unexpected prior in metadata")


def render_report(rows: list[dict[str, str]], metadata: dict) -> str:
    validate_summary(rows, metadata)
    by_key = {(float(r["true_h0"]), float(r["fractional_distance_error"])): r for r in rows}
    lines = [
        "# Known-host calibration surface: results report",
        "",
        "> **Scope:** synthetic, single-host, matched low-redshift distance-model diagnostic. This is not an observational gravitational-wave result or a measurement of the real Hubble constant.",
        "",
        f"- Pinned source simulator commit: `{metadata.get('original_commit', 'not recorded')}`",
        f"- Base seed: `{metadata.get('seed', 'not recorded')}`",
        f"- Replicates per truth/noise cell: **{metadata.get('replicates_per_cell', 'not recorded')}**",
        f"- Prior: **[{metadata['prior'][0]}, {metadata['prior'][1]}] km s⁻¹ Mpc⁻¹**",
        f"- Raw records represented: **{metadata['raw_row_count']}** across **36** cells.",
        "",
        "## Cell-level results",
        "",
        "Coverage is the empirical fraction of replicate credible intervals containing the fixed truth. Wilson intervals quantify binomial Monte Carlo uncertainty within each cell; they do not account for model misspecification or uncertainty about whether the simulation design represents real events.",
        "",
        "| True H₀ | Fractional distance error | n | 68% coverage [Wilson 95% CI] | 95% coverage [Wilson 95% CI] | Mean median error | Median RMSE | Mean 68% width | Mean 95% width | Near-edge median fraction |",
        "|---:|---:|---:|---:|---:|---:|---:|---:|---:|---:|",
    ]
    for truth in TRUTHS:
        for noise in NOISES:
            r = by_key[(truth, noise)]
            lines.append(
                f"| {truth:.0f} | {100*noise:.0f}% | {int(r['n'])} "
                f"| {float(r['coverage_68']):.3f} [{float(r['coverage_68_wilson_low']):.3f}, {float(r['coverage_68_wilson_high']):.3f}] "
                f"| {float(r['coverage_95']):.3f} [{float(r['coverage_95_wilson_low']):.3f}, {float(r['coverage_95_wilson_high']):.3f}] "
                f"| {float(r['mean_median_error']):+.3f} | {float(r['rmse_median']):.3f} "
                f"| {float(r['mean_width_68']):.3f} | {float(r['mean_width_95']):.3f} "
                f"| {float(r['median_near_prior_edge_fraction']):.3f} |"
            )
    lines.extend([
        "",
        "## Descriptive boundary-versus-interior comparison",
        "",
        "The two boundary truths (50 and 90) are summarized separately from the seven interior truths. Values below are unweighted means across truth cells at each noise level, not pooled-event estimates. Because the truth grid and prior-boundary status differ by group, these summaries are descriptive and should not be interpreted as a causal estimate of the effect of being near a prior boundary.",
        "",
        "| Fractional distance error | Truth group | Mean 68% coverage | Mean 95% coverage | Mean absolute cell bias | Mean cell RMSE | Mean 68% width | Mean 95% width |",
        "|---:|---|---:|---:|---:|---:|---:|---:|",
    ])
    for noise in NOISES:
        for label, truths in (("Boundary (50, 90)", (50.0, 90.0)), ("Interior (55–85)", (55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0))):
            group = [by_key[(truth, noise)] for truth in truths]
            avg = lambda field: sum(float(r[field]) for r in group) / len(group)
            abs_bias = sum(abs(float(r["mean_median_error"])) for r in group) / len(group)
            lines.append(
                f"| {100*noise:.0f}% | {label} | {avg('coverage_68'):.3f} | {avg('coverage_95'):.3f} "
                f"| {abs_bias:.3f} | {avg('rmse_median'):.3f} | {avg('mean_width_68'):.3f} | {avg('mean_width_95'):.3f} |"
            )
    lines.extend([
        "",
        "## Interpretation safeguards",
        "",
        "1. Inspect each truth/noise cell and its Monte Carlo interval before discussing nominal coverage.",
        "2. Treat boundary and interior truths separately; prior truncation can alter interval behavior.",
        "3. Increasing interval width with noise is not itself evidence of correct calibration; coverage and error must be considered together.",
        "4. This is a one-host, synthetic, matched-model experiment. It does not validate catalogue incompleteness corrections, host-galaxy weighting under model mismatch, or real LVK inference.",
        "5. Any claims about a real calibration surface must be based on the raw artifact after the independent auditor reports PASS.",
        "",
    ])
    return "\n".join(lines)


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("output_dir", nargs="?", default="results/known-host-calibration-surface")
    args = parser.parse_args()
    out = Path(args.output_dir)
    rows = read_summary(out / "cell_summary.csv")
    metadata = json.loads((out / "cell_summary.json").read_text(encoding="utf-8"))
    report = render_report(rows, metadata)
    path = out / "calibration_surface_report.md"
    path.write_text(report, encoding="utf-8")
    print(f"Wrote {path}")


if __name__ == "__main__":
    main()
