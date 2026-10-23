#!/usr/bin/env python3
"""Create publication-oriented diagnostic figures from paired factorial CSV outputs.

Example:
  python scripts/plot_paired_factorial_results.py \
    --input-dir results/paired-factorial-sensitivity \
    --output-dir results/paired-factorial-sensitivity/figures
"""
from __future__ import annotations

import argparse
import ast
from pathlib import Path

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import pandas as pd


SCENARIO_LABELS = {
    "complete_uniform_uniform": "Complete\nUniform",
    "random70_uniform_uniform": "Random 70%\nUniform",
    "redshift70_uniform_uniform": "Redshift 70%\nUniform",
    "complete_luminosity_luminosity": "Complete\nLum. matched",
    "complete_luminosity_uniform": "Complete\nLum. host / uniform infer.",
    "random70_luminosity_luminosity": "Random 70%\nLum. matched",
    "random70_luminosity_uniform": "Random 70%\nLum. host / uniform infer.",
    "redshift70_luminosity_luminosity": "Redshift 70%\nLum. matched",
    "redshift70_luminosity_uniform": "Redshift 70%\nLum. host / uniform infer.",
}


def parse_interval(value: str) -> tuple[float, float]:
    parsed = ast.literal_eval(value)
    if not isinstance(parsed, (list, tuple)) or len(parsed) != 2:
        raise ValueError(f"Expected two-element interval, got {value!r}")
    lo, hi = float(parsed[0]), float(parsed[1])
    if lo > hi:
        raise ValueError(f"Reversed interval: {value!r}")
    return lo, hi


def load_inputs(input_dir: Path) -> tuple[pd.DataFrame, pd.DataFrame]:
    summary_path = input_dir / "scenario_summary.csv"
    contrasts_path = input_dir / "paired_contrasts.csv"
    if not summary_path.is_file() or not contrasts_path.is_file():
        raise FileNotFoundError(
            f"Expected {summary_path} and {contrasts_path}; run the paired experiment first."
        )
    summary = pd.read_csv(summary_path)
    contrasts = pd.read_csv(contrasts_path)
    required_summary = {"scenario", "rmse_median", "coverage_68", "coverage_95", "mean_width_68", "host_observed_fraction"}
    required_contrasts = {"label", "metric", "mean_difference_B_minus_A", "bootstrap_95_interval"}
    if missing := required_summary - set(summary.columns):
        raise ValueError(f"Summary CSV missing columns: {sorted(missing)}")
    if missing := required_contrasts - set(contrasts.columns):
        raise ValueError(f"Contrast CSV missing columns: {sorted(missing)}")
    if summary["scenario"].duplicated().any():
        raise ValueError("Duplicate scenario rows in summary CSV")
    return summary, contrasts


def scenario_metrics(summary: pd.DataFrame, output_dir: Path) -> None:
    metrics = [
        ("rmse_median", "RMSE of posterior median (km s$^{-1}$ Mpc$^{-1}$)", "factorial_rmse.png"),
        ("coverage_68", "Empirical 68% credible-interval coverage", "factorial_coverage68.png"),
        ("mean_width_68", "Mean 68% credible-interval width (km s$^{-1}$ Mpc$^{-1}$)", "factorial_width68.png"),
        ("host_observed_fraction", "Fraction of runs with host in observed catalogue", "factorial_host_visibility.png"),
    ]
    frame = summary.copy()
    frame["label"] = frame["scenario"].map(SCENARIO_LABELS).fillna(frame["scenario"])
    for column, ylabel, filename in metrics:
        fig, ax = plt.subplots(figsize=(12, 6))
        ax.bar(frame["label"], frame[column])
        ax.set_ylabel(ylabel)
        ax.set_xlabel("Prespecified simulation scenario")
        ax.tick_params(axis="x", labelrotation=25)
        for tick in ax.get_xticklabels():
            tick.set_horizontalalignment("right")
        if column == "coverage_68":
            ax.axhline(0.68, linestyle="--", linewidth=1, label="Nominal 68%")
            ax.legend()
            ax.set_ylim(0, 1.05)
        if column == "host_observed_fraction":
            ax.set_ylim(0, 1.05)
        ax.grid(axis="y", alpha=0.25)
        fig.tight_layout()
        fig.savefig(output_dir / filename, dpi=180)
        plt.close(fig)


def contrast_forest(contrasts: pd.DataFrame, output_dir: Path) -> None:
    selected = contrasts[contrasts["metric"].isin(["mean_error", "median_error"])].copy()
    if selected.empty:
        raise ValueError("No mean_error or median_error paired contrasts found")
    parsed = selected["bootstrap_95_interval"].map(parse_interval)
    selected["ci_low"] = [x[0] for x in parsed]
    selected["ci_high"] = [x[1] for x in parsed]
    selected = selected.sort_values(["metric", "label"]).reset_index(drop=True)
    fig, ax = plt.subplots(figsize=(12, max(6, 0.42 * len(selected))))
    y = list(range(len(selected)))
    estimate = selected["mean_difference_B_minus_A"].astype(float)
    ax.errorbar(
        estimate, y,
        xerr=[estimate - selected["ci_low"], selected["ci_high"] - estimate],
        fmt="o", capsize=3,
    )
    labels = [f"{row.label} [{row.metric}]" for row in selected.itertuples()]
    ax.set_yticks(y, labels)
    ax.axvline(0, linestyle="--", linewidth=1)
    ax.set_xlabel("Paired difference, B − A (km s$^{-1}$ Mpc$^{-1}$)")
    ax.set_title("Paired posterior-error contrasts with bootstrap 95% intervals")
    ax.grid(axis="x", alpha=0.25)
    fig.tight_layout()
    fig.savefig(output_dir / "factorial_paired_error_contrasts.png", dpi=180)
    plt.close(fig)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input-dir", type=Path, required=True)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    summary, contrasts = load_inputs(args.input_dir)
    args.output_dir.mkdir(parents=True, exist_ok=True)
    scenario_metrics(summary, args.output_dir)
    contrast_forest(contrasts, args.output_dir)
    print(f"Wrote diagnostic figures to {args.output_dir}")


if __name__ == "__main__":
    main()
