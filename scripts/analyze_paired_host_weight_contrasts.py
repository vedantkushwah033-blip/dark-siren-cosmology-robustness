"""Paired, truth/noise-stratified contrasts for correct vs reversed host weights.

This analyzes the same simulated observation under both weighting assumptions.
It is a synthetic-model diagnostic, not a real GW analysis.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np

TRUTHS = (55.0, 65.0, 75.0, 85.0)
NOISES = (0.03, 0.10, 0.20)
SUPPORTS = ("narrow_50_90", "broad_20_220")
MEASURES = ("uniform_H0", "log_uniform_H0")
METRICS = ("bias_difference", "rmse_difference", "coverage68_difference",
           "coverage95_difference", "mean_width68_difference", "mean_width95_difference")


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as f:
        return list(csv.DictReader(f))


def verify_artifact(root):
    meta = json.loads((root / "run_metadata.json").read_text(encoding="utf-8"))
    audit = json.loads((root / "audit_report.json").read_text(encoding="utf-8"))
    if meta.get("audit_status") != "PASS" or audit.get("audit") != "PASS":
        raise ValueError("paired contrasts require an independently audited artifact")
    for name in ("per_run.csv", "summary.csv"):
        expected = meta.get("sha256", {}).get(name)
        actual = hashlib.sha256((root / name).read_bytes()).hexdigest()
        if not expected or expected != actual:
            raise ValueError(f"SHA-256 mismatch for {name}")
    rows = read_csv(root / "per_run.csv")
    reps = int(meta["replicates_per_truth_noise_cell"])
    if reps < 100 or len(rows) != len(TRUTHS) * len(NOISES) * reps:
        raise ValueError("raw observation count or replicate count does not match design")
    return rows, reps


def _paired_metrics(rows, indices, prefix, truth):
    c, r = prefix + "correct", prefix + "reversed"
    e0 = np.array([float(rows[i][c + "_median"]) - truth for i in indices])
    e1 = np.array([float(rows[i][r + "_median"]) - truth for i in indices])
    cov680 = np.array([int(rows[i][c + "_covered68"]) for i in indices])
    cov681 = np.array([int(rows[i][r + "_covered68"]) for i in indices])
    cov950 = np.array([int(rows[i][c + "_covered95"]) for i in indices])
    cov951 = np.array([int(rows[i][r + "_covered95"]) for i in indices])
    w680 = np.array([float(rows[i][c + "_upper68"]) - float(rows[i][c + "_lower68"]) for i in indices])
    w681 = np.array([float(rows[i][r + "_upper68"]) - float(rows[i][r + "_lower68"]) for i in indices])
    w950 = np.array([float(rows[i][c + "_upper95"]) - float(rows[i][c + "_lower95"]) for i in indices])
    w951 = np.array([float(rows[i][r + "_upper95"]) - float(rows[i][r + "_lower95"]) for i in indices])
    return np.array([
        np.mean(e1 - e0),
        np.sqrt(np.mean(e1**2)) - np.sqrt(np.mean(e0**2)),
        np.mean(cov681 - cov680),
        np.mean(cov951 - cov950),
        np.mean(w681 - w680),
        np.mean(w951 - w950),
    ], dtype=float)


def analyze(rows, seed=20261031, resamples=4000):
    if resamples < 100:
        raise ValueError("at least 100 bootstrap resamples are required")
    cells = {}
    for truth in TRUTHS:
        for noise in NOISES:
            cell = [r for r in rows if float(r["true_H0"]) == truth
                    and float(r["noise_fraction"]) == noise]
            cell.sort(key=lambda r: int(r["replicate"]))
            if not cell or len({int(r["replicate"]) for r in cell}) != len(cell):
                raise ValueError(f"missing or duplicate paired observations for {(truth, noise)}")
            cells[(truth, noise)] = cell
    rng = np.random.default_rng(seed)
    cell_results, aggregate_results = [], []
    for support in SUPPORTS:
        for measure in MEASURES:
            prefix = f"{support}_{measure}_"
            bootstrap_by_cell = {}
            point_by_cell = {}
            for (truth, noise), cell in cells.items():
                # Paired rows preserve the exact same simulated observation across arms.
                point = _paired_metrics(cell, np.arange(len(cell)), prefix, truth)
                boot = np.empty((resamples, len(METRICS)))
                for b in range(resamples):
                    ix = rng.integers(0, len(cell), len(cell))
                    boot[b] = _paired_metrics(cell, ix, prefix, truth)
                point_by_cell[(truth, noise)] = point
                bootstrap_by_cell[(truth, noise)] = boot
                row = {"true_H0": truth, "noise_fraction": noise,
                       "prior_support": support, "prior_measure": measure,
                       "replicates": len(cell)}
                for j, metric in enumerate(METRICS):
                    row[metric + "_estimate"] = float(point[j])
                    row[metric + "_ci95_low"] = float(np.quantile(boot[:, j], .025))
                    row[metric + "_ci95_high"] = float(np.quantile(boot[:, j], .975))
                cell_results.append(row)
            # Equal-cell macro-average: resample independently within every design cell.
            keys = list(cells)
            point_macro = np.mean([point_by_cell[k] for k in keys], axis=0)
            boot_macro = np.empty((resamples, len(METRICS)))
            for b in range(resamples):
                boot_macro[b] = np.mean([
                    bootstrap_by_cell[k][rng.integers(0, resamples)] for k in keys
                ], axis=0)
            for j, metric in enumerate(METRICS):
                aggregate_results.append({
                    "prior_support": support, "prior_measure": measure, "metric": metric,
                    "equal_cell_mean_difference": float(point_macro[j]),
                    "ci95_low": float(np.quantile(boot_macro[:, j], .025)),
                    "ci95_high": float(np.quantile(boot_macro[:, j], .975)),
                    "truth_noise_cells": len(keys), "replicates_per_cell": len(next(iter(cells.values()))),
                    "bootstrap_resamples": resamples,
                })
    return cell_results, aggregate_results


def write_csv(path, rows):
    if not rows:
        raise ValueError("refusing to write an empty result")
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def write_report(path, aggregate):
    lines = [
        "# Paired host-weight contrasts",
        "",
        "Differences are **reversed weights minus correct weights** on identical simulated observations.",
        "Coverage differences are proportions (multiply by 100 for percentage points).",
        "The summary is an equal-weight average of 12 truth/noise cells; intervals use a stratified paired bootstrap.",
        "These results are conditional on the synthetic two-host model, not real gravitational-wave data.",
        "",
        "| Prior support | Prior measure | Metric | Equal-cell mean difference | Paired bootstrap 95% interval |",
        "|---|---|---|---:|---:|",
    ]
    for r in aggregate:
        lines.append(f"| {r['prior_support']} | {r['prior_measure']} | {r['metric']} | "
                     f"{r['equal_cell_mean_difference']:.6g} | "
                     f"[{r['ci95_low']:.6g}, {r['ci95_high']:.6g}] |")
    path.write_text("\n".join(lines) + "\n", encoding="utf-8")


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    args = parser.parse_args()
    rows, _ = verify_artifact(args.results_dir)
    cells, aggregate = analyze(rows, args.seed, args.bootstrap_resamples)
    write_csv(args.results_dir / "paired_host_weight_cell_contrasts.csv", cells)
    write_csv(args.results_dir / "paired_host_weight_contrasts.csv", aggregate)
    write_report(args.results_dir / "paired_host_weight_contrasts.md", aggregate)
    print(f"PASS: wrote {len(cells)} paired cell contrasts and {len(aggregate)} aggregate metric rows")


if __name__ == "__main__":
    main()
