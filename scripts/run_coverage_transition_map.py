"""Extend fixed-truth coverage diagnostics across the prior-boundary transition.

The simulation generator and inference implementation are reused unchanged.
This is a synthetic diagnostic, not an observational H0 analysis.
"""
from __future__ import annotations
import argparse, csv, importlib.util
import numpy as np
from pathlib import Path

HERE = Path(__file__).resolve().parent
SPEC = importlib.util.spec_from_file_location("prior_support_source", HERE / "run_prior_support_sensitivity.py")
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Cannot load prior-support experiment")
prior = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(prior)

TRUTH_GRID = (52.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 88.0)
PRIOR_GRID = ((50.0, 90.0), (40.0, 100.0))

def write_csv(path, rows):
    if not rows:
        raise RuntimeError(f"Refusing to write empty output: {path}")
    path.parent.mkdir(parents=True, exist_ok=True)
    fields = list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    with path.open(newline="", encoding="utf-8") as f:
        reader = csv.DictReader(f)
        check = list(reader)
        if reader.fieldnames != fields or len(check) != len(rows):
            raise RuntimeError(f"CSV round-trip validation failed: {path}")
    return fields

def run(replicates=300, seed=20261025, grid_size=4001):
    if replicates < 100:
        raise ValueError("Use at least 100 replicates per truth")
    old_truths, old_priors = prior.TRUTHS, prior.PRIORS
    try:
        prior.TRUTHS = TRUTH_GRID
        prior.PRIORS = PRIOR_GRID
        rows, summaries = prior.run_audit(replicates=replicates, seed=seed, grid_size=grid_size)
    finally:
        prior.TRUTHS, prior.PRIORS = old_truths, old_priors
    expected_rows = len(TRUTH_GRID) * len(PRIOR_GRID) * replicates
    expected_summaries = len(TRUTH_GRID) * len(PRIOR_GRID) * len(prior.METHODS)
    if len(rows) != expected_rows or len(summaries) != expected_summaries:
        raise RuntimeError("Coverage-transition experiment is incomplete")
    return rows, summaries

def paired_method_contrasts(rows, seed=20261026, bootstrap_draws=5000):
    """Paired corrected-minus-uncorrected contrasts with percentile bootstrap CIs."""
    if bootstrap_draws < 1000:
        raise ValueError("Use at least 1000 bootstrap draws")
    rng = np.random.default_rng(seed)
    grouped = {}
    for row in rows:
        key = (row["true_H0"], row["prior_min"], row["prior_max"], row["replicate"])
        grouped.setdefault(key, row)
    cells = {}
    for row in rows:
        key = (row["true_H0"], row["prior_min"], row["prior_max"])
        cells.setdefault(key, []).append(row)
    output = []
    metrics = {
        "posterior_median": lambda r, m: r[f"{m}_median"],
        "absolute_error": lambda r, m: abs(r[f"{m}_median"] - r["true_H0"]),
        "squared_error": lambda r, m: (r[f"{m}_median"] - r["true_H0"]) ** 2,
        "width_68": lambda r, m: r[f"{m}_upper_68"] - r[f"{m}_lower_68"],
        "width_95": lambda r, m: r[f"{m}_upper_95"] - r[f"{m}_lower_95"],
        "coverage_68_indicator": lambda r, m: float(r[f"{m}_lower_68"] <= r["true_H0"] <= r[f"{m}_upper_68"]),
        "coverage_95_indicator": lambda r, m: float(r[f"{m}_lower_95"] <= r["true_H0"] <= r[f"{m}_upper_95"]),
    }
    for (truth, pmin, pmax), cell in sorted(cells.items()):
        if len(cell) < 2:
            raise RuntimeError("Expected paired corrected and uncorrected rows")
        cell = sorted(cell, key=lambda r: r["replicate"])
        for metric, fn in metrics.items():
            deltas = np.array([fn(r, "corrected") - fn(r, "uncorrected") for r in cell], dtype=float)
            if not np.all(np.isfinite(deltas)):
                raise RuntimeError(f"Non-finite paired contrast for {metric}")
            boot_idx = rng.integers(0, len(deltas), size=(bootstrap_draws, len(deltas)))
            boot_means = deltas[boot_idx].mean(axis=1)
            output.append({
                "true_H0": truth, "prior_min": pmin, "prior_max": pmax,
                "metric": metric, "n_pairs": len(deltas),
                "mean_corrected_minus_uncorrected": float(deltas.mean()),
                "median_pair_difference": float(np.median(deltas)),
                "bootstrap_95_low": float(np.quantile(boot_means, 0.025)),
                "bootstrap_95_high": float(np.quantile(boot_means, 0.975)),
                "bootstrap_draws": bootstrap_draws, "bootstrap_seed": seed,
            })
    expected = len(TRUTH_GRID) * len(PRIOR_GRID) * len(metrics)
    if len(output) != expected:
        raise RuntimeError(f"Expected {expected} paired contrast summaries, got {len(output)}")
    return output

def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=300)
    parser.add_argument("--seed", type=int, default=20261025)
    parser.add_argument("--grid-size", type=int, default=4001)
    parser.add_argument("--output-dir", default="results/coverage-transition")
    args = parser.parse_args()
    rows, summaries = run(args.replicates, args.seed, args.grid_size)
    out = Path(args.output_dir)
    write_csv(out / "per_run.csv", rows)
    write_csv(out / "summary.csv", summaries)
    print(f"PASS: {len(rows)} observation/prior rows; {len(summaries)} summaries")
    for row in summaries:
        print(row)

if __name__ == "__main__":
    main()
