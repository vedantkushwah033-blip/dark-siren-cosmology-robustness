"""Independent reconstruction audit for paired host-weight contrast artifacts.

This deliberately does not import the primary contrast-analysis script. It
recomputes paired point estimates and stratified bootstrap intervals directly
from audited per-run posterior summaries, then compares every output row.
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
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))

def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()

def _arm_arrays(rows, prefix, truth):
    e_correct = np.asarray([float(x[f"{prefix}correct_median"]) - truth for x in rows])
    e_reversed = np.asarray([float(x[f"{prefix}reversed_median"]) - truth for x in rows])
    c68_correct = np.asarray([int(x[f"{prefix}correct_covered68"]) for x in rows], dtype=float)
    c68_reversed = np.asarray([int(x[f"{prefix}reversed_covered68"]) for x in rows], dtype=float)
    c95_correct = np.asarray([int(x[f"{prefix}correct_covered95"]) for x in rows], dtype=float)
    c95_reversed = np.asarray([int(x[f"{prefix}reversed_covered95"]) for x in rows], dtype=float)
    w68_correct = np.asarray([float(x[f"{prefix}correct_upper68"]) - float(x[f"{prefix}correct_lower68"]) for x in rows])
    w68_reversed = np.asarray([float(x[f"{prefix}reversed_upper68"]) - float(x[f"{prefix}reversed_lower68"]) for x in rows])
    w95_correct = np.asarray([float(x[f"{prefix}correct_upper95"]) - float(x[f"{prefix}correct_lower95"]) for x in rows])
    w95_reversed = np.asarray([float(x[f"{prefix}reversed_upper95"]) - float(x[f"{prefix}reversed_lower95"]) for x in rows])
    return (e_correct, e_reversed, c68_correct, c68_reversed, c95_correct, c95_reversed,
            w68_correct, w68_reversed, w95_correct, w95_reversed)

def _contrasts(arrays, indices):
    ec, er, c68c, c68r, c95c, c95r, w68c, w68r, w95c, w95r = arrays
    ec, er = ec[indices], er[indices]
    c68c, c68r, c95c, c95r = c68c[indices], c68r[indices], c95c[indices], c95r[indices]
    w68c, w68r, w95c, w95r = w68c[indices], w68r[indices], w95c[indices], w95r[indices]
    return np.stack((
        np.mean(er - ec, axis=-1),
        np.sqrt(np.mean(er**2, axis=-1)) - np.sqrt(np.mean(ec**2, axis=-1)),
        np.mean(c68r - c68c, axis=-1),
        np.mean(c95r - c95c, axis=-1),
        np.mean(w68r - w68c, axis=-1),
        np.mean(w95r - w95c, axis=-1),
    ), axis=-1)

def audit(root: Path, seed=20261031, resamples=4000, atol=2e-10):
    meta = json.loads((root / "run_metadata.json").read_text(encoding="utf-8"))
    report = json.loads((root / "audit_report.json").read_text(encoding="utf-8"))
    if meta.get("audit_status") != "PASS" or report.get("audit") != "PASS":
        raise ValueError("paired contrast audit requires a posterior-audited artifact")
    for name in ("per_run.csv", "summary.csv"):
        expected_hash = meta.get("sha256", {}).get(name)
        if not expected_hash or expected_hash != sha256(root / name):
            raise ValueError(f"input SHA-256 mismatch: {name}")
    raw = read_csv(root / "per_run.csv")
    reps = int(meta["replicates_per_truth_noise_cell"])
    expected_cells = {(t, n) for t in TRUTHS for n in NOISES}
    grouped = {}
    for row in raw:
        key = (float(row["true_H0"]), float(row["noise_fraction"]))
        if key not in expected_cells:
            raise ValueError(f"unexpected truth/noise cell: {key}")
        grouped.setdefault(key, []).append(row)
    if set(grouped) != expected_cells or len(raw) != len(expected_cells) * reps:
        raise ValueError("raw data do not match the declared 12-cell design")
    for key, block in grouped.items():
        block.sort(key=lambda x: int(x["replicate"]))
        ids = [int(x["replicate"]) for x in block]
        if len(block) != reps or ids != list(range(reps)):
            raise ValueError(f"missing, duplicated, or unordered replicate IDs in {key}")
    rng = np.random.default_rng(seed)
    expected_cell_rows = []
    expected_aggregate_rows = []
    for support in SUPPORTS:
        for measure in MEASURES:
            prefix = f"{support}_{measure}_"
            cell_points = {}
            cell_boots = {}
            for truth in TRUTHS:
                for noise in NOISES:
                    block = grouped[(truth, noise)]
                    arrays = _arm_arrays(block, prefix, truth)
                    point = _contrasts(arrays, np.arange(reps))[0]
                    indices = rng.integers(0, reps, size=(resamples, reps))
                    boot = _contrasts(arrays, indices)
                    cell_points[(truth, noise)] = point
                    cell_boots[(truth, noise)] = boot
                    row = {"true_H0": truth, "noise_fraction": noise,
                           "prior_support": support, "prior_measure": measure,
                           "replicates": reps}
                    for j, metric in enumerate(METRICS):
                        row[metric + "_estimate"] = float(point[j])
                        row[metric + "_ci95_low"] = float(np.quantile(boot[:, j], .025))
                        row[metric + "_ci95_high"] = float(np.quantile(boot[:, j], .975))
                    expected_cell_rows.append(row)
            keys = list(grouped)
            point_macro = np.mean([cell_points[k] for k in keys], axis=0)
            macro_boot = np.mean([cell_boots[k] for k in keys], axis=0)
            for j, metric in enumerate(METRICS):
                expected_aggregate_rows.append({
                    "prior_support": support, "prior_measure": measure, "metric": metric,
                    "equal_cell_mean_difference": float(point_macro[j]),
                    "ci95_low": float(np.quantile(macro_boot[:, j], .025)),
                    "ci95_high": float(np.quantile(macro_boot[:, j], .975)),
                    "truth_noise_cells": len(keys), "replicates_per_cell": reps,
                    "bootstrap_resamples": resamples,
                })
    actual_cells = read_csv(root / "paired_host_weight_cell_contrasts.csv")
    actual_aggregate = read_csv(root / "paired_host_weight_contrasts.csv")
    if len(actual_cells) != 48 or len(actual_aggregate) != 24:
        raise ValueError(f"wrong output row counts: cells={len(actual_cells)}, aggregate={len(actual_aggregate)}")
    def compare(expected_rows, actual_rows, key_fields, numeric_fields, label):
        actual_index = {tuple(row[k] for k in key_fields): row for row in actual_rows}
        if len(actual_index) != len(actual_rows):
            raise ValueError(f"duplicate keys in {label}")
        for expected in expected_rows:
            key = tuple(str(expected[k]) for k in key_fields)
            if key not in actual_index:
                raise ValueError(f"missing {label} row {key}")
            actual = actual_index[key]
            for field in numeric_fields:
                delta = abs(float(actual[field]) - float(expected[field]))
                if not np.isfinite(delta) or delta > atol:
                    raise ValueError(f"{label} mismatch at {key}, {field}: delta={delta}")
        if len(actual_index) != len(expected_rows):
            raise ValueError(f"unexpected extra rows in {label}")
    compare(expected_cell_rows, actual_cells,
            ("true_H0", "noise_fraction", "prior_support", "prior_measure"),
            [f"{m}_{suffix}" for m in METRICS for suffix in ("estimate", "ci95_low", "ci95_high")],
            "paired cell contrasts")
    compare(expected_aggregate_rows, actual_aggregate,
            ("prior_support", "prior_measure", "metric"),
            ("equal_cell_mean_difference", "ci95_low", "ci95_high"),
            "paired aggregate contrasts")
    result = {
        "audit": "PASS", "raw_records": len(raw), "truth_noise_cells": len(expected_cells),
        "cell_contrast_rows_reconstructed": len(expected_cell_rows),
        "aggregate_contrast_rows_reconstructed": len(expected_aggregate_rows),
        "seed": seed, "bootstrap_resamples": resamples, "input_hashes_verified": True,
        "all_paired_point_estimates_and_bootstrap_intervals_reconstructed": True,
        "scope": "internal reproducibility under the synthetic two-host model only",
    }
    (root / "paired_contrast_audit.json").write_text(json.dumps(result, indent=2) + "\\n", encoding="utf-8")
    print(f"PASS: independently reconstructed {len(expected_cell_rows)} cell rows and {len(expected_aggregate_rows)} aggregate rows")
    return result

def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", type=Path)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--bootstrap-resamples", type=int, default=4000)
    args = parser.parse_args()
    audit(args.results_dir, args.seed, args.bootstrap_resamples)

if __name__ == "__main__":
    main()
