"""Independent audit of prior-measure sensitivity raw rows and summaries."""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np
from scipy.integrate import trapezoid
from scipy.stats import norm

C = 299792.458
H0 = 70.0
Z1 = 0.02
Z2S = (0.025, 0.04, 0.06)
NOISES = (0.03, 0.10, 0.20)
SUPPORTS = {"narrow_50_90": (50.0, 90.0), "broad_20_220": (20.0, 220.0)}
MEASURES = ("uniform_H0", "log_uniform_H0")
WEIGHTS = {"correct": (0.7, 0.3), "reversed": (0.3, 0.7)}
GRID = np.linspace(10.0, 250.0, 9601)


def read_csv(path):
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def host_masses(distance, sigma, z2, weights, support, measure):
    mask = (GRID >= support[0]) & (GRID <= support[1])
    h = GRID[mask]
    density = np.ones_like(h) if measure == "uniform_H0" else 1.0 / h
    raw = np.array([
        weights[i] * trapezoid(norm.pdf(distance, loc=C * z / h, scale=sigma) * density, h)
        for i, z in enumerate((Z1, z2))
    ])
    if not np.all(np.isfinite(raw)) or raw.sum() <= 0:
        raise ValueError("invalid reconstructed host component evidence")
    return raw / raw.sum()


def percentile_ci(values, seed, resamples):
    rng = np.random.default_rng(seed)
    values = np.asarray(values, dtype=float)
    draws = np.empty(resamples)
    for start in range(0, resamples, 250):
        stop = min(start + 250, resamples)
        idx = rng.integers(0, len(values), size=(stop - start, len(values)))
        draws[start:stop] = values[idx].mean(axis=1)
    return float(np.quantile(draws, 0.025)), float(np.quantile(draws, 0.975))


def audit(root):
    raw_path, summary_path, metadata_path = (
        root / "per_run.csv", root / "summary.csv", root / "run_metadata.json"
    )
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("audit_status") != "PENDING":
        raise ValueError("metadata must be PENDING before independent audit")
    for name, path in (("per_run.csv", raw_path), ("summary.csv", summary_path)):
        if metadata.get("sha256", {}).get(name) != sha256(path):
            raise ValueError(f"SHA-256 mismatch for {name}")
    reps = int(metadata["replicates_per_redshift_noise_cell"])
    seed = int(metadata["base_seed"])
    boots = int(metadata["bootstrap_resamples"])
    if reps < 200 or boots < 1000:
        raise ValueError("below protocol minimum")
    raw, summaries = read_csv(raw_path), read_csv(summary_path)
    if len(raw) != len(Z2S) * len(NOISES) * reps:
        raise ValueError("raw row count mismatch")
    if len(summaries) != len(Z2S) * len(NOISES) * len(SUPPORTS) * len(MEASURES):
        raise ValueError("summary row count mismatch")

    cells = {}
    expected_columns = set()
    for sname in SUPPORTS:
        for measure in MEASURES:
            for weight in WEIGHTS:
                expected_columns.update({
                    f"{sname}_{measure}_{weight}_host1_mass",
                    f"{sname}_{measure}_{weight}_host2_mass",
                })
    for i, row in enumerate(raw):
        z2, noise, rep = float(row["z2"]), float(row["noise_fraction"]), int(row["replicate"])
        if z2 not in Z2S or noise not in NOISES or not 0 <= rep < reps:
            raise ValueError(f"invalid design key at raw row {i}")
        if int(row["base_seed"]) != seed or float(row["true_H0"]) != H0:
            raise ValueError(f"provenance/truth mismatch at raw row {i}")
        host = int(row["true_host_index"])
        expected_z = Z1 if host == 0 else z2 if host == 1 else None
        if expected_z is None or not math.isclose(float(row["true_host_redshift"]), expected_z, abs_tol=1e-12):
            raise ValueError(f"host index mismatch at raw row {i}")
        sigma = noise * C * 0.03 / H0
        if not math.isclose(float(row["sigma_mpc"]), sigma, rel_tol=1e-12, abs_tol=1e-12):
            raise ValueError(f"sigma mismatch at raw row {i}")
        distance = float(row["observed_distance_mpc"])
        if not math.isfinite(distance) or distance <= 0:
            raise ValueError(f"invalid observed distance at raw row {i}")
        key = (z2, noise)
        cell = cells.setdefault(key, {})
        if rep in cell:
            raise ValueError(f"duplicate replicate {rep} in cell {key}")
        cell[rep] = row
        for support_name, support in SUPPORTS.items():
            for measure in MEASURES:
                for weight_name, weights in WEIGHTS.items():
                    expected = host_masses(distance, sigma, z2, weights, support, measure)
                    prefix = f"{support_name}_{measure}_{weight_name}"
                    for host_idx, name in enumerate(("host1_mass", "host2_mass")):
                        actual = float(row[f"{prefix}_{name}"])
                        if not math.isfinite(actual) or not math.isclose(actual, expected[host_idx], rel_tol=1e-9, abs_tol=1e-10):
                            raise ValueError(f"{prefix}_{name} reconstruction mismatch at raw row {i}")
                    mass1, mass2 = float(row[f"{prefix}_host1_mass"]), float(row[f"{prefix}_host2_mass"])
                    if not math.isclose(mass1 + mass2, 1.0, abs_tol=1e-10):
                        raise ValueError(f"host masses do not sum to one at raw row {i}")
    expected_cells = {(z, n) for z in Z2S for n in NOISES}
    if set(cells) != expected_cells or any(set(cell) != set(range(reps)) for cell in cells.values()):
        raise ValueError("incomplete or duplicate replicate design")

    expected = []
    cell_index = 0
    for z2 in Z2S:
        for noise in NOISES:
            block = [cells[(z2, noise)][r] for r in range(reps)]
            for support_index, (support_name, support) in enumerate(SUPPORTS.items()):
                for measure_index, measure in enumerate(MEASURES):
                    ckey = f"{support_name}_{measure}_correct_host2_mass"
                    rkey = f"{support_name}_{measure}_reversed_host2_mass"
                    diffs = np.array([float(r[rkey]) - float(r[ckey]) for r in block])
                    # Seed schedule matches the prespecified runner protocol.
                    boot_seed = seed + 104729 * cell_index + 1009 * support_index + 37 * measure_index
                    lo, hi = percentile_ci(diffs, boot_seed, boots)
                    expected.append({
                        "z1": Z1, "z2": z2, "noise_fraction": noise,
                        "prior_support": support_name, "prior_min": support[0], "prior_max": support[1],
                        "prior_measure": measure, "replicates": reps,
                        "correct_mean_host2_mass": float(np.mean([float(r[ckey]) for r in block])),
                        "reversed_mean_host2_mass": float(np.mean([float(r[rkey]) for r in block])),
                        "mean_paired_host2_mass_difference_reversed_minus_correct": float(diffs.mean()),
                        "ci95_low": lo, "ci95_high": hi,
                        "bootstrap_resamples": boots, "bootstrap_seed": boot_seed,
                    })
                    cell_index += 1
    if len(expected) != len(summaries):
        raise ValueError("summary count mismatch")
    numeric = (
        "z1", "z2", "noise_fraction", "prior_min", "prior_max",
        "correct_mean_host2_mass", "reversed_mean_host2_mass",
        "mean_paired_host2_mass_difference_reversed_minus_correct", "ci95_low", "ci95_high",
    )
    for i, (actual, wanted) in enumerate(zip(summaries, expected)):
        for field in ("prior_support", "prior_measure"):
            if actual[field] != wanted[field]:
                raise ValueError(f"{field} mismatch at summary row {i}")
        for field in numeric:
            a, b = float(actual[field]), float(wanted[field])
            if not math.isfinite(a) or not math.isclose(a, b, rel_tol=1e-9, abs_tol=1e-10):
                raise ValueError(f"{field} mismatch at summary row {i}")
        for field in ("replicates", "bootstrap_resamples", "bootstrap_seed"):
            if int(actual[field]) != int(wanted[field]):
                raise ValueError(f"{field} mismatch at summary row {i}")
    metadata["audit_status"] = "PASS"
    metadata["independent_audit"] = {
        "raw_records": len(raw), "summary_rows": len(summaries),
        "host_component_masses_reconstructed": True,
        "paired_summary_estimates_and_bootstrap_intervals_reconstructed": True,
        "hashes_verified": True,
    }
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    report = {
        "audit": "PASS", "raw_records": len(raw), "summary_rows": len(summaries),
        "host_component_masses_reconstructed": True,
        "paired_summary_estimates_and_bootstrap_intervals_reconstructed": True,
        "hashes_verified": True,
    }
    (root / "audit_report.json").write_text(json.dumps(report, indent=2) + "\n", encoding="utf-8")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", nargs="?", default="results/prior-measure-sensitivity")
    args = parser.parse_args()
    report = audit(Path(args.results_dir))
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
