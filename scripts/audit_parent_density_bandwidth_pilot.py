"""Independent raw-sample audit for parent-density KDE bandwidth pilots."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde


Z_MIN, Z_MAX = 0.01, 0.40


def _hash(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _target(grid: np.ndarray) -> tuple[np.ndarray, np.ndarray]:
    intensity = (
        1.0
        + 2.2 * np.exp(-0.5 * ((grid - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((grid - 0.275) / 0.040) ** 2)
    )
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    raw = intensity * merger
    return merger, raw / float(trapezoid(raw, grid))


def _independent_density(grid: np.ndarray, sample: np.ndarray, multiplier: float,
                         merger: np.ndarray) -> np.ndarray:
    augmented = np.concatenate((sample, 2.0 * grid[0] - sample, 2.0 * grid[-1] - sample))
    n = sample.size
    original_mean = float(np.mean(sample))
    original_variance = float(np.sum((sample-original_mean)**2) / (n-1))
    p_aug = np.full(augmented.size, 1.0 / augmented.size)
    augmented_mean = float(p_aug @ augmented)
    augmented_variance = float(
        np.sum(p_aug * (augmented-augmented_mean)**2) / (1.0 - float(p_aug @ p_aug))
    )
    if original_variance <= 0 or augmented_variance <= 0:
        raise AssertionError("undefined covariance in independent reflected KDE")
    covariance_adjustment = np.sqrt(original_variance / augmented_variance)
    factor = float(multiplier) * n ** (-1.0/5.0) * covariance_adjustment
    weights = np.full(augmented.size, 1.0 / augmented.size)
    kde = gaussian_kde(augmented, weights=weights, bw_method=lambda _kde: factor)
    raw = np.asarray(kde(grid), dtype=float) * merger
    area = float(trapezoid(raw, grid))
    if not np.isfinite(area) or area <= 0:
        raise AssertionError("invalid independent density normalization")
    return raw / area


def audit_bandwidth_pilot(root: Path) -> dict:
    root = root.resolve()
    manifest_path = root / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    if manifest.get("status") != "PASS":
        raise AssertionError("source manifest status is not PASS")
    expected_files = set(manifest["files"])
    actual_files = {
        p.name for p in root.iterdir()
        if p.is_file() and p.name != "manifest.json"
    }
    if expected_files != actual_files or len(actual_files) != int(manifest["file_count"]):
        raise AssertionError("source artifact file set mismatch")
    for record in manifest["files"].values():
        path = root / record["path"]
        if path.stat().st_size != int(record["bytes"]) or _hash(path) != record["sha256"]:
            raise AssertionError(f"source artifact checksum/size mismatch: {record['path']}")

    config = json.loads((root / "config.json").read_text(encoding="utf-8"))
    attempts = pd.read_csv(root / "attempt_metrics.csv")
    samples = pd.read_csv(root / "parent_samples.csv")
    truth_frame = pd.read_csv(root / "truth_density.csv")
    n_rep, n_parent = int(config["n_replicates"]), int(config["n_parent"])
    candidates = sorted(set(float(x) for x in config["bandwidth_multipliers"]))
    grid = truth_frame.redshift_grid.to_numpy(float)
    merger, truth = _target(grid)
    if len(grid) != int(config["grid_points"]) or not np.allclose(
        truth_frame.analytic_target_density.to_numpy(float), truth, rtol=1e-11, atol=1e-13
    ):
        raise AssertionError("analytic target density does not match independent DGP")
    if samples.duplicated(["replicate", "object_id"]).any():
        raise AssertionError("duplicate parent sample identifiers")
    if len(samples) != n_rep*n_parent:
        raise AssertionError("raw parent sample row count mismatch")
    if set(samples.replicate.astype(int)) != set(range(n_rep)):
        raise AssertionError("raw parent sample replicate IDs mismatch")
    if len(attempts) != n_rep*len(candidates):
        raise AssertionError("attempt grid row count mismatch")
    if attempts.duplicated(["replicate", "bandwidth_multiplier"]).any():
        raise AssertionError("duplicate bandwidth attempts")
    if set(attempts.bandwidth_multiplier.astype(float)) != set(candidates):
        raise AssertionError("candidate bandwidth set mismatch")
    maximum_error = 0.0
    for replicate in range(n_rep):
        group = samples.loc[samples.replicate == replicate].sort_values("object_id")
        sample = group.redshift.to_numpy(float)
        if len(sample) != n_parent or not np.all(np.isfinite(sample)) or np.any(
            (sample < Z_MIN) | (sample > Z_MAX)
        ):
            raise AssertionError(f"invalid raw parent sample in replicate {replicate}")
        rows = attempts.loc[attempts.replicate == replicate]
        for row in rows.itertuples(index=False):
            if row.status != "SUCCESS":
                raise AssertionError("non-success parent KDE attempt in bandwidth source")
            density = _independent_density(grid, sample, float(row.bandwidth_multiplier), merger)
            expected_iae = float(trapezoid(np.abs(density-truth), grid))
            difference = abs(expected_iae-float(row.iae))
            maximum_error = max(maximum_error, difference)
            if difference > 2e-10:
                raise AssertionError(
                    f"independent IAE mismatch at replicate={replicate}, multiplier={row.bandwidth_multiplier}"
                )
    return {
        "status": "PASS",
        "replicates_verified": n_rep,
        "parent_units_verified": len(samples),
        "bandwidth_attempts_reconstructed": len(attempts),
        "candidate_multipliers": candidates,
        "maximum_iae_absolute_difference": maximum_error,
        "scope": "Independent reconstruction of synthetic parent-density IAE; no H0 inference.",
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_dir", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    report = audit_bandwidth_pilot(args.source_dir)
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps(report, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
