"""Generate a lightweight parent-density KDE bandwidth training/holdout source artifact."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import sys
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid

from dark_siren_robustness.catalogue_density import estimate_host_density


Z_MIN, Z_MAX = 0.01, 0.40
DEFAULT_MULTIPLIERS = (0.05, 0.08, 0.10, 0.12, 0.15, 0.20, 0.25, 0.30, 0.40)


def parent_intensity(z: np.ndarray) -> np.ndarray:
    z = np.asarray(z, dtype=float)
    return (
        1.0
        + 2.2 * np.exp(-0.5 * ((z - 0.105) / 0.025) ** 2)
        + 1.6 * np.exp(-0.5 * ((z - 0.275) / 0.040) ** 2)
    )


def draw_parent(seed: int, replicate: int, n_parent: int) -> np.ndarray:
    grid = np.linspace(Z_MIN, Z_MAX, 20001)
    intensity = parent_intensity(grid)
    cdf = np.r_[0.0, cumulative_trapezoid(intensity, grid)]
    cdf /= cdf[-1]
    rng = np.random.default_rng(np.random.SeedSequence([int(seed), int(replicate), 1]))
    return np.interp(rng.random(n_parent), cdf, grid)


def run_bandwidth_pilot(
    *, n_replicates: int = 100, n_parent: int = 4000,
    seed: int = 20261114, grid_points: int = 801,
    bandwidth_multipliers: tuple[float, ...] = DEFAULT_MULTIPLIERS,
) -> tuple[pd.DataFrame, pd.DataFrame, pd.DataFrame, dict]:
    multipliers = tuple(sorted(set(float(x) for x in bandwidth_multipliers)))
    if n_replicates < 4 or n_parent < 50 or grid_points < 101:
        raise ValueError("require >=4 replicates, >=50 parent units, and >=101 grid points")
    if len(multipliers) < 2 or any(not np.isfinite(x) or x <= 0 for x in multipliers):
        raise ValueError("at least two distinct positive finite multipliers are required")
    grid = np.linspace(Z_MIN, Z_MAX, grid_points)
    merger = 1.0 + 1.5 * (grid - Z_MIN) / (Z_MAX - Z_MIN)
    truth_raw = parent_intensity(grid) * merger
    truth = truth_raw / float(trapezoid(truth_raw, grid))
    attempts, sample_rows = [], []
    for replicate in range(n_replicates):
        sample = draw_parent(seed, replicate, n_parent)
        for object_id, z in enumerate(sample):
            sample_rows.append({
                "replicate": replicate, "object_id": object_id, "redshift": float(z)
            })
        for multiplier in multipliers:
            density = estimate_host_density(
                grid, sample, merger_weight=merger,
                scott_multiplier=multiplier,
            )
            iae = float(trapezoid(np.abs(density - truth), grid))
            attempts.append({
                "replicate": replicate, "seed": int(seed),
                "gamma": 0.0, "target_retention": 0.25, "n_validation": 300,
                "degree": 0, "bandwidth_multiplier": multiplier,
                "arm": "parent_sample_kde", "status": "SUCCESS", "iae": iae,
                "n_parent": n_parent,
            })
    attempts_frame = pd.DataFrame(attempts)
    samples_frame = pd.DataFrame(sample_rows)
    truth_frame = pd.DataFrame({"redshift_grid": grid, "analytic_target_density": truth})
    config = {
        "schema_version": 1,
        "status": "synthetic_parent_density_bandwidth_training_holdout_source",
        "seed": int(seed), "n_replicates": int(n_replicates),
        "n_parent": int(n_parent), "grid_points": int(grid_points),
        "grid_bounds": [Z_MIN, Z_MAX],
        "bandwidth_multipliers": list(multipliers),
        "dgp": "parent intensity 1 + 2.2*N(0.105,0.025^2) + 1.6*N(0.275,0.040^2)",
        "target": "normalized analytic parent intensity times merger weight 1+1.5*(z-z_min)/(z_max-z_min)",
        "rng_protocol": "numpy SeedSequence([master_seed, replicate_id, 1])",
        "selection_metric": "integrated absolute error against analytic parent density",
        "claim_limit": "Synthetic density-recovery bandwidth selection only; no H0 inference or observational validation.",
    }
    return attempts_frame, samples_frame, truth_frame, config


def _sha256(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def write_artifacts(
    attempts: pd.DataFrame, samples: pd.DataFrame, truth: pd.DataFrame,
    config: dict, output_dir: Path,
) -> dict:
    output_dir.mkdir(parents=True, exist_ok=True)
    if any(output_dir.iterdir()):
        raise ValueError("output directory must be empty")
    payloads = {
        "attempt_metrics.csv": attempts.to_csv(index=False),
        "parent_samples.csv": samples.to_csv(index=False),
        "truth_density.csv": truth.to_csv(index=False),
        "config.json": json.dumps(config, indent=2, sort_keys=True) + "\n",
    }
    for name, content in payloads.items():
        (output_dir / name).write_text(content, encoding="utf-8")
    files = {}
    for name in sorted(payloads):
        path = output_dir / name
        files[name] = {"path": name, "bytes": path.stat().st_size, "sha256": _sha256(path)}
    manifest = {
        "status": "PASS", "file_count": len(files), "files": files,
        "runtime": {
            "python": sys.version, "platform": platform.platform(),
            "numpy": np.__version__, "pandas": pd.__version__,
        },
    }
    (output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=100)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--seed", type=int, default=20261114)
    parser.add_argument("--grid-points", type=int, default=801)
    parser.add_argument("--bandwidth-multipliers", type=lambda raw: tuple(
        float(value.strip()) for value in raw.split(",") if value.strip()
    ), default=DEFAULT_MULTIPLIERS)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()
    attempts, samples, truth, config = run_bandwidth_pilot(
        n_replicates=args.replicates, n_parent=args.n_parent,
        seed=args.seed, grid_points=args.grid_points,
        bandwidth_multipliers=args.bandwidth_multipliers,
    )
    manifest = write_artifacts(attempts, samples, truth, config, args.output_dir)
    print(json.dumps({"config": config, "manifest": manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
