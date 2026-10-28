"""Run the synthetic selection-function misspecification pilot."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import subprocess
import sys
import time
from importlib.metadata import PackageNotFoundError, version
from pathlib import Path

from dark_siren_robustness.selection_misspecification import run_selection_misspecification


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as handle:
        for block in iter(lambda: handle.read(1024 * 1024), b""):
            digest.update(block)
    return digest.hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replicates", type=int, default=30)
    parser.add_argument("--n-parent", type=int, default=6000)
    parser.add_argument("--n-detected", type=int, default=40)
    parser.add_argument("--seed", type=int, default=20261029)
    parser.add_argument("--h0-points", type=int, default=161)
    parser.add_argument("--z-points", type=int, default=801)
    parser.add_argument("--bandwidth-factor", type=float, default=0.8)
    parser.add_argument("--output-dir", type=Path, default=Path("artifacts/selection_misspecification"))
    args = parser.parse_args()
    started = time.perf_counter()
    rows, summary = run_selection_misspecification(
        n_replicates=args.replicates, n_parent=args.n_parent,
        n_detected=args.n_detected, seed=args.seed,
        h0_points=args.h0_points, z_points=args.z_points,
        bandwidth_factor=args.bandwidth_factor,
    )
    elapsed = time.perf_counter() - started
    args.output_dir.mkdir(parents=True, exist_ok=True)
    raw_path = args.output_dir / "replicate_metrics.csv"
    summary_path = args.output_dir / "summary.json"
    rows.to_csv(raw_path, index=False)
    summary["runtime_seconds"] = elapsed
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    try:
        commit = subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        commit = "unavailable"
    packages = {}
    for package in ("numpy", "scipy", "pandas"):
        try:
            packages[package] = version(package)
        except PackageNotFoundError:
            packages[package] = "unavailable"
    manifest = {
        "experiment": "selection_function_misspecification_v1",
        "git_commit": commit,
        "python": sys.version,
        "platform": platform.platform(),
        "packages": packages,
        "command": " ".join(sys.argv),
        "runtime_seconds": elapsed,
        "replicate_count": args.replicates,
        "configuration": {
            "n_parent": args.n_parent, "n_detected": args.n_detected,
            "seed_start": args.seed, "h0_points": args.h0_points,
            "z_points": args.z_points, "bandwidth_factor": args.bandwidth_factor,
        },
        "files": {
            raw_path.name: {"sha256": _sha256(raw_path), "bytes": raw_path.stat().st_size},
            summary_path.name: {"sha256": _sha256(summary_path), "bytes": summary_path.stat().st_size},
        },
    }
    (args.output_dir / "manifest.json").write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({"summary": summary, "manifest": manifest}, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
