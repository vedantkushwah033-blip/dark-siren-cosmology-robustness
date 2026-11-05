"""Configuration-driven wrapper for the paired nonlinear inclusion screen.

The underlying screen runner remains the single implementation of simulation and
independent auditing. This wrapper validates and hashes a versioned manifest,
passes its settings explicitly, and records the exact canonical configuration.
"""
from __future__ import annotations

import argparse
import hashlib
import importlib.metadata
import json
import platform
from pathlib import Path
import subprocess
import sys

SCHEMA_VERSION = 1
REQUIRED = {
    "schema_version", "replicates", "n_parent", "n_detected",
    "validation_sizes", "seed", "h0_points", "z_points",
    "curvature_levels", "output_dir",
}


def load_config(path: Path) -> tuple[dict, str]:
    try:
        config = json.loads(path.read_text(encoding="utf-8"))
    except (OSError, json.JSONDecodeError) as exc:
        raise ValueError(f"cannot read configuration JSON: {exc}") from exc
    if not isinstance(config, dict):
        raise ValueError("configuration root must be a JSON object")
    missing, unknown = REQUIRED - set(config), set(config) - REQUIRED
    if missing or unknown:
        raise ValueError(f"configuration keys invalid; missing={sorted(missing)}, unknown={sorted(unknown)}")
    if config["schema_version"] != SCHEMA_VERSION:
        raise ValueError(f"unsupported schema_version {config['schema_version']!r}")
    integer_fields = ("replicates", "n_parent", "n_detected", "seed", "h0_points", "z_points")
    if any(type(config[key]) is not int for key in integer_fields):
        raise ValueError("replicates, n_parent, n_detected, seed, h0_points and z_points must be integers")
    if config["replicates"] < 1 or config["n_parent"] < 100 or config["n_detected"] < 1:
        raise ValueError("replicates/n_detected must be positive and n_parent >= 100")
    if config["h0_points"] < 31 or config["z_points"] < 101:
        raise ValueError("h0_points >= 31 and z_points >= 101 are required")
    sizes = config["validation_sizes"]
    if not isinstance(sizes, list) or not sizes or any(type(x) is not int for x in sizes):
        raise ValueError("validation_sizes must be a nonempty list of integer counts")
    if len(set(sizes)) != len(sizes) or any(x < 10 or x > config["n_parent"] for x in sizes):
        raise ValueError("validation_sizes must be unique and lie in [10, n_parent]")
    curvatures = config["curvature_levels"]
    if not isinstance(curvatures, list) or not curvatures:
        raise ValueError("curvature_levels must be a nonempty list")
    if any(type(x) not in (int, float) for x in curvatures) or len(set(curvatures)) != len(curvatures):
        raise ValueError("curvature_levels must contain unique numeric values")
    if not isinstance(config["output_dir"], str) or not config["output_dir"].strip():
        raise ValueError("output_dir must be a nonempty string")
    canonical = json.dumps(config, sort_keys=True, separators=(",", ":"))
    digest = hashlib.sha256(canonical.encode("utf-8")).hexdigest()
    return config, digest


def collect_provenance(output_dir: Path) -> dict:
    try:
        revision = subprocess.run(
            ["git", "rev-parse", "HEAD"], check=True, capture_output=True, text=True
        ).stdout.strip()
    except (OSError, subprocess.CalledProcessError):
        revision = None
    packages = {}
    for package in ("numpy", "scipy", "pandas", "matplotlib"):
        try:
            packages[package] = importlib.metadata.version(package)
        except importlib.metadata.PackageNotFoundError:
            packages[package] = None
    artifact_hashes = {}
    for path in sorted(output_dir.rglob("*")):
        if path.is_file() and path.name not in {"screen_report.json", "provenance.json"}:
            relative = str(path.relative_to(output_dir))
            artifact_hashes[relative] = hashlib.sha256(path.read_bytes()).hexdigest()
    return {
        "source_revision": revision,
        "python_version": sys.version,
        "platform": platform.platform(),
        "package_versions": packages,
        "artifact_sha256": artifact_hashes,
    }


def main() -> None:
    parser = argparse.ArgumentParser(description="Run the nonlinear inclusion screen from a validated JSON manifest.")
    parser.add_argument("--config", type=Path, required=True)
    args = parser.parse_args()
    try:
        config, digest = load_config(args.config)
    except ValueError as exc:
        parser.error(str(exc))

    output_dir = Path(config["output_dir"])
    command = [
        sys.executable, "scripts/run_estimated_inclusion_curvature_screen.py",
        "--replicates", str(config["replicates"]),
        "--n-parent", str(config["n_parent"]),
        "--n-detected", str(config["n_detected"]),
        "--validation-sizes", ",".join(map(str, config["validation_sizes"])),
        "--seed", str(config["seed"]),
        "--h0-points", str(config["h0_points"]),
        "--z-points", str(config["z_points"]),
        "--output-dir", str(output_dir),
    ]
    subprocess.run(command, check=True)
    report_path = output_dir / "screen_report.json"
    report = json.loads(report_path.read_text(encoding="utf-8"))
    report["configuration_schema_version"] = SCHEMA_VERSION
    report["configuration"] = config
    report["configuration_sha256"] = digest
    provenance = collect_provenance(output_dir)
    report["provenance"] = provenance
    (output_dir / "provenance.json").write_text(
        json.dumps(provenance, indent=2, sort_keys=True) + "\\n", encoding="utf-8"
    )
    (output_dir / "experiment_config.canonical.json").write_text(
        json.dumps(config, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    report_path.write_text(json.dumps(report, indent=2, sort_keys=True) + "\n", encoding="utf-8")
    print(json.dumps({"status": report["status"], "configuration_sha256": digest, "output_dir": str(output_dir)}, indent=2))


if __name__ == "__main__":
    main()
