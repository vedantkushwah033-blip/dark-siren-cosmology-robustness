"""Independently verify the posterior-calibration raw rows against the declared RNG design.

This checks simulation provenance separately from posterior reconstruction: every
host draw and observed distance must be reproducible from the recorded seed and
prespecified data-generating process. It is not a check of astrophysical realism.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

import numpy as np

C = 299792.458
TRUTHS = (55.0, 65.0, 75.0, 85.0)
NOISES = (0.03, 0.10, 0.20)
TRUE_WEIGHTS = (0.7, 0.3)
Z_HOSTS = (0.02, 0.04)


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def audit_generation(root: Path) -> dict:
    metadata_path = root / "run_metadata.json"
    raw_path = root / "per_run.csv"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    if metadata.get("audit_status") != "PENDING":
        raise ValueError("generation audit must run before posterior audit")
    expected_hash = metadata.get("sha256", {}).get("per_run.csv")
    actual_hash = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    if not expected_hash or expected_hash != actual_hash:
        raise ValueError("SHA-256 mismatch for per_run.csv")

    seed = int(metadata["base_seed"])
    replicates = int(metadata["replicates_per_truth_noise_cell"])
    if replicates < 100:
        raise ValueError("replicate count below protocol minimum")
    rows = read_csv(raw_path)
    expected_count = len(TRUTHS) * len(NOISES) * replicates
    if len(rows) != expected_count or int(metadata.get("raw_records", -1)) != expected_count:
        raise ValueError("raw record count does not match declared design or metadata")

    rng = np.random.default_rng(seed)
    row_index = 0
    for truth in TRUTHS:
        for noise in NOISES:
            sigma = noise * C * 0.04 / truth
            for replicate in range(replicates):
                row = rows[row_index]
                expected_key = (truth, noise, replicate, seed)
                actual_key = (
                    float(row["true_H0"]), float(row["noise_fraction"]),
                    int(row["replicate"]), int(row["base_seed"]),
                )
                if actual_key != expected_key:
                    raise ValueError(f"design/order/seed mismatch at raw row {row_index}")

                host = int(rng.choice(2, p=TRUE_WEIGHTS))
                z_host = Z_HOSTS[host]
                distance_true = C * z_host / truth
                distance = float(rng.normal(distance_true, sigma))
                while distance <= 0:
                    distance = float(rng.normal(distance_true, sigma))

                comparisons = {
                    "true_host_index": (float(row["true_host_index"]), float(host)),
                    "true_host_redshift": (float(row["true_host_redshift"]), z_host),
                    "sigma_mpc": (float(row["sigma_mpc"]), sigma),
                    "observed_distance_mpc": (float(row["observed_distance_mpc"]), distance),
                }
                for field, (actual, expected) in comparisons.items():
                    if not math.isfinite(actual) or not math.isclose(
                        actual, expected, rel_tol=1e-12, abs_tol=1e-12
                    ):
                        raise ValueError(
                            f"generation mismatch for {field} at raw row {row_index}"
                        )
                row_index += 1

    report = {
        "generation_audit": "PASS",
        "raw_records": len(rows),
        "seed": seed,
        "replicates_per_truth_noise_cell": replicates,
        "rng_and_data_generating_process_reproduced": True,
        "input_hash_verified": True,
        "scope": "synthetic provenance only; not astrophysical validation",
    }
    (root / "generation_audit.json").write_text(
        json.dumps(report, indent=2) + "\n", encoding="utf-8"
    )
    return report


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("results_dir", nargs="?", default="results/posterior-calibration")
    args = parser.parse_args()
    report = audit_generation(Path(args.results_dir))
    print(f"Generation audit {report['generation_audit']}: {report['raw_records']} rows.")


if __name__ == "__main__":
    main()
