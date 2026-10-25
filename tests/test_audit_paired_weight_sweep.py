import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "paired_weight_sweep_audit", ROOT / "scripts" / "audit_paired_weight_sweep.py"
)
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def make_fixture(root, replicates=200, bootstrap=1000, seed=123):
    raw = []
    for z2 in AUDIT.Z2_VALUES:
        for noise in AUDIT.NOISE_FRACTIONS:
            sigma = noise * AUDIT.C_KM_S * 0.03 / AUDIT.H0_TRUE
            for rep in range(replicates):
                host_index = rep % 2
                host_z = AUDIT.Z1 if host_index == 0 else z2
                row = {
                    "z1": AUDIT.Z1, "z2": z2, "redshift_separation": z2 - AUDIT.Z1,
                    "noise_fraction": noise, "replicate": rep, "true_host_index": host_index,
                    "true_host_redshift": host_z, "true_H0": AUDIT.H0_TRUE, "base_seed": seed,
                    "observed_distance_mpc": 100.0 + rep / 1000, "sigma_mpc": sigma,
                }
                for support in AUDIT.PRIOR_SUPPORTS:
                    medians = {}
                    for label, offset in (("correct", 0.0), ("reversed", 0.25)):
                        median = 70.0 + ((rep % 7) - 3) * 0.1 + offset + noise * 2
                        half = 1.0 + noise * 10
                        lower, upper = median - half, median + half
                        row[f"{support}_{label}_median"] = median
                        row[f"{support}_{label}_median_error"] = median - 70.0
                        row[f"{support}_{label}_lower_68"] = lower
                        row[f"{support}_{label}_upper_68"] = upper
                        row[f"{support}_{label}_covered_68"] = int(lower <= 70.0 <= upper)
                        row[f"{support}_{label}_width_68"] = upper - lower
                        medians[label] = median
                    row[f"paired_abs_error_difference_{support}"] = (
                        abs(medians["reversed"] - 70.0) - abs(medians["correct"] - 70.0)
                    )
                raw.append(row)
    raw_path = root / "per_run.csv"
    write_csv(raw_path, raw)

    by_cell = {}
    for row in raw:
        by_cell.setdefault((float(row["z2"]), float(row["noise_fraction"])), []).append(row)
    summaries = []
    cell_index = 0
    for z2 in AUDIT.Z2_VALUES:
        for noise in AUDIT.NOISE_FRACTIONS:
            cell_rows = by_cell[(z2, noise)]
            for support_index, (support, bounds) in enumerate(AUDIT.PRIOR_SUPPORTS.items()):
                summaries.append(AUDIT.reconstruct_summary(
                    cell_rows, support, bounds, replicates, bootstrap, seed,
                    cell_index, support_index,
                ))
            cell_index += 1
    summary_path = root / "summary.csv"
    write_csv(summary_path, summaries)
    metadata = {
        "experiment": "test fixture",
        "audit_status": "PENDING",
        "base_seed": seed,
        "replicates_per_redshift_noise_cell": replicates,
        "bootstrap_resamples": bootstrap,
        "prior_supports": {k: list(v) for k, v in AUDIT.PRIOR_SUPPORTS.items()},
        "z2_values": list(AUDIT.Z2_VALUES),
        "noise_fractions": list(AUDIT.NOISE_FRACTIONS),
        "raw_records": len(raw),
        "summary_rows": len(summaries),
        "posterior_fits": 4 * len(raw),
        "sha256": {
            "per_run.csv": hashlib.sha256(raw_path.read_bytes()).hexdigest(),
            "summary.csv": hashlib.sha256(summary_path.read_bytes()).hexdigest(),
        },
    }
    (root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def test_independent_audit_reconstructs_summaries_and_hashes(tmp_path):
    make_fixture(tmp_path)
    report = AUDIT.audit_directory(tmp_path)
    assert report["audit"] == "PASS"
    assert report["raw_records"] == 9 * 200
    assert report["summary_rows"] == 18
    assert report["posterior_fits"] == 4 * 9 * 200
    metadata = json.loads((tmp_path / "run_metadata.json").read_text())
    assert metadata["audit_status"] == "PASS"
    assert metadata["independent_audit"]["all_summary_metrics_and_bootstrap_intervals_reconstructed"] is True


def test_independent_audit_rejects_summary_tampering(tmp_path):
    make_fixture(tmp_path)
    path = tmp_path / "summary.csv"
    summaries = AUDIT.read_csv(path)
    summaries[0]["correct_rmse"] = str(float(summaries[0]["correct_rmse"]) + 0.5)
    write_csv(path, summaries)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["sha256"]["summary.csv"] = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="correct_rmse mismatch"):
        AUDIT.audit_directory(tmp_path)


def test_bootstrap_reconstruction_is_deterministic():
    values = np.array([0.0, 1.0, 1.0, 0.0, -1.0])
    first = AUDIT.percentile_ci(values, np.random.default_rng(17), 1000)
    second = AUDIT.percentile_ci(values, np.random.default_rng(17), 1000)
    assert first == second
    assert first[0] <= first[1]
