import csv
import json
from pathlib import Path
import sys

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_distance_uncertainty_sensitivity import audit_artifacts
from dark_siren_robustness.uncertainty_sensitivity import ARMS, run_experiment


def _write_artifacts(folder, replicates=2, seed=123, grid_points=101):
    raw, summaries = run_experiment(replicates, seed, grid_points)
    folder.mkdir(parents=True)
    raw_path = folder / "replicates.csv"
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw[0]))
        writer.writeheader()
        writer.writerows(raw)
    metadata = {
        "status": "synthetic sensitivity diagnostic; not an observational GW analysis",
        "seed": seed,
        "replicates_per_cell": replicates,
        "grid_points": grid_points,
        "prior": [50.0, 90.0],
        "truths": [50, 55, 60, 65, 70, 75, 80, 85, 90],
        "fractional_errors": [0.05, 0.10, 0.15, 0.20],
        "arms": list(ARMS),
        "raw_record_count": len(raw),
        "summary_cell_count": len(summaries),
        "boundary_truths_must_be_interpreted_separately": True,
        "summary": summaries,
    }
    (folder / "summary.json").write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")


def test_independent_audit_reconstructs_raw_rows_and_summaries(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    report = audit_artifacts(folder)
    assert report["status"] == "PASS"
    assert report["raw_record_count"] == 9 * 4 * 2 * 3
    assert report["summary_cell_count"] == 9 * 4 * 3
    assert report["a_b_pairs_verified"] == 9 * 4 * 2
    assert report["posterior_rows_reconstructed"] == report["raw_record_count"]


def test_independent_audit_rejects_tampered_observation(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    raw_path = folder / "replicates.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["observed_distance_mpc"] = str(float(rows[0]["observed_distance_mpc"]) + 1.0)
    with raw_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="observed distance"):
        audit_artifacts(folder)


def test_independent_audit_rejects_tampered_summary(tmp_path):
    folder = tmp_path / "artifacts"
    _write_artifacts(folder)
    path = folder / "summary.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["summary"][0]["coverage_68"] = 1.0
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="coverage_68 mismatch"):
        audit_artifacts(folder)
