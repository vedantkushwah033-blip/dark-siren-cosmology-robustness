import csv
import json
import sys
from pathlib import Path

import pytest

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))
from audit_host_density_benchmark import audit_artifacts
from dark_siren_robustness.host_density_benchmark import run_benchmark


def _write_artifacts(folder):
    rows, summary = run_benchmark(
        n_replicates=2, n_parent=1200, n_detected=8, seed=55123,
        h0_points=51, z_points=201,
    )
    folder.mkdir(parents=True)
    rows.to_csv(folder / "replicate_metrics.csv", index=False)
    (folder / "summary.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")


def test_independent_audit_reconstructs_paired_rows_and_aggregates(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    report = audit_artifacts(folder)
    assert report["status"] == "PASS"
    assert report["raw_record_count"] == 6
    assert report["replicates_verified"] == 2
    assert report["posterior_normalizations_reconstructed"] == 6
    assert report["posterior_summaries_reconstructed"] == 6
    assert report["posterior_grid_independently_reconstructed"] is True
    assert "independently integrated" in report["audit_boundary"]
    assert report["aggregate_metrics_reconstructed"] == 3
    assert len(report["replicate_metrics_sha256"]) == 64


def test_independent_audit_rejects_tampered_interval_width(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    path = folder / "replicate_metrics.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["ci68_width"] = str(float(rows[0]["ci68_width"]) + 1.0)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="interval width"):
        audit_artifacts(folder)


def test_independent_audit_rejects_tampered_aggregate_summary(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    path = folder / "summary.json"
    metadata = json.loads(path.read_text(encoding="utf-8"))
    metadata["aggregate_metrics"][0]["rmse"] += 1.0
    path.write_text(json.dumps(metadata), encoding="utf-8")
    with pytest.raises(ValueError, match="aggregate rmse mismatch"):
        audit_artifacts(folder)


def test_independent_audit_rejects_missing_paired_arm(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    path = folder / "replicate_metrics.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows = [row for row in rows if not (row["replicate"] == "1" and row["arm"] == "oracle")]
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="raw record count"):
        audit_artifacts(folder)
def test_independent_audit_rejects_tampered_posterior_grid(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    path = folder / "replicate_metrics.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    # Alter the grid while preserving all summary columns: a scalar-only audit
    # would miss this inconsistency.
    density = json.loads(rows[0]["posterior_density_grid_json"])
    density[10] *= 1.25
    rows[0]["posterior_density_grid_json"] = json.dumps(density, separators=(",", ":"))
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="reconstructed"):
        audit_artifacts(folder)


def test_independent_audit_rejects_wrong_posterior_grid_length(tmp_path):
    folder = tmp_path / "benchmark"
    _write_artifacts(folder)
    path = folder / "replicate_metrics.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
        columns = list(rows[0])
    rows[0]["posterior_density_grid_json"] = json.dumps([0.0, 1.0])
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=columns)
        writer.writeheader()
        writer.writerows(rows)
    with pytest.raises(ValueError, match="grid length"):
        audit_artifacts(folder)
