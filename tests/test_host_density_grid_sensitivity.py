import sys
from pathlib import Path

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "scripts"))

import pytest

from audit_host_density_grid_sensitivity import audit_grid_sensitivity
from run_host_density_grid_sensitivity import run_grid_sensitivity


def test_paired_host_density_grid_sensitivity_is_independently_audited(tmp_path):
    folder = tmp_path / "grid_sensitivity"
    run_grid_sensitivity(
        output_dir=folder, grid_points=[51, 101, 201], n_replicates=1,
        n_parent=500, n_detected=5, seed=71123, z_points=101,
    )
    report = audit_grid_sensitivity(folder, tolerance=0.5)
    assert report["status"] == "PASS"
    assert report["finest_grid_points"] == 201
    assert report["replicate_arm_pairs"] == 3
    assert len(report["per_grid"]) == 3


def test_grid_sensitivity_auditor_rejects_tampered_raw_artifact(tmp_path):
    folder = tmp_path / "grid_sensitivity"
    run_grid_sensitivity(
        output_dir=folder, grid_points=[51, 101], n_replicates=1,
        n_parent=500, n_detected=5, seed=71124, z_points=101,
    )
    raw = folder / "grid_51" / "replicate_metrics.csv"
    raw.write_text(raw.read_text(encoding="utf-8") + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match="checksum"):
        audit_grid_sensitivity(folder, tolerance=0.5)
