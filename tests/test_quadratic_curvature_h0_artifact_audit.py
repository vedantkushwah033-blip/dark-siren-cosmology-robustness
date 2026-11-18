from __future__ import annotations

import runpy
from pathlib import Path
from types import SimpleNamespace

import pytest


def _smoke_artifact(tmp_path: Path) -> Path:
    runner = runpy.run_path("scripts/run_quadratic_curvature_h0_smoke.py")
    args = SimpleNamespace(
        seed=2026111911,
        validation_seed=2026111912,
        event_seed=2026111913,
        n_parent=500,
        n_validation_frame=300,
        n_validation=100,
        n_events=2,
        gamma=0.5,
        beta=-3.9,
        h0_true=70.0,
        sigma_mpc=140.0,
        threshold=250.0,
        z_min=0.01,
        z_max=0.40,
        h0_min=55.0,
        h0_max=85.0,
        h0_points=31,
        orders=(3, 5),
        output_dir=str(tmp_path / "artifact"),
    )
    runner["run"](args)
    return Path(args.output_dir)


def test_independent_h0_artifact_auditor_reconstructs_likelihood(tmp_path):
    root = _smoke_artifact(tmp_path)
    auditor = runpy.run_path("scripts/audit_quadratic_curvature_h0_artifact.py")
    report = auditor["audit"](root)
    assert report["status"] == "PASS"
    assert report["files_manifest_verified"] == 9
    assert report["detected_events"] == 2
    assert report["coefficient_nodes_reconstructed"] > 0
    assert report["max_abs_reconstruction_error"] < 2e-10


def test_independent_h0_artifact_auditor_rejects_manifest_tampering(tmp_path):
    root = _smoke_artifact(tmp_path)
    path = root / "observable/events.csv"
    path.write_text(path.read_text() + "\n")
    auditor = runpy.run_path("scripts/audit_quadratic_curvature_h0_artifact.py")
    with pytest.raises(AssertionError, match="manifest integrity failure"):
        auditor["audit"](root)
