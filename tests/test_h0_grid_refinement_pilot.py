from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))
from run_h0_grid_refinement_pilot import (\n    compare_grid_pair, verify_manifest, _assert_same_generated_data, DATA_FILES,\n)
from dark_siren_robustness.h0_posterior import normalize_h0_posterior


def test_identical_smooth_posterior_passes_grid_comparison():
    coarse_grid = np.linspace(20.0, 300.0, 281)
    fine_grid = np.linspace(20.0, 300.0, 561)
    coarse = normalize_h0_posterior(
        coarse_grid, -0.5 * ((coarse_grid - 70.0) / 8.0) ** 2,
        prior_min=20.0, prior_max=300.0,
    )
    fine = normalize_h0_posterior(
        fine_grid, -0.5 * ((fine_grid - 70.0) / 8.0) ** 2,
        prior_min=20.0, prior_max=300.0,
    )
    result = compare_grid_pair(coarse, fine)
    assert result["passed_predeclared_tolerances"] is True
    assert result["max_abs_quantile_difference"] < 0.05


def test_shifted_posterior_fails_grid_comparison():
    coarse_grid = np.linspace(20.0, 300.0, 141)
    fine_grid = np.linspace(20.0, 300.0, 561)
    coarse = normalize_h0_posterior(
        coarse_grid, -0.5 * ((coarse_grid - 70.0) / 8.0) ** 2,
        prior_min=20.0, prior_max=300.0,
    )
    shifted = normalize_h0_posterior(
        fine_grid, -0.5 * ((fine_grid - 72.0) / 8.0) ** 2,
        prior_min=20.0, prior_max=300.0,
    )
    result = compare_grid_pair(coarse, shifted)
    assert result["passed_predeclared_tolerances"] is False


def test_manifest_verification_rejects_tampered_file(tmp_path):
    target = tmp_path / "record.txt"
    target.write_text("original")
    digest = hashlib.sha256(target.read_bytes()).hexdigest()
    (tmp_path / "manifest.json").write_text(json.dumps({"record.txt": digest}))
    assert verify_manifest(tmp_path) == 1
    target.write_text("tampered")
    with pytest.raises(AssertionError, match="manifest integrity failures"):
        verify_manifest(tmp_path)


def _write_generated_data_pair(baseline: Path, refined: Path):
    for relative in DATA_FILES:
        left = baseline / relative
        right = refined / relative
        left.parent.mkdir(parents=True, exist_ok=True)
        right.parent.mkdir(parents=True, exist_ok=True)
        pd.DataFrame({
            "object_id": [10, 11],
            "redshift": [0.10, 0.20],
            "retained": [1, 0],
        }).to_csv(left, index=False)
        pd.DataFrame({
            "object_id": [10, 11],
            "redshift": [0.10 + 1e-15, 0.20 - 1e-15],
            "retained": [1, 0],
        }).to_csv(right, index=False)


def test_generated_data_audit_allows_only_floating_roundoff(tmp_path):
    baseline, refined = tmp_path / "baseline", tmp_path / "refined"
    _write_generated_data_pair(baseline, refined)
    _assert_same_generated_data(baseline, refined)


def test_generated_data_audit_rejects_changed_ids_and_material_float_drift(tmp_path):
    baseline, refined = tmp_path / "baseline", tmp_path / "refined"
    _write_generated_data_pair(baseline, refined)
    target = refined / DATA_FILES[0]
    changed = pd.read_csv(target)
    changed.loc[0, "object_id"] = 99
    changed.to_csv(target, index=False)
    with pytest.raises(AssertionError, match="exact mismatch"):
        _assert_same_generated_data(baseline, refined)

    _write_generated_data_pair(baseline, refined)
    target = refined / DATA_FILES[0]
    changed = pd.read_csv(target)
    changed.loc[0, "redshift"] += 1e-5
    changed.to_csv(target, index=False)
    with pytest.raises(AssertionError, match="numeric mismatch"):
        _assert_same_generated_data(baseline, refined)
