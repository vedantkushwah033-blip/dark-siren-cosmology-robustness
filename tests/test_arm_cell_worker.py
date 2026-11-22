from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest


SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "run_arm_cell.py"
SPEC = importlib.util.spec_from_file_location("run_arm_cell", SCRIPT)
WORKER = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(WORKER)


def _sha256(path: Path) -> str:
    digest = hashlib.sha256()
    with path.open("rb") as stream:
        for block in iter(lambda: stream.read(65536), b""):
            digest.update(block)
    return digest.hexdigest()


def _bundle(root: Path) -> None:
    (root / "observable").mkdir(parents=True)
    (root / "truth").mkdir()
    parent_z = np.linspace(0.0101, 0.3999, 600)
    parent = pd.DataFrame({
        "object_id": np.arange(parent_z.size),
        "redshift": parent_z,
        "inclusion_probability": np.full(parent_z.size, 0.5),
    })
    catalogue = parent.iloc[::2].copy()
    events = pd.DataFrame({"distance_mpc": [275.0, 340.0]})
    validation = pd.DataFrame({
        "validation_id": np.arange(20),
        "redshift": np.linspace(0.015, 0.395, 20),
        "retained": (np.arange(20) % 3 != 0).astype(int),
    })
    parent.to_csv(root / "truth" / "parent_truth.csv", index=False)
    catalogue.to_csv(root / "observable" / "catalogue_observed.csv", index=False)
    events.to_csv(root / "observable" / "gw_events_observed.csv", index=False)
    validation.to_csv(root / "observable" / "validation_n20.csv", index=False)
    config = {
        "catalogue_model": "quadratic_logit",
        "validation_sizes": [20],
        "h0_prior_support": [55.0, 85.0],
        "h0_grid_points": 41,
        "redshift_grid_points": 201,
        "curvature_gamma": 0.0,
        "h0_true": 70.0,
        "sigma_mpc": 140.0,
        "detection_threshold_mpc": 250.0,
        "scott_multiplier": 0.7,
        "quadrature_order": 3,
        "inclusion_model_degree": 1,
        "quadratic_calibration": {"x2_center": 0.2845005579191241},
    }
    (root / "config.json").write_text(json.dumps(config), encoding="utf-8")
    inputs = [
        "config.json",
        "observable/catalogue_observed.csv",
        "observable/gw_events_observed.csv",
        "truth/parent_truth.csv",
    ]
    manifest = {
        "status": "PASS",
        "files": {
            relative: {
                "bytes": (root / relative).stat().st_size,
                "sha256": _sha256(root / relative),
            }
            for relative in inputs
        },
    }
    (root / "manifest.json").write_text(json.dumps(manifest), encoding="utf-8")


@pytest.mark.parametrize("arm", ["A_catalogue_only", "B_parent_oracle"])
def test_worker_executes_one_arm_without_other_arm_inputs(tmp_path, arm):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    _bundle(bundle)
    if arm == "A_catalogue_only":
        (bundle / "truth" / "parent_truth.csv").unlink()
        manifest = json.loads((bundle / "manifest.json").read_text())
        manifest["files"].pop("truth/parent_truth.csv")
        (bundle / "manifest.json").write_text(json.dumps(manifest))
    if arm == "B_parent_oracle":
        (bundle / "observable" / "catalogue_observed.csv").unlink()
        manifest = json.loads((bundle / "manifest.json").read_text())
        manifest["files"].pop("observable/catalogue_observed.csv")
        (bundle / "manifest.json").write_text(json.dumps(manifest))

    output = tmp_path / "cell"
    result = WORKER.run(SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=3,
        validation_size=20, arm=arm,
    ))
    assert result["status"] == "PASS"
    summary = json.loads((output / "summary.json").read_text())
    assert summary["arm"] == arm
    assert summary["status"] == "SUCCESS"
    assert summary["replicate"] == 3
    assert (output / "arm_cell_arrays.npz").is_file()
    output_manifest = json.loads((output / "manifest.json").read_text())
    assert output_manifest["status"] == "PASS"
    for name, entry in output_manifest["files"].items():
        path = output / name
        assert path.stat().st_size == entry["bytes"]
        assert _sha256(path) == entry["sha256"]


def test_worker_rejects_tampered_shared_input(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    _bundle(bundle)
    with (bundle / "observable" / "catalogue_observed.csv").open("a") as stream:
        stream.write("\n")
    with pytest.raises(ValueError, match="byte count|SHA-256"):
        WORKER.run(SimpleNamespace(
            input_dir=bundle, output_dir=tmp_path / "cell", replicate=0,
            validation_size=20, arm="A_catalogue_only",
        ))


def test_worker_rejects_validation_size_not_in_frozen_bundle(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    _bundle(bundle)
    with pytest.raises(ValueError, match="not present"):
        WORKER.run(SimpleNamespace(
            input_dir=bundle, output_dir=tmp_path / "cell", replicate=0,
            validation_size=100, arm="A_catalogue_only",
        ))

AUDIT_SCRIPT = Path(__file__).resolve().parents[1] / "scripts" / "audit_arm_cell.py"
AUDIT_SPEC = importlib.util.spec_from_file_location("audit_arm_cell", AUDIT_SCRIPT)
AUDITOR = importlib.util.module_from_spec(AUDIT_SPEC)
assert AUDIT_SPEC.loader is not None
AUDIT_SPEC.loader.exec_module(AUDITOR)


@pytest.mark.parametrize("arm", ["A_catalogue_only", "B_parent_oracle"])
def test_independent_auditor_reconstructs_worker_outputs(tmp_path, arm):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    _bundle(bundle)
    if arm == "A_catalogue_only":
        (bundle / "truth" / "parent_truth.csv").unlink()
        manifest = json.loads((bundle / "manifest.json").read_text())
        manifest["files"].pop("truth/parent_truth.csv")
        (bundle / "manifest.json").write_text(json.dumps(manifest))
    if arm == "B_parent_oracle":
        (bundle / "observable" / "catalogue_observed.csv").unlink()
        manifest = json.loads((bundle / "manifest.json").read_text())
        manifest["files"].pop("observable/catalogue_observed.csv")
        (bundle / "manifest.json").write_text(json.dumps(manifest))
    output = tmp_path / "cell"
    WORKER.run(SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=3,
        validation_size=20, arm=arm,
    ))
    report = AUDITOR.audit_cell(bundle, output)
    assert report["status"] == "PASS"
    assert report["arm"] == arm
    assert report["maximum_log_likelihood_absolute_difference"] <= 2e-10


def test_independent_auditor_rejects_tampered_cell_array(tmp_path):
    bundle = tmp_path / "bundle"
    bundle.mkdir()
    _bundle(bundle)
    output = tmp_path / "cell"
    WORKER.run(SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=0,
        validation_size=20, arm="A_catalogue_only",
    ))
    with np.load(output / "arm_cell_arrays.npz") as archive:
        arrays = {key: archive[key] for key in archive.files}
    arrays["marginal_log_likelihood"] = arrays["marginal_log_likelihood"] + 0.1
    np.savez_compressed(output / "arm_cell_arrays.npz", **arrays)
    with pytest.raises(AssertionError, match="output manifest mismatch"):
        AUDITOR.audit_cell(bundle, output)
