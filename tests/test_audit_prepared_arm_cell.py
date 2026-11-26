import json
from pathlib import Path
import subprocess
import sys

import pytest


ROOT = Path(__file__).resolve().parents[1]
ARMS = [
    "A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
    "D_map_plugin_ipw", "E_uncertainty_integrated_ipw",
]


def _run(command):
    return subprocess.run(command, cwd=ROOT, capture_output=True, text=True, check=True)


def _prepare(tmp_path, degree):
    bundle = tmp_path / f"bundle_degree_{degree}"
    _run([
        sys.executable, str(ROOT / "scripts/prepare_shared_input_bundle.py"),
        "--seed", str(20261130 + degree), "--event-seed", str(20261150 + degree), "--validation-seed", str(20261140 + degree), "--validation-sample-seed", str(20261160 + degree),
        "--n-parent", "120", "--n-detected", "2", "--validation-frame", "40",
        "--validation-sizes", "20", "--h0-points", "41",
        "--quadrature-order", "3", "--inclusion-degree", str(degree),
        "--curvature-gamma", "0.5", "--output-dir", str(bundle),
    ])
    return bundle


def _execute_and_audit(bundle, tmp_path, arm, degree):
    cell = tmp_path / f"cell_d{degree}_{arm}"
    _run([
        sys.executable, str(ROOT / "scripts/run_arm_cell_from_bundle.py"),
        "--bundle-dir", str(bundle), "--output-dir", str(cell),
        "--replicate", "1", "--gamma", "0.5", "--validation-size", "20", "--arm", arm,
    ])
    report_path = tmp_path / f"audit_d{degree}_{arm}.json"
    completed = _run([
        sys.executable, str(ROOT / "scripts/audit_arm_cell_from_bundle.py"),
        "--bundle-dir", str(bundle), "--cell-dir", str(cell),
        "--output", str(report_path),
    ])
    report = json.loads(report_path.read_text())
    assert report["status"] == "PASS"
    assert report["audit_status"] == "PASS"
    assert report["cell_key"]["arm"] == arm
    assert report["cell_key"]["validation_size"] == 20
    assert "raw prepared inputs" in report["scope"]
    return completed


def test_all_five_arms_audit_directly_from_preparation_only_bundle(tmp_path):
    bundle = _prepare(tmp_path, degree=1)
    for arm in ARMS:
        _execute_and_audit(bundle, tmp_path, arm, degree=1)
    assert not (bundle / "results").exists()


@pytest.mark.parametrize("arm", ["D_map_plugin_ipw", "E_uncertainty_integrated_ipw"])
def test_degree_two_estimated_arms_audit_from_preparation_only_bundle(tmp_path, arm):
    bundle = _prepare(tmp_path, degree=2)
    _execute_and_audit(bundle, tmp_path, arm, degree=2)
    assert not (bundle / "results").exists()


def test_preparation_only_bundle_runs_through_audit_gated_serial_lifecycle(tmp_path):
    bundle = _prepare(tmp_path, degree=1)
    plan = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    cells = [{
        "replicate": 1, "gamma": 0.5, "validation_size": 20,
        "arm": arm, "bundle_dir": str(bundle),
        "output_dir": str(tmp_path / "coordinated" / arm),
    } for arm in ARMS]
    plan.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    completed = _run([
        sys.executable, str(ROOT / "scripts/run_arm_cell_plan.py"),
        "--plan", str(plan), "--checkpoint", str(checkpoint),
    ])
    report = json.loads(completed.stdout)
    assert report["status"] == "PASS"
    assert report["status_counts"] == {
        "SUCCESS": 5, "NON_EVALUABLE": 0, "FAILED": 0, "NOT_RUN": 0,
    }
    from dark_siren_robustness.cell_ledger_store import read_checkpoint
    saved = read_checkpoint(checkpoint)
    assert len(saved["cell_ledger"]) == 5
    for row in saved["cell_ledger"]:
        assert row["status"] == "SUCCESS"
        assert row["audit_status"] == "PASS"
        assert Path(row["audit_report_path"]).is_file()
