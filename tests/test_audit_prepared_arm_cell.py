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
        "--seed", str(20261130 + degree), "--validation-seed", str(20261140 + degree),
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
