import ast
import json
from pathlib import Path
import subprocess
import sys

from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest


REPO_ROOT = Path(__file__).resolve().parents[1]


def test_preparation_script_has_no_inference_imports():
    source = (REPO_ROOT / "scripts/prepare_shared_input_bundle.py").read_text()
    tree = ast.parse(source)
    imported = set()
    for node in ast.walk(tree):
        if isinstance(node, ast.Import):
            imported.update(alias.name for alias in node.names)
        elif isinstance(node, ast.ImportFrom) and node.module:
            imported.add(node.module)
    assert not any(
        name.endswith(("five_arm_comparison", "five_arm_kde_likelihood",
                       "arm_cell_execution", "inclusion_posterior",
                       "h0_posterior", "catalogue_density"))
        for name in imported
    )


def test_preparation_only_bundle_is_frozen_and_contains_no_inference_outputs(tmp_path):
    output = tmp_path / "prepared_bundle"
    command = [
        sys.executable, str(REPO_ROOT / "scripts/prepare_shared_input_bundle.py"),
        "--seed", "20261122", "--event-seed", "20261124", "--validation-seed", "20261123", "--validation-sample-seed", "20261125",
        "--n-parent", "120", "--n-detected", "2",
        "--validation-frame", "60", "--validation-sizes", "20,40",
        "--curvature-gamma", "0.5", "--h0-points", "41",
        "--output-dir", str(output),
    ]
    completed = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    result = json.loads(completed.stdout)
    assert result["status"] == "PASS"
    assert result["preparation_only"] is True
    assert result["inference_executed"] is False
    frozen = load_shared_input_manifest(
        output, expected_validation_sizes=[20, 40],
        expected_provenance={"seed": 20261122, "curvature_gamma": 0.5},
    )
    assert frozen["status"] == "FROZEN"
    assert (output / "observable/validation_n20.csv").is_file()
    assert (output / "observable/validation_n40.csv").is_file()
    assert not (output / "results").exists()
    assert not list(output.rglob("*.npz"))
    config = json.loads((output / "config.json").read_text())
    assert config["status"] == "shared_inputs_only_no_inference"
    assert config["validation_sizes"] == [20, 40]
    audited = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts/audit_prepared_shared_input_bundle.py"), str(output)],
        cwd=REPO_ROOT, capture_output=True, text=True, check=True,
    )
    report = json.loads(audited.stdout)
    assert report["status"] == "PASS"
    assert report["inference_executed"] is False
    assert report["validation_sizes_verified"] == [20, 40]


def test_preparation_refuses_to_overwrite_existing_bundle(tmp_path):
    output = tmp_path / "already_exists"
    output.mkdir()
    command = [
        sys.executable, str(REPO_ROOT / "scripts/prepare_shared_input_bundle.py"),
        "--seed", "200", "--event-seed", "201", "--validation-seed", "202", "--validation-sample-seed", "203", "--n-parent", "120", "--n-detected", "1", "--validation-frame", "20",
        "--validation-sizes", "10", "--output-dir", str(output),
    ]
    completed = subprocess.run(command, cwd=REPO_ROOT, capture_output=True, text=True)
    assert completed.returncode != 0
    assert "refusing to overwrite existing bundle" in completed.stderr


def test_preparation_auditor_rejects_tampered_raw_input(tmp_path):
    output = tmp_path / "tampered_bundle"
    prepare = [
        sys.executable, str(REPO_ROOT / "scripts/prepare_shared_input_bundle.py"),
        "--seed", "300", "--event-seed", "301", "--validation-seed", "302", "--validation-sample-seed", "303", "--n-parent", "120", "--n-detected", "2", "--validation-frame", "30",
        "--validation-sizes", "20", "--h0-points", "41", "--output-dir", str(output),
    ]
    subprocess.run(prepare, cwd=REPO_ROOT, capture_output=True, text=True, check=True)
    with (output / "truth/parent_truth.csv").open("a") as stream:
        stream.write("tampered\\n")
    audited = subprocess.run(
        [sys.executable, str(REPO_ROOT / "scripts/audit_prepared_shared_input_bundle.py"), str(output)],
        cwd=REPO_ROOT, capture_output=True, text=True,
    )
    assert audited.returncode != 0
    assert "hash/size mismatch" in audited.stderr or "hash/size mismatch" in audited.stdout
