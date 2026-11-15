import csv
import hashlib
import importlib.util
import json
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]


def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


RUNNER = load("prior_measure_runner", "run_prior_measure_sensitivity.py")
AUDIT = load("prior_measure_auditor", "audit_prior_measure_sensitivity.py")


def write_rows(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def build_fixture(root):
    raw, summaries = RUNNER.run(replicates=200, seed=811, bootstrap_resamples=1000)
    write_rows(root / "per_run.csv", raw)
    write_rows(root / "summary.csv", summaries)
    metadata = {
        "experiment": "test fixture",
        "audit_status": "PENDING",
        "base_seed": 811,
        "replicates_per_redshift_noise_cell": 200,
        "bootstrap_resamples": 1000,
        "raw_records": len(raw),
        "summary_rows": len(summaries),
        "sha256": {
            name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in ("per_run.csv", "summary.csv")
        },
    }
    (root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")


def test_component_masses_normalize_and_measures_can_differ():
    uniform = RUNNER.component_masses(150.0, 10.0, [0.02, 0.06], [0.7, 0.3], (20, 220), "uniform_H0")
    log_uniform = RUNNER.component_masses(150.0, 10.0, [0.02, 0.06], [0.7, 0.3], (20, 220), "log_uniform_H0")
    assert sum(uniform) == pytest.approx(1.0)
    assert sum(log_uniform) == pytest.approx(1.0)
    assert abs(uniform[1] - log_uniform[1]) > 1e-5


def test_independent_audit_reconstructs_masses_and_bootstrap_summaries(tmp_path):
    build_fixture(tmp_path)
    report = AUDIT.audit(tmp_path)
    assert report["audit"] == "PASS"
    assert report["raw_records"] == 9 * 200
    assert report["summary_rows"] == 36
    metadata = json.loads((tmp_path / "run_metadata.json").read_text())
    assert metadata["audit_status"] == "PASS"
    assert metadata["independent_audit"]["host_component_masses_reconstructed"] is True


def test_independent_audit_rejects_tampered_summary_even_after_hash_refresh(tmp_path):
    build_fixture(tmp_path)
    path = tmp_path / "summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        summaries = list(csv.DictReader(handle))
    summaries[0]["correct_mean_host2_mass"] = str(float(summaries[0]["correct_mean_host2_mass"]) + 0.01)
    write_rows(path, summaries)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["sha256"]["summary.csv"] = hashlib.sha256(path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    with pytest.raises(ValueError, match="correct_mean_host2_mass mismatch"):
        AUDIT.audit(tmp_path)
