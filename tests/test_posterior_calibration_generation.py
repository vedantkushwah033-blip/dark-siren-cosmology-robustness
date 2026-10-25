import csv, hashlib, importlib.util, json
from pathlib import Path
import pytest

ROOT = Path(__file__).resolve().parents[1]

def load(name, path):
    spec = importlib.util.spec_from_file_location(name, ROOT / "scripts" / path)
    mod = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(mod)
    return mod

RUN = load("generation_runner", "run_posterior_calibration.py")
GEN = load("generation_auditor", "audit_posterior_calibration_generation.py")

def write_csv(path, rows):
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)

def make_fixture(root):
    raw, summaries = RUN.run(replicates=100, seed=99127)
    write_csv(root / "per_run.csv", raw)
    write_csv(root / "summary.csv", summaries)
    metadata = {
        "audit_status": "PENDING",
        "base_seed": 99127,
        "replicates_per_truth_noise_cell": 100,
        "raw_records": len(raw),
        "summary_rows": len(summaries),
        "sha256": {
            name: hashlib.sha256((root / name).read_bytes()).hexdigest()
            for name in ("per_run.csv", "summary.csv")
        },
    }
    (root / "run_metadata.json").write_text(json.dumps(metadata, indent=2) + "\n")

def test_generation_audit_reproduces_all_seeded_raw_rows(tmp_path):
    make_fixture(tmp_path)
    report = GEN.audit_generation(tmp_path)
    assert report["generation_audit"] == "PASS"
    assert report["raw_records"] == 1200
    assert report["rng_and_data_generating_process_reproduced"] is True

def test_generation_audit_rejects_tampered_observation_even_after_hash_refresh(tmp_path):
    make_fixture(tmp_path)
    raw_path = tmp_path / "per_run.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[17]["observed_distance_mpc"] = str(float(rows[17]["observed_distance_mpc"]) + 0.01)
    write_csv(raw_path, rows)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["sha256"]["per_run.csv"] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    with pytest.raises(ValueError, match="generation mismatch"):
        GEN.audit_generation(tmp_path)

def test_generation_audit_rejects_wrong_row_order(tmp_path):
    make_fixture(tmp_path)
    raw_path = tmp_path / "per_run.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0], rows[1] = rows[1], rows[0]
    write_csv(raw_path, rows)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text())
    metadata["sha256"]["per_run.csv"] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n")
    with pytest.raises(ValueError, match="design/order/seed mismatch"):
        GEN.audit_generation(tmp_path)


def _rewrite_raw_and_refresh_hash(root, mutate):
    raw_path = root / "per_run.csv"
    with raw_path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    mutate(rows)
    write_csv(raw_path, rows)
    metadata_path = root / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["sha256"]["per_run.csv"] = hashlib.sha256(raw_path.read_bytes()).hexdigest()
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\\n", encoding="utf-8")


def test_generation_audit_rejects_modified_host_identity_even_after_hash_refresh(tmp_path):
    make_fixture(tmp_path)

    def mutate(rows):
        rows[25]["true_host_index"] = "1" if rows[25]["true_host_index"] == "0" else "0"

    _rewrite_raw_and_refresh_hash(tmp_path, mutate)
    with pytest.raises(ValueError, match="generation mismatch for true_host_index"):
        GEN.audit_generation(tmp_path)


def test_generation_audit_rejects_modified_uncertainty_even_after_hash_refresh(tmp_path):
    make_fixture(tmp_path)

    def mutate(rows):
        rows[31]["sigma_mpc"] = str(float(rows[31]["sigma_mpc"]) * 1.01)

    _rewrite_raw_and_refresh_hash(tmp_path, mutate)
    with pytest.raises(ValueError, match="generation mismatch for sigma_mpc"):
        GEN.audit_generation(tmp_path)


def test_generation_audit_rejects_missing_or_inconsistent_hash(tmp_path):
    make_fixture(tmp_path)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata["sha256"]["per_run.csv"] = "0" * 64
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\\n", encoding="utf-8")
    with pytest.raises(ValueError, match="SHA-256 mismatch"):
        GEN.audit_generation(tmp_path)
