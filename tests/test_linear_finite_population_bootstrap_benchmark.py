from __future__ import annotations

import json
import sys

import numpy as np
import pytest

from scripts.run_linear_finite_population_bootstrap_benchmark import main, run_benchmark
from scripts.audit_linear_finite_population_bootstrap_benchmark import audit_artifact


def test_linear_benchmark_records_all_methods_and_exact_sample_indices():
    summaries, replicates, frames = run_benchmark(
        scenarios=1,
        repetitions=2,
        population_size=40,
        sample_sizes=(8,),
        bootstrap_replicates=20,
        seed=20261031,
    )
    assert len(summaries) == 3
    assert len(replicates) == 2 * 3
    assert {row["method"] for row in summaries} == {
        "studentized_pseudopop_srswor", "wald_srswor",
        "ordinary_percentile_bootstrap",
    }
    for row in replicates:
        indices = np.asarray(__import__("json").loads(row["sample_indices"]))
        assert indices.shape == (8,)
        assert len(np.unique(indices)) == 8
        assert np.all((indices >= 0) & (indices < 40))
        assert row["interval_lower"] <= row["interval_upper"] or row["failure_code"]


def test_linear_benchmark_rejects_invalid_design_inputs():
    with pytest.raises(ValueError):
        run_benchmark(scenarios=0, repetitions=2, population_size=40,
                      sample_sizes=(8,), bootstrap_replicates=20)
    with pytest.raises(ValueError):
        run_benchmark(scenarios=1, repetitions=2, population_size=40,
                      sample_sizes=(40,), bootstrap_replicates=20)
    with pytest.raises(ValueError):
        run_benchmark(scenarios=1, repetitions=2, population_size=40,
                      sample_sizes=(8, 8), bootstrap_replicates=20)


def test_linear_benchmark_artifact_passes_independent_replay_audit(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_linear_finite_population_bootstrap_benchmark.py",
        "--scenarios", "1",
        "--repetitions", "2",
        "--population-size", "40",
        "--sample-sizes", "8",
        "--bootstrap-replicates", "20",
        "--seed", "20261031",
        "--output-dir", str(tmp_path),
    ])
    main()
    report = audit_artifact(tmp_path)
    assert report["status"] == "PASS"
    assert report["replicates_independently_reconstructed"] == 6
    assert report["finite_frames_verified"] == 1


def test_linear_benchmark_auditor_rejects_tampered_interval_even_with_updated_checksum(tmp_path, monkeypatch):
    from pathlib import Path
    import csv
    import hashlib

    monkeypatch.setattr(sys, "argv", [
        "run_linear_finite_population_bootstrap_benchmark.py",
        "--scenarios", "1", "--repetitions", "2", "--population-size", "40",
        "--sample-sizes", "8", "--bootstrap-replicates", "20",
        "--seed", "20261031", "--output-dir", str(tmp_path),
    ])
    main()
    replicate_path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with replicate_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fieldnames = reader.fieldnames
        rows = list(reader)
    rows[0]["interval_lower"] = str(float(rows[0]["interval_lower"]) + 0.25)
    with replicate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fieldnames)
        writer.writeheader()
        writer.writerows(rows)
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["replicate_csv_sha256"] = hashlib.sha256(replicate_path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="endpoint does not reconstruct"):
        audit_artifact(Path(tmp_path))



def _write_small_linear_artifact(tmp_path, monkeypatch):
    monkeypatch.setattr(sys, "argv", [
        "run_linear_finite_population_bootstrap_benchmark.py",
        "--scenarios", "1", "--repetitions", "2", "--population-size", "40",
        "--sample-sizes", "8", "--bootstrap-replicates", "20",
        "--seed", "20261031", "--output-dir", str(tmp_path),
    ])
    main()


def _refresh_manifest_checksum(tmp_path, key, filename):
    import hashlib
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest[key] = hashlib.sha256((tmp_path / filename).read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")


def test_linear_auditor_rejects_missing_attempt_even_if_checksum_is_updated(tmp_path, monkeypatch):
    import csv

    _write_small_linear_artifact(tmp_path, monkeypatch)
    path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    rows.pop()
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _refresh_manifest_checksum(
        tmp_path, "replicate_csv_sha256", "linear_finite_population_bootstrap_replicates.csv"
    )
    with pytest.raises(ValueError, match="attempt grid is incomplete"):
        audit_artifact(tmp_path)


def test_linear_auditor_rejects_modified_sample_seed_even_if_checksum_is_updated(tmp_path, monkeypatch):
    import csv

    _write_small_linear_artifact(tmp_path, monkeypatch)
    path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    rows[0]["sample_seed"] = str(int(rows[0]["sample_seed"]) + 1)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _refresh_manifest_checksum(
        tmp_path, "replicate_csv_sha256", "linear_finite_population_bootstrap_replicates.csv"
    )
    with pytest.raises(ValueError, match="sample seed does not match"):
        audit_artifact(tmp_path)


def test_linear_auditor_rejects_tampered_inner_attempt_accounting_even_with_updated_checksum(tmp_path, monkeypatch):
    import csv

    _write_small_linear_artifact(tmp_path, monkeypatch)
    path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    studentized = next(row for row in rows if row["method"] == "studentized_pseudopop_srswor")
    studentized["inner_failures"] = str(int(studentized["inner_failures"]) + 1)
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _refresh_manifest_checksum(
        tmp_path, "replicate_csv_sha256", "linear_finite_population_bootstrap_replicates.csv"
    )
    with pytest.raises(ValueError, match="inner attempt provenance mismatch"):
        audit_artifact(tmp_path)


def test_linear_auditor_requires_explicit_na_for_nonstudentized_attempt_fields(tmp_path, monkeypatch):
    import csv

    _write_small_linear_artifact(tmp_path, monkeypatch)
    path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    rows[0]["inner_attempts"] = "0"
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _refresh_manifest_checksum(
        tmp_path, "replicate_csv_sha256", "linear_finite_population_bootstrap_replicates.csv"
    )
    with pytest.raises(ValueError, match="non-studentized method must mark inner attempt fields NA"):
        audit_artifact(tmp_path)


def test_linear_auditor_preserves_v1_artifact_compatibility(tmp_path, monkeypatch):
    import csv
    import hashlib

    _write_small_linear_artifact(tmp_path, monkeypatch)
    replicate_path = tmp_path / "linear_finite_population_bootstrap_replicates.csv"
    with replicate_path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    legacy_fields = [field for field in fields if not field.startswith("inner_")]
    with replicate_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=legacy_fields)
        writer.writeheader()
        writer.writerows([{key: row[key] for key in legacy_fields} for row in rows])
    manifest_path = tmp_path / "manifest.json"
    manifest = json.loads(manifest_path.read_text(encoding="utf-8"))
    manifest["schema"] = "linear_finite_population_bootstrap_benchmark_v1"
    manifest["replicate_csv_sha256"] = hashlib.sha256(replicate_path.read_bytes()).hexdigest()
    manifest_path.write_text(json.dumps(manifest), encoding="utf-8")
    report = audit_artifact(tmp_path)
    assert report["status"] == "PASS"
    assert report["artifact_schema"] == "linear_finite_population_bootstrap_benchmark_v1"
    assert report["inner_attempt_accounting"] == "NOT_AVAILABLE_IN_HISTORICAL_V1"
    assert report["summary_estimands"] == ["conditional_coverage", "conditional_mcse"]


def test_linear_auditor_rejects_tampered_unconditional_coverage_even_with_updated_checksum(tmp_path, monkeypatch):
    import csv

    _write_small_linear_artifact(tmp_path, monkeypatch)
    path = tmp_path / "linear_finite_population_bootstrap_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        reader = csv.DictReader(handle)
        fields, rows = reader.fieldnames, list(reader)
    rows[0]["coverage_unconditional_success_and_cover"] = str(
        float(rows[0]["coverage_unconditional_success_and_cover"]) + 0.1
    )
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)
    _refresh_manifest_checksum(tmp_path, "summary_csv_sha256", "linear_finite_population_bootstrap_summary.csv")
    with pytest.raises(ValueError, match="summary unconditional success-and-coverage mismatch"):
        audit_artifact(tmp_path)
