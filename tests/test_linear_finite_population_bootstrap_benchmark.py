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
