from __future__ import annotations

import csv
import json

from scripts.run_multiframe_v3_pilot import run_multiframe_v3
from dark_siren_robustness.multiframe_design import MANIFEST_SCHEMA


def pilot_manifest():
    return {
        "schema": MANIFEST_SCHEMA, "root_seed": 20261031,
        "population_size": 80, "frames_per_scenario": 2,
        "outer_repetitions": 3, "bootstrap_replicates": 25,
        "nominal_level": 0.95, "sample_sizes": [10, 30],
        "methods": [
            {"method_id": "wald_srswor"},
            {"method_id": "ordinary_percentile_bootstrap"},
            {"method_id": "studentized_pseudopop_srswor"},
        ],
        "scenarios": [
            {"scenario_id": 3, "name": "Normal", "family": "normal",
             "params": {"loc": 0.0, "scale": 1.0}},
        ],
    }


def test_v3_runner_writes_complete_paired_pilot(tmp_path):
    artifact = run_multiframe_v3(pilot_manifest(), tmp_path)
    assert artifact["status"] == "PILOT_NOT_CONFIRMATORY"
    assert artifact["frame_rows"] == 160
    assert artifact["outer_sample_keys"] == 12
    assert artifact["replicate_rows"] == 36
    assert len(artifact["replicate_csv_sha256"]) == 64

    with (tmp_path / "replicates.csv").open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    assert len(rows) == 36
    grouped = {}
    for row in rows:
        key = (row["scenario_id"], row["frame_id"], row["sample_size"], row["replicate_id"])
        grouped.setdefault(key, set()).add(row["sample_indices"])
    assert len(grouped) == 12
    assert all(len(indices) == 1 for indices in grouped.values())
    assert all(len(methods) == 3 for methods in grouped.values())


def test_v3_runner_is_deterministic_and_method_order_invariant(tmp_path):
    first = tmp_path / "first"
    second = tmp_path / "second"
    manifest = pilot_manifest()
    run_multiframe_v3(manifest, first)
    manifest["methods"] = list(reversed(manifest["methods"]))
    run_multiframe_v3(manifest, second)
    a = (first / "replicates.csv").read_bytes()
    b = (second / "replicates.csv").read_bytes()
    # Row order follows declared methods, so compare keyed numerical payloads.
    def keyed(path):
        with path.open(newline="", encoding="utf-8") as f:
            rows = list(csv.DictReader(f))
        return {
            (r["scenario_id"], r["frame_id"], r["sample_size"], r["replicate_id"], r["method_id"]):
            (r["method_seed"], r["sample_seed"], r["sample_indices"], r["interval_lower"], r["interval_upper"])
            for r in rows
        }
    assert keyed(first / "replicates.csv") == keyed(second / "replicates.csv")
    assert json.loads((first / "manifest.json").read_text())["status"] == "PILOT_NOT_CONFIRMATORY"
