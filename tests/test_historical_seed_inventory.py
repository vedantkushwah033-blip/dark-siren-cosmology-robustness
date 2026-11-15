from __future__ import annotations

import importlib.util
import json
import sys
from pathlib import Path

import pytest

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load():
    path = SCRIPTS / "audit_historical_seed_inventory.py"
    spec = importlib.util.spec_from_file_location("historical_seed_inventory", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _schedule(path: Path, training=None, holdout=None):
    path.write_text(json.dumps({
        "training": training or [
            {"population_seed": 1000, "validation_seed": 1001},
            {"population_seed": 1010, "validation_seed": 1011},
        ],
        "holdout": holdout or [
            {"population_seed": 2000, "validation_seed": 2001},
            {"population_seed": 2010, "validation_seed": 2011},
        ],
    }))
    return path


def test_inventory_extracts_nested_json_and_csv_seed_fields_and_hashes_sources(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "manifest.json").write_text(json.dumps({
        "run": {"population_seed": 11, "validation_seed": 12},
        "children": [{"seed": 13}, {"random_seed": 14}],
    }))
    (root / "ledger.csv").write_text(
        "replicate,population_seed,validation_seed\n0,21,22\n1,31,32\n"
    )
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "PASS"
    assert report["source_file_count"] == 2
    assert report["historical_seed_count"] == 8
    assert all(len(item["sha256"]) == 64 for item in report["source_files"])
    assert report["training_historical_overlap"] == []
    assert report["holdout_historical_overlap"] == []


def test_inventory_fails_closed_on_historical_training_overlap(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({"population_seed": 1000}))
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    assert report["training_historical_overlap"] == [1000]


def test_inventory_fails_closed_on_unsupported_yaml_artifact(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "child_config.yml").write_text("population_seed: 7\n")
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    assert any("unsupported YAML" in item for item in report["inventory_errors"])


def test_schedule_rejects_population_validation_collision(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    schedule = _schedule(
        tmp_path / "schedule.json",
        training=[{"population_seed": 5, "validation_seed": 5}],
    )
    with pytest.raises(ValueError, match="reuses population and validation seed"):
        module.audit([root], schedule)


def test_schedule_rejects_training_holdout_seed_reuse(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    schedule = _schedule(
        tmp_path / "schedule.json",
        training=[{"population_seed": 100, "validation_seed": 101}],
        holdout=[{"population_seed": 101, "validation_seed": 102}],
    )
    with pytest.raises(ValueError, match="training and holdout schedules share"):
        module.audit([root], schedule)


def test_inventory_rejects_non_integer_seed_value(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "bad.json").write_text(json.dumps({"population_seed": "unknown"}))
    with pytest.raises(ValueError, match="non-integer seed"):
        module.audit([root], _schedule(tmp_path / "schedule.json"))


def test_schedule_checks_additional_random_streams(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    schedule = _schedule(
        tmp_path / "schedule.json",
        training=[{"population_seed": 10, "validation_seed": 11, "bootstrap_seed": 12}],
        holdout=[{"population_seed": 20, "validation_seed": 21, "bootstrap_seed": 11}],
    )
    with pytest.raises(ValueError, match="training and holdout schedules share"):
        module.audit([root], schedule)


def test_inventory_does_not_pass_when_sources_contain_no_seed_fields(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "metrics.json").write_text(json.dumps({"iae": 0.2, "replicates": 3}))
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    assert report["historical_seed_count"] == 0


def test_seed_stride_parameter_is_not_misclassified_as_a_random_seed(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({
        "master_seed": 500,
        "seed_stride": 10,
        "population_seed": 501,
        "validation_seed": 502,
    }))
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["historical_seed_values"] == [500, 501, 502]


def test_inventory_ignores_descriptive_seed_schedule_metadata_but_keeps_integer_seeds(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({
        "master_seed": 3000,
        "seed_stride": 10,
        "seed_schedule": "population_seed=master_seed+replicate_id*seed_stride; validation_seed=population_seed+1",
    }))
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "PASS"
    assert report["historical_seed_values"] == [3000]



def test_inventory_marks_missing_artifact_root_not_evaluable(tmp_path):
    module = _load()
    missing = tmp_path / "does-not-exist"
    report = module.audit([missing], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    assert any("does not exist" in item for item in report["inventory_errors"])


def test_inventory_rejects_malformed_json_artifact(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "broken.json").write_text('{"population_seed": ')
    with pytest.raises(json.JSONDecodeError):
        module.audit([root], _schedule(tmp_path / "schedule.json"))


def test_schedule_rejects_missing_required_seed_stream(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    schedule = _schedule(
        tmp_path / "schedule.json",
        training=[{"population_seed": 1000}],
    )
    with pytest.raises(ValueError, match="population_seed and validation_seed"):
        module.audit([root], schedule)


def test_cli_writes_parseable_json_with_terminal_newline(tmp_path):
    import subprocess

    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({"population_seed": 3000}))
    schedule = _schedule(tmp_path / "schedule.json")
    output = tmp_path / "out" / "report.json"
    result = subprocess.run(
        [
            sys.executable,
            str(SCRIPTS / "audit_historical_seed_inventory.py"),
            "--artifact-root", str(root),
            "--schedule", str(schedule),
            "--output", str(output),
        ],
        check=False,
        capture_output=True,
        text=True,
    )
    assert result.returncode == 0
    payload = output.read_bytes()
    assert payload.endswith(b"\n")
    assert json.loads(payload.decode("utf-8"))["status"] == "PASS"


def test_inventory_reads_integer_seed_arrays_from_npz_without_pickle(tmp_path):
    import numpy as np

    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({"master_seed": 3000}))
    np.savez(root / "streams.npz",
             population_seed=np.array([4000, 4010], dtype=np.int64),
             validation_seed=np.array([4001, 4011], dtype=np.int64),
             posterior_nodes=np.array([0.1, 0.2]))
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "PASS"
    assert report["historical_seed_values"] == [3000, 4000, 4001, 4010, 4011]
    npz_records = [row for row in report["source_files"] if row["path"] == "streams.npz"]
    assert len(npz_records) == 1
    assert len(npz_records[0]["seeds"]) == 4


def test_inventory_rejects_object_seed_arrays_without_enabling_pickle(tmp_path):
    import numpy as np

    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    np.savez(root / "unsafe.npz", population_seed=np.array([3000], dtype=object))
    with pytest.raises(ValueError, match="Object arrays cannot be loaded when allow_pickle=False"):
        module.audit([root], _schedule(tmp_path / "schedule.json"))


def test_inventory_ignores_across_seeds_summary_metrics(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "config.json").write_text(json.dumps({"population_seed": 3000}))
    (root / "summary.csv").write_text(
        "mean_pit_mean_across_seeds,mean_coverage_across_seeds\\n"
        "0.5138380959669132,0.681\\n"
    )
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "PASS"
    assert report["historical_seed_values"] == [3000]



def test_ledger_diagnostic_allows_common_random_numbers_across_gamma(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "replicate_ledger.csv").write_text(
        "gamma,replicate,population_seed,validation_seed\n"
        "-1,0,10,11\n0,0,10,11\n"
        "-1,1,20,21\n0,1,20,21\n"
    )
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "PASS"
    assert report["historical_seed_warnings"] == []


def test_ledger_diagnostic_warns_when_validation_seed_is_reused_across_replicates(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "replicate_ledger.csv").write_text(
        "gamma,replicate,population_seed,validation_seed\n"
        "-1,0,10,99\n0,0,10,99\n"
        "-1,1,20,99\n0,1,20,99\n"
    )
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    warnings = report["historical_seed_warnings"]
    assert len(warnings) == 1
    assert warnings[0]["kind"] == "seed_reused_across_replicates"
    assert "seed=99" in warnings[0]["detail"]



def test_ledger_detects_cross_stream_seed_reuse_between_replicates(tmp_path):
    module = _load()
    root = tmp_path / "history"
    root.mkdir()
    (root / "replicate_ledger.csv").write_text(
        "gamma,replicate,population_seed,validation_seed\\n"
        "-1,0,10,11\\n0,0,10,11\\n"
        "-1,1,20,10\\n0,1,20,10\\n"
    )
    report = module.audit([root], _schedule(tmp_path / "schedule.json"))
    assert report["status"] == "NOT_EVALUABLE"
    warnings = report["historical_seed_warnings"]
    assert len(warnings) == 1
    assert warnings[0]["kind"] == "seed_reused_across_replicates"
    assert "seed=10" in warnings[0]["detail"]
