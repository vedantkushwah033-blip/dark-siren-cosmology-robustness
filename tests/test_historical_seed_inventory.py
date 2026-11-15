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
