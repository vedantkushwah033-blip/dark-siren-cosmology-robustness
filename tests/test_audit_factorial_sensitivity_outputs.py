import csv
import json
from pathlib import Path

import pytest

from scripts.audit_factorial_sensitivity_outputs import EXPECTED, audit


def _write_csv(path: Path, rows: list[dict]) -> None:
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)


def _fixture(root: Path) -> None:
    summary = []
    raw = []
    for scenario_index, (scenario, (target, missingness, host, weighting)) in enumerate(EXPECTED.items()):
        seed = 2026 + scenario_index * 1_000_000
        summary.append({
            "scenario": scenario,
            "target_completeness": target,
            "missingness": missingness,
            "host_generation": host,
            "inference_weighting": weighting,
            "n": 2,
            "seed": seed,
            "median_bias": 0.0,
            "mean_bias": 0.0,
            "rmse_median": 0.0,
            "rmse_mean": 0.0,
            "mean_width_68": 2.0,
            "mean_width_95": 4.0,
            "coverage_68": 1.0,
            "coverage_95": 1.0,
            "host_observed_fraction": 1.0 if target == 1.0 else 0.5,
        })
        for run in (0, 1):
            raw.append({
                "scenario": scenario,
                "run": run,
                "missingness": missingness,
                "host_generation": host,
                "inference_weighting": weighting,
                "target_completeness": target,
                "realized_completeness": target,
                "host_observed": "True" if target == 1.0 or run == 0 else "False",
                "median": 70.0,
                "mean": 70.0,
                "lower_68": 69.0,
                "upper_68": 71.0,
                "lower_95": 68.0,
                "upper_95": 72.0,
                "width_68": 2.0,
                "width_95": 4.0,
            })
    _write_csv(root / "scenario_summary.csv", summary)
    _write_csv(root / "per_run_results.csv", raw)
    (root / "scenario_summary.json").write_text(json.dumps({
        "scenario_count": 9,
        "requested_replicates_per_scenario": 2,
        "base_seed": 2026,
        "scenarios": [
            {**row, "requested_replicates": 2}
            for row in summary
        ],
    }), encoding="utf-8")


def test_audit_accepts_complete_valid_design(tmp_path):
    _fixture(tmp_path)
    result = audit(tmp_path)
    assert result["status"] == "PASS"
    assert result["raw_rows_checked"] == 18


def test_audit_rejects_wrong_completeness(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "per_run_results.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    target_row = next(row for row in rows if row["scenario"] == "random70_uniform_matched")
    target_row["realized_completeness"] = "0.71"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="exact 70% completeness"):
        audit(tmp_path)


def test_audit_rejects_incomplete_complete_catalogue(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "per_run_results.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    target_row = next(row for row in rows if row["scenario"] == "complete_uniform_matched")
    target_row["realized_completeness"] = "0.99"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="complete catalogue \(100%\)"):
        audit(tmp_path)


def test_audit_rejects_nonfinite_posterior(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "per_run_results.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["median"] = "nan"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="non-finite"):
        audit(tmp_path)


def test_audit_rejects_summary_metric_not_reproduced_by_raw_rows(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["median_bias"] = "4.0"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="summary metric median_bias disagrees"):
        audit(tmp_path)


def test_audit_rejects_inconsistent_interval_width(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "per_run_results.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[0]["width_68"] = "3.0"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="width_68 inconsistent"):
        audit(tmp_path)

def test_audit_rejects_json_manifest_design_mismatch(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["scenarios"][0]["inference_weighting"] = "luminosity"
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="JSON manifest design mismatch"):
        audit(tmp_path)


def test_audit_rejects_json_csv_design_mismatch(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["scenarios"][0]["target_completeness"] = 0.7
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="JSON manifest design mismatch"):
        audit(tmp_path)

def test_audit_rejects_json_csv_metric_mismatch(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["scenarios"][0]["coverage_68"] = 0.5
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="JSON/CSV mismatch for summary metric coverage_68"):
        audit(tmp_path)


def test_audit_rejects_manifest_replicate_count_mismatch(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.json"
    manifest = json.loads(path.read_text(encoding="utf-8"))
    manifest["scenarios"][0]["requested_replicates"] = 3
    path.write_text(json.dumps(manifest), encoding="utf-8")
    with pytest.raises(ValueError, match="JSON manifest replicate count mismatch"):
        audit(tmp_path)

def test_audit_rejects_wrong_scenario_seed(tmp_path):
    _fixture(tmp_path)
    path = tmp_path / "scenario_summary.csv"
    with path.open(newline="", encoding="utf-8") as handle:
        rows = list(csv.DictReader(handle))
    rows[1]["seed"] = "2027"
    _write_csv(path, rows)
    with pytest.raises(ValueError, match="seed mismatch"):
        audit(tmp_path)
