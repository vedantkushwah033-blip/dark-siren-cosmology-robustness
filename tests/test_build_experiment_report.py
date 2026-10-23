import csv
import importlib.util
import json
from pathlib import Path

import pytest

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "build_experiment_report.py"
SPEC = importlib.util.spec_from_file_location("build_experiment_report", MODULE_PATH)
report = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(report)


def test_format_handles_numbers_and_invalid_values():
    assert report.fmt(1.23456) == "1.235"
    assert report.fmt(float("nan")) == "INVALID"
    assert report.fmt(None) == "None"


def test_load_json_and_csv_reject_missing_or_empty_files(tmp_path):
    with pytest.raises(SystemExit, match="Missing or empty JSON"):
        report.load_json(tmp_path / "missing.json")
    empty = tmp_path / "empty.csv"
    empty.write_text("", encoding="utf-8")
    with pytest.raises(SystemExit, match="Missing or empty CSV"):
        report.load_csv(empty)


def test_report_builds_from_validated_synthetic_fixture(tmp_path, monkeypatch):
    results = tmp_path / "results"
    factorial_dir = results / "factorial-sensitivity"
    factorial_dir.mkdir(parents=True)
    scenario_names = [f"scenario_{i}" for i in range(9)]
    metadata = {
        "original_commit": "abc123",
        "base_seed": 2026,
        "requested_replicates_per_scenario": 2,
        "scenarios": [{"scenario": name} for name in scenario_names],
    }
    (factorial_dir / "scenario_summary.json").write_text(json.dumps(metadata), encoding="utf-8")

    summary_fields = [
        "scenario", "n", "coverage_68", "coverage_68_mc_se",
        "coverage_95", "coverage_95_mc_se", "median_bias",
        "rmse_median", "mean_width_68", "mean_width_95",
    ]
    with (factorial_dir / "scenario_summary.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=summary_fields)
        writer.writeheader()
        for name in scenario_names:
            writer.writerow({
                "scenario": name, "n": 2, "coverage_68": 0.5,
                "coverage_68_mc_se": 0.35, "coverage_95": 1.0,
                "coverage_95_mc_se": 0.0, "median_bias": 1.25,
                "rmse_median": 2.5, "mean_width_68": 4.0, "mean_width_95": 8.0,
            })

    with (factorial_dir / "per_run_results.csv").open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=["scenario", "run"])
        writer.writeheader()
        for name in scenario_names:
            for run in range(2):
                writer.writerow({"scenario": name, "run": run})

    contrast_data = {
        "bootstrap_replicates": 100,
        "contrasts": [{
            "label": "test contrast", "metric": "median_error",
            "difference_B_minus_A": 0.75,
            "bootstrap_95_interval": [0.1, 1.2],
        }],
    }
    (factorial_dir / "contrast_analysis.json").write_text(json.dumps(contrast_data), encoding="utf-8")

    output = tmp_path / "report.md"
    monkeypatch.setattr(
        "sys.argv",
        ["build_experiment_report.py", "--results-root", str(results), "--output", str(output)],
    )
    report.main()
    text = output.read_text(encoding="utf-8")
    assert "Main per-run records: **18**" in text
    assert "1.250" in text
    assert "test contrast" in text
    assert "synthetic simulation experiment" in text


def test_report_rejects_incomplete_scenario_count(tmp_path, monkeypatch):
    results = tmp_path / "results"
    factorial_dir = results / "factorial-sensitivity"
    factorial_dir.mkdir(parents=True)
    (factorial_dir / "scenario_summary.json").write_text(
        json.dumps({"requested_replicates_per_scenario": 2, "scenarios": []}),
        encoding="utf-8",
    )
    (factorial_dir / "scenario_summary.csv").write_text("scenario\n", encoding="utf-8")
    (factorial_dir / "per_run_results.csv").write_text("scenario\n", encoding="utf-8")
    monkeypatch.setattr(
        "sys.argv",
        ["build_experiment_report.py", "--results-root", str(results)],
    )
    with pytest.raises(SystemExit, match="Expected exactly nine"):
        report.main()
