import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "analyze_factorial_contrasts.py"
SPEC = importlib.util.spec_from_file_location("analyze_factorial_contrasts", MODULE_PATH)
analysis = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(analysis)


def test_value_preserves_valid_zero_pull_median():
    row = {
        "median": "72.0",
        "mean": "72.0",
        "lower_68": "69.0",
        "upper_68": "73.0",
        "lower_95": "67.0",
        "upper_95": "75.0",
        "posterior_sd": "2.0",
        "pull_median": "0",
    }
    assert analysis.value(row, "pull_median") == 0.0


def test_value_calculates_pull_when_field_is_empty():
    row = {
        "median": "72.0",
        "mean": "72.0",
        "lower_68": "69.0",
        "upper_68": "73.0",
        "lower_95": "67.0",
        "upper_95": "75.0",
        "posterior_sd": "2.0",
        "pull_median": "",
    }
    assert analysis.value(row, "pull_median") == 1.0

def test_contrast_input_validation_rejects_mismatched_replicate_counts(tmp_path, monkeypatch):
    import csv
    import json
    import pytest

    names = sorted({s for a, b, _ in analysis.SCENARIO_CONTRASTS for s in (a, b)})
    input_path = tmp_path / "per_run_results.csv"
    fields = ["scenario", "median", "mean", "lower_68", "upper_68", "lower_95", "upper_95", "posterior_sd", "pull_median"]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, name in enumerate(names):
            count = 2 if index == 0 else 3
            for _ in range(count):
                writer.writerow({"scenario": name, "median": 70, "mean": 70, "lower_68": 68, "upper_68": 72, "lower_95": 66, "upper_95": 74, "posterior_sd": 1, "pull_median": 0})
    output_path = tmp_path / "contrast.json"
    monkeypatch.setattr("sys.argv", ["analyze_factorial_contrasts.py", "--input", str(input_path), "--output", str(output_path), "--bootstrap", "100"])
    with pytest.raises(SystemExit, match="Inconsistent replicate counts"):
        analysis.main()


def test_contrast_input_validation_rejects_reversed_interval(tmp_path, monkeypatch):
    import csv
    import pytest

    names = sorted({s for a, b, _ in analysis.SCENARIO_CONTRASTS for s in (a, b)})
    input_path = tmp_path / "per_run_results.csv"
    fields = ["scenario", "median", "mean", "lower_68", "upper_68", "lower_95", "upper_95", "posterior_sd", "pull_median"]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for name in names:
            writer.writerow({"scenario": name, "median": 70, "mean": 70, "lower_68": 73, "upper_68": 72, "lower_95": 66, "upper_95": 74, "posterior_sd": 1, "pull_median": 0})
    output_path = tmp_path / "contrast.json"
    monkeypatch.setattr("sys.argv", ["analyze_factorial_contrasts.py", "--input", str(input_path), "--output", str(output_path), "--bootstrap", "100"])
    with pytest.raises(SystemExit, match="Reversed credible interval"):
        analysis.main()

def test_contrast_input_validation_rejects_nonfinite_posterior_summary(tmp_path, monkeypatch):
    import csv
    import pytest

    names = sorted({s for a, b, _ in analysis.SCENARIO_CONTRASTS for s in (a, b)})
    input_path = tmp_path / "per_run_results.csv"
    fields = ["scenario", "median", "mean", "lower_68", "upper_68", "lower_95", "upper_95", "posterior_sd", "pull_median"]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, name in enumerate(names):
            writer.writerow({"scenario": name, "median": "nan" if index == 0 else 70, "mean": 70, "lower_68": 68, "upper_68": 72, "lower_95": 66, "upper_95": 74, "posterior_sd": 1, "pull_median": 0})
    output_path = tmp_path / "contrast.json"
    monkeypatch.setattr("sys.argv", ["analyze_factorial_contrasts.py", "--input", str(input_path), "--output", str(output_path), "--bootstrap", "100"])
    with pytest.raises(SystemExit, match="Non-finite posterior summary"):
        analysis.main()


def test_contrast_input_validation_rejects_nonpositive_posterior_sd(tmp_path, monkeypatch):
    import csv
    import pytest

    names = sorted({s for a, b, _ in analysis.SCENARIO_CONTRASTS for s in (a, b)})
    input_path = tmp_path / "per_run_results.csv"
    fields = ["scenario", "median", "mean", "lower_68", "upper_68", "lower_95", "upper_95", "posterior_sd", "pull_median"]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, name in enumerate(names):
            writer.writerow({"scenario": name, "median": 70, "mean": 70, "lower_68": 68, "upper_68": 72, "lower_95": 66, "upper_95": 74, "posterior_sd": 0 if index == 0 else 1, "pull_median": 0})
    output_path = tmp_path / "contrast.json"
    monkeypatch.setattr("sys.argv", ["analyze_factorial_contrasts.py", "--input", str(input_path), "--output", str(output_path), "--bootstrap", "100"])
    with pytest.raises(SystemExit, match="Non-positive posterior_sd"):
        analysis.main()

def test_contrast_input_validation_rejects_non_nested_intervals(tmp_path, monkeypatch):
    import csv
    import pytest

    names = sorted({s for a, b, _ in analysis.SCENARIO_CONTRASTS for s in (a, b)})
    input_path = tmp_path / "per_run_results.csv"
    fields = ["scenario", "median", "mean", "lower_68", "upper_68", "lower_95", "upper_95", "posterior_sd", "pull_median"]
    with input_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        for index, name in enumerate(names):
            writer.writerow({"scenario": name, "median": 70, "mean": 70, "lower_68": 68, "upper_68": 72, "lower_95": 66 if index else 69, "upper_95": 74, "posterior_sd": 1, "pull_median": 0})
    output_path = tmp_path / "contrast.json"
    monkeypatch.setattr("sys.argv", ["analyze_factorial_contrasts.py", "--input", str(input_path), "--output", str(output_path), "--bootstrap", "100"])
    with pytest.raises(SystemExit, match="not nested inside 95% interval"):
        analysis.main()