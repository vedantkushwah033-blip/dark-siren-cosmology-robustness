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