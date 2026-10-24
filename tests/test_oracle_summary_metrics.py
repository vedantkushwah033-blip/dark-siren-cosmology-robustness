import importlib.util
from pathlib import Path

import pytest


ROOT = Path(__file__).resolve().parents[1] / "scripts"


def load_script(filename):
    path = ROOT / filename
    spec = importlib.util.spec_from_file_location(filename.removesuffix(".py"), path)
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    return module


@pytest.mark.parametrize(
    "filename",
    ["run_oracle_missing_host_control.py", "run_full_population_oracle.py"],
)
def test_oracle_median_bias_is_median_of_per_run_median_errors(filename):
    module = load_script(filename)
    rows = [
        {"median": 69.0},
        {"median": 71.0},
        {"median": 80.0},
    ]
    # Errors are -1, +1, +10: their median is +1, not their mean (+10/3).
    assert module.median_bias_from_rows(rows, true_h0=70.0) == 1.0


@pytest.mark.parametrize(
    "filename",
    ["run_oracle_missing_host_control.py", "run_full_population_oracle.py"],
)
def test_oracle_median_bias_rejects_empty_rows(filename):
    module = load_script(filename)
    with pytest.raises(ValueError, match="must not be empty"):
        module.median_bias_from_rows([], true_h0=70.0)
