import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_factorial_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("run_factorial_sensitivity", MODULE_PATH)
factorial = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(factorial)


def test_median_bias_uses_median_not_mean():
    # Errors are [0, 0, 0, 100]; mean is 25, median is 0.
    assert factorial.median_bias([70, 70, 70, 170], 70) == 0.0


def test_median_bias_even_sample_averages_middle_errors():
    assert factorial.median_bias([68, 72], 70) == 0.0


def test_median_bias_rejects_empty_input():
    import pytest
    with pytest.raises(ValueError, match="must not be empty"):
        factorial.median_bias([], 70)