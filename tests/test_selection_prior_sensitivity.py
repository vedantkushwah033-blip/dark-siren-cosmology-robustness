import importlib.util
from pathlib import Path
import numpy as np
import pytest

ROOT = Path(__file__).resolve().parents[1]
PATH = ROOT / "scripts" / "run_selection_prior_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("selection_prior_sensitivity", PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_prior_sensitivity_reproducible_and_has_expected_cells():
    rows1, summaries1 = module.run_audit(replicates_per_truth=100, seed=45, grid_size=1001)
    rows2, summaries2 = module.run_audit(replicates_per_truth=100, seed=45, grid_size=1001)
    assert rows1 == rows2
    assert summaries1 == summaries2
    assert len(rows1) == 900
    assert len(summaries1) == 27
    assert all(0 < r["observed_distance_mpc"] < module.DETECTION_DISTANCE_MAX for r in rows1)


def test_prior_sensitivity_records_paired_uncertainty():
    _, summaries = module.run_audit(replicates_per_truth=100, seed=46, grid_size=1001)
    paired = [r for r in summaries if r["method"] == "paired_corrected_minus_uncorrected_absolute_error"]
    assert len(paired) == 9
    assert all(np.isfinite(r["paired_change_mcse"]) and r["paired_change_mcse"] >= 0 for r in paired)


def test_prior_sensitivity_rejects_small_replicate_count():
    with pytest.raises(ValueError, match="at least 100"):
        module.run_audit(replicates_per_truth=20)
