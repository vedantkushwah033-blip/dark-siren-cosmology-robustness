from __future__ import annotations

import hashlib
import json
import sys
from pathlib import Path

import numpy as np
import pandas as pd
import pytest

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT / "scripts"))

from evaluate_heldout_weight_diagnostic_prediction import (
    _fit_predict,
    _holm_adjust,
    _one_cell,
    _verify_manifest,
)


def test_fit_predict_recovers_linear_signal_using_training_scaling():
    x_train = np.linspace(0.0, 1.0, 50)
    y_train = 2.0 + 3.0 * x_train
    x_test = np.array([0.1, 0.5, 0.9])
    pred, baseline, slope = _fit_predict(x_train, y_train, x_test)
    assert pred == pytest.approx(2.0 + 3.0 * x_test, abs=1e-12)
    assert baseline == pytest.approx(float(y_train.mean()))
    assert slope == pytest.approx(3.0, abs=1e-12)


def test_fit_predict_falls_back_to_training_mean_for_constant_diagnostic():
    pred, baseline, slope = _fit_predict(
        np.ones(20), np.linspace(1.0, 2.0, 20), np.array([1.0, 1.0, 1.0])
    )
    assert pred == pytest.approx(np.full(3, baseline))
    assert slope == 0.0


def test_holm_adjustment_is_monotone_in_sorted_p_values():
    raw = [0.01, 0.04, 0.03, 0.20]
    adjusted = _holm_adjust(raw)
    assert adjusted == pytest.approx([0.04, 0.09, 0.09, 0.20])
    assert all(0.0 <= value <= 1.0 for value in adjusted)


def test_one_cell_bootstrap_is_deterministic_and_uses_holdout():
    rng = np.random.default_rng(123)
    x_train = rng.uniform(0.0, 1.0, 200)
    y_train = 1.0 + 4.0 * x_train + rng.normal(0.0, 0.1, 200)
    x_holdout = rng.uniform(0.0, 1.0, 100)
    y_holdout = 1.0 + 4.0 * x_holdout + rng.normal(0.0, 0.1, 100)
    train = pd.DataFrame({"x": x_train, "absolute_median_error": y_train})
    holdout = pd.DataFrame({"x": x_holdout, "absolute_median_error": y_holdout})
    a = _one_cell(train, holdout, x_column="x", resamples=1000, seed=77)
    b = _one_cell(train, holdout, x_column="x", resamples=1000, seed=77)
    assert a == b
    assert a["holdout_delta_mae_model_minus_baseline"] < 0.0
    assert a["delta_mae_bootstrap_high"] < 0.0


def test_one_cell_rejects_wrong_training_denominator():
    train = pd.DataFrame({"x": np.arange(20), "absolute_median_error": np.arange(20)})
    holdout = pd.DataFrame({"x": np.arange(100), "absolute_median_error": np.arange(100)})
    with pytest.raises(AssertionError, match="expected 200 training"):
        _one_cell(train, holdout, x_column="x", resamples=1000, seed=1)



def test_manifest_verifier_checks_declared_file_hash_and_size(tmp_path: Path):
    payload = tmp_path / "replicate_ledger.csv"
    payload.write_text("cell_id,status\\n1,SUCCESS\\n", encoding="utf-8")
    digest = hashlib.sha256(payload.read_bytes()).hexdigest()
    (tmp_path / "manifest.json").write_text(json.dumps({
        "file_count": 1,
        "files": {"replicate_ledger.csv": {
            "bytes": payload.stat().st_size,
            "sha256": digest,
        }},
    }), encoding="utf-8")
    assert _verify_manifest(tmp_path) == 1

    payload.write_text("cell_id,status\\n1,FAILED\\n", encoding="utf-8")
    with pytest.raises(AssertionError, match="manifest checksum/size mismatch"):
        _verify_manifest(tmp_path)


def test_manifest_verifier_rejects_path_traversal(tmp_path: Path):
    (tmp_path / "manifest.json").write_text(json.dumps({
        "file_count": 1,
        "files": {"../outside.csv": {"bytes": 0, "sha256": "0" * 64}},
    }), encoding="utf-8")
    with pytest.raises(AssertionError, match="unsafe manifest path"):
        _verify_manifest(tmp_path)
