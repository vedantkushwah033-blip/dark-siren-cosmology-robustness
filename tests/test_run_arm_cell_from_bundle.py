from pathlib import Path
import json
from types import SimpleNamespace

import numpy as np
import pandas as pd
import pytest

from dark_siren_robustness.shared_input_bundle import write_shared_input_manifest
from scripts import run_arm_cell_from_bundle as runner


def _bundle(root: Path) -> Path:
    rows = {
        "observable/catalogue_observed.csv": pd.DataFrame({
            "object_id": [1, 2], "redshift": [0.1, 0.2],
        }),
        "observable/gw_events_observed.csv": pd.DataFrame({
            "distance_mpc": [300.0],
        }),
        "observable/validation_n5.csv": pd.DataFrame({
            "redshift": [0.1, 0.15, 0.2, 0.25, 0.3],
            "retained": [1, 0, 1, 0, 1],
        }),
        "truth/parent_truth.csv": pd.DataFrame({
            "object_id": [1, 2, 3], "redshift": [0.1, 0.2, 0.3],
            "inclusion_probability": [0.8, 0.7, 0.6],
        }),
        "truth/gw_event_truth.csv": pd.DataFrame({"event_id": [1]}),
        "truth/gw_generation_trace_truth.csv": pd.DataFrame({"attempt": [1]}),
        "truth/validation_frame_truth.csv": pd.DataFrame({"object_id": [1]}),
    }
    for relative, frame in rows.items():
        path = root / relative
        path.parent.mkdir(parents=True, exist_ok=True)
        frame.to_csv(path, index=False)
    config = {
        "seed": 11, "validation_sizes": [5], "sigma_mpc": 140.0,
        "h0_prior_support": [50.0, 90.0], "h0_grid_points": 5,
        "redshift_grid_points": 7, "inclusion_model_degree": 1,
        "h0_true": 70.0, "detection_threshold_mpc": 250.0,
        "scott_multiplier": 0.7, "quadrature_order": 3,
        "fitted_quadratic_basis_center": 0.0,
    }
    (root / "config.json").write_text(json.dumps(config), encoding="utf-8")
    write_shared_input_manifest(
        root, validation_sizes=[5],
        provenance={"seed": 11, "curvature_gamma": 0.5},
    )
    return root


def _fake_result():
    return SimpleNamespace(
        raw_weights=np.array([1.0, 1.0]),
        component_weights=np.array([[1.0, 1.0]]),
        mixture_weights=np.array([1.0]),
        component_density=np.ones((1, 7)) / 0.39,
        component_log_likelihood=np.zeros((1, 5)),
        marginal_log_likelihood=np.zeros(5),
        component_log_detection_efficiency=np.zeros((1, 5)),
        posterior_summary={"posterior_median": 70.0, "q025": 55.0, "q975": 85.0},
        weight_diagnostics={"kish_ess": 2.0},
        realized_bandwidths=np.array([0.03]),
        inclusion_mode=None, inclusion_nodes=None,
        inclusion_quadrature_weights=None,
    )


def test_runner_executes_one_cell_and_marks_output_unaudited(tmp_path, monkeypatch):
    bundle = _bundle(tmp_path / "bundle")
    captured = {}

    def fake_evaluate(arm, **kwargs):
        captured["arm"] = arm
        captured["kwargs"] = kwargs
        return _fake_result()

    monkeypatch.setattr(runner, "evaluate_arm_cell", fake_evaluate)
    output = tmp_path / "cells" / "one"
    result = runner.run_cell(
        bundle, output, replicate=3, gamma=0.5,
        validation_size=5, arm="A_catalogue_only",
    )
    assert result["status"] == "COMPUTED_NOT_AUDITED"
    assert captured["arm"] == "A_catalogue_only"
    assert captured["kwargs"]["validation_redshift"] is None
    assert captured["kwargs"]["validation_retained"] is None
    metadata = json.loads((output / "cell_result.json").read_text())
    assert metadata["status"] == "COMPUTED_NOT_AUDITED"
    assert metadata["cell_key"]["replicate"] == 3
    manifest = json.loads((output / "manifest.json").read_text())
    assert set(manifest["files"]) == {"cell_arrays.npz", "cell_result.json"}


def test_runner_rejects_gamma_or_validation_size_outside_frozen_bundle(tmp_path, monkeypatch):
    bundle = _bundle(tmp_path / "bundle")
    monkeypatch.setattr(runner, "evaluate_arm_cell", lambda *args, **kwargs: _fake_result())
    with pytest.raises(ValueError, match="gamma"):
        runner.run_cell(bundle, tmp_path / "wrong_gamma", replicate=1,
                        gamma=-0.5, validation_size=5, arm="A_catalogue_only")
    with pytest.raises(ValueError, match="validation_size"):
        runner.run_cell(bundle, tmp_path / "wrong_size", replicate=1,
                        gamma=0.5, validation_size=10, arm="A_catalogue_only")


def test_runner_never_overwrites_existing_cell_output(tmp_path, monkeypatch):
    bundle = _bundle(tmp_path / "bundle")
    output = tmp_path / "existing"
    output.mkdir()
    monkeypatch.setattr(runner, "evaluate_arm_cell", lambda *args, **kwargs: _fake_result())
    with pytest.raises(FileExistsError):
        runner.run_cell(bundle, output, replicate=1, gamma=0.5,
                        validation_size=5, arm="A_catalogue_only")
