from __future__ import annotations

import importlib.util
import json
from pathlib import Path

import numpy as np
import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "audit_inclusion_posterior_calibration",
    ROOT / "scripts" / "audit_inclusion_posterior_calibration.py",
)
AUDIT = importlib.util.module_from_spec(SPEC)
assert SPEC.loader is not None
SPEC.loader.exec_module(AUDIT)


def _write_pair(root: Path, n: int, ids: list[int]) -> Path:
    path = root / f"n{n}"
    (path / "observable").mkdir(parents=True)
    (path / "truth").mkdir()
    config = {
        "analysis_generator_config": {"seed": 10, "n_parent": 20},
        "validation_generator_config": {"seed": 11, "n_validation": n},
    }
    (path / "config.json").write_text(json.dumps(config))
    pd.DataFrame({"object_id": [0, 1], "redshift": [0.1, 0.2]}).to_csv(
        path / "observable" / "analysis_catalogue.csv", index=False
    )
    pd.DataFrame({"event_id": [0], "distance_mpc": [300.0]}).to_csv(
        path / "observable" / "gw_events_observed.csv", index=False
    )
    frame = pd.DataFrame({
        "frame_id": [0, 1, 2],
        "redshift": [0.1, 0.2, 0.3],
        "inclusion_probability": [0.8, 0.6, 0.4],
        "retained": [1, 0, 1],
        "selected_for_validation": [int(i in ids) for i in range(3)],
    })
    frame.to_csv(path / "truth" / "validation_frame_truth.csv", index=False)
    frame.set_index("frame_id").loc[ids].reset_index().rename(
        columns={"frame_id": "validation_id"}
    )[["validation_id", "redshift", "retained"]].to_csv(
        path / "observable" / "validation_observed.csv", index=False
    )
    return path


def test_independent_pit_reconstruction_matches_trapezoid_cdf():
    frame = pd.DataFrame({"H0": [0.0, 1.0, 2.0], "density": [0.0, 1.0, 0.0]})
    assert AUDIT._pit(frame, 1.0, "density") == pytest.approx(0.5)


def test_nested_validation_audit_checks_prefixes_and_fixed_analysis(tmp_path):
    first = _write_pair(tmp_path, 2, [2, 0])
    second = _write_pair(tmp_path, 3, [2, 0, 1])
    result = AUDIT._assert_nested_validation([first, second])
    assert result["paired_sizes_checked"] == 2
    assert result["analysis_files_checked"] == 2


def test_nested_validation_audit_rejects_changed_analysis_data(tmp_path):
    first = _write_pair(tmp_path, 2, [2, 0])
    second = _write_pair(tmp_path, 3, [2, 0, 1])
    pd.DataFrame({"object_id": [0], "redshift": [0.9]}).to_csv(
        second / "observable" / "analysis_catalogue.csv", index=False
    )
    with pytest.raises(AssertionError, match="analysis_catalogue"):
        AUDIT._assert_nested_validation([first, second])


def test_nested_validation_audit_rejects_non_nested_validation_ids(tmp_path):
    first = _write_pair(tmp_path, 2, [2, 0])
    second = _write_pair(tmp_path, 3, [0, 2, 1])
    with pytest.raises(AssertionError, match="nested permutation prefixes"):
        AUDIT._assert_nested_validation([first, second])
