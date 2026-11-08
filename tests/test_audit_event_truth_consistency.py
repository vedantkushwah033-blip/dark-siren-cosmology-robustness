from __future__ import annotations

import importlib.util
from pathlib import Path

import pandas as pd
import pytest


ROOT = Path(__file__).resolve().parents[1]
SPEC = importlib.util.spec_from_file_location(
    "audit_independent_validation_inclusion_h0_smoke",
    ROOT / "scripts" / "audit_independent_validation_inclusion_h0_smoke.py",
)
assert SPEC is not None and SPEC.loader is not None
AUDIT = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(AUDIT)


def _tables():
    parent = pd.DataFrame({
        "object_id": [0, 1, 2],
        "redshift": [0.05, 0.20, 0.35],
    })
    event_truth = pd.DataFrame({
        "event_id": [0, 1],
        "host_object_id": [0, 2],
        "host_redshift": [0.05, 0.35],
        "latent_mean_distance_mpc": [
            299_792.458 * 0.05 / 70.0,
            299_792.458 * 0.35 / 70.0,
        ],
    })
    return event_truth, parent


def test_event_truth_reconstruction_accepts_consistent_parent_and_h0():
    event_truth, parent = _tables()
    AUDIT._validate_event_truth_against_parent(
        event_truth, parent, h0_true=70.0
    )


def test_event_truth_reconstruction_rejects_wrong_host_redshift():
    event_truth, parent = _tables()
    event_truth.loc[0, "host_redshift"] += 0.01
    with pytest.raises(AssertionError, match="host redshifts"):
        AUDIT._validate_event_truth_against_parent(
            event_truth, parent, h0_true=70.0
        )


def test_event_truth_reconstruction_rejects_wrong_generating_h0():
    event_truth, parent = _tables()
    with pytest.raises(AssertionError, match="latent mean distances"):
        AUDIT._validate_event_truth_against_parent(
            event_truth, parent, h0_true=71.0
        )


def test_event_truth_reconstruction_rejects_unknown_parent_host():
    event_truth, parent = _tables()
    event_truth.loc[0, "host_object_id"] = 99
    with pytest.raises(AssertionError, match="unknown parent hosts"):
        AUDIT._validate_event_truth_against_parent(
            event_truth, parent, h0_true=70.0
        )
