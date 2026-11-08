from pathlib import Path

import pandas as pd
import pytest

from dark_siren_robustness.estimated_inclusion_benchmark import (
    run_probability_validation_benchmark,
)
from scripts.audit_finite_validation_probability_benchmark import audit


def test_independent_probability_auditor_reconstructs_raw_records(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500,
        validation_sizes=(100,),
        curvature_values=(0.0, 1.2),
        degrees=(1, 2),
        replicates=2,
        seed=91,
        raw_output=raw_path,
    )
    summary.to_csv(summary_path, index=False)

    report = audit(summary_path, raw_path)
    assert report["status"] == "PASS"
    assert report["experiments_reconstructed"] == len(summary)
    raw = pd.read_csv(raw_path)
    assert raw["record_type"].eq("heldout_parent").all()
    assert not raw.duplicated([
        "scenario_seed", "replicate", "gamma", "n_validation",
        "model_degree", "object_id",
    ]).any()
    assert (raw.groupby([
        "scenario_seed", "replicate", "gamma", "n_validation",
        "model_degree",
    ]).size() == 400).all()


def test_independent_probability_auditor_rejects_tampered_summary(tmp_path: Path):
    raw_path = tmp_path / "raw_unit_records.csv"
    summary_path = tmp_path / "replicates.csv"
    summary = run_probability_validation_benchmark(
        n_parent=500,
        validation_sizes=(100,),
        curvature_values=(0.0,),
        degrees=(1,),
        replicates=1,
        seed=92,
        raw_output=raw_path,
    )
    summary.loc[0, "heldout_brier_score"] += 0.01
    summary.to_csv(summary_path, index=False)
    with pytest.raises(ValueError, match="mismatch"):
        audit(summary_path, raw_path)
