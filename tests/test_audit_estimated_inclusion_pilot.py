import json

import pandas as pd
import pytest

from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot
from scripts.audit_estimated_inclusion_pilot import audit_pilot_artifacts


def test_independent_pilot_auditor_reconstructs_output(tmp_path):
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100, 200),
        seed=20261031,
        h0_points=41,
        z_points=201,
    )
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    report = audit_pilot_artifacts(tmp_path)
    assert report["status"] == "PASS"
    assert report["n_rows"] == 4


def test_independent_pilot_auditor_detects_corrupted_summary(tmp_path):
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100,),
        seed=20261032,
        h0_points=41,
        z_points=201,
    )
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    if summary["arm_summaries"]:
        summary["arm_summaries"][0]["n_attempted"] += 1
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(AssertionError, match="n_attempted"):
        audit_pilot_artifacts(tmp_path)
