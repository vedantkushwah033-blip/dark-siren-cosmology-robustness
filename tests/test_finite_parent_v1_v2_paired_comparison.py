from __future__ import annotations

import json
import subprocess
import sys

import pandas as pd

from scripts.audit_finite_parent_validation_benchmark import audit


def test_paired_v1_v2_runner_uses_identical_inputs_and_audits_both_arms(tmp_path):
    output = tmp_path / "paired"
    subprocess.run(
        [
            sys.executable,
            "scripts/run_finite_parent_v1_v2_paired_comparison.py",
            "--scenarios", "1",
            "--draws", "2",
            "--n-parent", "1200",
            "--n-detected", "8",
            "--validation-sizes", "100",
            "--seed", "20261107",
            "--output-dir", str(output),
        ],
        check=True,
        capture_output=True,
        text=True,
    )

    v1_audit = audit(output / "v1")
    v2_audit = audit(output / "v2")
    assert v1_audit["status"] == "INDEPENDENT_NUMERICAL_AUDIT_PASS"
    assert v2_audit["status"] == "INDEPENDENT_NUMERICAL_AUDIT_PASS"

    summary = json.loads((output / "paired_summary.json").read_text())
    records = pd.read_csv(output / "paired_records.csv")
    assert summary["schema"] == "finite_parent_v1_v2_paired_comparison_v1"
    assert summary["status"] == "SMOKE_ONLY_NOT_CONFIRMATORY"
    assert summary["n_attempted"] == 2
    assert (
        summary["n_both_success"]
        + summary["n_v1_only_success"]
        + summary["n_v2_only_success"]
        + summary["n_both_failed"]
        == 2
    )
    assert len(records) == 2
    assert set(records["_merge"]) == {"both"}
