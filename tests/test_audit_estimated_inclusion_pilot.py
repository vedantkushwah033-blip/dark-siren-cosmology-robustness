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
    assert report["n_replicate_arm_pairs_expected"] == 4
    assert report["n_replicate_arm_pairs_observed"] == 4
    assert set(report["per_arm_counts"]) == {"catalogue_only", "oracle", "estimated_100", "estimated_200"}
    assert "n_boundary_limited_successful_fits" in report


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


def test_independent_pilot_auditor_rejects_unexpected_contrast(tmp_path):
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100,),
        seed=20261033,
        h0_points=41,
        z_points=201,
    )
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    summary["contrasts"].append({
        "comparison_arm": "fabricated_arm",
        "reference_arm": "catalogue_only",
        "n_paired_replicates": 1,
        "mean_absolute_error_reduction": 0.0,
        "standard_error": None,
        "normal_95_low": None,
        "normal_95_high": None,
    })
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(AssertionError, match="unexpected stored contrasts"):
        audit_pilot_artifacts(tmp_path)



def _write_small_pilot(tmp_path, seed=20261040):
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=1,
        n_parent=1000,
        n_detected=5,
        validation_sizes=(100,),
        seed=seed,
        h0_points=41,
        z_points=201,
    )
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    return frame, summary


def test_independent_pilot_auditor_rejects_coverage_flag_inconsistent_with_endpoints(tmp_path):
    frame, summary = _write_small_pilot(tmp_path)
    successful = frame.index[frame["status"] == "success"]
    assert len(successful) > 0
    index = successful[0]
    frame.loc[index, "covered_68"] = not bool(frame.loc[index, "covered_68"])
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(AssertionError, match="covered_68 flags disagree"):
        audit_pilot_artifacts(tmp_path)


def test_independent_pilot_auditor_requires_explicit_h0_truth(tmp_path):
    _frame, summary = _write_small_pilot(tmp_path, seed=20261041)
    del summary["h0_true"]
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(AssertionError, match="summary.h0_true"):
        audit_pilot_artifacts(tmp_path)


def test_independent_pilot_auditor_rejects_missing_replicate_arm_pair(tmp_path):
    frame, summary = _write_small_pilot(tmp_path, seed=20261042)
    frame = frame.iloc[:-1]
    frame.to_csv(tmp_path / "replicates.csv", index=False)
    (tmp_path / "summary.json").write_text(json.dumps(summary), encoding="utf-8")
    with pytest.raises(AssertionError, match="incomplete replicate-arm matrix"):
        audit_pilot_artifacts(tmp_path)


def test_strict_boolean_parser_does_not_treat_string_false_as_true():
    from scripts.audit_estimated_inclusion_pilot import _strict_bool

    assert _strict_bool("False", "fixture") is False
    assert _strict_bool("True", "fixture") is True
    with pytest.raises(AssertionError, match="explicit Boolean"):
        _strict_bool("not-a-bool", "fixture")
