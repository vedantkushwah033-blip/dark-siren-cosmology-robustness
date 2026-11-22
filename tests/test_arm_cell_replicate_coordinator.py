from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

import pandas as pd
import pytest

from test_shared_input_generator import GENERATOR


def _load_script(name: str):
    path = Path(__file__).resolve().parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _bundle(root: Path) -> None:
    generator = _load_script("generate_shared_inputs.py")
    generator.run(SimpleNamespace(
        output_dir=root,
        seed=15001,
        event_seed=15002,
        validation_seed=15003,
        validation_sample_seed=15004,
        n_parent=600,
        n_detected=2,
        validation_frame=30,
        validation_sizes=[20],
        inclusion_degree=1,
        curvature_gamma=0.0,
        curvature_beta=-3.9,
        h0_true=70.0,
        h0_prior_min=50.0,
        h0_prior_max=90.0,
        h0_points=41,
        scott_multiplier=0.7,
        quadrature_order=3,
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
    ))


def test_coordinator_executes_and_independently_audits_cells_then_resumes_idempotently(tmp_path):
    bundle = tmp_path / "shared"
    output = tmp_path / "run"
    _bundle(bundle)
    coordinator = _load_script("run_arm_cell_replicate.py")
    aggregate_auditor = _load_script("audit_arm_cell_replicate.py")
    args = SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=4,
        arms=["A_catalogue_only", "B_parent_oracle"],
        max_workers=1, resume=False,
    )
    first = coordinator.run(args)
    assert first["status"] == "PASS"
    assert first["planned_cells"] == 2
    assert first["status_counts"] == {
        "FAILED": 0, "NON_EVALUABLE": 0, "NOT_RUN": 0, "SUCCESS": 2,
    }
    report = aggregate_auditor.audit(output, bundle)
    assert report["status"] == "PASS"
    assert report["successful_cells_reaudited"] == 2

    args.resume = True
    second = coordinator.run(args)
    assert second["status"] == "PASS"
    assert second["completed_this_invocation"] == 0
    assert second["successful_this_invocation"] == 0
    report_after_resume = aggregate_auditor.audit(output, bundle)
    assert report_after_resume["status"] == "PASS"
    checkpoint = json.loads((output / "cell_lifecycle_checkpoint.json").read_text())
    assert len(checkpoint["cell_ledger"]) == 2
    assert all(row["status"] == "SUCCESS" for row in checkpoint["cell_ledger"])



def _successful_two_cell_run(tmp_path: Path):
    bundle = tmp_path / "shared"
    output = tmp_path / "run"
    _bundle(bundle)
    coordinator = _load_script("run_arm_cell_replicate.py")
    aggregate_auditor = _load_script("audit_arm_cell_replicate.py")
    result = coordinator.run(SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=4,
        arms=["A_catalogue_only", "B_parent_oracle"],
        max_workers=1, resume=False,
    ))
    assert result["status"] == "PASS"
    assert aggregate_auditor.audit(output, bundle)["status"] == "PASS"
    return output, bundle, aggregate_auditor


def test_aggregate_auditor_rejects_edited_ledger_csv_values_even_if_row_count_is_unchanged(tmp_path):
    output, bundle, auditor = _successful_two_cell_run(tmp_path)
    ledger_path = output / "cell_ledger.csv"
    text = ledger_path.read_text(encoding="utf-8")
    assert "SUCCESS" in text
    ledger_path.write_text(text.replace("SUCCESS", "FAILED", 1), encoding="utf-8")
    with pytest.raises(AssertionError, match="CSV ledger export differs"):
        auditor.audit(output, bundle)


def test_aggregate_auditor_rejects_edited_summary_csv_values_even_if_row_count_is_unchanged(tmp_path):
    output, bundle, auditor = _successful_two_cell_run(tmp_path)
    summary_path = output / "cell_summaries.csv"
    summaries = pd.read_csv(summary_path)
    assert len(summaries) == 2
    summaries.loc[0, "posterior_median"] += 1.0
    summaries.to_csv(summary_path, index=False)
    with pytest.raises(AssertionError, match="CSV cell summaries differ"):
        auditor.audit(output, bundle)



def test_interrupted_attempt_resume_retains_immutable_first_start_timestamp(tmp_path, monkeypatch):
    bundle = tmp_path / "shared"
    output = tmp_path / "run"
    _bundle(bundle)
    coordinator = _load_script("run_arm_cell_replicate.py")
    aggregate_auditor = _load_script("audit_arm_cell_replicate.py")
    args = SimpleNamespace(
        input_dir=bundle, output_dir=output, replicate=9,
        arms=["A_catalogue_only"], max_workers=1, resume=False,
    )

    real_run = coordinator.subprocess.run

    def interrupt_before_worker(*_args, **_kwargs):
        raise KeyboardInterrupt("simulated coordinator interruption")

    monkeypatch.setattr(coordinator.subprocess, "run", interrupt_before_worker)
    with pytest.raises(KeyboardInterrupt, match="simulated coordinator interruption"):
        coordinator.run(args)

    cell_root = next((output / "cells").iterdir())
    first_attempt = json.loads(
        (cell_root / "attempt_0001" / "attempt_started.json").read_text(encoding="utf-8")
    )
    checkpoint_before_resume = json.loads(
        (output / "cell_lifecycle_checkpoint.json").read_text(encoding="utf-8")
    )
    original_first_start = checkpoint_before_resume["cell_ledger"][0]["started_at"]
    assert original_first_start == first_attempt["first_cell_started_at"]

    monkeypatch.setattr(coordinator.subprocess, "run", real_run)
    args.resume = True
    resumed = coordinator.run(args)
    assert resumed["status"] == "PASS"

    second_attempt = json.loads(
        (cell_root / "attempt_0002" / "attempt_started.json").read_text(encoding="utf-8")
    )
    assert second_attempt["attempt_started_at"] != original_first_start
    assert second_attempt["first_cell_started_at"] == original_first_start
    assert aggregate_auditor.audit(output, bundle)["status"] == "PASS"



def test_global_multi_bundle_scheduler_regression_suite(tmp_path, monkeypatch):
    # Keep this integration suite in the focused worker workflow's existing
    # test target so the stacked PR receives the complete Python matrix.
    from test_global_arm_cell_coordinator import (
        test_global_auditor_rejects_a_cell_bound_to_the_wrong_bundle,
        test_global_scheduler_isolates_one_cell_failure_and_audits_remaining_cells,
        test_global_scheduler_recovers_interrupted_cell_without_losing_first_start,
        test_global_scheduler_rejects_unpaired_exogenous_streams_across_curvature,
        test_global_scheduler_runs_multiple_replicates_and_curvatures_with_one_checkpoint,
    )

    cases = [
        test_global_scheduler_runs_multiple_replicates_and_curvatures_with_one_checkpoint,
        test_global_auditor_rejects_a_cell_bound_to_the_wrong_bundle,
        test_global_scheduler_rejects_unpaired_exogenous_streams_across_curvature,
        test_global_scheduler_isolates_one_cell_failure_and_audits_remaining_cells,
        test_global_scheduler_recovers_interrupted_cell_without_losing_first_start,
    ]
    for index, case in enumerate(cases):
        case_root = tmp_path / f"global_case_{index}"
        case_root.mkdir()
        if "monkeypatch" in case.__code__.co_varnames:
            case(case_root, monkeypatch)
        else:
            case(case_root)
