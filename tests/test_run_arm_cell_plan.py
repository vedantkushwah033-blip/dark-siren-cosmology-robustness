import json
from pathlib import Path

import pytest

from dark_siren_robustness.cell_ledger_store import (
    initialize_checkpoint,
    mark_cell_started,
    read_checkpoint,
    record_cell_outcome,
)
from scripts import run_arm_cell_plan as coordinator


def _write_plan(path: Path, *, output_suffix=""):
    cells = []
    for arm in ("A_catalogue_only", "B_parent_oracle"):
        cells.append({
            "replicate": 1, "gamma": 0.0, "validation_size": 100,
            "arm": arm, "bundle_dir": "bundle",
            "output_dir": f"cells/{arm}{output_suffix}",
        })
    path.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    return cells


def _complete(checkpoint: Path, kwargs: dict, *, status="SUCCESS"):
    key = (kwargs["replicate"], kwargs["gamma"], kwargs["validation_size"], kwargs["arm"])
    mark_cell_started(checkpoint, key, "2026-11-22T00:00:00+00:00")
    if status == "SUCCESS":
        record_cell_outcome(
            checkpoint, key, status="SUCCESS", audit_status="PASS",
            reason_code="", completed_at="2026-11-22T00:00:01+00:00",
            audit_report_path="/tmp/audit.json",
        )
        return {"status": "SUCCESS", "cell_key": list(key)}
    record_cell_outcome(
        checkpoint, key, status="FAILED", audit_status="NOT_RUN",
        reason_code="CHILD_PROCESS_ERROR", reason_detail="injected",
        completed_at="2026-11-22T00:00:01+00:00",
    )
    return {"status": "FAILED", "cell_key": list(key)}


def test_plan_initializes_all_cells_before_first_execution(tmp_path, monkeypatch):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    _write_plan(plan_path)
    observed = []

    def fake_runner(checkpoint_path, bundle, output, **kwargs):
        state = read_checkpoint(checkpoint_path)
        observed.append((len(state["planned_cells"]), [r["status"] for r in state["cell_ledger"]]))
        return _complete(Path(checkpoint_path), kwargs)

    monkeypatch.setattr(coordinator, "run_audited_cell", fake_runner)
    result = coordinator.run_plan(plan_path, checkpoint)

    assert observed == [(2, ["NOT_RUN", "NOT_RUN"]), (2, ["SUCCESS", "NOT_RUN"])]
    assert result["status"] == "PASS"
    assert result["status_counts"]["SUCCESS"] == 2
    assert result["execution_mode"] == "SERIAL"


def test_coordinator_continues_after_one_cell_failure(tmp_path, monkeypatch):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    _write_plan(plan_path)
    calls = []

    def fake_runner(checkpoint_path, bundle, output, **kwargs):
        calls.append(kwargs["arm"])
        return _complete(
            Path(checkpoint_path), kwargs,
            status="FAILED" if kwargs["arm"] == "A_catalogue_only" else "SUCCESS",
        )

    monkeypatch.setattr(coordinator, "run_audited_cell", fake_runner)
    result = coordinator.run_plan(plan_path, checkpoint)

    assert calls == ["A_catalogue_only", "B_parent_oracle"]
    assert result["status"] == "COMPLETED_WITH_NON_SUCCESS"
    assert result["ledger_validation_status"] == "PASS"
    assert result["status_counts"]["FAILED"] == 1
    assert result["status_counts"]["SUCCESS"] == 1


def test_resume_finalizes_interrupted_cells_and_skips_terminal_cells(tmp_path, monkeypatch):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    _write_plan(plan_path)
    planned = coordinator._normalize_plan(plan_path)
    initialize_checkpoint(checkpoint, planned)
    first_key = (1, 0.0, 100, "A_catalogue_only")
    mark_cell_started(checkpoint, first_key, "2026-11-22T00:00:00+00:00")
    _complete(checkpoint, {
        "replicate": 1, "gamma": 0.0, "validation_size": 100,
        "arm": "B_parent_oracle",
    })

    def fake_runner(checkpoint_path, bundle, output, **kwargs):
        state = read_checkpoint(checkpoint_path)
        row = next(r for r in state["cell_ledger"] if r["arm"] == kwargs["arm"])
        assert row["reason_code"] == "RUN_INTERRUPTED"
        assert row["started_at"] == "2026-11-22T00:00:00+00:00"
        return _complete(Path(checkpoint_path), kwargs)

    monkeypatch.setattr(coordinator, "run_audited_cell", fake_runner)
    result = coordinator.run_plan(plan_path, checkpoint, recover_interrupted=True)

    assert result["status_counts"]["SUCCESS"] == 2
    assert len(result["cell_results"]) == 2
    assert result["cell_results"][1]["skipped"] is True
    final = read_checkpoint(checkpoint)
    assert final["cell_ledger"][0]["started_at"] == "2026-11-22T00:00:00+00:00"


def test_explicit_recovery_requeues_untouched_cells_and_preserves_terminal_cells(tmp_path, monkeypatch):
    plan_path = tmp_path / "recovery-plan.json"
    checkpoint = tmp_path / "recovery-checkpoint.json"
    cells = [
        {
            "replicate": 1, "gamma": 0.0, "validation_size": 100,
            "arm": arm, "bundle_dir": "bundle", "output_dir": f"cells/{arm}",
        }
        for arm in ("A_catalogue_only", "B_parent_oracle", "E_uncertainty_integrated_ipw")
    ]
    plan_path.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    planned = coordinator._normalize_plan(plan_path)
    initialize_checkpoint(checkpoint, planned)

    interrupted_key = (1, 0.0, 100, "A_catalogue_only")
    terminal_key = (1, 0.0, 100, "B_parent_oracle")
    untouched_key = (1, 0.0, 100, "E_estimated_ipw")
    mark_cell_started(checkpoint, interrupted_key, "2026-11-22T00:00:00+00:00")
    _complete(checkpoint, {
        "replicate": 1, "gamma": 0.0, "validation_size": 100,
        "arm": "B_parent_oracle",
    })

    calls = []

    def fake_runner(checkpoint_path, bundle, output, **kwargs):
        calls.append(kwargs["arm"])
        state = read_checkpoint(checkpoint_path)
        row = next(r for r in state["cell_ledger"] if r["arm"] == kwargs["arm"])
        if kwargs["arm"] == "A_catalogue_only":
            assert row["reason_code"] == "RUN_INTERRUPTED"
            assert row["started_at"] == "2026-11-22T00:00:00+00:00"
        else:
            assert kwargs["arm"] == "E_estimated_ipw"
            assert row["reason_code"] == "RUN_NOT_STARTED"
            assert row["started_at"] is None
        return _complete(Path(checkpoint_path), kwargs)

    monkeypatch.setattr(coordinator, "run_audited_cell", fake_runner)
    result = coordinator.run_plan(plan_path, checkpoint, recover_interrupted=True)

    assert calls == ["A_catalogue_only", "E_estimated_ipw"]
    assert result["status"] == "PASS"
    assert result["status_counts"]["SUCCESS"] == 3
    assert len(result["cell_results"]) == 3
    assert result["cell_results"][1]["skipped"] is True
    final = read_checkpoint(checkpoint)
    rows = {(
        row["replicate"], row["gamma"], row["validation_size"], row["arm"]
    ): row for row in final["cell_ledger"]}
    assert rows[interrupted_key]["status"] == "SUCCESS"
    assert rows[interrupted_key]["started_at"] == "2026-11-22T00:00:00+00:00"
    assert rows[terminal_key]["status"] == "SUCCESS"
    assert rows[untouched_key]["status"] == "SUCCESS"


def test_coordinator_rejects_changed_immutable_plan(tmp_path, monkeypatch):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    _write_plan(plan_path)
    initialize_checkpoint(checkpoint, coordinator._normalize_plan(plan_path))
    _write_plan(plan_path, output_suffix="_changed")
    with pytest.raises(ValueError, match="immutable plan"):
        coordinator.run_plan(plan_path, checkpoint)


def test_plan_rejects_shared_output_destination(tmp_path):
    plan_path = tmp_path / "plan.json"
    cells = _write_plan(plan_path)
    cells[1]["output_dir"] = cells[0]["output_dir"]
    plan_path.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    with pytest.raises(ValueError, match="share an output_dir"):
        coordinator._normalize_plan(plan_path)


def test_plan_rejects_unknown_arm_before_checkpoint_creation(tmp_path):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    cells = _write_plan(plan_path)
    cells[0]["arm"] = "not_an_arm"
    plan_path.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    with pytest.raises(ValueError, match="unknown arm"):
        coordinator.run_plan(plan_path, checkpoint)
    assert not checkpoint.exists()
