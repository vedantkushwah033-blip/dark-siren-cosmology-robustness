from pathlib import Path

import pytest

from dark_siren_robustness.cell_ledger_store import (
    finalize_interrupted_cells,
    initialize_checkpoint,
    mark_cell_started,
    read_checkpoint,
)
from scripts import run_audited_arm_cell as runner


KEY = (1, 0.0, 100, "A_catalogue_only")


def _checkpoint(path: Path):
    initialize_checkpoint(path, [{
        "replicate": KEY[0], "gamma": KEY[1],
        "validation_size": KEY[2], "arm": KEY[3],
    }])


def _fake_run(_bundle, output, **_kwargs):
    Path(output).mkdir(parents=True)


def _pass_audit(_bundle, _output):
    return {"status": "PASS", "audit_status": "PASS", "cell_key": {
        "replicate": 1, "gamma": 0.0, "validation_size": 100,
        "arm": "A_catalogue_only",
    }}


def test_audited_cell_records_success_only_after_independent_pass(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    output = tmp_path / "cell"
    audit_path = tmp_path / "cell.audit.json"
    _checkpoint(checkpoint)
    monkeypatch.setattr(runner, "run_cell", _fake_run)
    monkeypatch.setattr(runner, "audit_cell", _pass_audit)

    result = runner.run_audited_cell(
        checkpoint, tmp_path / "bundle", output,
        replicate=1, gamma=0.0, validation_size=100,
        arm="A_catalogue_only", audit_report_path=audit_path,
    )

    assert result["status"] == "SUCCESS"
    assert audit_path.is_file()
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["status"] == "SUCCESS"
    assert row["audit_status"] == "PASS"
    assert row["audit_report_path"] == str(audit_path.resolve())


def test_execution_exception_records_failed_not_non_evaluable(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    _checkpoint(checkpoint)

    def fail_run(*_args, **_kwargs):
        raise RuntimeError("injected evaluator failure")

    monkeypatch.setattr(runner, "run_cell", fail_run)
    result = runner.run_audited_cell(
        checkpoint, tmp_path / "bundle", tmp_path / "cell",
        replicate=1, gamma=0.0, validation_size=100,
        arm="A_catalogue_only",
    )

    assert result["status"] == "FAILED"
    assert result["reason_code"] == "CHILD_PROCESS_ERROR"
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["status"] == "FAILED"
    assert row["audit_status"] == "NOT_RUN"


def test_audit_exception_records_child_audit_failed(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    _checkpoint(checkpoint)
    monkeypatch.setattr(runner, "run_cell", _fake_run)

    def fail_audit(*_args, **_kwargs):
        raise ValueError("injected independent-audit mismatch")

    monkeypatch.setattr(runner, "audit_cell", fail_audit)
    result = runner.run_audited_cell(
        checkpoint, tmp_path / "bundle", tmp_path / "cell",
        replicate=1, gamma=0.0, validation_size=100,
        arm="A_catalogue_only",
    )

    assert result["status"] == "FAILED"
    assert result["reason_code"] == "CHILD_AUDIT_FAILED"
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["status"] == "FAILED"
    assert row["audit_status"] == "NOT_PASS"


def test_active_cell_requires_explicit_interruption_recovery(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    _checkpoint(checkpoint)
    mark_cell_started(checkpoint, KEY, "2026-11-22T00:00:00+00:00")
    monkeypatch.setattr(runner, "run_cell", _fake_run)

    with pytest.raises(RuntimeError, match="stop all workers"):
        runner.run_audited_cell(
            checkpoint, tmp_path / "bundle", tmp_path / "cell",
            replicate=1, gamma=0.0, validation_size=100,
            arm="A_catalogue_only",
        )
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["status"] == "NOT_RUN"
    assert row["started_at"] == "2026-11-22T00:00:00+00:00"


def test_interrupted_cell_restarts_and_preserves_first_start_timestamp(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    output = tmp_path / "cell"
    _checkpoint(checkpoint)
    first_start = "2026-11-22T00:00:00+00:00"
    mark_cell_started(checkpoint, KEY, first_start)
    finalize_interrupted_cells(checkpoint)

    calls = []
    def resumed_run(_bundle, destination, **_kwargs):
        calls.append("run")
        Path(destination).mkdir(parents=True)

    monkeypatch.setattr(runner, "run_cell", resumed_run)
    monkeypatch.setattr(runner, "audit_cell", _pass_audit)
    result = runner.run_audited_cell(
        checkpoint, tmp_path / "bundle", output,
        replicate=1, gamma=0.0, validation_size=100,
        arm="A_catalogue_only",
    )

    assert result["status"] == "SUCCESS"
    assert calls == ["run"]
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["started_at"] == first_start
    assert row["status"] == "SUCCESS"


def test_published_artifact_is_reaudited_after_interrupted_checkpoint(tmp_path, monkeypatch):
    checkpoint = tmp_path / "checkpoint.json"
    output = tmp_path / "cell"
    _checkpoint(checkpoint)
    first_start = "2026-11-22T00:00:00+00:00"
    mark_cell_started(checkpoint, KEY, first_start)
    output.mkdir()
    finalize_interrupted_cells(checkpoint)

    def must_not_run(*_args, **_kwargs):
        raise AssertionError("published cell artifact should be re-audited, not recomputed")

    monkeypatch.setattr(runner, "run_cell", must_not_run)
    monkeypatch.setattr(runner, "audit_cell", _pass_audit)
    result = runner.run_audited_cell(
        checkpoint, tmp_path / "bundle", output,
        replicate=1, gamma=0.0, validation_size=100,
        arm="A_catalogue_only",
    )

    assert result["status"] == "SUCCESS"
    row = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert row["started_at"] == first_start
    assert row["status"] == "SUCCESS"
