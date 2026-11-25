import json
import os
from pathlib import Path
import subprocess
import sys
import time

import pytest

from dark_siren_robustness.cell_ledger_store import read_checkpoint
from scripts import run_arm_cell_plan as coordinator


@pytest.mark.skipif(os.name != "posix", reason="hard-kill test requires POSIX process semantics")
def test_real_process_kill_is_recovered_and_first_start_is_preserved(tmp_path, monkeypatch):
    plan_path = tmp_path / "plan.json"
    checkpoint = tmp_path / "checkpoint.json"
    marker = tmp_path / "child-started.marker"
    cells = [{
        "replicate": 1, "gamma": 0.0, "validation_size": 100,
        "arm": "A_catalogue_only", "bundle_dir": "bundle",
        "output_dir": "cells/A_catalogue_only",
    }]
    plan_path.write_text(json.dumps({"schema_version": 1, "cells": cells}))
    planned = coordinator._normalize_plan(plan_path)
    from dark_siren_robustness.cell_ledger_store import initialize_checkpoint
    initialize_checkpoint(checkpoint, planned)

    child_code = """
import sys, time
from pathlib import Path
from dark_siren_robustness.cell_ledger_store import mark_cell_started
checkpoint, marker = sys.argv[1], sys.argv[2]
mark_cell_started(checkpoint, (1, 0.0, 100, "A_catalogue_only"),
                  "2026-11-22T03:00:00+00:00")
Path(marker).write_text("persisted")
time.sleep(60)
"""
    process = subprocess.Popen([sys.executable, "-c", child_code, str(checkpoint), str(marker)])
    try:
        deadline = time.monotonic() + 10
        while not marker.exists() and time.monotonic() < deadline:
            if process.poll() is not None:
                raise AssertionError(f"worker exited early with code {process.returncode}")
            time.sleep(0.02)
        assert marker.is_file(), "worker did not persist its start marker in time"
    finally:
        process.kill()
        process.wait(timeout=10)

    interrupted = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert interrupted["status"] == "NOT_RUN"
    assert interrupted["reason_code"] == "RUN_NOT_STARTED"
    assert interrupted["started_at"] == "2026-11-22T03:00:00+00:00"

    def resume_runner(checkpoint_path, bundle, output, **kwargs):
        from dark_siren_robustness.cell_ledger_store import mark_cell_started, record_cell_outcome
        state = read_checkpoint(checkpoint_path)
        row = state["cell_ledger"][0]
        assert row["reason_code"] == "RUN_INTERRUPTED"
        assert row["started_at"] == "2026-11-22T03:00:00+00:00"
        key = (1, 0.0, 100, "A_catalogue_only")
        mark_cell_started(checkpoint_path, key, "2026-11-22T03:00:10+00:00")
        record_cell_outcome(
            checkpoint_path, key, status="SUCCESS", audit_status="PASS",
            reason_code="", reason_detail="hard-kill recovery test",
            completed_at="2026-11-22T03:00:11+00:00",
            audit_report_path=str(tmp_path / "synthetic-audit-report.json"),
        )
        return {"status": "SUCCESS", "cell_key": list(key)}

    monkeypatch.setattr(coordinator, "run_audited_cell", resume_runner)
    result = coordinator.run_plan(plan_path, checkpoint, recover_interrupted=True)
    final = read_checkpoint(checkpoint)["cell_ledger"][0]
    assert result["status"] == "PASS"
    assert final["status"] == "SUCCESS"
    assert final["started_at"] == "2026-11-22T03:00:00+00:00"
