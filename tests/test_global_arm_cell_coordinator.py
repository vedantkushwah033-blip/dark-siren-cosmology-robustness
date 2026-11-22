from __future__ import annotations

import hashlib
import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

from test_arm_cell_replicate_coordinator import _load_script


def _load_global(name: str):
    path = Path(__file__).resolve().parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _make_bundle(root: Path, *, seed: int, gamma: float, h0: float, event_seed: int | None = None) -> None:
    generator = _load_script("generate_shared_inputs.py")
    generator.run(SimpleNamespace(
        output_dir=root,
        seed=seed,
        event_seed=seed + 1 if event_seed is None else event_seed,
        validation_seed=seed + 2,
        validation_sample_seed=seed + 3,
        n_parent=600,
        n_detected=2,
        validation_frame=30,
        validation_sizes=[20],
        inclusion_degree=1,
        curvature_gamma=gamma,
        curvature_beta=-3.9,
        h0_true=h0,
        h0_prior_min=50.0,
        h0_prior_max=90.0,
        h0_points=41,
        scott_multiplier=0.7,
        quadrature_order=3,
        sigma_mpc=140.0,
        detection_threshold_mpc=250.0,
    ))


def _plan(tmp_path: Path) -> Path:
    bundles = []
    for replicate, seed, h0 in ((1, 25001, 68.0), (2, 25011, 72.0)):
        for gamma, label in ((-1.0, "gm1"), (1.0, "gp1")):
            root = tmp_path / f"bundle_r{replicate}_{label}"
            _make_bundle(root, seed=seed, gamma=gamma, h0=h0)
            bundles.append({"replicate": replicate, "input_dir": str(root)})
    path = tmp_path / "bundle_plan.json"
    path.write_text(json.dumps({"bundles": bundles}, indent=2) + "\n", encoding="utf-8")
    return path



def test_global_scheduler_runs_multiple_replicates_and_curvatures_with_one_checkpoint(tmp_path):
    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    args = SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only", "B_parent_oracle"],
        max_workers=2, resume=False,
    )
    first = scheduler.run(args)
    assert first["status"] == "PASS"
    assert first["bundle_count"] == 4
    assert first["replicate_count"] == 2
    assert first["curvature_count"] == 2
    assert first["planned_cells"] == 8
    assert first["status_counts"] == {
        "FAILED": 0, "NON_EVALUABLE": 0, "NOT_RUN": 0, "SUCCESS": 8,
    }
    report = auditor.audit(output)
    assert report["status"] == "PASS"
    assert report["successful_cells_reaudited"] == 8

    checkpoint = json.loads((output / "cell_lifecycle_checkpoint.json").read_text())
    assert len(checkpoint["cell_ledger"]) == 8
    assert {(row["replicate"], row["gamma"]) for row in checkpoint["planned_cells"]} == {
        (1, -1.0), (1, 1.0), (2, -1.0), (2, 1.0),
    }

    args.resume = True
    second = scheduler.run(args)
    assert second["status"] == "PASS"
    assert second["completed_this_invocation"] == 0
    assert auditor.audit(output)["status"] == "PASS"


def test_global_auditor_rejects_a_cell_bound_to_the_wrong_bundle(tmp_path):
    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    scheduler.run(SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only"], max_workers=1, resume=False,
    ))
    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text())
    checkpoint["planned_cells"][0]["input_dir"] = checkpoint["planned_cells"][1]["input_dir"]
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    with __import__("pytest").raises(AssertionError, match="immutable global plan hash"):
        auditor.audit(output)



def test_global_scheduler_rejects_unpaired_exogenous_streams_across_curvature(tmp_path):
    first, second = tmp_path / "bundle_gm1", tmp_path / "bundle_gp1"
    _make_bundle(first, seed=26001, gamma=-1.0, h0=70.0)
    _make_bundle(second, seed=26009, gamma=1.0, h0=70.0)
    bundle_plan = tmp_path / "unpaired_plan.json"
    bundle_plan.write_text(json.dumps({"bundles": [
        {"replicate": 3, "input_dir": str(first)},
        {"replicate": 3, "input_dir": str(second)},
    ]}), encoding="utf-8")
    scheduler = _load_global("run_global_arm_cells.py")
    with __import__("pytest").raises(ValueError, match="not paired across curvature"):
        scheduler.run(SimpleNamespace(
            bundle_plan=bundle_plan, output_dir=tmp_path / "run",
            arms=["A_catalogue_only"], max_workers=1, resume=False,
        ))



def test_global_scheduler_isolates_one_cell_failure_and_audits_remaining_cells(tmp_path, monkeypatch):
    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    import run_arm_cell_replicate as cell_coordinator

    real_run = cell_coordinator.subprocess.run
    calls = {"count": 0}

    def fail_first_worker(command, *args, **kwargs):
        if str(command[1]).endswith("run_arm_cell.py") and calls["count"] == 0:
            calls["count"] += 1
            return SimpleNamespace(returncode=2, stdout="", stderr="intentional test failure")
        return real_run(command, *args, **kwargs)

    monkeypatch.setattr(cell_coordinator.subprocess, "run", fail_first_worker)
    result = scheduler.run(SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only"], max_workers=1, resume=False,
    ))
    assert result["status"] == "INCOMPLETE"
    assert result["status_counts"]["FAILED"] == 1
    assert result["status_counts"]["SUCCESS"] == 3
    report = auditor.audit(output)
    assert report["status"] == "PASS"
    assert report["status_counts"]["FAILED"] == 1
    assert report["successful_cells_reaudited"] == 3



def test_global_scheduler_recovers_interrupted_cell_without_losing_first_start(tmp_path, monkeypatch):
    import pytest
    import run_arm_cell_replicate as cell_coordinator

    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    real_run = cell_coordinator.subprocess.run
    calls = {"count": 0}

    def interrupt_first_worker(command, *args, **kwargs):
        if str(command[1]).endswith("run_arm_cell.py") and calls["count"] == 0:
            calls["count"] += 1
            raise KeyboardInterrupt("simulated global scheduler interruption")
        return real_run(command, *args, **kwargs)

    monkeypatch.setattr(cell_coordinator.subprocess, "run", interrupt_first_worker)
    args = SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only"], max_workers=1, resume=False,
    )
    with pytest.raises(KeyboardInterrupt, match="simulated global scheduler interruption"):
        scheduler.run(args)

    checkpoint = json.loads((output / "cell_lifecycle_checkpoint.json").read_text())
    interrupted = next(row for row in checkpoint["cell_ledger"] if row["started_at"])
    first_started = interrupted["started_at"]
    cell_root = output / "cells" / interrupted["cell_id"]
    attempt1 = json.loads(
        (cell_root / "attempt_0001" / "attempt_started.json").read_text()
    )
    assert attempt1["first_cell_started_at"] == first_started

    monkeypatch.setattr(cell_coordinator.subprocess, "run", real_run)
    args.resume = True
    resumed = scheduler.run(args)
    assert resumed["status"] == "PASS"
    attempt2 = json.loads(
        (cell_root / "attempt_0002" / "attempt_started.json").read_text()
    )
    assert attempt2["first_cell_started_at"] == first_started
    assert attempt2["attempt_started_at"] != first_started
    assert auditor.audit(output)["status"] == "PASS"



def test_global_scheduler_rejects_reused_independent_event_seed_stream(tmp_path):
    first, second = tmp_path / "bundle_r1", tmp_path / "bundle_r2"
    _make_bundle(first, seed=27001, gamma=0.0, h0=68.0, event_seed=27999)
    _make_bundle(second, seed=27011, gamma=0.0, h0=72.0, event_seed=27999)
    bundle_plan = tmp_path / "reused_stream_plan.json"
    bundle_plan.write_text(json.dumps({"bundles": [
        {"replicate": 1, "input_dir": str(first)},
        {"replicate": 2, "input_dir": str(second)},
    ]}), encoding="utf-8")
    scheduler = _load_global("run_global_arm_cells.py")
    with __import__("pytest").raises(ValueError, match="reuse the event_seed stream"):
        scheduler.run(SimpleNamespace(
            bundle_plan=bundle_plan, output_dir=tmp_path / "run",
            arms=["A_catalogue_only"], max_workers=1, resume=False,
        ))



def test_global_auditor_rejects_seed_provenance_tampering_even_with_rehashed_plan(tmp_path):
    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    scheduler.run(SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only"], max_workers=1, resume=False,
    ))
    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text())
    planned = checkpoint["planned_cells"][0]
    key = (planned["replicate"], planned["gamma"], planned["validation_size"], planned["arm"])
    planned["event_seed"] = int(planned["event_seed"]) + 123
    for row in checkpoint["cell_ledger"]:
        if (row["replicate"], row["gamma"], row["validation_size"], row["arm"]) == key:
            row["event_seed"] = planned["event_seed"]
    payload = json.dumps(checkpoint["planned_cells"], sort_keys=True, separators=(",", ":"), allow_nan=False)
    checkpoint["plan_sha256"] = hashlib.sha256(payload.encode("utf-8")).hexdigest()
    checkpoint_path.write_text(json.dumps(checkpoint, indent=2, sort_keys=True) + "\n")
    with __import__("pytest").raises(AssertionError, match="event_seed provenance differs"):
        auditor.audit(output)


def test_global_scheduler_recovers_after_real_coordinator_process_termination(tmp_path):
    """Kill the scheduler OS process after a cell claim, then resume and audit."""
    import subprocess
    import sys
    import time

    bundle = tmp_path / "bundle"
    _make_bundle(bundle, seed=31001, gamma=0.5, h0=70.0)
    bundle_plan = tmp_path / "bundle_plan.json"
    bundle_plan.write_text(json.dumps({"bundles": [
        {"replicate": 1, "input_dir": str(bundle)},
    ]}), encoding="utf-8")
    output = tmp_path / "global_run"
    marker = tmp_path / "worker_dispatch_started"

    repo_root = Path(__file__).resolve().parents[1]
    child_code = r'''
import sys
import time
from pathlib import Path
from types import SimpleNamespace

repo_root, plan_path, output_path, marker_path = sys.argv[1:]
sys.path.insert(0, str(Path(repo_root) / "src"))
sys.path.insert(0, str(Path(repo_root) / "scripts"))
import run_global_arm_cells as scheduler
import run_arm_cell_replicate as coordinator

def block_worker_dispatch(*args, **kwargs):
    Path(marker_path).write_text("dispatch reached", encoding="utf-8")
    time.sleep(120)

coordinator.subprocess.run = block_worker_dispatch
scheduler.run(SimpleNamespace(
    bundle_plan=Path(plan_path),
    output_dir=Path(output_path),
    arms=["A_catalogue_only"],
    max_workers=1,
    resume=False,
))
'''
    process = subprocess.Popen(
        [sys.executable, "-c", child_code, str(repo_root), str(bundle_plan),
         str(output), str(marker)],
        cwd=repo_root,
    )
    try:
        deadline = time.monotonic() + 20.0
        while time.monotonic() < deadline and process.poll() is None and not marker.exists():
            time.sleep(0.05)
        assert marker.exists(), "child did not reach worker dispatch before timeout"
        process.terminate()
        process.wait(timeout=10)
        assert process.returncode != 0, "terminated coordinator unexpectedly exited successfully"
    finally:
        if process.poll() is None:
            process.kill()
            process.wait(timeout=10)

    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    checkpoint = json.loads(checkpoint_path.read_text(encoding="utf-8"))
    assert len(checkpoint["cell_ledger"]) == 1
    interrupted = checkpoint["cell_ledger"][0]
    assert interrupted["status"] == "NOT_RUN"
    assert interrupted["reason_code"] == "RUN_NOT_STARTED"
    first_started = interrupted["started_at"]
    assert first_started

    cell_root = output / "cells" / interrupted["cell_id"]
    attempt1 = json.loads(
        (cell_root / "attempt_0001" / "attempt_started.json").read_text(encoding="utf-8")
    )
    assert attempt1["first_cell_started_at"] == first_started

    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    resumed = scheduler.run(SimpleNamespace(
        bundle_plan=bundle_plan,
        output_dir=output,
        arms=["A_catalogue_only"],
        max_workers=1,
        resume=True,
    ))
    assert resumed["status"] == "PASS"
    assert resumed["status_counts"]["SUCCESS"] == 1
    attempt2 = json.loads(
        (cell_root / "attempt_0002" / "attempt_started.json").read_text(encoding="utf-8")
    )
    assert attempt2["first_cell_started_at"] == first_started
    assert attempt2["attempt_started_at"] != first_started
    report = auditor.audit(output)
    assert report["status"] == "PASS"
    assert report["successful_cells_reaudited"] == 1
