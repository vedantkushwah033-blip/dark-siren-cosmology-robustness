from __future__ import annotations

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


def _make_bundle(root: Path, *, seed: int, gamma: float, h0: float) -> None:
    generator = _load_script("generate_shared_inputs.py")
    generator.run(SimpleNamespace(
        output_dir=root,
        seed=seed,
        event_seed=seed + 1,
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
