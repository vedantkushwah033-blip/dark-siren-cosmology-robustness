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
    first, second = tmp_path / "bundle_r1_gm1", tmp_path / "bundle_r2_gp1"
    _make_bundle(first, seed=25001, gamma=-1.0, h0=68.0)
    _make_bundle(second, seed=25011, gamma=1.0, h0=72.0)
    path = tmp_path / "bundle_plan.json"
    path.write_text(json.dumps({"bundles": [
        {"replicate": 1, "input_dir": str(first)},
        {"replicate": 2, "input_dir": str(second)},
    ]}, indent=2) + "\n", encoding="utf-8")
    return path


def test_global_scheduler_runs_multiple_replicates_and_curvatures_with_one_checkpoint(tmp_path):
    bundle_plan = _plan(tmp_path)
    output = tmp_path / "global_run"
    scheduler = _load_global("run_global_arm_cells.py")
    auditor = _load_global("audit_global_arm_cells.py")
    args = SimpleNamespace(
        bundle_plan=bundle_plan, output_dir=output,
        arms=["A_catalogue_only", "B_parent_oracle"],
        max_workers=1, resume=False,
    )
    first = scheduler.run(args)
    assert first["status"] == "PASS"
    assert first["bundle_count"] == 2
    assert first["replicate_count"] == 2
    assert first["curvature_count"] == 2
    assert first["planned_cells"] == 4
    assert first["status_counts"] == {
        "FAILED": 0, "NON_EVALUABLE": 0, "NOT_RUN": 0, "SUCCESS": 4,
    }
    report = auditor.audit(output)
    assert report["status"] == "PASS"
    assert report["successful_cells_reaudited"] == 4

    checkpoint = json.loads((output / "cell_lifecycle_checkpoint.json").read_text())
    assert len(checkpoint["cell_ledger"]) == 4
    assert {(row["replicate"], row["gamma"]) for row in checkpoint["planned_cells"]} == {
        (1, -1.0), (2, 1.0),
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
