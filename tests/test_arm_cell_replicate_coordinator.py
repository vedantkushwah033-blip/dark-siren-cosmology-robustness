from __future__ import annotations

import importlib.util
import json
from pathlib import Path
from types import SimpleNamespace

from test_shared_input_generator import GENERATOR


def _load_script(name: str):
    path = Path(__file__).resolve().parents[1] / "scripts" / name
    spec = importlib.util.spec_from_file_location(path.stem, path)
    module = importlib.util.module_from_spec(spec)
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def _bundle(root: Path) -> None:
    generator = _load_script("generate_shared_inputs.py")\n    generator.run(SimpleNamespace(
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
