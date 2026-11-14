from __future__ import annotations

import importlib.util
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
SCRIPTS = ROOT / "scripts"
sys.path.insert(0, str(SCRIPTS))


def _load(name: str, filename: str):
    spec = importlib.util.spec_from_file_location(name, SCRIPTS / filename)
    module = importlib.util.module_from_spec(spec)
    sys.modules[name] = module
    assert spec.loader is not None
    spec.loader.exec_module(module)
    return module


def test_all_gamma_zero_reference_cells_failed_remains_auditable(tmp_path, monkeypatch):
    from argparse import Namespace

    runner = _load(
        "curvature_runner_missing_reference",
        "run_quadratic_logit_curvature_density_replicates.py",
    )
    auditor = _load(
        "curvature_auditor_missing_reference",
        "audit_quadratic_logit_curvature_density_replicates.py",
    )
    original = runner.run_one

    def fail_reference_cells(args):
        if float(args.gamma) == 0.0:
            raise RuntimeError("injected failure of gamma-zero reference cell")
        return original(args)

    monkeypatch.setattr(runner, "run_one", fail_reference_cells)
    output = tmp_path / "curvature_missing_reference"
    args = Namespace(
        output_dir=output,
        master_seed=20261220,
        seed_stride=10,
        replicates=2,
        gammas=(-1.0, 0.0, 1.0),
        n_parent=350,
        validation_frame=120,
        validation_sizes=(40, 80),
        beta=-3.9,
        quadrature_order=3,
        grid_points=101,
        scott_multiplier=0.7,
        absolute_widths=(0.012, 0.024),
    )

    result = runner.run(args)
    assert result["status"] == "COMPLETED_WITH_FAILURES"
    assert result["failed_cells"] == 2
    report = auditor.audit(output)
    assert report["status"] == "INCOMPLETE"
    assert report["failed_cells"] == 2
    assert report["successful_cells"] == 4
    assert report["paired_contrast_rows_reconstructed"] == 0
    assert report["paired_contrast_summary_rows_reconstructed"] == 0
    import pandas as pd
    assert pd.read_csv(output / "paired_contrasts_all_cells.csv").empty
    assert pd.read_csv(output / "paired_contrast_summary.csv").empty
