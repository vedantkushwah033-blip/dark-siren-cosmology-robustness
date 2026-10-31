from __future__ import annotations

import json
from pathlib import Path

import pandas as pd

from scripts import run_estimated_inclusion_curvature_screen as screen


def test_curvature_screen_runs_all_conditions_and_audits(tmp_path, monkeypatch):
    output = tmp_path / "curvature_screen"
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_estimated_inclusion_curvature_screen.py",
            "--replicates", "1",
            "--n-parent", "1000",
            "--n-detected", "5",
            "--validation-sizes", "100,300,500",
            "--seed", "20261031",
            "--h0-points", "31",
            "--z-points", "101",
            "--output-dir", str(output),
        ],
    )

    screen.main()

    report = json.loads((output / "screen_report.json").read_text(encoding="utf-8"))
    assert report["status"] == "SCREEN_NOT_CONFIRMATORY"
    assert report["curvature_levels"] == [-2.0, -1.0, 0.0, 1.0, 2.0]
    assert report["common_random_numbers_across_curvature"] is True
    assert len(report["conditions"]) == 5
    assert all(item["audit_status"] == "PASS" for item in report["conditions"])

    combined = pd.read_csv(output / "all_conditions_replicates.csv")
    assert len(combined) == 25
    assert set(combined["true_quadratic"]) == {-2.0, -1.0, 0.0, 1.0, 2.0}
    for item in report["conditions"]:
        condition_dir = output / item["condition"]
        assert (condition_dir / "replicates.csv").is_file()
        assert (condition_dir / "summary.json").is_file()
        assert (condition_dir / "independent_audit.json").is_file()
