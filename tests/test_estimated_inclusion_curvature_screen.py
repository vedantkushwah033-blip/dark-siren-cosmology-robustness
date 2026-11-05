from __future__ import annotations

import json
from pathlib import Path

import pandas as pd
import pytest

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
    assert len(combined) == 40  # 5 curvatures x (catalogue + oracle + 6 estimated arms)
    assert set(combined["true_quadratic"]) == {-2.0, -1.0, 0.0, 1.0, 2.0}
    expected_arms = {
        "catalogue_only", "oracle",
        "estimated_linear_100", "estimated_linear_300", "estimated_linear_500",
        "estimated_quadratic_100", "estimated_quadratic_300", "estimated_quadratic_500",
    }
    assert set(combined["arm"]) == expected_arms
    for item in report["conditions"]:
        condition_dir = output / item["condition"]
        assert (condition_dir / "replicates.csv").is_file()
        assert (condition_dir / "summary.json").is_file()
        assert (condition_dir / "independent_audit.json").is_file()


def test_curvature_screen_refuses_to_overwrite_existing_artifacts(tmp_path, monkeypatch):
    output = tmp_path / "existing"
    output.mkdir()
    sentinel = output / "keep.txt"
    sentinel.write_text("preserve prior artifact", encoding="utf-8")
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_estimated_inclusion_curvature_screen.py",
            "--replicates", "1",
            "--n-parent", "1000",
            "--n-detected", "5",
            "--validation-sizes", "100",
            "--output-dir", str(output),
        ],
    )
    with pytest.raises(SystemExit):
        screen.main()
    assert sentinel.read_text(encoding="utf-8") == "preserve prior artifact"



@pytest.mark.parametrize(
    "curvature_levels",
    ["nan,1", "inf,1", "0.01,0.02"],
)
def test_curvature_screen_rejects_nonfinite_or_colliding_levels(
    tmp_path, monkeypatch, curvature_levels
):
    output = tmp_path / "invalid_curvature_screen"
    monkeypatch.setattr(
        "sys.argv",
        [
            "run_estimated_inclusion_curvature_screen.py",
            "--replicates", "1",
            "--n-parent", "1000",
            "--n-detected", "5",
            "--validation-sizes", "100",
            "--curvature-levels", curvature_levels,
            "--output-dir", str(output),
        ],
    )

    with pytest.raises(SystemExit):
        screen.main()

    # Validation must fail before creating an artifact directory.
    assert not output.exists()
