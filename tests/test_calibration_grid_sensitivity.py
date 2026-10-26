"""Tests for the independently audited PIT grid-sensitivity experiment."""
import csv
import json

import pytest

from scripts.run_calibration_grid_sensitivity import run_sensitivity


def test_grid_sensitivity_replays_and_audits_each_resolution(tmp_path):
    results = run_sensitivity(
        grid_points=[1201, 2401],
        n_replicates=100,
        seed=86420,
        output_dir=tmp_path,
    )
    assert [row["grid_points"] for row in results] == [1201, 2401]
    assert all(row["independent_audit"] == "PASS" for row in results)
    assert all(row["n_replicates"] == 100 for row in results)
    for points in (1201, 2401):
        case_dir = tmp_path / f"grid_{points}"
        assert (case_dir / "pit_replicates.csv").exists()
        assert (case_dir / "pit_summary.json").exists()
        report = json.loads((case_dir / "audit_report.json").read_text(encoding="utf-8"))
        assert report["status"] == "PASS"
        assert report["n_replicates_audited"] == 100
    with (tmp_path / "grid_sensitivity.csv").open(newline="", encoding="utf-8") as handle:
        table = list(csv.DictReader(handle))
    assert len(table) == 2
    assert [int(row["grid_points"]) for row in table] == [1201, 2401]


@pytest.mark.parametrize(
    ("grid_points", "message"),
    [
        ([2401], "at least two grid resolutions"),
        ([1201, 1201], "must be unique"),
        ([1201, 1], "at least two points"),
    ],
)
def test_grid_sensitivity_rejects_invalid_grid_sets(tmp_path, grid_points, message):
    with pytest.raises(ValueError, match=message):
        run_sensitivity(grid_points, n_replicates=100, seed=1, output_dir=tmp_path)
