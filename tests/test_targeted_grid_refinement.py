import csv
import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "run_targeted_grid_refinement.py"
SPEC = importlib.util.spec_from_file_location("run_targeted_grid_refinement", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_targeted_refinement_uses_flagged_cases_and_all_four_grids(tmp_path, monkeypatch):
    flagged_path = tmp_path / "flagged.csv"
    fields = list(module.IDENTITY) + ["summary", "flag_for_targeted_refinement"]
    key = ["70.0", "0.04", "0.03", "2", "correct"]
    with flagged_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        writer.writerow(dict(zip(module.IDENTITY, key), summary="upper_95",
                             flag_for_targeted_refinement="True"))
        writer.writerow(dict(zip(module.IDENTITY, key), summary="lower_95",
                             flag_for_targeted_refinement="False"))

    per_run_path = tmp_path / "per_run.csv"
    per_fields = list(module.IDENTITY) + ["grid_size", "observed_distance_mpc"]
    with per_run_path.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=per_fields)
        writer.writeheader()
        writer.writerow(dict(zip(module.IDENTITY, key), grid_size=4001, observed_distance_mpc=500.0))

    monkeypatch.setattr(module, "posterior_h0", lambda grid, *args, **kwargs: grid)
    monkeypatch.setattr(module, "posterior_summary", lambda grid, density: {
        "median": 70.0, "lower_68": 65.0, "upper_68": 75.0,
        "lower_95": 60.0, "upper_95": 80.0,
    })
    rows = module.run(flagged_path, per_run_path)
    assert len(rows) == 4
    assert {row["grid_size"] for row in rows} == {8001, 16001, 32001, 64001}
    assert all(row["observed_distance_mpc"] == 500.0 for row in rows)
    assert all(row["inference_weighting"] == "correct" for row in rows)
