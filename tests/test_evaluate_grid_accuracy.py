import csv
import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "evaluate_grid_accuracy.py"
SPEC = importlib.util.spec_from_file_location("evaluate_grid_accuracy", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def write_fixture(path, coarse_upper_shift=0.02):
    fields = [
        "H0_true", "z2", "noise_fraction", "replicate", "inference_weighting",
        "grid_size", "observed_distance_mpc", "posterior_area_on_prior_support",
        "posterior_median", "lower_68", "upper_68", "lower_95", "upper_95",
    ]
    rows = []
    for rep in range(3):
        for grid, delta in ((1001, coarse_upper_shift), (2001, coarse_upper_shift / 4), (4001, 0.0)):
            rows.append({
                "H0_true": 70.0, "z2": 0.04, "noise_fraction": 0.03,
                "replicate": rep, "inference_weighting": "correct",
                "grid_size": grid, "observed_distance_mpc": 100.0 + rep,
                "posterior_area_on_prior_support": 1.0,
                "posterior_median": 70.0 + delta / 10,
                "lower_68": 65.0 - delta / 10, "upper_68": 75.0 + delta / 10,
                "lower_95": 60.0 - delta / 10, "upper_95": 80.0 + delta,
            })
    with path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=fields)
        writer.writeheader()
        writer.writerows(rows)


def test_accuracy_evaluator_flags_outliers_and_writes_all_case_summaries(tmp_path):
    source = tmp_path / "per_run.csv"
    write_fixture(source, coarse_upper_shift=0.08)
    assessment, cells, cases = module.evaluate(source)
    assert assessment["matched_cases"] == 3
    assert assessment["design_cells"] == 1
    assert assessment["status"] == "FAIL"
    assert assessment["failed_cell_summary_tests"] >= 1
    upper95 = [r for r in cases if r["summary"] == "upper_95"]
    assert len(upper95) == 3
    assert all(r["flag_for_targeted_refinement"] for r in upper95)
    assert all(r["targeted_refinement_status"] == "required" for r in upper95)


def test_accuracy_evaluator_passes_small_grid_differences(tmp_path):
    source = tmp_path / "per_run.csv"
    write_fixture(source, coarse_upper_shift=0.001)
    assessment, cells, cases = module.evaluate(source)
    assert assessment["status"] == "PASS_FOR_GRID_DIFFERENCES_ONLY"
    assert assessment["failed_cell_summary_tests"] == 0
    assert assessment["flagged_case_summary_pairs"] == 0


def test_accuracy_evaluator_rejects_missing_resolution(tmp_path):
    source = tmp_path / "per_run.csv"
    write_fixture(source)
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    rows = [row for row in rows if not (row["replicate"] == "1" and row["grid_size"] == "2001")]
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    try:
        module.evaluate(source)
    except ValueError as exc:
        assert "Incomplete resolution triplet" in str(exc)
    else:
        raise AssertionError("Expected incomplete grid triplet to be rejected")



def test_accuracy_evaluator_fails_if_independent_area_is_not_normalized(tmp_path):
    source = tmp_path / "per_run.csv"
    write_fixture(source, coarse_upper_shift=0.001)
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    for row in rows:
        row["posterior_area_on_prior_support"] = "1.0001"
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    assessment, cells, cases = module.evaluate(source)
    assert assessment["normalization_checked_rows"] == 9
    assert assessment["normalization_failures"] == 9
    assert assessment["status"] == "FAIL"



def test_accuracy_evaluator_rejects_missing_normalization_diagnostic(tmp_path):
    source = tmp_path / "per_run.csv"
    write_fixture(source, coarse_upper_shift=0.001)
    rows = list(csv.DictReader(source.open(encoding="utf-8")))
    for row in rows:
        row.pop("posterior_area_on_prior_support")
    with source.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(rows[0]))
        writer.writeheader()
        writer.writerows(rows)
    try:
        module.evaluate(source)
    except ValueError as exc:
        assert "Missing posterior_area_on_prior_support" in str(exc)
    else:
        raise AssertionError("Expected missing normalization diagnostics to be rejected")
