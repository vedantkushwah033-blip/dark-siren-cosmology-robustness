import csv
import importlib.util
from pathlib import Path

MODULE_PATH = Path(__file__).resolve().parents[1] / "scripts" / "evaluate_grid_accuracy.py"
SPEC = importlib.util.spec_from_file_location("evaluate_grid_accuracy", MODULE_PATH)
module = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(module)


def test_targeted_refinement_resolves_flagged_pair_but_does_not_hide_cell_failure(tmp_path):
    per_run = tmp_path / "per_run.csv"
    fields = [
        "H0_true", "z2", "noise_fraction", "replicate", "inference_weighting",
        "grid_size", "observed_distance_mpc", "posterior_area_on_prior_support", "posterior_median", "lower_68",
        "upper_68", "lower_95", "upper_95",
    ]
    with per_run.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for rep in range(2):
            for grid, delta in ((1001, 0.08), (2001, 0.02), (4001, 0.0)):
                writer.writerow({
                    "H0_true": 70, "z2": 0.04, "noise_fraction": 0.03,
                    "replicate": rep, "inference_weighting": "correct",
                    "grid_size": grid, "observed_distance_mpc": 500 + rep,
                    "posterior_median": 70, "lower_68": 65, "upper_68": 75,
                    "lower_95": 60, "upper_95": 80 + delta,
                })
    targeted = tmp_path / "targeted.csv"
    target_fields = list(module.IDENTITY) + [
        "grid_size", "observed_distance_mpc", "posterior_median",
        "lower_68", "upper_68", "lower_95", "upper_95",
    ]
    with targeted.open("w", newline="", encoding="utf-8") as f:
        writer = csv.DictWriter(f, fieldnames=target_fields)
        writer.writeheader()
        for rep in range(2):
            for grid, delta in ((8001, 0.003), (16001, 0.0015), (32001, 0.0005), (64001, 0.0)):
                writer.writerow({
                    "H0_true": 70, "z2": 0.04, "noise_fraction": 0.03,
                    "replicate": rep, "inference_weighting": "correct",
                    "grid_size": grid, "observed_distance_mpc": 500 + rep,
                    "posterior_median": 70, "lower_68": 65, "upper_68": 75,
                    "lower_95": 60, "upper_95": 80 + delta,
                })
    assessment, cells, cases = module.evaluate(per_run, targeted)
    flagged = [row for row in cases if row["flag_for_targeted_refinement"]]
    assert len(flagged) == 2
    assert all(row["targeted_refinement_status"] == "resolved_within_tolerance" for row in flagged)
    assert assessment["targeted_refinement_resolved_pairs"] == 2
    assert assessment["targeted_refinement_unresolved_pairs"] == 0
    assert assessment["status"] == "FAIL"  # cell p95 still fails; targeted results cannot mask it
