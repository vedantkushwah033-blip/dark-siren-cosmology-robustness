import csv
import json
import subprocess
import sys
from pathlib import Path

from dark_siren_robustness.uncertainty_sensitivity import run_experiment


def test_stratified_paired_summary_checks_pairs_and_writes_report(tmp_path):
    raw, _ = run_experiment(
        replicates=2,
        seed=123,
        grid_points=101,
        truths=(55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0),
        fractional_errors=(0.1,),
    )
    csv_path = tmp_path / "replicates.csv"
    with csv_path.open("w", newline="", encoding="utf-8") as handle:
        writer = csv.DictWriter(handle, fieldnames=list(raw[0]))
        writer.writeheader()
        writer.writerows(raw)

    script = Path(__file__).resolve().parents[1] / "scripts" / "summarize_distance_uncertainty_paired_contrasts.py"
    output = tmp_path / "paired_summary.json"
    subprocess.run(
        [sys.executable, str(script), str(csv_path), "--seed", "11",
         "--bootstrap-resamples", "100", "--output", str(output)],
        check=True,
        capture_output=True,
        text=True,
    )
    report = json.loads(output.read_text(encoding="utf-8"))
    assert report["status"] == "PASS"
    assert report["pairing_verified"] is True
    assert report["all_a_b_pairs_checked"] == 14
    assert report["interior_a_b_pairs_checked"] == 14
    assert len(report["results"]) == 5
    assert all(item["replicates_per_cell"] == 2 for item in report["results"])
    assert all(item["stratified_bootstrap_95_low"] <= item["stratified_bootstrap_95_high"]
               for item in report["results"])
