"""Fail-fast validation for generated calibration artifacts."""
from __future__ import annotations

import csv
import json
from pathlib import Path
import math
import sys


def main() -> None:
    root = Path(sys.argv[1] if len(sys.argv) > 1 else "results/original-audit")
    summary_path = root / "summary.json"
    curve_path = root / "coverage_curve.csv"
    unknown_path = root / "per_run_unknown_host.csv"
    known_path = root / "per_run_known_host.csv"
    for path in (summary_path, curve_path, unknown_path, known_path):
        if not path.is_file() or path.stat().st_size == 0:
            raise SystemExit(f"Missing or empty required artifact: {path}")

    summary = json.loads(summary_path.read_text(encoding="utf-8"))
    n = int(summary["requested_replicates"])
    for scenario in ("unknown_host", "known_host"):
        metrics = summary[scenario]
        if int(metrics["n"]) != n:
            raise SystemExit(f"{scenario}: expected {n} runs, got {metrics['n']}")
        for key in ("coverage_68", "coverage_95", "coverage_68_mc_se", "coverage_95_mc_se"):
            value = float(metrics[key])
            if not math.isfinite(value) or not 0.0 <= value <= 1.0:
                raise SystemExit(f"{scenario}.{key} is invalid: {value}")

    with curve_path.open(newline="", encoding="utf-8") as f:
        rows = list(csv.DictReader(f))
    expected_levels = {0.50, 0.68, 0.80, 0.90, 0.95}
    if len(rows) != 10:
        raise SystemExit(f"Expected 10 coverage-curve rows, got {len(rows)}")
    observed = set()
    widths_by_scenario = {"unknown_host": {}, "known_host": {}}
    for row in rows:
        scenario = row["scenario"]
        level = float(row["nominal_level"])
        replicates = int(row["replicates"])
        coverage = float(row["coverage"])
        width = float(row["mean_interval_width"])
        mcse = float(row["coverage_mc_se"])
        if scenario not in {"unknown_host", "known_host"} or level not in expected_levels:
            raise SystemExit(f"Unexpected scenario/level: {row}")
        if replicates != n:
            raise SystemExit(f"{scenario} at level {level}: expected {n} replicates, got {replicates}")
        if not all(math.isfinite(x) for x in (coverage, width, mcse)):
            raise SystemExit(f"Non-finite coverage curve value: {row}")
        if not (0.0 <= coverage <= 1.0 and width > 0.0 and 0.0 <= mcse <= 0.5):
            raise SystemExit(f"Invalid coverage curve values: {row}")
        if not math.isclose(mcse, math.sqrt(coverage * (1.0 - coverage) / n), rel_tol=1e-8, abs_tol=1e-12):
            raise SystemExit(f"Coverage MC standard error inconsistent with coverage and n: {row}")
        widths_by_scenario[scenario][level] = width
        observed.add((scenario, level))
    expected = {(s, level) for s in ("unknown_host", "known_host") for level in expected_levels}
    if observed != expected:
        raise SystemExit("Coverage curve does not contain every scenario × level combination")
    for scenario, by_level in widths_by_scenario.items():
        ordered = [by_level[level] for level in sorted(expected_levels)]
        if any(right < left for left, right in zip(ordered, ordered[1:])):
            raise SystemExit(f"{scenario}: mean interval width decreases as nominal level increases")

    for path in (unknown_path, known_path):
        with path.open(newline="", encoding="utf-8") as f:
            count = sum(1 for _ in csv.DictReader(f))
        if count != n:
            raise SystemExit(f"{path}: expected {n} rows, found {count}")

    print(f"Validated calibration artifacts in {root} (n={n} per scenario).")


if __name__ == "__main__":
    main()
