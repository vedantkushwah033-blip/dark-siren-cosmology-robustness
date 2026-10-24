"""Independent structural audit for factorial sensitivity artifacts.

This validator checks experiment design, row counts, finite outputs, exact
catalogue masks, posterior interval ordering, and scenario/manifest agreement.
It does not independently validate the simulator's scientific assumptions.
"""
from __future__ import annotations

import argparse
import csv
import json
import math
from pathlib import Path

EXPECTED = {
    "complete_uniform_matched": (1.0, "complete", "uniform", "uniform"),
    "complete_luminosity_matched": (1.0, "complete", "luminosity", "luminosity"),
    "complete_luminosity_host_uniform_inference": (1.0, "complete", "luminosity", "uniform"),
    "random70_uniform_matched": (0.7, "random", "uniform", "uniform"),
    "random70_luminosity_matched": (0.7, "random", "luminosity", "luminosity"),
    "random70_luminosity_host_uniform_inference": (0.7, "random", "luminosity", "uniform"),
    "redshift70_uniform_matched": (0.7, "redshift", "uniform", "uniform"),
    "redshift70_luminosity_matched": (0.7, "redshift", "luminosity", "luminosity"),
    "redshift70_luminosity_host_uniform_inference": (0.7, "redshift", "luminosity", "uniform"),
}
POSTERIOR_FIELDS = ("median", "lower_68", "upper_68", "lower_95", "upper_95", "width_68", "width_95")


def read_csv(path: Path) -> list[dict[str, str]]:
    with path.open(newline="", encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def finite(row: dict[str, str], key: str, label: str) -> float:
    try:
        value = float(row[key])
    except (KeyError, TypeError, ValueError) as exc:
        raise ValueError(f"{label}: missing/non-numeric {key}") from exc
    if not math.isfinite(value):
        raise ValueError(f"{label}: non-finite {key}={value}")
    return value


def audit(root: Path) -> dict:
    summary_json = json.loads((root / "scenario_summary.json").read_text(encoding="utf-8"))
    summary = read_csv(root / "scenario_summary.csv")
    runs = read_csv(root / "per_run_results.csv")
    n = int(summary_json["requested_replicates_per_scenario"])
    if n < 2:
        raise ValueError("replicate count must be at least 2")
    if int(summary_json["scenario_count"]) != len(EXPECTED) or len(summary) != len(EXPECTED):
        raise ValueError("expected exactly nine factorial scenarios")
    if len(runs) != len(EXPECTED) * n:
        raise ValueError(f"expected {len(EXPECTED)*n} per-run rows, found {len(runs)}")

    summary_by_name = {row["scenario"]: row for row in summary}
    if set(summary_by_name) != set(EXPECTED):
        raise ValueError("summary scenario names do not match prespecified design")
    grouped: dict[str, list[dict[str, str]]] = {name: [] for name in EXPECTED}
    for row in runs:
        name = row.get("scenario", "")
        if name not in EXPECTED:
            raise ValueError(f"unexpected scenario in raw output: {name!r}")
        grouped[name].append(row)

    checked = 0
    for name, design in EXPECTED.items():
        target, missingness, host_generation, inference_weighting = design
        srow = summary_by_name[name]
        for field, expected in (
            ("target_completeness", target),
            ("missingness", missingness),
            ("host_generation", host_generation),
            ("inference_weighting", inference_weighting),
        ):
            if field not in srow or (field in ("missingness", "host_generation", "inference_weighting")
                                     and srow[field] != expected):
                raise ValueError(f"{name}: summary design mismatch for {field}")
            if field == "target_completeness" and not math.isclose(float(srow[field]), expected, abs_tol=1e-12):
                raise ValueError(f"{name}: target completeness mismatch")
        if int(srow["n"]) != n or len(grouped[name]) != n:
            raise ValueError(f"{name}: replicate count mismatch")
        indices = set()
        for row in grouped[name]:
            label = f"{name} run={row.get('run', '?')}"
            indices.add(row.get("run"))
            target_value = finite(row, "target_completeness", label)
            realized = finite(row, "realized_completeness", label)
            if not math.isclose(target_value, target, abs_tol=1e-12):
                raise ValueError(f"{label}: per-run target mismatch")
            if target == 0.7 and not math.isclose(realized, 0.7, abs_tol=1e-12):
                raise ValueError(f"{label}: expected exact 70% completeness, got {realized}")
            vals = {key: finite(row, key, label) for key in POSTERIOR_FIELDS}
            if not (vals["lower_95"] <= vals["lower_68"] <= vals["median"]
                    <= vals["upper_68"] <= vals["upper_95"]):
                raise ValueError(f"{label}: posterior interval ordering invalid")
            if not math.isclose(vals["width_68"], vals["upper_68"] - vals["lower_68"], rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"{label}: width_68 inconsistent with interval endpoints")
            if not math.isclose(vals["width_95"], vals["upper_95"] - vals["lower_95"], rel_tol=1e-9, abs_tol=1e-9):
                raise ValueError(f"{label}: width_95 inconsistent with interval endpoints")
            checked += 1
        if len(indices) != n:
            raise ValueError(f"{name}: duplicate/missing run identifiers")
    return {"status": "PASS", "scenarios": len(EXPECTED), "replicates_per_scenario": n,
            "raw_rows_checked": checked, "summary_rows_checked": len(summary),
            "scientific_scope": "synthetic structural audit only"}


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("results_dir", nargs="?", default="results/factorial-sensitivity")
    args = parser.parse_args()
    print(json.dumps(audit(Path(args.results_dir)), indent=2))


if __name__ == "__main__":
    main()
