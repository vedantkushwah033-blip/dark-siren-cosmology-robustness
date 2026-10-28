"""Independently summarize fixed-truth coverage heterogeneity from raw posterior rows.

This script reads archived raw records, not the experiment runner or its summary
JSON. It is descriptive: it does not assert nominal calibration or astrophysical validity.
"""
from __future__ import annotations

import argparse
import csv
import hashlib
import json
import math
from pathlib import Path

REQUIRED = {
    "true_h0", "fractional_error", "replicate", "arm",
    "lower_68", "upper_68", "lower_95", "upper_95",
}
EXPECTED_ARMS = (
    "A_fixed_realized_sigma",
    "B_fractional_sigma_H0",
    "C_fixed_reference_sigma",
)


def wilson(successes: int, n: int, z: float = 1.959963984540054) -> list[float]:
    if n <= 0 or not 0 <= successes <= n:
        raise ValueError("Wilson interval requires n > 0 and 0 <= successes <= n")
    p = successes / n
    denom = 1.0 + z * z / n
    center = (p + z * z / (2.0 * n)) / denom
    half = z * math.sqrt(p * (1.0 - p) / n + z * z / (4.0 * n * n)) / denom
    low = 0.0 if successes == 0 else max(0.0, center - half)
    high = 1.0 if successes == n else min(1.0, center + half)
    return [low, high]


def audit(path: Path, *, require_full_design: bool = False) -> dict:
    raw = path.read_bytes()
    digest = hashlib.sha256(raw).hexdigest()
    rows = []
    with path.open("r", encoding="utf-8-sig", newline="") as stream:
        reader = csv.DictReader(stream)
        missing = REQUIRED - set(reader.fieldnames or [])
        if missing:
            raise ValueError(f"missing required columns: {sorted(missing)}")
        for line_number, row in enumerate(reader, start=2):
            try:
                truth = float(row["true_h0"])
                frac = float(row["fractional_error"])
                replicate = int(row["replicate"])
                arm = row["arm"]
                bounds = [float(row[k]) for k in ("lower_68", "upper_68", "lower_95", "upper_95")]
            except (TypeError, ValueError) as exc:
                raise ValueError(f"invalid numeric value at CSV line {line_number}") from exc
            if not all(math.isfinite(x) for x in [truth, frac, *bounds]):
                raise ValueError(f"non-finite value at CSV line {line_number}")
            if bounds[0] > bounds[1] or bounds[2] > bounds[3]:
                raise ValueError(f"reversed interval at CSV line {line_number}")
            if bounds[2] > bounds[0] or bounds[3] < bounds[1]:
                raise ValueError(f"95% interval does not contain 68% interval at CSV line {line_number}")
            if replicate < 0 or frac <= 0:
                raise ValueError(f"negative replicate index or nonpositive fractional error at CSV line {line_number}")
            rows.append({
                "true_h0": truth,
                "fractional_error": frac,
                "replicate": replicate,
                "arm": arm,
                "covered68": bounds[0] <= truth <= bounds[1],
                "covered95": bounds[2] <= truth <= bounds[3],
            })
    if not rows:
        raise ValueError("input CSV has no data rows")
    keys = [(r["true_h0"], r["fractional_error"], r["replicate"], r["arm"]) for r in rows]
    if len(set(keys)) != len(keys):
        raise ValueError("duplicate (truth, fractional_error, replicate, arm) keys")
    arms = sorted({r["arm"] for r in rows})
    if set(arms) != set(EXPECTED_ARMS):
        raise ValueError(f"unexpected arm set: {arms}")

    design_complete = None
    if require_full_design:
        expected_truths = {50.0, 55.0, 60.0, 65.0, 70.0, 75.0, 80.0, 85.0, 90.0}
        expected_fractions = {0.05, 0.10, 0.15, 0.20}
        expected_replicates = set(range(300))
        if {r["true_h0"] for r in rows} != expected_truths:
            raise ValueError("full design has unexpected or missing true_h0 cells")
        if {r["fractional_error"] for r in rows} != expected_fractions:
            raise ValueError("full design has unexpected or missing fractional_error cells")
        for truth in expected_truths:
            for frac in expected_fractions:
                for arm in EXPECTED_ARMS:
                    cell = [r for r in rows if r["true_h0"] == truth and r["fractional_error"] == frac and r["arm"] == arm]
                    if len(cell) != 300 or {r["replicate"] for r in cell} != expected_replicates:
                        raise ValueError(f"incomplete replicate cell: true_h0={truth}, fractional_error={frac}, arm={arm}")
        if len(rows) != 32400:
            raise ValueError(f"full design expected 32400 rows, found {len(rows)}")
        design_complete = True

    def summarize(group: list[dict]) -> dict:
        n = len(group)
        out = {"n": n}
        for metric in ("covered68", "covered95"):
            successes = sum(r[metric] for r in group)
            out[metric.replace("covered", "coverage_")] = {
                "successes": successes,
                "n": n,
                "estimate": successes / n,
                "wilson_95": wilson(successes, n),
            }
        return out

    all_by_arm = {}
    interior_by_arm = {}
    boundary_by_arm = {}
    per_truth = []
    per_truth_fraction = []
    for arm in EXPECTED_ARMS:
        group = [r for r in rows if r["arm"] == arm]
        all_by_arm[arm] = summarize(group)
        interior = [r for r in group if 50.0 < r["true_h0"] < 90.0]
        boundary = [r for r in group if r["true_h0"] in (50.0, 90.0)]
        interior_by_arm[arm] = summarize(interior)
        boundary_by_arm[arm] = summarize(boundary)

    truths = sorted({r["true_h0"] for r in rows})
    fractions = sorted({r["fractional_error"] for r in rows})
    for truth in truths:
        for arm in EXPECTED_ARMS:
            group = [r for r in rows if r["true_h0"] == truth and r["arm"] == arm]
            per_truth.append({"true_h0": truth, "arm": arm, **summarize(group)})
            for frac in fractions:
                cell = [r for r in group if r["fractional_error"] == frac]
                if cell:
                    per_truth_fraction.append({
                        "true_h0": truth, "fractional_error": frac, "arm": arm, **summarize(cell)
                    })

    return {
        "status": "PASS",
        "scope": "independent descriptive audit of archived synthetic raw posterior records",
        "not_real_gw_data": True,
        "source_file": path.name,
        "source_sha256": digest,
        "raw_record_count": len(rows),
        "unique_record_key_count": len(set(keys)),
        "full_design_required": require_full_design,
        "full_design_complete": design_complete,
        "coverage_rule": "closed interval: lower <= true_h0 <= upper",
        "boundary_truths": [50.0, 90.0],
        "interior_truth_rule": "50 < true_h0 < 90; this experiment's interior truths are 55,60,65,70,75,80,85",
        "overall_by_arm": all_by_arm,
        "interior_by_arm": interior_by_arm,
        "boundary_by_arm": boundary_by_arm,
        "by_truth_and_arm": per_truth,
        "by_truth_fraction_and_arm": per_truth_fraction,
        "interpretation_warning": (
            "Overall coverage pools distinct fixed-truth cells and is not a single-parameter "
            "calibration estimand. Exact prior-boundary truths have zero coverage under the "
            "stored central equal-tail intervals; report those cells separately from interior truths."
        ),
    }


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input_csv", type=Path)
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--require-full-design", action="store_true", help="Require all 9 truths × 4 error fractions × 300 replicates × 3 arms")
    args = parser.parse_args()
    report = audit(args.input_csv, require_full_design=args.require_full_design)
    text = json.dumps(report, indent=2, sort_keys=True) + "\n"
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(text, encoding="utf-8")
    print(json.dumps({
        "status": report["status"],
        "raw_record_count": report["raw_record_count"],
        "source_sha256": report["source_sha256"],
        "overall_68": {k: v["coverage_68"]["estimate"] for k, v in report["overall_by_arm"].items()},
        "interior_68": {k: v["coverage_68"]["estimate"] for k, v in report["interior_by_arm"].items()},
    }, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
