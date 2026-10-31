"""Run the prespecified nonlinear inclusion-curvature screening campaign.

This is a pipeline/overlap screen, not a confirmatory analysis. Each curvature
condition reuses the same deterministic replicate seed schedule so parent
populations, GW events, and catalogue-uniform draws are paired across curvature.
Every condition is independently audited before it is included in the aggregate.
"""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot
from scripts.audit_estimated_inclusion_pilot import audit_pilot_artifacts


CURVATURES = (-2.0, -1.0, 0.0, 1.0, 2.0)


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run paired nonlinear inclusion-curvature smoke screen (not confirmatory)."
    )
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-sizes", default="100,300,1000")
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--h0-points", type=int, default=41)
    parser.add_argument("--z-points", type=int, default=201)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.replicates < 1:
        parser.error("--replicates must be positive")
    sizes = tuple(int(v.strip()) for v in args.validation_sizes.split(",") if v.strip())
    if not sizes or len(set(sizes)) != len(sizes) or any(v < 10 for v in sizes):
        parser.error("--validation-sizes must be distinct integer counts >= 10")

    args.output_dir.mkdir(parents=True, exist_ok=True)
    condition_frames: list[pd.DataFrame] = []
    condition_summaries: list[dict] = []
    for curvature in CURVATURES:
        label = f"g_{curvature:+.1f}".replace("+", "plus").replace("-", "minus").replace(".", "p")
        condition_dir = args.output_dir / label
        condition_dir.mkdir(parents=True, exist_ok=True)

        # The same root seed across g gives common random numbers across
        # conditions; only the retention probability's curvature changes.
        frame, summary = run_estimated_inclusion_pilot(
            n_replicates=args.replicates,
            n_parent=args.n_parent,
            n_detected=args.n_detected,
            validation_sizes=sizes,
            seed=args.seed,
            true_quadratic=curvature,
            h0_points=args.h0_points,
            z_points=args.z_points,
        )
        frame = frame.copy()
        frame["true_quadratic"] = float(curvature)
        frame.to_csv(condition_dir / "replicates.csv", index=False)
        (condition_dir / "summary.json").write_text(
            json.dumps(summary, indent=2, sort_keys=True) + "\\n", encoding="utf-8"
        )

        audit = audit_pilot_artifacts(condition_dir)
        (condition_dir / "independent_audit.json").write_text(
            json.dumps(audit, indent=2, sort_keys=True) + "\\n", encoding="utf-8"
        )
        if audit.get("status") != "PASS":
            raise RuntimeError(f"independent audit failed for curvature {curvature}: {audit}")

        condition_frames.append(frame)
        condition_summaries.append({
            "true_quadratic": float(curvature),
            "condition": label,
            "audit_status": audit["status"],
            "n_rows": int(len(frame)),
            "n_failed_rows": int((frame["status"] != "success").sum()),
            "n_boundary_limited_successful_fits": int(audit["n_boundary_limited_successful_fits"]),
            "n_catalogue_rows": int((frame["arm"] == "catalogue_only").sum()),
            "probability_min_min": float(frame["true_probability_min"].min()),
            "probability_max_max": float(frame["true_probability_max"].max()),
            "fraction_below_005_max": float(frame["true_probability_below_005"].max()),
            "fraction_below_010_max": float(frame["true_probability_below_010"].max()),
            "max_true_inverse_probability": float(frame["true_max_inverse_probability"].max()),
            "audit": audit,
            "pilot_summary": summary,
        })

    combined = pd.concat(condition_frames, ignore_index=True)
    combined.to_csv(args.output_dir / "all_conditions_replicates.csv", index=False)
    report = {
        "study": "paired_nonlinear_inclusion_curvature_screen",
        "status": "SCREEN_NOT_CONFIRMATORY",
        "seed": int(args.seed),
        "replicates_per_curvature": int(args.replicates),
        "n_parent": int(args.n_parent),
        "n_detected": int(args.n_detected),
        "validation_sizes": list(sizes),
        "curvature_levels": list(CURVATURES),
        "common_random_numbers_across_curvature": True,
        "h0_grid_points": int(args.h0_points),
        "redshift_grid_points": int(args.z_points),
        "conditions": condition_summaries,
        "aggregate_rows": int(len(combined)),
        "aggregate_failed_rows": int((combined["status"] != "success").sum()),
        "interpretation_limit": (
            "Pipeline and overlap screening only. The five curvature levels and "
            "validation sizes are exploratory as a family; no confirmatory claim, "
            "method winner, or real-observation inference is supported."
        ),
    }
    (args.output_dir / "screen_report.json").write_text(
        json.dumps(report, indent=2, sort_keys=True) + "\\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": report["status"],
        "conditions": len(condition_summaries),
        "replicates_per_curvature": args.replicates,
        "aggregate_rows": len(combined),
        "aggregate_failed_rows": report["aggregate_failed_rows"],
        "audit_statuses": [item["audit_status"] for item in condition_summaries],
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
