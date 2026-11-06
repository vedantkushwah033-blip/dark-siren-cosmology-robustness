from __future__ import annotations

import argparse
import json
from pathlib import Path

import pandas as pd

from dark_siren_robustness.finite_parent_sampling import (
    run_finite_parent_benchmark_checkpointed,
)
from dark_siren_robustness.scenario_lifecycle_audit import audit_scenario_lifecycle


def main() -> None:
    parser = argparse.ArgumentParser(
        description="Run the finite-parent SRSWOR benchmark with scenario-level checkpoints."
    )
    parser.add_argument("--scenarios", type=int, default=10)
    parser.add_argument("--draws", type=int, default=20)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-sizes", default="100,300")
    parser.add_argument("--seed", type=int, default=20261032)
    parser.add_argument("--solver-version", choices=("v1", "v2"), default="v1")
    parser.add_argument("--h0-points", type=int, default=81)
    parser.add_argument("--z-points", type=int, default=401)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    sizes = tuple(int(value.strip()) for value in args.validation_sizes.split(",") if value.strip())
    args.output_dir.mkdir(parents=True, exist_ok=True)
    records, ledger, state = run_finite_parent_benchmark_checkpointed(
        checkpoint_dir=args.output_dir,
        n_scenarios=args.scenarios,
        n_validation_draws=args.draws,
        n_parent=args.n_parent,
        n_detected=args.n_detected,
        validation_sizes=sizes,
        seed=args.seed,
        solver_version=args.solver_version,
        h0_points=args.h0_points,
        z_points=args.z_points,
    )
    records.to_csv(args.output_dir / "records.csv", index=False)
    audit = audit_scenario_lifecycle(args.output_dir)
    successful_records = (
        int(records["status"].eq("success").sum()) if "status" in records else 0
    )
    failed_records = (
        int(records["status"].eq("failed").sum()) if "status" in records else 0
    )
    summary = {
        "schema": "finite_parent_srswor_lifecycle_v1",
        "status": "SMOKE_ONLY_NOT_CONFIRMATORY",
        "solver_version": args.solver_version,
        "seed": args.seed,
        "validation_sizes": list(sizes),
        "scenarios": args.scenarios,
        "draws_per_scenario_and_size": args.draws,
        "n_planned_attempts": int(state["attempt_counts"]["n_planned"]),
        "attempt_counts": state["attempt_counts"],
        "n_returned_scientific_records": int(len(records)),
        "n_successful_scientific_records": successful_records,
        "n_failed_draw_records": failed_records,
        "scenario_states": state["scenario_states"],
        "lifecycle_audit": audit,
        "limitations": [
            "Synthetic conditional experiment only; not real LVK data or an observed H0 measurement.",
            "The attempt ledger is the complete planned denominator; records.csv contains returned scientific draw-level outputs only.",
            "The lifecycle auditor validates accounting and checkpoint integrity, not numerical KDE, likelihood, or posterior reconstruction.",
            "A durable last checkpoint is not a resume implementation; interrupted runs must be restarted in a fresh output directory.",
            "Smoke execution is not confirmatory evidence and makes no estimator-superiority claim.",
        ],
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": summary["status"],
        "n_planned_attempts": summary["n_planned_attempts"],
        "n_returned_scientific_records": summary["n_returned_scientific_records"],
        "attempt_counts": summary["attempt_counts"],
        "output_dir": str(args.output_dir),
    }, indent=2))


if __name__ == "__main__":
    main()
