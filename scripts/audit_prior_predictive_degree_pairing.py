"""Independent cross-degree pairing audit for Stage 0 prior-predictive H0 runs."""
from __future__ import annotations

import argparse
from pathlib import Path

import pandas as pd


def _equal_csv(left: Path, right: Path, label: str) -> None:
    a = pd.read_csv(left)
    b = pd.read_csv(right)
    try:
        pd.testing.assert_frame_equal(
            a, b, check_dtype=False, check_exact=True, check_like=False
        )
    except AssertionError as exc:
        raise AssertionError(f"degree-pair input mismatch for {label}: {exc}") from exc


def audit(degree_one: Path, degree_two: Path) -> dict:
    d1, d2 = degree_one.resolve(), degree_two.resolve()
    c1 = __import__("json").loads((d1 / "config.json").read_text())
    c2 = __import__("json").loads((d2 / "config.json").read_text())
    if int(c1["inclusion_model_degree"]) != 1 or int(c2["inclusion_model_degree"]) != 2:
        raise AssertionError("input directories are not degree-one and degree-two runs")
    paired_keys = (
        "master_seed", "replicates", "gammas", "n_parent", "n_detected",
        "validation_frame", "validation_sizes", "h0_prior_support",
        "h0_truth_mode", "h0_truth_schedule", "seed_schedule",
    )
    for key in paired_keys:
        if c1[key] != c2[key]:
            raise AssertionError(f"degree-one/degree-two top-level config mismatch: {key}")
    if c1["h0_truth_mode"] != "prior_predictive_uniform":
        raise AssertionError("paired runs are not prior-predictive")
    l1 = pd.read_csv(d1 / "replicate_ledger.csv")
    l2 = pd.read_csv(d2 / "replicate_ledger.csv")
    if not (l1["status"] == "SUCCESS").all() or not (l2["status"] == "SUCCESS").all():
        raise AssertionError("one or more planned cells did not succeed in both degrees")
    if not (l1["audit_status"] == "PASS").all() or not (l2["audit_status"] == "PASS").all():
        raise AssertionError("one or more child audits did not pass in both degrees")
    fields = [
        "cell_id", "replicate", "gamma", "validation_size", "arm",
        "analysis_seed", "event_seed", "validation_seed",
        "validation_sample_seed", "h0_truth_seed", "h0_true",
    ]
    a = l1[fields].sort_values("cell_id").reset_index(drop=True)
    b = l2[fields].sort_values("cell_id").reset_index(drop=True)
    pd.testing.assert_frame_equal(a, b, check_dtype=False, check_exact=True)
    reps = int(c1["replicates"])
    gammas = [float(x) for x in c1["gammas"]]
    sizes = [int(x) for x in c1["validation_sizes"]]
    common_files = [
        "truth/parent_truth.csv",
        "truth/gw_event_truth.csv",
        "truth/gw_generation_trace_truth.csv",
        "truth/validation_frame_truth.csv",
        "observable/catalogue_observed.csv",
        "observable/gw_events_observed.csv",
        *[f"observable/validation_n{n}.csv" for n in sizes],
    ]
    compared = 0
    for rep in range(reps):
        for gamma in gammas:
            cell = f"rep{rep:03d}_g{gamma:+g}".replace("+", "p").replace("-", "m")
            for relative in common_files:
                _equal_csv(
                    d1 / "replicates" / cell / relative,
                    d2 / "replicates" / cell / relative,
                    f"replicate={rep}, gamma={gamma}, file={relative}",
                )
                compared += 1
    return {
        "status": "PASS",
        "replicates": reps,
        "curvatures": len(gammas),
        "validation_sizes": sizes,
        "planned_cells_per_degree": len(l1),
        "common_input_files_compared": compared,
        "interpretation": "Independent cross-degree common-input and seed/truth pairing audit; not calibration evidence.",
    }


def main() -> None:
    parser = argparse.ArgumentParser()
    parser.add_argument("degree_one_dir", type=Path)
    parser.add_argument("degree_two_dir", type=Path)
    args = parser.parse_args()
    import json
    print(json.dumps(audit(args.degree_one_dir, args.degree_two_dir), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
