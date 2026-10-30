"""Persist validation-bootstrap sensitivity draws from one pilot replicate/arm."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.estimated_inclusion_study import (
    run_validation_bootstrap_posterior_sensitivity,
)


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--pilot-csv", type=Path, required=True)
    parser.add_argument("--replicate", type=int, required=True)
    parser.add_argument("--validation-size", type=int, required=True)
    parser.add_argument("--n-bootstrap", type=int, default=500)
    parser.add_argument("--seed", type=int, default=20261031)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    if args.replicate < 0 or args.validation_size < 10:
        parser.error("replicate must be nonnegative and validation-size at least 10")
    frame = pd.read_csv(args.pilot_csv)
    arm = f"estimated_{args.validation_size}"
    selected = frame.loc[
        (frame["replicate"] == args.replicate)
        & (frame["arm"] == arm)
        & (frame["status"] == "success")
    ]
    if len(selected) != 1:
        raise SystemExit(
            f"Expected exactly one successful row for replicate={args.replicate}, "
            f"arm={arm}; found {len(selected)}"
        )
    row = selected.iloc[0]
    payload = json.loads(row["inference_audit_payload_json"])
    if payload.get("fitted_inclusion_model") is None:
        raise SystemExit("Selected row has no fitted inclusion model")
    if payload.get("validation_redshift") is None or payload.get("validation_included") is None:
        raise SystemExit("Selected row does not persist validation records")

    result = run_validation_bootstrap_posterior_sensitivity(
        redshift_grid=np.asarray(payload["redshift_grid"], dtype=float),
        observed_redshift=np.asarray(payload["catalogue_redshift"], dtype=float),
        validation_redshift=np.asarray(payload["validation_redshift"], dtype=float),
        validation_included=np.asarray(payload["validation_included"], dtype=int),
        h0_grid=np.asarray(payload["h0_grid"], dtype=float),
        observed_gw_distances_mpc=np.asarray(payload["observed_distances_mpc"], dtype=float),
        sigma_mpc=float(payload["sigma_mpc"]),
        detection_threshold_mpc=float(payload["detection_threshold_mpc"]),
        n_bootstrap=args.n_bootstrap,
        seed=args.seed,
        merger_weight=np.asarray(payload["merger_weight_grid"], dtype=float),
        bandwidth_factor=payload.get("bandwidth_factor"),
    )

    args.output_dir.mkdir(parents=True, exist_ok=True)
    array_names = (
        "density_draws", "intercept_draws", "slope_draws",
        "coefficient_boundary_limited_draws", "inverse_weight_ess_draws",
        "inverse_weight_maximum_draws", "inverse_weight_p99_draws",
        "posterior_mean_draws", "posterior_median_draws",
        "posterior_ci68_low_draws", "posterior_ci68_high_draws",
    )
    np.savez_compressed(
        args.output_dir / "bootstrap_draws.npz",
        **{name: np.asarray(result[name]) for name in array_names},
    )
    # Only inference inputs are copied. Scoring-only parent truth is deliberately excluded.
    inputs = {
        "schema": "validation_bootstrap_inputs_v1",
        "replicate": int(args.replicate),
        "arm": arm,
        "validation_size": int(args.validation_size),
        "inference_payload": payload,
    }
    (args.output_dir / "inputs.json").write_text(
        json.dumps(inputs, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    summary = {
        "schema": "validation_bootstrap_artifact_v1",
        "replicate": int(args.replicate),
        "arm": arm,
        "validation_size": int(args.validation_size),
        "n_attempted": int(result["n_attempted"]),
        "n_successful": int(result["n_successful"]),
        "n_failed": int(result["n_failed"]),
        "success_fraction": float(result["n_successful"] / result["n_attempted"]),
        "seed": int(result["seed"]),
        "bootstrap_method": "paired validation-record resampling with replacement",
        "finite_population_correction_applied": False,
        "uncertainty_scope": result["interpretation"],
        "posterior_mean_bootstrap_percentile_95": result[
            "posterior_mean_bootstrap_percentile_95"
        ],
        "failure_reasons": list(result["failure_reasons"]),
        "array_file": "bootstrap_draws.npz",
        "input_file": "inputs.json",
        "warning": (
            "This is validation-sample sensitivity conditional on the fixed observed "
            "catalogue, fixed GW events, logistic model family, and toy likelihood. "
            "It is not an H0 posterior credible interval or full uncertainty interval."
        ),
    }
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "status": "WRITTEN",
        "output_dir": str(args.output_dir),
        "n_attempted": summary["n_attempted"],
        "n_successful": summary["n_successful"],
        "n_failed": summary["n_failed"],
    }, indent=2))


if __name__ == "__main__":
    main()
