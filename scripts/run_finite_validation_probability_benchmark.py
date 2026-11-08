"""Run the finite-validation inclusion-probability benchmark."""
from __future__ import annotations

import argparse
import hashlib
import json
import platform
import time
from pathlib import Path

import numpy as np
import pandas as pd
import scipy

from dark_siren_robustness.estimated_inclusion_benchmark import (
    run_probability_validation_benchmark,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--raw-output", type=Path, required=True)
    parser.add_argument("--config", type=Path, default=None)
    parser.add_argument("--n-parent", type=int, default=None)
    parser.add_argument("--replicates", type=int, default=None)
    parser.add_argument("--seed", type=int, default=None)
    parser.add_argument("--validation-sizes", type=str, default=None)
    parser.add_argument("--curvatures", type=str, default=None)
    parser.add_argument("--degrees", type=str, default=None)
    args = parser.parse_args()

    settings = {
        "n_parent": 4000,
        "replicates": 20,
        "seed": 20261108,
        "validation_sizes": (100, 300, 1000),
        "curvature_values": (-2.0, -1.2, 0.0, 1.2, 2.0),
        "degrees": (1, 2),
    }
    if args.config is not None:
        config = json.loads(args.config.read_text(encoding="utf-8"))
        if config.get("schema_version") != 1:
            raise ValueError("unsupported finite-validation config schema_version")
        settings.update({
            "n_parent": config["n_parent"],
            "replicates": config["replicates"],
            "seed": config["seed"],
            "validation_sizes": tuple(config["validation_sizes"]),
            "curvature_values": tuple(config["curvature_values"]),
            "degrees": tuple(config["model_degrees"]),
        })
    for name in ("n_parent", "replicates", "seed"):
        value = getattr(args, name)
        if value is not None:
            settings[name] = value
    if args.validation_sizes is not None:
        settings["validation_sizes"] = tuple(int(v) for v in args.validation_sizes.split(","))
    if args.curvatures is not None:
        settings["curvature_values"] = tuple(float(v) for v in args.curvatures.split(","))
    if args.degrees is not None:
        settings["degrees"] = tuple(int(v) for v in args.degrees.split(","))

    started = time.perf_counter()
    result = run_probability_validation_benchmark(
        **settings,
        raw_output=args.raw_output,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.raw_output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    runtime_seconds = time.perf_counter() - started
    config_bytes = (
        args.config.read_bytes() if args.config is not None
        else json.dumps(settings, sort_keys=True).encode("utf-8")
    )
    manifest = {
        "artifact_schema_version": 1,
        "experiment": "finite_validation_probability_stage",
        "config": settings,
        "config_sha256": hashlib.sha256(config_bytes).hexdigest(),
        "python_version": platform.python_version(),
        "numpy_version": np.__version__,
        "pandas_version": pd.__version__,
        "scipy_version": scipy.__version__,
        "runtime_seconds": runtime_seconds,
        "summary_rows": int(len(result)),
        "raw_unit_records": int(sum(settings["n_parent"] for _ in range(len(result)))),
        "fit_status_counts": {
            str(key): int(value)
            for key, value in result["fit_status"].value_counts(dropna=False).items()
        },
        "claim_limit": "Probability-stage pilot only; no host-KDE or H0 posterior claim.",
        "summary_file": args.output.name,
        "raw_file": args.raw_output.name,
    }
    manifest_path = args.output.parent / "run_manifest.json"
    manifest_path.write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\\n", encoding="utf-8")
    print(f"wrote {len(result)} attempt rows to {args.output}")
    print(f"wrote unit-level audit records to {args.raw_output}")
    print(result.groupby(
        ["gamma", "n_validation", "model_degree", "fit_status"],
        dropna=False,
    ).size().to_string())


if __name__ == "__main__":
    main()
