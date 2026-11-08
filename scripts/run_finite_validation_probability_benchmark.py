"""Run the finite-validation inclusion-probability benchmark."""
from __future__ import annotations
import argparse
from pathlib import Path

from dark_siren_robustness.estimated_inclusion_benchmark import (
    run_probability_validation_benchmark,
)


def main():
    parser = argparse.ArgumentParser()
    parser.add_argument("--output", type=Path, required=True)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--replicates", type=int, default=20)
    parser.add_argument("--seed", type=int, default=20261108)
    parser.add_argument("--raw-output", type=Path, required=True)
    args = parser.parse_args()
    result = run_probability_validation_benchmark(
        n_parent=args.n_parent, replicates=args.replicates, seed=args.seed,
        raw_output=args.raw_output,
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    result.to_csv(args.output, index=False)
    print(f"wrote {len(result)} attempt rows to {args.output}")
    print(f"wrote unit-level audit records to {args.raw_output}")
    print(result.groupby(
        ["gamma", "n_validation", "model_degree", "fit_status"],
        dropna=False,
    ).size().to_string())


if __name__ == "__main__":
    main()
