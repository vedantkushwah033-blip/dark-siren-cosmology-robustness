from __future__ import annotations

import argparse
import json
from pathlib import Path

from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a non-confirmatory estimated-inclusion pilot.")
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-sizes", default="100,300,1000")
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--h0-points", type=int, default=81)
    parser.add_argument("--z-points", type=int, default=401)
    parser.add_argument("--output-dir", type=Path, required=True)
    args = parser.parse_args()

    sizes = tuple(int(value.strip()) for value in args.validation_sizes.split(",") if value.strip())
    frame, summary = run_estimated_inclusion_pilot(
        n_replicates=args.replicates,
        n_parent=args.n_parent,
        n_detected=args.n_detected,
        validation_sizes=sizes,
        seed=args.seed,
        h0_points=args.h0_points,
        z_points=args.z_points,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    frame.to_csv(args.output_dir / "replicates.csv", index=False)
    (args.output_dir / "summary.json").write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "rows": int(len(frame)),
        "failed_rows": int((frame["status"] != "success").sum()),
        "status": summary["status"],
    }, indent=2))


if __name__ == "__main__":
    main()
