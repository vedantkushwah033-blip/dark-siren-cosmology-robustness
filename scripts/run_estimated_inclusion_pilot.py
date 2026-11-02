from __future__ import annotations

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform

from dark_siren_robustness.estimated_inclusion_study import run_estimated_inclusion_pilot


def _file_record(path: Path, *, rows: int | None = None) -> dict[str, object]:
    payload = path.read_bytes()
    return {
        "bytes": len(payload),
        "sha256": hashlib.sha256(payload).hexdigest(),
        "rows": rows,
    }


def write_pilot_manifest(
    output_dir: str | Path,
    *,
    frame,
    summary: dict,
    arguments: dict[str, object],
    source_commit: str | None = None,
) -> dict[str, object]:
    """Write a deterministic, per-file provenance manifest for a pilot output.

    The manifest covers the two primary raw outputs, not itself or the later
    independent-audit reports. A GitHub Actions run supplies GITHUB_SHA;
    local runs may legitimately record null when no source revision is known.
    """
    root = Path(output_dir)
    csv_path = root / "replicates.csv"
    summary_path = root / "summary.json"
    for path in (csv_path, summary_path):
        if not path.is_file():
            raise FileNotFoundError(f"cannot manifest missing pilot output: {path.name}")

    if source_commit is not None:
        source_commit = str(source_commit).strip().lower()
        if len(source_commit) != 40 or any(ch not in "0123456789abcdef" for ch in source_commit):
            raise ValueError("source_commit must be a full 40-character hexadecimal SHA")

    manifest = {
        "schema": "estimated_inclusion_pilot_manifest_v1",
        "study": "estimated_catalogue_inclusion_probability_pilot",
        "status": str(summary.get("status", "UNKNOWN")),
        "source_commit": source_commit,
        "source_commit_source": "GITHUB_SHA" if source_commit else "not_available",
        "arguments": arguments,
        "environment": {
            "python": platform.python_version(),
            "platform": platform.platform(),
            "software_versions": summary.get("software_versions", {}),
        },
        "expected_rows": int(len(frame)),
        "files": {
            "replicates.csv": _file_record(csv_path, rows=int(len(frame))),
            "summary.json": _file_record(summary_path),
        },
    }
    manifest_path = root / "manifest.json"
    manifest_path.write_text(
        json.dumps(manifest, indent=2, sort_keys=True) + "\n",
        encoding="utf-8",
    )
    return manifest


def main() -> None:
    parser = argparse.ArgumentParser(description="Run a non-confirmatory estimated-inclusion pilot.")
    parser.add_argument("--replicates", type=int, default=10)
    parser.add_argument("--n-parent", type=int, default=4000)
    parser.add_argument("--n-detected", type=int, default=30)
    parser.add_argument("--validation-sizes", default="100,300,1000")
    parser.add_argument("--seed", type=int, default=20261030)
    parser.add_argument("--true-quadratic", type=float, default=0.0,
                        help="Quadratic logit coefficient for nonlinear selection stress; zero preserves linear selection.")
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
        true_quadratic=args.true_quadratic,
        h0_points=args.h0_points,
        z_points=args.z_points,
    )
    args.output_dir.mkdir(parents=True, exist_ok=True)
    csv_path = args.output_dir / "replicates.csv"
    summary_path = args.output_dir / "summary.json"
    frame.to_csv(csv_path, index=False)
    summary_path.write_text(
        json.dumps(summary, indent=2, sort_keys=True) + "\n", encoding="utf-8"
    )
    manifest = write_pilot_manifest(
        args.output_dir,
        frame=frame,
        summary=summary,
        arguments={
            "replicates": args.replicates,
            "n_parent": args.n_parent,
            "n_detected": args.n_detected,
            "validation_sizes": list(sizes),
            "seed": args.seed,
            "true_quadratic": args.true_quadratic,
            "h0_points": args.h0_points,
            "z_points": args.z_points,
        },
        source_commit=os.environ.get("GITHUB_SHA"),
    )
    print(json.dumps({
        "output_dir": str(args.output_dir),
        "rows": int(len(frame)),
        "failed_rows": int((frame["status"] != "success").sum()),
        "manifest_schema": manifest["schema"],
        "source_commit": manifest["source_commit"],
        "status": summary["status"],
    }, indent=2))


if __name__ == "__main__":
    main()
