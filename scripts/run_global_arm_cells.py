"""Global restartable scheduler for arm × validation-size × replicate × curvature cells.

Input preparation remains a separate, inference-free step. Supply a JSON plan
whose bundles each point to one frozen shared-input bundle. One global checkpoint
then owns every requested cell across all bundles.
"""
from __future__ import annotations

import argparse
from concurrent.futures import ProcessPoolExecutor, as_completed
import hashlib
import json
from pathlib import Path
import sys

import pandas as pd

sys.path.insert(0, str(Path(__file__).resolve().parent))

from dark_siren_robustness.arm_cell_execution import ARM_NAMES
from dark_siren_robustness.cell_ledger_store import (
    finalize_interrupted_cells, initialize_checkpoint, read_checkpoint,
)
from dark_siren_robustness.shared_input_bundle import load_shared_input_manifest

from run_arm_cell_replicate import _cell_id, _execute_one, _sha256


def _read_bundle_plan(path: Path) -> list[dict]:
    raw = json.loads(path.read_text(encoding="utf-8"))
    bundles = raw.get("bundles") if isinstance(raw, dict) else None
    if not isinstance(bundles, list) or not bundles:
        raise ValueError("bundle plan must contain a nonempty 'bundles' list")
    normalized = []
    seen = set()
    for item in bundles:
        if not isinstance(item, dict) or not {"replicate", "input_dir"} <= set(item):
            raise ValueError("each bundle requires replicate and input_dir")
        replicate = int(item["replicate"])
        root = Path(item["input_dir"]).expanduser().resolve()
        if replicate < 0 or not root.is_dir():
            raise ValueError("replicate must be nonnegative and input_dir must exist")
        config = json.loads((root / "config.json").read_text(encoding="utf-8"))
        if config.get("status") != "shared_inputs_only_no_inference":
            raise ValueError(f"bundle is not inference-free: {root}")
        gamma = float(config["curvature_gamma"])
        key = (replicate, gamma)
        if key in seen:
            raise ValueError(f"duplicate replicate/curvature bundle: {key}")
        seen.add(key)
        shared = load_shared_input_manifest(root, expected_validation_sizes=config["validation_sizes"])
        normalized.append({
            "replicate": replicate, "gamma": gamma, "input_dir": root,
            "config": config, "shared": shared,
        })
    normalized.sort(key=lambda row: (row["replicate"], row["gamma"]))
    return normalized


def run(args: argparse.Namespace) -> dict:
    bundles = _read_bundle_plan(args.bundle_plan.resolve())
    if args.max_workers < 1:
        raise ValueError("max_workers must be positive")
    arms = list(args.arms)
    if not arms or len(set(arms)) != len(arms) or any(arm not in ARM_NAMES for arm in arms):
        raise ValueError("arms must be a nonempty unique subset of declared arms")

    output = args.output_dir.resolve()
    output.mkdir(parents=True, exist_ok=True)
    checkpoint_path = output / "cell_lifecycle_checkpoint.json"
    planned = []
    for bundle in bundles:
        root, config, shared = bundle["input_dir"], bundle["config"], bundle["shared"]
        manifest_hash = _sha256(root / "shared_input_manifest.json")
        config_hash = _sha256(root / "config.json")
        provenance = shared["provenance"]
        for nval in config["validation_sizes"]:
            for arm in arms:
                row = {
                    "cell_id": None,
                    "replicate": bundle["replicate"],
                    "gamma": bundle["gamma"],
                    "validation_size": int(nval),
                    "arm": arm,
                    "input_dir": str(root),
                    "shared_input_manifest_sha256": manifest_hash,
                    "config_sha256": config_hash,
                    "analysis_seed": int(provenance["seed"]),
                    "event_seed": provenance.get("event_seed"),
                    "validation_seed": int(provenance["validation_seed"]),
                    "validation_sample_seed": provenance.get("validation_sample_seed"),
                    "h0_true": float(config["h0_true"]),
                    "inclusion_model_degree": int(config["inclusion_model_degree"]),
                }
                row["cell_id"] = _cell_id(row)
                planned.append(row)
    keys = [(r["replicate"], r["gamma"], r["validation_size"], r["arm"]) for r in planned]
    if len(keys) != len(set(keys)):
        raise ValueError("global plan contains duplicate cell identities")

    if args.resume:
        if not checkpoint_path.is_file():
            raise ValueError("--resume requires an existing global checkpoint")
        existing = read_checkpoint(checkpoint_path)
        if existing["planned_cells"] != planned:
            raise ValueError("bundle plan or execution plan differs from immutable global checkpoint")
        finalize_interrupted_cells(checkpoint_path, "Global scheduler restart recovered unfinished cells.")
    else:
        if any(path.name != ".gitkeep" for path in output.iterdir()):
            raise ValueError("output directory must be empty unless --resume is specified")
        initialize_checkpoint(checkpoint_path, planned)

    checkpoint = read_checkpoint(checkpoint_path)
    terminal = {"SUCCESS", "NON_EVALUABLE", "FAILED"}
    pending = {
        (r["replicate"], r["gamma"], r["validation_size"], r["arm"])
        for r in checkpoint["cell_ledger"] if r["status"] not in terminal
    }
    tasks = [{
        "input_dir": row["input_dir"], "output_dir": str(output),
        "checkpoint_path": str(checkpoint_path), "plan_row": row,
    } for row in planned
        if (row["replicate"], row["gamma"], row["validation_size"], row["arm"]) in pending]

    outcomes = []
    if args.max_workers == 1:
        outcomes = [_execute_one(task) for task in tasks]
    elif tasks:
        with ProcessPoolExecutor(max_workers=args.max_workers) as pool:
            futures = [pool.submit(_execute_one, task) for task in tasks]
            for future in as_completed(futures):
                outcomes.append(future.result())

    checkpoint = read_checkpoint(checkpoint_path)
    ledger = checkpoint["cell_ledger"]
    ledger_frame = pd.DataFrame(ledger)
    summary_columns = [
        "replicate", "gamma", "validation_size", "arm", "status", "h0_true",
        "median_error", "posterior_median", "q025", "q05", "q16", "q84", "q95", "q975",
        "weight_diagnostics", "realized_bandwidths", "input_sha256",
    ]
    summaries = []
    for row in ledger:
        if row["status"] == "SUCCESS":
            report_path = output / row["audit_report_path"]
            summaries.append(json.loads((report_path.parent / "result" / "summary.json").read_text()))
    summary_frame = pd.DataFrame(summaries, columns=summary_columns)
    ledger_frame.to_csv(output / "cell_ledger.csv", index=False)
    summary_frame.to_csv(output / "cell_summaries.csv", index=False)
    status_counts = {s: int((ledger_frame["status"] == s).sum())
                     for s in ("FAILED", "NON_EVALUABLE", "NOT_RUN", "SUCCESS")}
    status = "PASS" if len(ledger) and (ledger_frame["status"] == "SUCCESS").all() and (ledger_frame["audit_status"] == "PASS").all() else "INCOMPLETE"
    run_summary = {
        "status": status, "planned_cells": len(planned), "status_counts": status_counts,
        "replicate_count": len({r["replicate"] for r in planned}),
        "curvature_count": len({r["gamma"] for r in planned}),
        "bundle_count": len(bundles),
        "completed_this_invocation": len(outcomes),
        "successful_this_invocation": sum(x["status"] == "SUCCESS" for x in outcomes),
        "failed_this_invocation": sum(x["status"] == "FAILED" for x in outcomes),
        "max_workers": args.max_workers,
        "claim_limit": "Global per-cell execution lifecycle; not confirmatory calibration.",
    }
    (output / "run_summary.json").write_text(json.dumps(run_summary, indent=2, sort_keys=True) + "\n")
    files = sorted(p for p in output.rglob("*") if p.is_file()
                   and p.name != "run_manifest.json" and not p.name.endswith(".lock"))
    manifest = {"status": status, "file_count": len(files), "files": {
        str(p.relative_to(output)): {"bytes": p.stat().st_size, "sha256": _sha256(p)}
        for p in files}}
    (output / "run_manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True) + "\n")
    return run_summary


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bundle-plan", type=Path, required=True,
                        help="JSON object with bundles:[{replicate,input_dir}, ...]")
    parser.add_argument("--output-dir", type=Path, required=True)
    parser.add_argument("--arms", type=lambda x: [v.strip() for v in x.split(",")],
                        default=list(ARM_NAMES))
    parser.add_argument("--max-workers", type=int, default=1)
    parser.add_argument("--resume", action="store_true")
    args = parser.parse_args()
    print(json.dumps(run(args), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
