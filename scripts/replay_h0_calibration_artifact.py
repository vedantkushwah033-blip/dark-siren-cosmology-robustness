"""Replay the corrected independent auditor on a historical complete H0 pilot artifact."""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import pandas as pd

from audit_quadratic_logit_h0_calibration_pilot import audit


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main() -> None:
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, required=True)
    parser.add_argument("--output", type=Path, required=True)
    args = parser.parse_args()
    root = args.root.resolve()

    original_manifest = json.loads((root / "manifest.json").read_text())
    original_bad = [
        name for name, digest in original_manifest.items()
        if not (root / name).is_file() or sha(root / name) != digest
    ]
    if original_bad:
        raise AssertionError(f"source artifact manifest failures: {original_bad}")

    planned = pd.read_csv(root / "planned_cells.csv")
    ledger = pd.read_csv(root / "replicate_ledger.csv")
    if len(planned) != 90 or len(ledger) != 90:
        raise AssertionError("expected the complete 90-cell pilot artifact")
    if set(planned.cell_id) != set(ledger.cell_id):
        raise AssertionError("planned cells and ledger cells differ")
    if not (ledger.status == "SUCCESS").all():
        raise AssertionError("historical artifact is not 90/90 SUCCESS")

    # The source ledger records runner-local absolute paths. Rebase only this
    # path metadata after first verifying the original manifest.
    ledger["artifact_dir"] = ledger["cell_id"].map(
        lambda cell: str(root / "replicates" / cell)
    )
    ledger.to_csv(root / "replicate_ledger.csv", index=False)

    # Re-sign the replay copy after the deterministic path rebasing. Child
    # manifests and original artifact files remain unchanged.
    replay_manifest = {
        str(path.relative_to(root)): sha(path)
        for path in sorted(root.rglob("*"))
        if path.is_file() and path.name != "manifest.json"
    }
    (root / "manifest.json").write_text(
        json.dumps(replay_manifest, indent=2, sort_keys=True) + "\n"
    )

    report = audit(root)
    report["source_manifest_entries_verified_before_path_rebase"] = len(original_manifest)
    report["source_manifest_failures"] = 0
    report["replay_manifest_entries_verified"] = len(replay_manifest)
    report["historical_summary_missing_not_run_column"] = (
        "not_run" not in pd.read_csv(root / "calibration_summary.csv").columns
    )
    report["historical_not_run_count_inferred_from_complete_ledger"] = int(
        (ledger.status == "NOT_RUN").sum()
    )
    report["interpretation"] = (
        "Independent artifact, pairing and arithmetic audit only; not SBC, "
        "frequentist calibration confirmation, or real-data validation."
    )
    args.output.parent.mkdir(parents=True, exist_ok=True)
    args.output.write_text(json.dumps(report, indent=2) + "\n")
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    main()
