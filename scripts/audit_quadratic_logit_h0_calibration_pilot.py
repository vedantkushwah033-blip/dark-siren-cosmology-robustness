"""Independent orchestration audit for the quadratic-logit H0 calibration pilot."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def audit(root: Path) -> dict:
    root=root.resolve()
    planned=pd.read_csv(root/"planned_cells.csv")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    summary=pd.read_csv(root/"calibration_summary.csv")
    manifest=json.loads((root/"manifest.json").read_text())
    bad=[name for name,digest in manifest.items()
         if not (root/name).is_file() or sha(root/name)!=digest]
    if bad: raise AssertionError(f"aggregate manifest failures: {bad}")
    if planned.cell_id.duplicated().any() or ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate cell identifiers")
    if set(planned.cell_id)!=set(ledger.cell_id):
        raise AssertionError("planned cells and ledger do not match")
    allowed={"SUCCESS","NON_EVALUABLE","FAILED","NOT_RUN"}
    if not set(ledger.status).issubset(allowed):
        raise AssertionError("unknown failure-ledger status")
    if (ledger.status=="NOT_RUN").any():
        raise AssertionError("pilot ended with cells not attempted")
    audited=0
    dirs={}
    for row in ledger.itertuples(index=False):
        d=Path(row.artifact_dir)
        if row.status in {"SUCCESS","NON_EVALUABLE"}:
            report=json.loads((d/"audit_report.json").read_text())
            if report.get("status")!="PASS": raise AssertionError(f"child audit failed: {row.cell_id}")
            cfg=json.loads((d/"config.json").read_text())
            if not np.isclose(cfg["inference"]["curvature_gamma"],row.gamma):
                raise AssertionError(f"gamma mismatch: {row.cell_id}")
            if cfg["analysis_generator_config"]["n_catalogue_observed"]!=len(pd.read_csv(d/"observable/analysis_catalogue.csv")):
                raise AssertionError(f"catalogue count mismatch: {row.cell_id}")
            dirs[row.cell_id]=d; audited+=1
        elif row.status=="FAILED":
            # Failed attempts are retained in the ledger; partial raw files remain available.
            if not row.error: raise AssertionError(f"failed cell has no error: {row.cell_id}")
    paired=0
    for (rep,gamma),block in ledger.groupby(["replicate","gamma"]):
        usable=block[block.cell_id.isin(dirs)]
        if len(usable)<2: continue
        usable=usable.sort_values("n_validation")
        first=dirs[usable.iloc[0].cell_id]
        ref_frame=pd.read_csv(first/"truth/validation_frame_truth.csv")
        ref_catalogue=pd.read_csv(first/"observable/analysis_catalogue.csv")
        ref_events=pd.read_csv(first/"observable/gw_events_observed.csv")
        prior_ids=None
        for row in usable.itertuples(index=False):
            d=dirs[row.cell_id]
            pd.testing.assert_frame_equal(pd.read_csv(d/"observable/analysis_catalogue.csv"),ref_catalogue)
            pd.testing.assert_frame_equal(pd.read_csv(d/"observable/gw_events_observed.csv"),ref_events)
            frame=pd.read_csv(d/"truth/validation_frame_truth.csv")
            cols=["frame_id","redshift","selected_for_validation"]
            pd.testing.assert_frame_equal(frame[cols],ref_frame[cols])
            ids=pd.read_csv(d/"observable/validation_observed.csv").validation_id.to_numpy()
            if prior_ids is not None and not np.array_equal(prior_ids,ids[:len(prior_ids)]):
                raise AssertionError("validation SRS samples are not nested prefixes")
            prior_ids=ids
        paired+=1
    # GW events and parent redshift realization must be invariant across gamma
    for rep,block in ledger.groupby("replicate"):
        usable=block[block.cell_id.isin(dirs)]
        if usable.empty: continue
        reference=None
        for row in usable.itertuples(index=False):
            d=dirs[row.cell_id]
            parent=pd.read_csv(d/"truth/analysis_parent_truth.csv")[["object_id","redshift"]]
            events=pd.read_csv(d/"observable/gw_events_observed.csv")
            current=(parent,events)
            if reference is None: reference=current
            else:
                pd.testing.assert_frame_equal(parent,reference[0])
                pd.testing.assert_frame_equal(events,reference[1])
    counts=ledger.status.value_counts().to_dict()
    for row in summary.itertuples(index=False):
        group=ledger[(ledger.gamma==row.gamma)&(ledger.n_validation==row.n_validation)]
        if len(group)!=row.planned or int((group.status=="SUCCESS").sum())!=row.success:
            raise AssertionError("summary table does not reconcile with replicate ledger")
    return {"status":"PASS","planned_cells":len(planned),"status_counts":counts,
            "child_audits_verified":audited,"paired_nested_blocks_checked":paired,
            "aggregate_manifest_entries_verified":len(manifest),
            "interpretation":"Arithmetic, pairing, nested design and artifact integrity only; not calibration evidence."}


def main():
    p=argparse.ArgumentParser(); p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
    a=p.parse_args(); result=audit(a.artifact_dir)
    payload=json.dumps(result,indent=2)
    if a.output: a.output.write_text(payload+"\n")
    print(payload)

if __name__=="__main__": main()
