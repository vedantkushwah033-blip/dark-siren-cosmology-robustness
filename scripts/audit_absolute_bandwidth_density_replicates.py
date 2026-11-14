"""Independent orchestration auditor for multi-replicate density pilot.

Each successful replicate is independently reconstructed by the separate
single-replicate raw-input auditor. This layer audits the planned ledger,
seed schedule, outer manifest, and exact aggregation across replicates.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_absolute_bandwidth_density_smoke import audit as audit_one

def _hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def audit(root:Path):
    root=root.resolve()
    manifest=json.loads((root/"manifest.json").read_text())
    checked=0
    for rel,meta in manifest["files"].items():
        path=root/rel
        if not path.is_file() or path.stat().st_size!=int(meta["bytes"]) or _hash(path)!=meta["sha256"]:
            raise AssertionError(f"outer manifest verification failed: {rel}")
        checked+=1
    if checked!=int(manifest["file_count"]):
        raise AssertionError("outer manifest file count mismatch")
    cfg=json.loads((root/"config.json").read_text())
    if cfg.get("status")!="DEVELOPMENT_PILOT_NOT_CONFIRMATORY_NOT_BANDWIDTH_SELECTION":
        raise AssertionError("unexpected experiment claim status")
    planned=int(cfg["replicates"])
    if planned<2: raise AssertionError("not a multi-replicate experiment")
    ledger=pd.read_csv(root/"replicate_ledger.csv").fillna("")
    if len(ledger)!=planned or ledger.replicate_id.tolist()!=list(range(planned)):
        raise AssertionError("replicate ledger does not enumerate the complete planned design")
    if set(ledger.status)-{"SUCCESS","FAILED"}:
        raise AssertionError("ledger contains invalid lifecycle status")
    expected_seed=[int(cfg["master_seed"])+i*int(cfg["seed_stride"]) for i in range(planned)]
    if ledger.population_seed.astype(int).tolist()!=expected_seed:
        raise AssertionError("population seed schedule mismatch")
    if ledger.validation_seed.astype(int).tolist()!=[x+1 for x in expected_seed]:
        raise AssertionError("validation seed schedule mismatch")
    all_seeds=ledger.population_seed.astype(int).tolist()+ledger.validation_seed.astype(int).tolist()
    if len(set(all_seeds))!=len(all_seeds):
        raise AssertionError("seed streams overlap")
    if int(cfg["attempted"])!=planned or int(cfg["successful"])!=int((ledger.status=="SUCCESS").sum()) or int(cfg["failed"])!=int((ledger.status=="FAILED").sum()):
        raise AssertionError("config lifecycle totals disagree with ledger")
    expected_frames=[]
    child_reports=[]
    for row in ledger.to_dict("records"):
        rid=int(row["replicate_id"])
        child=root/"replicates"/f"replicate_{rid:04d}"
        if row["status"]=="SUCCESS":
            child_cfg=json.loads((child/"config.json").read_text())
            expected_seed=int(row["population_seed"])
            if int(child_cfg["seed"])!=expected_seed or int(child_cfg["validation_seed"])!=int(row["validation_seed"]):
                raise AssertionError(f"replicate {rid} child seed/config mismatch")
            if int(child_cfg["n_parent"])!=int(cfg["n_parent"]):
                raise AssertionError(f"replicate {rid} parent-size/config mismatch")
            if list(child_cfg["validation_sizes"])!=list(cfg["validation_sizes"]):
                raise AssertionError(f"replicate {rid} validation-size/config mismatch")
            if list(child_cfg["absolute_bandwidth_candidates"])!=list(cfg["absolute_widths"]):
                raise AssertionError(f"replicate {rid} bandwidth-grid/config mismatch")
            if int(child_cfg["grid_points"])!=int(cfg["grid_points"]) or int(child_cfg["quadrature_order"])!=int(cfg["quadrature_order"]):
                raise AssertionError(f"replicate {rid} numerical-resolution/config mismatch")
            if not np.isclose(float(child_cfg["catalogue_floor"]),float(cfg["catalogue_floor"])):
                raise AssertionError(f"replicate {rid} catalogue-floor/config mismatch")
            report=audit_one(child)
            child_reports.append({"replicate_id":rid,**report})
            frame=pd.read_csv(child/"results"/"density_metrics.csv")
            if len(frame)!=int(row["metric_rows"]):
                raise AssertionError(f"replicate {rid} metric-row ledger mismatch")
            if not (frame.status=="SUCCESS").all():
                raise AssertionError(f"replicate {rid} has non-success metric rows despite SUCCESS status")
            frame["replicate"]=rid
            expected_frames.append(frame)
        else:
            if int(row["metric_rows"])!=0:
                raise AssertionError(f"failed replicate {rid} has nonzero metric count")
            if not row["error"]:
                raise AssertionError(f"failed replicate {rid} has no recorded error")
            if not child.exists():
                raise AssertionError(f"failed replicate {rid} has no preserved attempt directory")
    expected=pd.concat(expected_frames,ignore_index=True) if expected_frames else pd.DataFrame()
    actual=pd.read_csv(root/"density_metrics_all_replicates.csv")
    if len(actual)!=len(expected):
        raise AssertionError("aggregate metric row count mismatch")
    key=["replicate","validation_size","absolute_bandwidth_sd","arm","node"]
    if len(expected) and (expected.duplicated(key).any() or actual.duplicated(key).any()):
        raise AssertionError("duplicate aggregate metric key")
    expected=expected.sort_values(key).reset_index(drop=True)
    actual=actual.sort_values(key).reset_index(drop=True)
    for col in key+["status"]:
        if actual[col].tolist()!=expected[col].tolist():
            raise AssertionError(f"aggregate metric key/status mismatch: {col}")
    for col in ["iae","ise","boundary_iae","kish_ess"]:
        if col in expected and not np.allclose(actual[col].to_numpy(float),expected[col].to_numpy(float),rtol=2e-10,atol=2e-10,equal_nan=False):
            raise AssertionError(f"aggregate metric values mismatch: {col}")
    return {"status":"PASS","manifest_files_verified":checked,
      "planned_replicates":planned,"successful_replicates":int((ledger.status=="SUCCESS").sum()),
      "failed_replicates":int((ledger.status=="FAILED").sum()),
      "child_raw_audits_passed":len(child_reports),"aggregate_metric_rows_reconstructed":len(expected),
      "scope":"Outer seed/ledger/manifest/aggregation audit plus independent raw-input reconstruction of each successful replicate. No bandwidth selection or H0 inference."}

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
    args=p.parse_args()
    report=audit(args.artifact_dir)
    payload=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.parent.mkdir(parents=True,exist_ok=True)
        args.output.write_text(payload)
    print(payload)
if __name__=="__main__": main()
