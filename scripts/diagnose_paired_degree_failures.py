#!/usr/bin/env python3
"""Inspect failed degree-two child runs in a frozen paired-degree artifact."""
from __future__ import annotations
import argparse, json
from pathlib import Path
import pandas as pd

def label(rep: int, gamma: float) -> str:
    return f"rep{rep:03d}_g{gamma:+g}".replace("+","p").replace("-","m")

def main():
    p=argparse.ArgumentParser()
    p.add_argument("artifact_dir",type=Path)
    p.add_argument("--output",type=Path,default=Path("degree2_failure_diagnostics.json"))
    a=p.parse_args()
    root=a.artifact_dir
    d2=root/"degree2"
    ledger=pd.read_csv(d2/"replicate_ledger.csv")
    failed=ledger.loc[ledger.status!="SUCCESS"].copy()
    if failed.empty:
        report={"status":"NO_FAILURES","failed_cells":0}
    else:
        groups=failed.groupby(["replicate","gamma"],sort=True)
        children=[]
        for (rep,gamma),g in groups:
            cell=label(int(rep),float(gamma))
            parent=d2/"replicates"
            stderr=parent/(cell+".stderr.txt")
            audit_stderr=parent/(cell+".audit.stderr.txt")
            audit_stdout=parent/(cell+".audit.stdout.txt")
            def read(path):
                return path.read_text(errors="replace")[-12000:] if path.exists() else None
            children.append({
                "replicate":int(rep),"gamma":float(gamma),"failed_cells":int(len(g)),
                "status_values":sorted(g.status.astype(str).unique().tolist()),
                "audit_status_values":sorted(g.audit_status.astype(str).unique().tolist()),
                "ledger_error":str(g.error.iloc[0])[:5000],
                "runner_stderr_tail":read(stderr),
                "child_audit_stderr_tail":read(audit_stderr),
                "child_audit_stdout_tail":read(audit_stdout),
            })
        counts=failed.groupby(["gamma","error"],dropna=False).size().reset_index(name="failed_cells")
        report={
            "status":"DIAGNOSTIC_ONLY","failed_cells":int(len(failed)),
            "planned_cells":int(len(ledger)),
            "affected_replicate_gamma_runs":len(children),
            "failed_cells_by_gamma":{str(float(k)):int(v) for k,v in failed.groupby("gamma").size().items()},
            "error_groups":[{"gamma":float(r.gamma),"failed_cells":int(r.failed_cells),
                             "error":str(r.error)[:5000]} for r in counts.itertuples()],
            "children":children,
            "interpretation":"Frozen artifact diagnosis only; no incomplete result is promoted to a scientific result."
        }
    a.output.parent.mkdir(parents=True,exist_ok=True)
    a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    print(json.dumps(report,indent=2,sort_keys=True))

if __name__=="__main__": main()
