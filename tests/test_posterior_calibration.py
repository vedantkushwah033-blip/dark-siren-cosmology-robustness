import csv, hashlib, importlib.util, json
from pathlib import Path
import pytest

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
RUN=load("posterior_runner","run_posterior_calibration.py")
AUDIT=load("posterior_auditor","audit_posterior_calibration.py")

def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def fixture(root):
    raw,summaries=RUN.run(replicates=100,seed=1771)
    write_csv(root/"per_run.csv",raw); write_csv(root/"summary.csv",summaries)
    meta={"audit_status":"PENDING","base_seed":1771,"replicates_per_truth_noise_cell":100,
          "raw_records":len(raw),"summary_rows":len(summaries),
          "sha256":{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ("per_run.csv","summary.csv")}}
    (root/"run_metadata.json").write_text(json.dumps(meta,indent=2)+"\n")

def test_posterior_calibration_design_and_summary_metrics():
    raw,summaries=RUN.run(replicates=100,seed=1771)
    assert len(raw)==4*3*100
    assert len(summaries)==4*3*2*2*2
    assert {r["true_H0"] for r in raw}==set(RUN.TRUTHS)
    assert all(0<=int(r["covered68"])<=1 and 0<=int(r["covered95"])<=1 for r in raw)
    assert all(0<=r["coverage68"]<=1 and 0<=r["coverage95"]<=1 for r in summaries)

def test_independent_audit_reconstructs_posterior_outputs(tmp_path):
    fixture(tmp_path)
    report=AUDIT.audit(tmp_path)
    assert report["audit"]=="PASS"
    assert report["raw_records"]==1200 and report["summary_rows"]==96

def test_audit_rejects_tampered_summary_after_hash_refresh(tmp_path):
    fixture(tmp_path)
    path=tmp_path/"summary.csv"
    with path.open(newline="",encoding="utf-8") as f: rows=list(csv.DictReader(f))
    rows[0]["rmse_median"]=str(float(rows[0]["rmse_median"])+0.5)
    write_csv(path,rows)
    mp=tmp_path/"run_metadata.json"; meta=json.loads(mp.read_text())
    meta["sha256"]["summary.csv"]=hashlib.sha256(path.read_bytes()).hexdigest()
    mp.write_text(json.dumps(meta,indent=2)+"\n")
    with pytest.raises(ValueError,match="rmse_median mismatch"): AUDIT.audit(tmp_path)
