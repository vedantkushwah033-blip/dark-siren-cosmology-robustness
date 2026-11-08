import hashlib,json
import numpy as np
import pandas as pd
import pytest
from dark_siren_robustness.fixed_analysis_catalogue_benchmark import (
 run_fixed_catalogue_benchmark,write_fixed_catalogue_artifacts)
from dark_siren_robustness.fixed_analysis_catalogue_audit import audit_fixed_catalogue_artifact

def _smoke(tmp_path):
 a,u,s=run_fixed_catalogue_benchmark(n_replicates=2,n_analysis=1200,
  n_validation_max=1000,validation_sizes=(100,300,1000),gamma_values=(0.0,),
  retention_values=(0.5,),degrees=(1,),bandwidth_multipliers=(0.15,),
  grid_points=201,seed=1911)
 write_fixed_catalogue_artifacts(a,u,s,tmp_path)
 return a,u,s

def _refresh_attempt_hash(path):
 mp=path/"manifest.json";m=json.loads(mp.read_text())
 payload=(path/"attempt_metrics.csv").read_bytes()
 m["files"]["attempts"]["bytes"]=len(payload)
 m["files"]["attempts"]["sha256"]=hashlib.sha256(payload).hexdigest()
 mp.write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")

def test_fixed_catalogue_smoke_and_independent_audit(tmp_path):
 a,u,s=_smoke(tmp_path)
 assert len(a)==2*1*1*3*1*(3+1)
 assert s["design"]=="fixed_analysis_catalogue_independent_validation_frame"
 assert set(u.frame)=={"ANALYSIS","VALIDATION"}
 assert audit_fixed_catalogue_artifact(tmp_path)["status"]=="PASS"
 validation=u[(u.replicate==0)&(u.gamma==0.0)&(u.frame=="VALIDATION")].sort_values("validation_rank")
 assert abs(validation.redshift.corr(validation.retention_uniform)) < 0.12
 assert abs(validation.redshift.corr(validation.validation_rank)) < 0.12
 analysis=u[(u.replicate==0)&(u.gamma==0.0)&(u.frame=="ANALYSIS")]
 assert abs(analysis.redshift.corr(analysis.retention_uniform)) < 0.12

def test_analysis_catalogue_is_frozen_across_validation_sizes(tmp_path):
 a,u,s=_smoke(tmp_path)
 for _,g in u.groupby(["replicate","gamma","target_retention"]):
  analysis=g[g.frame=="ANALYSIS"]
  assert len(analysis)==1200
  assert analysis.object_id.nunique()==1200
 # Each attempt must report the same retained analysis-catalogue count for n=100 and n=300.
 for _,g in a.groupby(["replicate","gamma","target_retention","arm","degree","bandwidth_multiplier"]):
  assert g.n_analysis_catalogue.nunique()==1

def test_auditor_rejects_tampered_coefficient_with_refreshed_hash(tmp_path):
 a,u,s=_smoke(tmp_path)
 path=tmp_path/"attempt_metrics.csv";f=pd.read_csv(path)
 ok=f.index[(f.arm=="estimated_degree_1_ipw")&(f.status=="SUCCESS")]
 assert len(ok)>0
 i=ok[0];b=np.asarray(json.loads(f.loc[i,"coefficients_json"]),float);b[0]+=.1
 f.loc[i,"coefficients_json"]=json.dumps(b.tolist())
 f.to_csv(path,index=False);_refresh_attempt_hash(tmp_path)
 with pytest.raises(AssertionError,match="Fitted logistic score mismatch|Fitted logistic objective mismatch"):
  audit_fixed_catalogue_artifact(tmp_path)

def test_auditor_rejects_changed_analysis_membership(tmp_path):
 a,u,s=_smoke(tmp_path)
 path=tmp_path/"frame_units.csv";f=pd.read_csv(path)
 ix=f.index[(f.frame=="ANALYSIS")&(f.retained==0)][0]
 f.loc[ix,"retained"]=1
 f.to_csv(path,index=False)
 mp=tmp_path/"manifest.json";m=json.loads(mp.read_text());payload=path.read_bytes()
 m["files"]["units"]["bytes"]=len(payload);m["files"]["units"]["sha256"]=hashlib.sha256(payload).hexdigest()
 mp.write_text(json.dumps(m,indent=2,sort_keys=True)+"\n")
 with pytest.raises(AssertionError,match="Retention label does not follow"):
  audit_fixed_catalogue_artifact(tmp_path)
