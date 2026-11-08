"""Independent audit for fixed-analysis-catalogue artifacts; no production runner imports."""
from __future__ import annotations
import argparse,hashlib,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid,trapezoid
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import gaussian_kde

Z0,Z1=0.01,0.40
def _lam(z):
 z=np.asarray(z,float)
 return 1+2.2*np.exp(-.5*((z-.105)/.025)**2)+1.6*np.exp(-.5*((z-.275)/.040)**2)
def _center():
 z=np.linspace(Z0,Z1,20001);x=(z-Z0)/(Z1-Z0)
 return float(trapezoid(_lam(z)*x*x,z)/trapezoid(_lam(z),z))
def _alpha(g,r,c):
 z=np.linspace(Z0,Z1,20001);x=(z-Z0)/(Z1-Z0);l=_lam(z);den=trapezoid(l,z)
 return float(brentq(lambda a:trapezoid(l*expit(a-x+g*(x*x-c)),z)/den-r,-30,30))
def _target(z):
 m=1+1.5*(z-Z0)/(Z1-Z0);v=_lam(z)*m
 return v/trapezoid(v,z)
def _kde(grid,sample,prob,h,merger):
 sample=np.asarray(sample,float)
 if sample.size<2: raise AssertionError("Insufficient sample for successful KDE")
 w=np.ones(sample.size) if prob is None else 1/np.asarray(prob,float)
 w=w/w.max();neff=w.sum()**2/(w@w)
 aug=np.r_[sample,2*grid[0]-sample,2*grid[-1]-sample]
 aw=np.tile(w/3,3)
 factor=h*neff**(-.2)
 raw=gaussian_kde(aug,weights=aw,bw_method=lambda k:factor)(grid)*merger
 return raw/trapezoid(raw,grid)
def audit_fixed_catalogue_artifact(path):
 path=Path(path)
 a=pd.read_csv(path/"attempt_metrics.csv");u=pd.read_csv(path/"frame_units.csv")
 s=json.loads((path/"summary.json").read_text());c=json.loads((path/"config.json").read_text())
 for k,v in c.items():
  if s.get(k)!=v: raise AssertionError(f"Config-summary mismatch: {k}")
 m=json.loads((path/"manifest.json").read_text())
 for rec in m["files"].values():
  p=path/rec["path"];b=p.read_bytes()
  if len(b)!=rec["bytes"]: raise AssertionError("Manifest byte count mismatch")
  if hashlib.sha256(b).hexdigest()!=rec["sha256"]: raise AssertionError("Manifest checksum mismatch")
 if s["design"]!="fixed_analysis_catalogue_independent_validation_frame":
  raise AssertionError("Wrong design identifier")
 key=["replicate","gamma","target_retention","frame","object_id"]
 if u.duplicated(key).any(): raise AssertionError("Duplicate frame object key")
 attempt_key=["replicate","gamma","target_retention","n_validation","bandwidth_multiplier","degree","arm"]
 if a.duplicated(attempt_key).any(): raise AssertionError("Duplicate attempt key")
 center=_center()
 if not np.isclose(center,s["curvature_center"],rtol=1e-12,atol=1e-12): raise AssertionError("Curvature center mismatch")
 grid=np.linspace(*s["grid"]);truth=_target(grid);merger=1+1.5*(grid-Z0)/(Z1-Z0)
 if not np.allclose(truth,s["analytic_target_density"],rtol=1e-10,atol=1e-12): raise AssertionError("Analytic target mismatch")
 groups={}
 for keyval,g in u.groupby(["replicate","gamma","target_retention"],sort=False):
  af=g[g.frame=="ANALYSIS"].sort_values("object_id")
  vf=g[g.frame=="VALIDATION"].sort_values("object_id")
  if len(af)!=s["n_analysis_parent"] or len(vf)!=s["n_validation_max"]: raise AssertionError("Frame size mismatch")
  if af.object_id.nunique()!=len(af) or vf.object_id.nunique()!=len(vf): raise AssertionError("Frame IDs not unique")
  # object_id is frame-local; the frame column provides the namespace.\n  for frame in (af,vf):
   z=frame.redshift.to_numpy(float);y=frame.retained.to_numpy(int);uu=frame.retention_uniform.to_numpy(float)
   if np.any((uu<0)|(uu>=1)) or np.any((y!=0)&(y!=1)): raise AssertionError("Invalid retention records")
   gamma=float(keyval[1]);ret=float(keyval[2]);alpha=_alpha(gamma,ret,center)
   x=(z-Z0)/(Z1-Z0);pi=expit(alpha-x+gamma*(x*x-center))
   if not np.allclose(pi,frame.true_inclusion_probability,rtol=1e-11,atol=1e-12): raise AssertionError("DGP probability mismatch")
   if not np.array_equal((uu<pi).astype(int),y): raise AssertionError("Retention label does not follow stored uniform/DGP")
  if vf.validation_rank.min()<0 or vf.validation_rank.max()>=len(vf) or vf.validation_rank.nunique()!=len(vf): raise AssertionError("Invalid validation rank permutation")
  groups[keyval]=(af,vf)
 # Exact planned factorial
 expected=set()
 for rep in range(s["n_replicates"]):
  for g in s["gamma_values"]:
   for ret in s["retention_values"]:
    for n in s["validation_sizes"]:
     for h in s["bandwidth_multipliers"]:
      for arm in ("catalogue_only","known_inclusion_ipw","parent_sample_kde"):
       expected.add((rep,float(g),float(ret),int(n),float(h),0,arm))
      for d in s["degrees"]:
       expected.add((rep,float(g),float(ret),int(n),float(h),int(d),f"estimated_degree_{d}_ipw"))
 actual={(int(x.replicate),float(x.gamma),float(x.target_retention),int(x.n_validation),
          float(x.bandwidth_multiplier),int(x.degree),str(x.arm)) for x in a.itertuples()}
 if actual!=expected: raise AssertionError(f"Attempt design mismatch: missing={len(expected-actual)}, extra={len(actual-expected)}")
 checked=0
 for row in a.itertuples():
  k=(int(row.replicate),float(row.gamma),float(row.target_retention))
  af,vf=groups[k];az=af.redshift.to_numpy(float);ay=af.retained.to_numpy(int)
  ac=af[af.retained==1];cat=ac.redshift.to_numpy(float);api=ac.true_inclusion_probability.to_numpy(float)
  # Same fixed catalogue is implied by one immutable analysis frame, irrespective of n.
  vv=vf[vf.validation_rank<int(row.n_validation)].sort_values("validation_rank")
  if len(vv)!=row.n_validation: raise AssertionError("Nested validation sample size mismatch")
  sample=None;prob=None
  if row.arm=="catalogue_only": sample=cat
  elif row.arm=="known_inclusion_ipw": sample=cat;prob=api
  elif row.arm=="parent_sample_kde": sample=az
  elif row.arm.startswith("estimated_degree_"):
   sample=cat;d=int(row.degree);beta=np.asarray(json.loads(row.coefficients_json),float)
   if row.fit_status!="SUCCESS" or beta.shape!=(d+1,): raise AssertionError("Successful fit record malformed")
   vz=vv.redshift.to_numpy(float);vy=vv.retained.to_numpy(float)
   cfit=float(row.fit_center);sfit=float(row.fit_scale)
   if not np.isclose(cfit,vz.mean(),rtol=1e-12,atol=1e-12): raise AssertionError("Fit center mismatch")
   sd=float(vz.std());sd=sd if sd>0 else 1.
   if not np.isclose(sfit,sd,rtol=1e-12,atol=1e-12): raise AssertionError("Fit scale mismatch")
   x=(vz-cfit)/sfit;X=np.column_stack([x**j for j in range(d+1)])
   p=expit(X@beta);score=X.T@(p-vy)
   if np.linalg.norm(score,np.inf)>1e-6 or not np.isclose(np.linalg.norm(score,np.inf),row.fit_gradient_inf_norm,rtol=2e-5,atol=2e-8):
    raise AssertionError("Fitted logistic score mismatch")
   obj=float(np.sum(np.logaddexp(0,X@beta)-vy*(X@beta)))
   if not np.isclose(obj,row.fit_objective,rtol=2e-8,atol=2e-8): raise AssertionError("Fitted logistic objective mismatch")
   xa=(cat-cfit)/sfit;prob=expit(np.column_stack([xa**j for j in range(d+1)])@beta)
  else: raise AssertionError("Unknown arm")
  if row.status!="SUCCESS":
   if row.density_json and isinstance(row.density_json,str): raise AssertionError("Failed attempt contains density")
   continue
  density=np.asarray(json.loads(row.density_json),float)
  rebuilt=_kde(grid,sample,prob,float(row.bandwidth_multiplier),merger)
  if not np.allclose(density,rebuilt,rtol=2e-8,atol=2e-8): raise AssertionError("KDE reconstruction mismatch")
  iae=float(trapezoid(abs(density-truth),grid))
  if not np.isclose(iae,row.iae,rtol=2e-8,atol=2e-8): raise AssertionError("IAE mismatch")
  if int(row.n_analysis_catalogue)!=len(cat): raise AssertionError("Analysis catalogue count changed across n")
  checked+=1
 counts={str(k):int(v) for k,v in a.status.value_counts().items()}
 if counts!=s["status_counts"]: raise AssertionError("Status ledger mismatch")
 if not set(a.status).issubset({"SUCCESS","FAILED","NON_EVALUABLE","NOT_RUN"}):
  raise AssertionError("Unknown attempt status")
 expected_ledger=[]
 for key,g in a.groupby(["gamma","target_retention","n_validation","bandwidth_multiplier","arm"],sort=True):
  expected_ledger.append({"gamma":float(key[0]),"target_retention":float(key[1]),
   "n_validation":int(key[2]),"bandwidth_multiplier":float(key[3]),"arm":str(key[4]),
   "n_planned":int(len(g)),"n_success":int((g.status=="SUCCESS").sum()),
   "n_failed":int((g.status=="FAILED").sum()),
   "n_non_evaluable":int((g.status=="NON_EVALUABLE").sum()),
   "n_not_run":int((g.status=="NOT_RUN").sum()),
   "success_rate":float((g.status=="SUCCESS").mean()),
   "conditional_success_mean_iae":float(g.loc[g.status=="SUCCESS","iae"].mean())
     if (g.status=="SUCCESS").any() else None})
 if len(expected_ledger)!=len(s["condition_status_ledger"]):
  raise AssertionError("Condition failure-ledger size mismatch")
 for x,y in zip(expected_ledger,s["condition_status_ledger"]):
  for k,v in x.items():
   w=y[k]
   if v is None:
    if w is not None: raise AssertionError("Condition failure-ledger null mismatch")
   elif isinstance(v,float):
    if not np.isclose(v,w,rtol=1e-10,atol=1e-12): raise AssertionError("Condition failure-ledger value mismatch")
   elif v!=w: raise AssertionError("Condition failure-ledger value mismatch")
 expected_contrasts=[]
 for key,g in a.groupby(["gamma","target_retention","n_validation","bandwidth_multiplier"],sort=True):
  wide=g.pivot(index="replicate",columns="arm",values="iae")
  if "catalogue_only" not in wide: continue
  for degree in s["degrees"]:
   arm=f"estimated_degree_{degree}_ipw"
   if arm not in wide: continue
   pair=wide[["catalogue_only",arm]].dropna()
   diff=(pair[arm]-pair["catalogue_only"]).to_numpy(float)
   expected_contrasts.append({"gamma":float(key[0]),"target_retention":float(key[1]),
    "n_validation":int(key[2]),"bandwidth_multiplier":float(key[3]),"degree":int(degree),
    "n_pairs":int(len(diff)),"n_replicates_planned":int(s["n_replicates"]),
    "mean_iae_difference":float(diff.mean()) if len(diff) else None,
    "paired_mc_se":float(diff.std(ddof=1)/np.sqrt(len(diff))) if len(diff)>1 else None})
 if len(expected_contrasts)!=len(s["paired_ipw_minus_catalogue_only"]):
  raise AssertionError("Paired contrast count mismatch")
 for x,y in zip(expected_contrasts,s["paired_ipw_minus_catalogue_only"]):
  for k,v in x.items():
   w=y[k]
   if v is None:
    if w is not None: raise AssertionError("Paired contrast null mismatch")
   elif isinstance(v,float):
    if not np.isclose(v,w,rtol=1e-10,atol=1e-12): raise AssertionError("Paired contrast mismatch")
   elif v!=w: raise AssertionError("Paired contrast mismatch")
 rates={str(k):float((v.status=="SUCCESS").mean()) for k,v in a.groupby("arm")}
 if rates.keys()!=s["success_rate_by_arm"].keys(): raise AssertionError("Success-rate keys mismatch")
 for k,v in rates.items():
  if not np.isclose(v,s["success_rate_by_arm"][k]): raise AssertionError("Success rate mismatch")
 return {"status":"PASS","design":s["design"],"checked_frame_scenarios":len(groups),
  "checked_density_estimates":checked,"attempt_count":len(a),"status_counts":counts,
  "limitations":["Synthetic density recovery only; no H0 posterior calibration.",
   "Shares SciPy gaussian_kde primitive with production KDE implementation.",
   "Independent validation frame assumes complete binary retention-label ascertainment."]}
def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--artifact-dir",type=Path,required=True)
 p.add_argument("--output",type=Path);x=p.parse_args();r=audit_fixed_catalogue_artifact(x.artifact_dir)
 out=json.dumps(r,indent=2,sort_keys=True)+"\n"
 if x.output:x.output.write_text(out)
 print(out)
if __name__=="__main__":main()
