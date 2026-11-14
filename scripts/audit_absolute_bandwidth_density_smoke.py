"""Independent raw-input auditor for absolute-bandwidth density smoke."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde
from scipy.special import expit, logsumexp, roots_hermitenorm
from scipy.optimize import minimize
from itertools import product
import math

Z_MIN,Z_MAX=0.01,0.40
BOUNDARY_LOWER,BOUNDARY_UPPER=0.049,0.361

def _hash(p):
 h=hashlib.sha256()
 with Path(p).open("rb") as f:
  for b in iter(lambda:f.read(1<<20),b""): h.update(b)
 return h.hexdigest()

def _close(a,b,label,atol=2e-10):
 if not np.allclose(a,b,rtol=2e-10,atol=atol,equal_nan=False):
  raise AssertionError(f"{label} mismatch; max abs diff={np.max(np.abs(np.asarray(a)-np.asarray(b)))}")

def _metric(d,t,z):
 x=np.asarray(d)-t; left=z<=BOUNDARY_LOWER; right=z>=BOUNDARY_UPPER
 return {"iae":float(trapezoid(np.abs(x),z)),
  "ise":float(trapezoid(x*x,z)),
  "boundary_iae":float(trapezoid(np.abs(x[left]),z[left])+trapezoid(np.abs(x[right]),z[right]))}

def _kish_ess(weights):
 w=np.asarray(weights,float)
 if w.ndim!=1 or w.size==0 or np.any(w<=0) or not np.all(np.isfinite(w)):
  raise AssertionError("invalid weights for independent Kish ESS")
 scaled=w/np.max(w)
 if np.any(scaled==0):
  raise AssertionError("positive weights underflowed during independent Kish ESS")
 p=scaled/np.sum(scaled)
 return float(1.0/np.sum(p*p))


def _independent_weight_diagnostics(weights):
 w=np.asarray(weights,float)
 if w.ndim!=1 or w.size==0 or not np.all(np.isfinite(w)) or np.any(w<=0):
  raise AssertionError("invalid raw weights for independent concentration diagnostics")
 scaled=w/np.max(w)
 if np.any(scaled==0):
  raise AssertionError("positive mass underflowed in independent concentration diagnostics")
 normalized=scaled/np.sum(scaled)
 ordered=np.sort(normalized)[::-1]
 k1=max(1,int(np.ceil(0.01*w.size)))
 k5=max(1,int(np.ceil(0.05*w.size)))
 ess=float(1.0/np.sum(normalized**2))
 positive=normalized[normalized>0]
 entropy_ess=float(np.exp(-np.sum(positive*np.log(positive))))
 return {
  "kish_ess":ess,
  "kish_ess_fraction":ess/w.size,
  "entropy_ess":entropy_ess,
  "maximum_normalized_weight":float(ordered[0]),
  "top_1_percent_mass":float(np.sum(ordered[:k1])),
  "top_5_percent_mass":float(np.sum(ordered[:k5])),
 }

def _independent_kde(z,sample,weights,rate,width):
 w=np.asarray(weights,float)
 if w.shape!=sample.shape or np.any(w<=0) or not np.all(np.isfinite(w)): raise AssertionError("invalid raw weights")
 scaled=w/np.max(w)
 if np.any(scaled==0): raise AssertionError("positive weights underflowed")
 reflected=np.r_[sample,2*z[0]-sample,2*z[-1]-sample]
 augw=np.tile(scaled/3.0,3)
 p=augw/np.sum(augw)
 den=1.0-float(p@p)
 mu=float(p@reflected)
 var=float(np.sum(p*(reflected-mu)**2)/den)
 if not np.isfinite(var) or var<=0: raise AssertionError("reflected weighted covariance invalid")
 factor=float(width)/np.sqrt(var)
 # Independent direct check that SciPy's covariance SD equals the requested width.
 kde=gaussian_kde(reflected,weights=augw,bw_method=factor)
 realized=float(np.sqrt(kde.covariance[0,0]))
 if not np.isclose(realized,width,rtol=2e-10,atol=2e-12):
  raise AssertionError("realized physical kernel SD differs from requested width")
 raw=np.asarray(kde(z))*rate
 area=float(trapezoid(raw,z))
 if not np.isfinite(area) or area<=0: raise AssertionError("invalid KDE area")
 d=raw/area
 if not np.isclose(trapezoid(d,z),1.0,rtol=1e-10,atol=1e-10): raise AssertionError("KDE density does not integrate to one")
 return d

def _independent_scott_width(sample,weights,multiplier):
 w=np.asarray(weights,float)
 sample=np.asarray(sample,float)
 if sample.ndim!=1 or w.shape!=sample.shape or np.any(w<=0) or not np.all(np.isfinite(w)):
  raise AssertionError("invalid raw inputs for independent Scott width")
 scaled=w/np.max(w)
 if np.any(scaled==0): raise AssertionError("Scott width weights underflowed")
 probability=scaled/np.sum(scaled)
 denominator=1.0-float(probability@probability)
 mean=float(probability@sample)
 variance=float(np.sum(probability*(sample-mean)**2)/denominator)
 neff=float(np.sum(scaled)**2/np.sum(scaled**2))
 width=float(multiplier)*neff**(-1.0/5.0)*np.sqrt(variance)
 if not np.isfinite(width) or width<=0: raise AssertionError("independent Scott width is invalid")
 return float(width)

def _independent_logistic_quadrature(validation_redshift, retained, cfg):
 """Rebuild the proper degree-one logistic posterior without production code.

 The reference uses BFGS plus an independently coded Newton polish, then
 reconstructs the Laplace-proposal tensor Gauss-Hermite importance weights.
 """
 z=np.asarray(validation_redshift,float); y=np.asarray(retained,float)
 if z.ndim!=1 or y.shape!=z.shape or z.size<5 or not np.all(np.isfinite(z)):
  raise AssertionError("invalid validation inputs for independent posterior")
 if np.any((y!=0)&(y!=1)): raise AssertionError("validation labels are not binary")
 degree=1
 center=float(cfg["coordinate_center"]); scale=float(cfg["coordinate_scale"])
 if not np.isfinite(center) or not np.isfinite(scale) or scale<=0:
  raise AssertionError("invalid coefficient coordinate system")
 x=(z-center)/scale
 design=np.column_stack((np.ones(z.size),x))
 prior=cfg["coefficient_prior"]
 if prior.get("family")!="independent_gaussian":
  raise AssertionError("unsupported coefficient prior family")
 mu=np.asarray(prior["mean"],float); sd=np.asarray(prior["sd"],float)
 if mu.shape!=(2,) or sd.shape!=(2,) or not np.all(np.isfinite(mu)) or np.any(sd<=0):
  raise AssertionError("invalid persisted coefficient prior")
 precision=1.0/(sd*sd)
 def objective(beta):
  eta=design@beta
  loglik=np.sum(y*eta-np.logaddexp(0.0,eta))
  standardized=(beta-mu)/sd
  logprior=-0.5*float(standardized@standardized)-float(np.sum(np.log(sd)))-math.log(2.0*math.pi)
  return float(-(loglik+logprior))
 def gradient(beta):
  return design.T@(expit(design@beta)-y)+(beta-mu)*precision
 def hessian(beta):
  p=expit(design@beta)
  return design.T@((p*(1-p))[:,None]*design)+np.diag(precision)
 fit=minimize(objective,mu.copy(),jac=gradient,method="BFGS",
              options={"gtol":1e-10,"maxiter":4000})
 mode=np.asarray(fit.x,float)
 if not np.all(np.isfinite(mode)): raise AssertionError("independent posterior mode is non-finite")
 # Deterministic damped Newton polishing, independent of the production trust-exact optimizer.
 for _ in range(100):
  score=gradient(mode)
  if np.linalg.norm(score,ord=np.inf)<=1e-6: break
  info=hessian(mode)
  try: direction=np.linalg.solve(info,score)
  except np.linalg.LinAlgError as exc: raise AssertionError("independent posterior Hessian solve failed") from exc
  decrement=float(score@direction)
  if not np.isfinite(decrement) or decrement<=0: raise AssertionError("invalid independent Newton decrement")
  current=objective(mode); step=1.0; accepted=False
  while step>=2.0**-35:
   candidate=mode-step*direction
   candidate_score=float(np.linalg.norm(gradient(candidate),ord=np.inf))
   candidate_objective=objective(candidate)
   # Near floating-point objective resolution, Armijo decrease can disappear
   # while the posterior score still improves. Accept either a valid Armijo
   # step or a strict reduction in the independently evaluated score.
   armijo=candidate_objective<=current-1e-4*step*decrement
   score_improved=np.isfinite(candidate_score) and candidate_score < float(np.linalg.norm(score,ord=np.inf))
   if np.isfinite(candidate_score) and (armijo or score_improved):
    mode=candidate; accepted=True; break
   step*=0.5
  if not accepted:
   if np.linalg.norm(gradient(mode),ord=np.inf)<=1e-6: break
   raise AssertionError("independent posterior Newton polish stalled")
 else:
  if np.linalg.norm(gradient(mode),ord=np.inf)>1e-6:
   raise AssertionError("independent posterior mode did not converge")
 information=hessian(mode)
 eigen=np.linalg.eigvalsh(information)
 if np.min(eigen)<=0 or not np.all(np.isfinite(eigen)):
  raise AssertionError("independent posterior information is not positive definite")
 covariance=np.linalg.inv(information)
 chol=np.linalg.cholesky(covariance)
 order=int(cfg["quadrature_order"])
 if order<3: raise AssertionError("invalid posterior quadrature order")
 gh_nodes,gh_weights=roots_hermitenorm(order)
 tuples=list(product(range(order),repeat=degree+1))
 standard=np.asarray([[gh_nodes[i] for i in ix] for ix in tuples],float)
 base_log=np.asarray([sum(math.log(gh_weights[i]) for i in ix)
                      -0.5*(degree+1)*math.log(2.0*math.pi) for ix in tuples])
 nodes=mode[None,:]+standard@chol.T
 log_proposal=(-0.5*(degree+1)*math.log(2.0*math.pi)
               -0.5*np.sum(standard*standard,axis=1)
               -float(np.sum(np.log(np.diag(chol)))))
 eta=design@nodes.T
 log_likelihood=np.sum(y[:,None]*eta-np.logaddexp(0.0,eta),axis=0)
 standardized_prior=(nodes-mu[None,:])/sd[None,:]
 log_prior=(-0.5*np.sum(standardized_prior**2,axis=1)
            -float(np.sum(np.log(sd)))-(degree+1)*0.5*math.log(2.0*math.pi))
 log_importance=base_log+log_likelihood+log_prior-log_proposal
 log_norm=float(logsumexp(log_importance))
 q=np.exp(log_importance-log_norm)
 q/=q.sum()
 return mode,nodes,q


def audit(root:Path):
 root=root.resolve()
 manifest=json.loads((root/"manifest.json").read_text())
 checked=0
 for rel,m in manifest["files"].items():
  p=root/rel
  if not p.is_file() or p.stat().st_size!=int(m["bytes"]) or _hash(p)!=m["sha256"]:
   raise AssertionError(f"manifest verification failed: {rel}")
  checked+=1
 if checked!=int(manifest["file_count"]): raise AssertionError("manifest file count mismatch")
 cfg=json.loads((root/"config.json").read_text())
 if cfg.get("status")!="ENGINEERING_SMOKE_NOT_SELECTION_NOT_CONFIRMATORY": raise AssertionError("unexpected claim status")
 parent=pd.read_csv(root/"truth"/"parent_truth.csv")
 cat=pd.read_csv(root/"observable"/"catalogue_observed.csv")
 vf=pd.read_csv(root/"truth"/"validation_frame_truth.csv")
 vfull=pd.read_csv(root/"observable"/"validation_max.csv")
 if len(parent)!=int(cfg["n_parent"]) or len(cat)!=int(cfg["n_catalogue"]): raise AssertionError("population counts mismatch")
 if set(cat.object_id.astype(int))!=set(parent.loc[parent.retained==1,"object_id"].astype(int)): raise AssertionError("catalogue does not match parent retention")
 zparent=parent.redshift.to_numpy(float); zcat=cat.redshift.to_numpy(float)
 _close(zcat,parent.set_index("object_id").loc[cat.object_id,"redshift"].to_numpy(float),"catalogue-parent redshift join")
 floor=float(cfg["catalogue_floor"]); midpoint=float(cfg["catalogue_midpoint"]); width=float(cfg["catalogue_width"])
 ptrue=floor+(1-floor)/(1+np.exp((zcat-midpoint)/width))
 _close(parent.inclusion_probability.to_numpy(float),floor+(1-floor)/(1+np.exp((zparent-midpoint)/width)),"parent inclusion probabilities")
 bp=1+1.5*(zparent-Z_MIN)/(Z_MAX-Z_MIN); bc=1+1.5*(zcat-Z_MIN)/(Z_MAX-Z_MIN)
 z=np.load(root/"results"/"density_reconstruction_inputs.npz")
 grid=z["redshift_grid"]; rate=z["merger_weight_grid"]; truth=z["truth_density"]
 _close(z["parent_redshift"],zparent,"persisted parent sample"); _close(z["catalogue_redshift"],zcat,"persisted catalogue sample")
 _close(truth, np.asarray(1+2.2*np.exp(-.5*((grid-.105)/.025)**2)+1.6*np.exp(-.5*((grid-.275)/.04)**2))*rate /
   trapezoid((1+2.2*np.exp(-.5*((grid-.105)/.025)**2)+1.6*np.exp(-.5*((grid-.275)/.04)**2))*rate,grid),"analytic target")
 widths=[float(x) for x in cfg["absolute_bandwidth_candidates"]]
 sizes=[int(x) for x in cfg["validation_sizes"]]
 scott_multiplier=float(cfg["scott_multiplier"])
 if not np.isfinite(scott_multiplier) or scott_multiplier<=0: raise AssertionError("invalid Scott multiplier")
 metrics=pd.read_csv(root/"results"/"density_metrics.csv")
 expected=[]; expected_scott=[]; expected_diagnostics=[]; maxdiff=0.0
 for n in sizes:
  v=pd.read_csv(root/"observable"/f"validation_n{n}.csv")
  if "validation_id" not in v or "validation_id" not in vfull or "frame_id" not in vf:
   raise AssertionError("validation identifiers are missing")
  ids=v.validation_id.to_numpy()
  full_ids=vfull.validation_id.to_numpy()[:n]
  if len(np.unique(ids))!=n or not np.array_equal(ids,full_ids):
   raise AssertionError(f"validation IDs are not the declared nested prefix n={n}")
  if vf.frame_id.duplicated().any():
   raise AssertionError("truth validation-frame IDs are not unique")
  _close(v.redshift.to_numpy(float),vfull.iloc[:n].redshift.to_numpy(float),f"nested validation redshift n={n}")
  _close(v.retained.to_numpy(int),vfull.iloc[:n].retained.to_numpy(int),f"nested validation outcomes n={n}",atol=0)
  truth_join=v.merge(vf[["frame_id","redshift","retained"]],left_on="validation_id",right_on="frame_id",suffixes=("_obs","_truth"),validate="one_to_one")
  if len(truth_join)!=n or truth_join.frame_id.isna().any():
   raise AssertionError(f"validation ID-to-truth join is incomplete n={n}")
  _close(truth_join.redshift_obs.to_numpy(float),truth_join.redshift_truth.to_numpy(float),f"validation frame join n={n}")
  _close(truth_join.retained_obs.to_numpy(int),truth_join.retained_truth.to_numpy(int),f"validation retention join n={n}",atol=0)
  nodes=z[f"posterior_nodes_n{n}"]; q=z[f"posterior_weights_n{n}"]; mode=z[f"posterior_mode_n{n}"]
  if nodes.ndim!=2 or nodes.shape[1]!=2 or q.shape!=(len(nodes),) or np.any(q<0) or not np.isclose(q.sum(),1): raise AssertionError("invalid persisted coefficient quadrature")
  ref_mode,ref_nodes,ref_q=_independent_logistic_quadrature(v.redshift.to_numpy(float),v.retained.to_numpy(int),cfg)
  _close(mode,ref_mode,"independently reconstructed inclusion posterior mode",atol=1e-7)
  _close(nodes,ref_nodes,"independently reconstructed inclusion posterior nodes",atol=1e-7)
  _close(q,ref_q,"independently reconstructed inclusion posterior weights",atol=1e-7)
  xv=(v.redshift.to_numpy(float)-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])
  design=np.column_stack([np.ones(n),xv])
  pmap=expit(np.column_stack([np.ones(len(zcat)),(zcat-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])])@mode)
  pnodes=expit(np.column_stack([np.ones(len(zcat)),(zcat-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])])@nodes.T).T
  _close(z[f"weights_A_catalogue_only_n{n}"],bc,"A base weights")
  _close(z[f"weights_B_parent_oracle_n{n}"],bp,"B base weights")
  _close(z[f"weights_C_known_probability_ipw_n{n}"],bc/ptrue,"C raw weights")
  _close(z[f"weights_D_map_plugin_ipw_n{n}"],bc/pmap,"D raw weights")
  _close(z[f"weights_E_nodes_n{n}"],bc[None,:]/pnodes,"E raw weights")
  diagnostic_inputs={
   "A_catalogue_only":bc,
   "B_parent_oracle":bp,
   "C_known_probability_ipw":bc/ptrue,
   "D_map_plugin_ipw":bc/pmap,
  }
  for diagnostic_arm,diagnostic_weights in diagnostic_inputs.items():
   expected_diagnostics.append({
    "validation_size":n,"arm":diagnostic_arm,"node":-1,
    **_independent_weight_diagnostics(diagnostic_weights),
   })
  e_node_diagnostics=[
   _independent_weight_diagnostics(row) for row in (bc[None,:]/pnodes)
  ]
  expected_diagnostics.append({
   "validation_size":n,"arm":"E_uncertainty_integrated_ipw","node":-1,
   "kish_ess":min(x["kish_ess"] for x in e_node_diagnostics),
   "kish_ess_fraction":min(x["kish_ess_fraction"] for x in e_node_diagnostics),
   "entropy_ess":min(x["entropy_ess"] for x in e_node_diagnostics),
   "maximum_normalized_weight":max(x["maximum_normalized_weight"] for x in e_node_diagnostics),
   "top_1_percent_mass":max(x["top_1_percent_mass"] for x in e_node_diagnostics),
   "top_5_percent_mass":max(x["top_5_percent_mass"] for x in e_node_diagnostics),
  })
  _close(z[f"weights_E_quadrature_n{n}"],q,"E quadrature weights")
  inputs={"A_catalogue_only":(zcat,bc),"B_parent_oracle":(zparent,bp),
    "C_known_probability_ipw":(zcat,bc/ptrue),"D_map_plugin_ipw":(zcat,bc/pmap)}
  for arm,(sample,raw_weights) in inputs.items():
   expected_scott.append({"validation_size":n,"arm":arm,"node":-1,
     "realized_bandwidth_redshift":_independent_scott_width(sample,raw_weights,scott_multiplier)})
  for node,raw_weights in enumerate(bc[None,:]/pnodes):
   expected_scott.append({"validation_size":n,"arm":"E_uncertainty_integrated_ipw","node":node,
     "realized_bandwidth_redshift":_independent_scott_width(zcat,raw_weights,scott_multiplier)})
  for bw in widths:
   node_ds=[]
   for arm,(sample,w) in inputs.items():
    d=_independent_kde(grid,sample,w,rate,bw)
    stored=z[f"density_{arm}_n{n}_w{bw:.6f}"]
    _close(stored,d,f"{arm} density n={n} width={bw}")
    met=_metric(d,truth,grid)
    expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":arm,"node":-1,**met,"status":"SUCCESS","kish_ess":_kish_ess(w)})
   for j,w in enumerate(bc[None,:]/pnodes):
    d=_independent_kde(grid,zcat,w,rate,bw); node_ds.append(d)
    _close(z[f"density_E_nodes_n{n}_w{bw:.6f}"][j],d,f"E node density n={n} width={bw} node={j}")
    expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":"E_uncertainty_integrated_ipw","node":j,**_metric(d,truth,grid),"status":"SUCCESS","kish_ess":_kish_ess(w)})
   node_ds=np.asarray(node_ds)
   qmetric={k:float(np.dot(q,[_metric(d,truth,grid)[k] for d in node_ds])) for k in ("iae","ise","boundary_iae")}
   _close(z[f"density_E_posterior_mean_n{n}_w{bw:.6f}"],np.tensordot(q,node_ds,axes=(0,0)),f"E posterior mean density n={n} width={bw}")
   mean_density=np.tensordot(q,node_ds,axes=(0,0))
   mean_metric=_metric(mean_density,truth,grid)
   expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":"E_posterior_mean_density","node":-1,**mean_metric,"status":"SUCCESS","kish_ess":float("nan")})
   expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":"E_posterior_weighted_node_metric","node":-1,**qmetric,"status":"SUCCESS","kish_ess":float("nan")})
 expected_df=pd.DataFrame(expected)
 key=["replicate","validation_size","absolute_bandwidth_sd","arm","node"]
 if metrics.duplicated(key).any() or expected_df.duplicated(key).any(): raise AssertionError("duplicate metric key")
 got=metrics.sort_values(key).reset_index(drop=True); want=expected_df.sort_values(key).reset_index(drop=True)
 if len(got)!=len(want): raise AssertionError("metric row count mismatch")
 for col in key+["status"]:
  if not got[col].equals(want[col]): raise AssertionError(f"metric key/status mismatch: {col}")
 for col in ["iae","ise","boundary_iae"]:
  _close(got[col].to_numpy(float),want[col].to_numpy(float),f"metric {col}")
 if "kish_ess" not in got or not np.array_equal(got.kish_ess.isna().to_numpy(),want.kish_ess.isna().to_numpy()):
  raise AssertionError("metric Kish ESS missingness mismatch")
 ess_got=got.kish_ess.to_numpy(float)
 ess_want=want.kish_ess.to_numpy(float)
 finite=np.isfinite(ess_want)
 _close(ess_got[finite],ess_want[finite],"metric Kish ESS")
 scott=pd.read_csv(root/"results"/"scott_widths.csv")
 expected_scott_df=pd.DataFrame(expected_scott)
 scott_key=["validation_size","arm","node"]
 if scott.duplicated(scott_key).any() or expected_scott_df.duplicated(scott_key).any():
  raise AssertionError("duplicate realized Scott width key")
 scott=scott.sort_values(scott_key).reset_index(drop=True)
 expected_scott_df=expected_scott_df.sort_values(scott_key).reset_index(drop=True)
 if len(scott)!=len(expected_scott_df): raise AssertionError("realized Scott width row count mismatch")
 for col in scott_key:
  if scott[col].tolist()!=expected_scott_df[col].tolist(): raise AssertionError(f"realized Scott width key mismatch: {col}")
 _close(scott.realized_bandwidth_redshift.to_numpy(float),
        expected_scott_df.realized_bandwidth_redshift.to_numpy(float),
        "independently reconstructed realized Scott widths")
 diagnostics=pd.read_csv(root/"results"/"weight_diagnostics.csv")
 expected_diagnostics_df=pd.DataFrame(expected_diagnostics)
 diagnostic_key=["validation_size","arm","node"]
 diagnostic_columns=[
  "kish_ess","kish_ess_fraction","entropy_ess","maximum_normalized_weight",
  "top_1_percent_mass","top_5_percent_mass",
 ]
 if diagnostics.duplicated(diagnostic_key).any() or expected_diagnostics_df.duplicated(diagnostic_key).any():
  raise AssertionError("duplicate weight-diagnostic key")
 diagnostics=diagnostics.sort_values(diagnostic_key).reset_index(drop=True)
 expected_diagnostics_df=expected_diagnostics_df.sort_values(diagnostic_key).reset_index(drop=True)
 if len(diagnostics)!=len(expected_diagnostics_df):
  raise AssertionError("weight-diagnostic row count mismatch")
 for col in diagnostic_key:
  if diagnostics[col].tolist()!=expected_diagnostics_df[col].tolist():
   raise AssertionError(f"weight-diagnostic key mismatch: {col}")
 for col in diagnostic_columns:
  _close(diagnostics[col].to_numpy(float),expected_diagnostics_df[col].to_numpy(float),
         f"independently reconstructed weight diagnostic {col}")
 return {"status":"PASS","manifest_files_verified":checked,"metric_rows_reconstructed":len(want),
  "scott_width_rows_reconstructed":len(expected_scott_df),
  "weight_diagnostic_rows_reconstructed":len(expected_diagnostics_df),
  "realized_scott_width_min":float(expected_scott_df.realized_bandwidth_redshift.min()),
  "realized_scott_width_max":float(expected_scott_df.realized_bandwidth_redshift.max()),
  "validation_sizes":sizes,"absolute_widths":widths,
  "scope":"Independent raw-input reconstruction of five-arm KDE density curves, IAE/ISE/boundary IAE, realized common-Scott physical widths, and raw-weight concentration diagnostics; stored validation-posterior nodes are consumed, not independently refitted. No H0 inference."}

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
 a=p.parse_args(); report=audit(a.artifact_dir); out=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if a.output: a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(out)
 print(out)
if __name__=="__main__": main()
