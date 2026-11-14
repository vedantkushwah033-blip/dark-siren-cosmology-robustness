"""Independent raw-input auditor for absolute-bandwidth density smoke."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid
from scipy.stats import gaussian_kde
from scipy.special import expit

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
 expected=[]; expected_scott=[]; maxdiff=0.0
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
  xv=(v.redshift.to_numpy(float)-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])
  design=np.column_stack([np.ones(n),xv])
  pmap=expit(np.column_stack([np.ones(len(zcat)),(zcat-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])])@mode)
  pnodes=expit(np.column_stack([np.ones(len(zcat)),(zcat-float(cfg["coordinate_center"]))/float(cfg["coordinate_scale"])])@nodes.T).T
  _close(z[f"weights_A_catalogue_only_n{n}"],bc,"A base weights")
  _close(z[f"weights_B_parent_oracle_n{n}"],bp,"B base weights")
  _close(z[f"weights_C_known_probability_ipw_n{n}"],bc/ptrue,"C raw weights")
  _close(z[f"weights_D_map_plugin_ipw_n{n}"],bc/pmap,"D raw weights")
  _close(z[f"weights_E_nodes_n{n}"],bc[None,:]/pnodes,"E raw weights")
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
    expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":arm,"node":-1,**met,"status":"SUCCESS"})
   for j,w in enumerate(bc[None,:]/pnodes):
    d=_independent_kde(grid,zcat,w,rate,bw); node_ds.append(d)
    _close(z[f"density_E_nodes_n{n}_w{bw:.6f}"][j],d,f"E node density n={n} width={bw} node={j}")
    expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":"E_uncertainty_integrated_ipw","node":j,**_metric(d,truth,grid),"status":"SUCCESS"})
   node_ds=np.asarray(node_ds)
   qmetric={k:float(np.dot(q,[_metric(d,truth,grid)[k] for d in node_ds])) for k in ("iae","ise","boundary_iae")}
   _close(z[f"density_E_posterior_mean_n{n}_w{bw:.6f}"],np.tensordot(q,node_ds,axes=(0,0)),f"E posterior mean density n={n} width={bw}")
   expected.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":bw,"arm":"E_posterior_weighted_node_metric","node":-1,**qmetric,"status":"SUCCESS"})
 expected_df=pd.DataFrame(expected)
 key=["replicate","validation_size","absolute_bandwidth_sd","arm","node"]
 if metrics.duplicated(key).any() or expected_df.duplicated(key).any(): raise AssertionError("duplicate metric key")
 got=metrics.sort_values(key).reset_index(drop=True); want=expected_df.sort_values(key).reset_index(drop=True)
 if len(got)!=len(want): raise AssertionError("metric row count mismatch")
 for col in key+["status"]:
  if not got[col].equals(want[col]): raise AssertionError(f"metric key/status mismatch: {col}")
 for col in ["iae","ise","boundary_iae"]:
  _close(got[col].to_numpy(float),want[col].to_numpy(float),f"metric {col}")
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
 return {"status":"PASS","manifest_files_verified":checked,"metric_rows_reconstructed":len(want),
  "scott_width_rows_reconstructed":len(expected_scott_df),
  "realized_scott_width_min":float(expected_scott_df.realized_bandwidth_redshift.min()),
  "realized_scott_width_max":float(expected_scott_df.realized_bandwidth_redshift.max()),
  "validation_sizes":sizes,"absolute_widths":widths,
  "scope":"Independent raw-input reconstruction of five-arm KDE density curves, IAE/ISE/boundary IAE, and realized common-Scott physical widths; stored validation-posterior nodes are consumed, not independently refitted. No H0 inference."}

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
 a=p.parse_args(); report=audit(a.artifact_dir); out=json.dumps(report,indent=2,sort_keys=True)+"\n"
 if a.output: a.output.parent.mkdir(parents=True,exist_ok=True); a.output.write_text(out)
 print(out)
if __name__=="__main__": main()
