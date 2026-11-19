"""Independent seed and posterior reconstruction audit for full-parent H0 SBC."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm

from dark_siren_robustness.generative_selection import generate_synthetic_population
C_KM_S=299792.458

def _seed(master,rep,stream):
    return int(np.random.SeedSequence([int(master),int(rep),int(stream)]).generate_state(1,dtype=np.uint32)[0])

def _rebuild_posterior(parent_z,distances,h,host_w,sigma,threshold):
    z=np.asarray(parent_z,float); w=np.asarray(host_w,float); d=np.asarray(distances,float)
    logw=np.log(w); logw-=np.logaddexp.reduce(logw)
    mu=C_KM_S*z[None,:]/h[:,None]
    logbeta=np.logaddexp.reduce(norm.logcdf((mu-threshold)/sigma)+logw[None,:],axis=1)
    lognum=np.zeros(len(h))
    for x in d:
        lognum+=np.logaddexp.reduce(norm.logpdf(x,loc=mu,scale=sigma)+logw[None,:],axis=1)
    lp=lognum-len(d)*logbeta
    density=np.exp(lp-np.max(lp)); density/=trapezoid(density,h)
    cdf=np.r_[0.0,cumulative_trapezoid(density,h)]; cdf/=cdf[-1]
    return cdf

def audit(root:Path,output:Path|None=None,expected_replicates=200):
    root=root.resolve()
    manifest=json.loads((root/"manifest.json").read_text())
    bad=[name for name,digest in manifest.items() if not (root/name).is_file() or hashlib.sha256((root/name).read_bytes()).hexdigest()!=digest]
    if bad: raise AssertionError(f"manifest mismatch: {bad[:10]}")
    design=json.loads((root/"design.json").read_text())
    frame=pd.read_csv(root/"sbc_replicates.csv")
    n=int(design["replicates"]); master=int(design["master_seed"])
    if n!=expected_replicates or len(frame)!=n or frame.replicate.duplicated().any():
        raise AssertionError("planned SBC replicate inventory mismatch")
    if design["status"]!="COMPLETE" or design["generation"]!="generate_synthetic_population; full parent and fixed-count detected events":
        raise AssertionError("unexpected SBC design metadata")
    lo,hi=map(float,design["h0_prior"]); h=np.linspace(lo,hi,int(design["h0_grid_points"]))
    rows=[]
    for rep in range(n):
        row=frame[frame.replicate==rep]
        if len(row)!=1: raise AssertionError(f"missing/duplicate replicate {rep}")
        row=row.iloc[0]
        ts=_seed(master,rep,0); ss=_seed(master,rep,1)
        if int(row.truth_seed)!=ts or int(row.simulation_seed)!=ss:
            raise AssertionError(f"seed provenance mismatch replicate {rep}")
        truth=float(np.random.default_rng(ts).uniform(lo,hi))
        if not np.isclose(row.h0_true,truth,rtol=0,atol=1e-12):
            raise AssertionError(f"prior draw mismatch replicate {rep}")
        sim=generate_synthetic_population(seed=ss,n_parent=int(design["n_parent"]),
          n_detected=int(design["n_detected"]),h0_true=truth)
        parent=sim["parent_truth"]; events=sim["gw_events_observed"]; cfg=sim["config"]
        z=parent.redshift.to_numpy(float); d=events.distance_mpc.to_numpy(float)
        zmin,zmax=map(float,cfg["redshift_domain"])
        w=1.0+1.5*(z-zmin)/(zmax-zmin)
        cdf=_rebuild_posterior(z,d,h,w,float(cfg["sigma_mpc"]),float(cfg["detection_threshold_mpc"]))
        rank=float(np.interp(truth,h,cdf))
        q=np.interp([.025,.05,.10,.16,.25,.50,.75,.84,.90,.95,.975],cdf,h)
        expected={"posterior_rank":rank,"h0_median":q[5],"median_bias":q[5]-truth,
          "covered50":int(q[4]<=truth<=q[6]),"covered68":int(q[3]<=truth<=q[7]),
          "covered80":int(q[2]<=truth<=q[8]),"covered90":int(q[1]<=truth<=q[9]),
          "covered95":int(q[0]<=truth<=q[10]),"n_parent":len(z),"n_detected":len(d),
          "n_catalogue":len(sim["catalogue_observed"])}
        for key,value in expected.items():
            if not np.isclose(float(row[key]),float(value),rtol=0,atol=2e-8):
                raise AssertionError(f"{key} mismatch at replicate {rep}: archived={row[key]}, rebuilt={value}")
        rows.append({"replicate":rep,"posterior_rank":rank,"median_bias":q[5]-truth})
    rebuilt=pd.DataFrame(rows)
    np.testing.assert_allclose(rebuilt.posterior_rank,frame.sort_values("replicate").posterior_rank,rtol=0,atol=2e-8)
    report={"status":"PASS","manifest_entries_verified":len(manifest),"replicates_regenerated":n,
      "parent_populations_regenerated":n,"detected_event_sequences_regenerated":n,
      "posterior_rank_and_quantile_rows_reconstructed":n,
      "scope":"Independent raw-seed reconstruction of a synthetic full-parent prior-predictive SBC experiment; not an inclusion-arm result or observational inference."}
    if output:
        output=Path(output); output.parent.mkdir(parents=True,exist_ok=True)
        output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
    return report

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_dir",type=Path)
    p.add_argument("--expected-replicates",type=int,default=200); p.add_argument("--output",type=Path,required=True)
    a=p.parse_args(); print(json.dumps(audit(a.artifact_dir,a.output,a.expected_replicates),indent=2))
if __name__=="__main__": main()
