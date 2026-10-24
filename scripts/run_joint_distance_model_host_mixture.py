"""Paired synthetic distance-model × host-mixture interaction experiment.

Synthetic diagnostic only; not an observational dark-siren analysis.
"""
from __future__ import annotations
import argparse, csv, hashlib, json
from pathlib import Path
import numpy as np
from dark_siren_robustness.cosmology import luminosity_distance_flat_lcdm
from dark_siren_robustness.cosmological_inference import distance_grid_for_model, posterior_h0_from_distance_grid

TRUTHS=(60.0,70.0,80.0)
HOSTS=((0.02,0.04),(0.02,0.10),(0.10,0.40))
WEIGHTS=((0.5,0.5),(0.7,0.3),(0.9,0.1))
PRIOR=(40.0,100.0)
SEED=20261025

def stable_seed(seed, config):
    payload=json.dumps([int(seed),config],sort_keys=True,separators=(",",":"))
    return int.from_bytes(hashlib.sha256(payload.encode()).digest()[:8],"big")%(2**63-1)

def summary(grid,p):
    # Normalize and integrate on supported prior points only.
    mask=(grid>=PRIOR[0])&(grid<=PRIOR[1])
    x=grid[mask]; y=np.asarray(p)[mask]
    y=y/np.trapezoid(y,x)
    c=np.zeros_like(x)
    c[1:]=np.cumsum((y[1:]+y[:-1])*0.5*np.diff(x))
    c=c/c[-1]
    def q(prob): return float(np.interp(prob,c,x))
    return {"median":q(.5),"lower_68":q(.16),"upper_68":q(.84),"lower_95":q(.025),"upper_95":q(.975)}

def run(replicates=300, seed=SEED, grid_size=4001, bootstrap_draws=5000):
    if replicates<1 or grid_size<1001 or bootstrap_draws<1000: raise ValueError("invalid run sizes")
    hgrid=np.linspace(25,125,grid_size)
    rows=[]
    arms=("lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed")
    for truth in TRUTHS:
      for zpair in HOSTS:
       for wpair in WEIGHTS:
        cfg={"truth":truth,"z":zpair,"w":wpair}
        rng=np.random.default_rng(stable_seed(seed,cfg))
        z=np.asarray(zpair); w=np.asarray(wpair)
        true_means=luminosity_distance_flat_lcdm(z,truth,omega_m=.3,integration_steps=4000)
        sigma=max(20.0,.05*float(np.dot(w,true_means)))
        zeff=float(np.dot(w,z))
        # Precompute distance grids for each model and host representation.
        grids={
          "lcdm_two_host":distance_grid_for_model(hgrid,z,"flat_lcdm",omega_m=.3,integration_steps=2000),
          "lowz_two_host":distance_grid_for_model(hgrid,z,"low_z"),
          "lcdm_collapsed":distance_grid_for_model(hgrid,[zeff],"flat_lcdm",omega_m=.3,integration_steps=2000),
          "lowz_collapsed":distance_grid_for_model(hgrid,[zeff],"low_z"),
        }
        for rep in range(replicates):
          host=int(rng.choice(2,p=w))
          obs=float(true_means[host]+rng.normal(0,sigma))
          base={"true_H0":truth,"host_z1":z[0],"host_z2":z[1],"weight1":w[0],"weight2":w[1],
                "effective_redshift":zeff,"replicate":rep,"observed_distance_mpc":obs,"sigma_mpc":sigma,
                "seed":stable_seed(seed,cfg),"generated_host_index":host}
          for arm in arms:
            p=posterior_h0_from_distance_grid(hgrid,grids[arm],obs,sigma,
                 host_weights=(w if "two_host" in arm else [1.0]),prior_min=PRIOR[0],prior_max=PRIOR[1])
            s=summary(hgrid,p)
            row=dict(base,inference_arm=arm)
            for k,v in s.items(): row[k]=v
            row["error_H0"]=s["median"]-truth
            row["absolute_error_H0"]=abs(row["error_H0"])
            row["covered_68"]=int(s["lower_68"]<=truth<=s["upper_68"])
            row["covered_95"]=int(s["lower_95"]<=truth<=s["upper_95"])
            row["width_68"]=s["upper_68"]-s["lower_68"]
            row["width_95"]=s["upper_95"]-s["lower_95"]
            rows.append(row)
    return rows

def write_results(rows,out,bootstrap_draws=5000,seed=SEED,grid_size=4001):
    out=Path(out);out.mkdir(parents=True,exist_ok=True)
    with (out/"per_run.csv").open("w",newline="",encoding="utf-8") as f:
      w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    keys=("true_H0","host_z1","host_z2","weight1","weight2")
    groups={}
    for r in rows: groups.setdefault(tuple(r[k] for k in keys),{}).setdefault(r["replicate"],{})[r["inference_arm"]]=r
    summary_rows=[]; contrast_rows=[]
    for cfg,reps in groups.items():
      for arm in ("lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed"):
        rr=[v[arm] for v in reps.values()]
        e=np.array([x["error_H0"] for x in rr]);n=len(e)
        summary_rows.append(dict(zip(keys,cfg),arm=arm,n=n,bias=float(e.mean()),mae=float(np.abs(e).mean()),rmse=float(np.sqrt(np.mean(e*e))),
           coverage_68=float(np.mean([x["covered_68"] for x in rr])),coverage_95=float(np.mean([x["covered_95"] for x in rr])),
           mean_width_68=float(np.mean([x["width_68"] for x in rr])),mean_width_95=float(np.mean([x["width_95"] for x in rr]))))
      pairs=[
        ("distance_mismatch_two_host","lowz_two_host","lcdm_two_host"),
        ("distance_mismatch_collapsed","lowz_collapsed","lcdm_collapsed"),
        ("host_collapse_lcdm","lcdm_collapsed","lcdm_two_host"),
        ("host_collapse_lowz","lowz_collapsed","lowz_two_host")]
      deltas={}
      for name,a,b in pairs:
        d=np.array([v[a]["absolute_error_H0"]-v[b]["absolute_error_H0"] for v in reps.values()])
        deltas[name]=d
        brng=np.random.default_rng(stable_seed(seed+1,[cfg,name]))
        boot=np.array([d[brng.integers(0,len(d),len(d))].mean() for _ in range(bootstrap_draws)])
        contrast_rows.append(dict(zip(keys,cfg),contrast=name,n=len(d),mean_delta=float(d.mean()),
          bootstrap_95_low=float(np.quantile(boot,.025)),bootstrap_95_high=float(np.quantile(boot,.975)),bootstrap_draws=bootstrap_draws))
      interaction=deltas["distance_mismatch_collapsed"]-deltas["distance_mismatch_two_host"]
      brng=np.random.default_rng(stable_seed(seed+1,[cfg,"interaction"]))
      boot=np.array([interaction[brng.integers(0,len(interaction),len(interaction))].mean() for _ in range(bootstrap_draws)])
      contrast_rows.append(dict(zip(keys,cfg),contrast="distance_mismatch_interaction",n=len(interaction),mean_delta=float(interaction.mean()),
        bootstrap_95_low=float(np.quantile(boot,.025)),bootstrap_95_high=float(np.quantile(boot,.975)),bootstrap_draws=bootstrap_draws))
    for name,data in (("summary.csv",summary_rows),("paired_contrasts.csv",contrast_rows)):
      with (out/name).open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
    manifest={"experiment":"joint_distance_model_host_mixture","status":"synthetic diagnostic only","seed":seed,
      "replicates_per_cell":len(rows)//(len(TRUTHS)*len(HOSTS)*len(WEIGHTS)*4),"cell_count":len(TRUTHS)*len(HOSTS)*len(WEIGHTS),
      "observation_count":len(rows)//4,"raw_inference_rows":len(rows),"grid_size":grid_size,"prior":list(PRIOR),
      "truths":list(TRUTHS),"host_pairs":[list(x) for x in HOSTS],"weight_pairs":[list(x) for x in WEIGHTS],
      "bootstrap_draws":bootstrap_draws,"arms":list(("lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed"))}
    (out/"experiment_manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
    print(f"PASS: {len(rows)} inference rows; {manifest['observation_count']} observations; {manifest['cell_count']} cells")
def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("--replicates",type=int,default=300);p.add_argument("--seed",type=int,default=SEED)
    p.add_argument("--grid-size",type=int,default=4001);p.add_argument("--bootstrap-draws",type=int,default=5000);p.add_argument("--output-dir",default="results/joint-distance-host-mixture")
    a=p.parse_args();rows=run(a.replicates,a.seed,a.grid_size,a.bootstrap_draws);write_results(rows,a.output_dir,a.bootstrap_draws,a.seed,a.grid_size)
if __name__=="__main__": main()
