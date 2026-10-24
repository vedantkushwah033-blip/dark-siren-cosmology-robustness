"""Posterior-grid resolution audit for the joint distance/host-mixture diagnostic.

Uses the same deterministic synthetic observations across grid resolutions. This
is a numerical sensitivity check, not a new scientific data-generating model.
"""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
from dark_siren_robustness.cosmology import luminosity_distance_flat_lcdm
from dark_siren_robustness.cosmological_inference import distance_grid_for_model, posterior_h0_from_distance_grid
from scripts.run_joint_distance_model_host_mixture import TRUTHS, HOSTS, WEIGHTS, PRIOR, SEED, stable_seed, summary

ARMS=("lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed")
GRIDS=(2001,4001,8001)

def run(replicates=30, seed=SEED, grids=GRIDS):
    if replicates < 1 or len(grids)<2 or any(g<1001 for g in grids):
        raise ValueError("need positive replicates and at least two grids of >=1001 points")
    rows=[]
    for truth in TRUTHS:
      for zpair in HOSTS:
       for wpair in WEIGHTS:
        cfg={"truth":truth,"z":zpair,"w":wpair}
        rng=np.random.default_rng(stable_seed(seed,cfg))
        z=np.asarray(zpair); w=np.asarray(wpair)
        true_d=luminosity_distance_flat_lcdm(z,truth,omega_m=.3,integration_steps=4000)
        sigma=max(20.,.05*float(np.dot(w,true_d)))
        zeff=float(np.dot(w,z))
        for rep in range(replicates):
          host=int(rng.choice(2,p=w))
          obs=float(true_d[host]+rng.normal(0,sigma))
          obskey=dict(true_H0=truth,host_z1=z[0],host_z2=z[1],weight1=w[0],weight2=w[1],replicate=rep,
                      observed_distance_mpc=obs,sigma_mpc=sigma)
          for gs in grids:
            hgrid=np.linspace(25,125,gs)
            models={
             "lcdm_two_host":(z,"flat_lcdm",w),
             "lowz_two_host":(z,"low_z",w),
             "lcdm_collapsed":([zeff],"flat_lcdm",[1.]),
             "lowz_collapsed":([zeff],"low_z",[1.])}
            for arm,(zs,model,weights) in models.items():
              pred=distance_grid_for_model(hgrid,zs,model,omega_m=.3,integration_steps=2000)
              p=posterior_h0_from_distance_grid(hgrid,pred,obs,sigma,host_weights=weights,prior_min=PRIOR[0],prior_max=PRIOR[1])
              sm=summary(hgrid,p)
              rows.append(dict(**obskey,grid_size=gs,arm=arm,**sm))
    return rows

def summarize(rows, grids=GRIDS):
    keyed={}
    for r in rows:
      k=(r["true_H0"],r["host_z1"],r["host_z2"],r["weight1"],r["weight2"],r["replicate"],r["arm"])
      keyed.setdefault(k,{})[int(r["grid_size"])]=r
    out=[]
    ref=max(grids)
    for k,vals in keyed.items():
      if len(vals)!=len(grids): raise ValueError(f"missing grid comparison for {k}")
      for gs in grids:
        r=vals[gs]
        out.append(dict(true_H0=k[0],host_z1=k[1],host_z2=k[2],weight1=k[3],weight2=k[4],replicate=k[5],arm=k[6],
          grid_size=gs,reference_grid=ref,median_abs_diff_vs_reference=abs(r["median"]-vals[ref]["median"]),
          lower68_abs_diff_vs_reference=abs(r["lower_68"]-vals[ref]["lower_68"]),
          upper68_abs_diff_vs_reference=abs(r["upper_68"]-vals[ref]["upper_68"]),
          lower95_abs_diff_vs_reference=abs(r["lower_95"]-vals[ref]["lower_95"]),
          upper95_abs_diff_vs_reference=abs(r["upper_95"]-vals[ref]["upper_95"])))
    return out

def main():
 p=argparse.ArgumentParser();p.add_argument("--replicates",type=int,default=30);p.add_argument("--seed",type=int,default=SEED);p.add_argument("--output-dir",default="results/grid-resolution-audit");a=p.parse_args()
 root=Path(a.output_dir);root.mkdir(parents=True,exist_ok=True)
 rows=run(a.replicates,a.seed); diffs=summarize(rows)
 for name,data in (("grid_summaries.csv",rows),("grid_differences.csv",diffs)):
  with (root/name).open("w",newline="",encoding="utf-8") as f:
   w=csv.DictWriter(f,fieldnames=list(data[0]));w.writeheader();w.writerows(data)
 fields=("median_abs_diff_vs_reference","lower68_abs_diff_vs_reference","upper68_abs_diff_vs_reference","lower95_abs_diff_vs_reference","upper95_abs_diff_vs_reference")
 audit={field:max(r[field] for r in diffs if r["grid_size"]<max(GRIDS)) for field in fields}
 manifest={"experiment":"joint_grid_resolution_audit","status":"synthetic numerical sensitivity check","seed":a.seed,"replicates_per_cell":a.replicates,
   "cell_count":27,"observations":27*a.replicates,"arms":list(ARMS),"grid_sizes":list(GRIDS),"prior":list(PRIOR),
   "raw_posterior_summaries":len(rows),"difference_rows":len(diffs),"maximum_difference_from_8001_grid":audit}
 (root/"manifest.json").write_text(json.dumps(manifest,indent=2)+"\n")
 print(json.dumps(manifest,indent=2))
if __name__=="__main__": main()
