"""Paired multi-replicate runner for quadratic-logit curvature density stress tests.

The same population/validation seed pair is reused across gamma cells within a
replicate, enabling paired contrasts under common random numbers. This remains
a density-only synthetic benchmark.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
from types import SimpleNamespace

import numpy as np
import pandas as pd

from run_quadratic_logit_curvature_density_smoke import run as run_one

DEFAULT_GAMMAS=(-2.0,-1.0,-0.5,0.0,0.5,1.0,2.0)
DEFAULT_WIDTHS=(0.008,0.012,0.016,0.024,0.032,0.040,0.050,0.060,0.080)

def _hash(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def run(args):
    root=Path(args.output_dir).resolve()
    root.mkdir(parents=True,exist_ok=True)
    if any(root.iterdir()): raise ValueError("output directory must be empty")
    gammas=tuple(float(x) for x in args.gammas)
    widths=tuple(sorted(set(float(x) for x in args.absolute_widths)))
    if len(gammas)<2 or len(set(gammas))!=len(gammas): raise ValueError("gamma grid must contain at least two unique values")
    if any(not __import__("math").isfinite(x) for x in gammas): raise ValueError("gamma values must be finite")
    if args.replicates<2 or args.seed_stride<2: raise ValueError("at least two replicates and seed_stride >= 2 are required")
    if len(widths)<2 or any(not __import__("math").isfinite(x) or x<=0 for x in widths): raise ValueError("invalid absolute bandwidth grid")
    if not any(x == 0.0 for x in gammas): raise ValueError("gamma grid must include zero for paired reference contrasts")
    if args.master_seed < 0: raise ValueError("master seed must be non-negative")
    if args.n_parent < 100 or args.validation_frame < 5: raise ValueError("invalid parent or validation-frame size")
    if not args.validation_sizes or sorted(set(args.validation_sizes)) != list(args.validation_sizes): raise ValueError("validation sizes must be strictly increasing")
    if args.validation_sizes[0] < 5 or args.validation_sizes[-1] > args.validation_frame: raise ValueError("validation sizes exceed the validation frame")
    if args.quadrature_order < 3 or args.grid_points < 51: raise ValueError("invalid quadrature order or redshift grid")
    if not __import__("math").isfinite(args.beta): raise ValueError("beta must be finite")
    if not __import__("math").isfinite(args.scott_multiplier) or args.scott_multiplier <= 0: raise ValueError("Scott multiplier must be finite and positive")
    cfg={
      "schema_version":1,
      "experiment":"paired_quadratic_logit_curvature_density_pilot",
      "status":"DEVELOPMENT_PILOT_NOT_CONFIRMATORY_NOT_BANDWIDTH_SELECTION",
      "master_seed":int(args.master_seed),"seed_stride":int(args.seed_stride),
      "replicates":int(args.replicates),"gammas":list(gammas),
      "n_parent":int(args.n_parent),"validation_frame":int(args.validation_frame),
      "validation_sizes":list(args.validation_sizes),"absolute_widths":list(widths),
      "beta":float(args.beta),"quadrature_order":int(args.quadrature_order),
      "grid_points":int(args.grid_points),"scott_multiplier":float(args.scott_multiplier),
      "seed_design":"Within each replicate, all gamma cells share population_seed and validation_seed (common random numbers); replicate seed pairs are disjoint.",
      "analysis_model":"degree_one_logistic",
      "arm_C":"true quadratic inclusion probabilities; oracle diagnostic only",
      "arm_D_E":"degree-one logistic fitted to independent validation outcomes",
      "arm_E_primary":"loss of posterior-mean density",
      "limitations":["synthetic density recovery only","no bandwidth selection or holdout validation",
        "no H0 likelihood, posterior, PIT or coverage","no observational or public GW data analysis",
        "Arm C is a truth-informed oracle diagnostic, not an implementable estimator"]
    }
    (root/"config.json").write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    (root/"cells").mkdir()
    ledger=[]
    for rid in range(args.replicates):
      seed=args.master_seed+rid*args.seed_stride
      validation_seed=seed+1
      for gamma in gammas:
        label=f"gamma_{gamma:+.2f}".replace("+","p").replace("-","m").replace(".","_")
        dest=root/"cells"/label/f"replicate_{rid:04d}"
        dest.mkdir(parents=True)
        row={"replicate":rid,"gamma":gamma,"population_seed":seed,"validation_seed":validation_seed,
             "status":"NOT_RUN","audit_status":"NOT_RUN","metric_rows":0,"error":""}
        try:
          result=run_one(SimpleNamespace(output_dir=dest,seed=seed,validation_seed=validation_seed,
            n_parent=args.n_parent,validation_frame=args.validation_frame,
            validation_sizes=tuple(args.validation_sizes),gamma=gamma,beta=args.beta,
            quadrature_order=args.quadrature_order,grid_points=args.grid_points,
            scott_multiplier=args.scott_multiplier,absolute_widths=widths))
          row.update(status="SUCCESS",metric_rows=int(result["metric_rows"]),
                     weight_diagnostic_rows=int(result["weight_diagnostic_rows"]),
                     scott_width_rows=int(result["scott_width_rows"]))
        except Exception as exc:
          row.update(status="FAILED",error=f"{type(exc).__name__}: {exc}")
        ledger.append(row)
        pd.DataFrame(ledger).to_csv(root/"cell_ledger.csv",index=False)
    metric_frames=[]; diagnostic_frames=[]; scott_frames=[]
    for row in ledger:
      if row["status"]!="SUCCESS": continue
      label=f"gamma_{row['gamma']:+.2f}".replace("+","p").replace("-","m").replace(".","_")
      child=root/"cells"/label/f"replicate_{int(row['replicate']):04d}"
      for filename,frames in [("density_metrics.csv",metric_frames),("weight_diagnostics.csv",diagnostic_frames),("scott_widths.csv",scott_frames)]:
        frame=pd.read_csv(child/"results"/filename)
        frame.insert(0,"replicate",int(row["replicate"]))
        frame.insert(1,"gamma",float(row["gamma"]))
        frames.append(frame)
    def combine(frames,name,columns):
      combined=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=columns)
      combined.to_csv(root/name,index=False)
    combine(metric_frames,"density_metrics_all_cells.csv",[
      "replicate","gamma","validation_size","absolute_bandwidth_sd","arm","node",
      "iae","ise","boundary_iae","status","kish_ess"])
    combine(diagnostic_frames,"weight_diagnostics_all_cells.csv",[
      "replicate","gamma","validation_size","arm","node","kish_ess","kish_ess_fraction",
      "maximum_normalized_weight","top_1_percent_mass","top_5_percent_mass"])
    combine(scott_frames,"scott_widths_all_cells.csv",[
      "replicate","gamma","validation_size","arm","node","realized_bandwidth_redshift"])
    # Paired gamma and arm contrasts use common random numbers.
    metric_columns = ["gamma", "replicate", "validation_size", "absolute_bandwidth_sd", "arm", "node", "iae", "ise", "boundary_iae", "status", "kish_ess"]
    metric_all = pd.concat(metric_frames, ignore_index=True) if metric_frames else pd.DataFrame(columns=metric_columns)
    primary_arms = ("A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
                    "D_map_plugin_ipw", "E_posterior_mean_density")
    metric_all = metric_all[metric_all["arm"].isin(primary_arms) & (metric_all["node"] == -1)].copy()
    metric_names = ("iae", "ise", "boundary_iae")
    contrast_rows = []
    pair_keys = ["replicate", "validation_size", "absolute_bandwidth_sd"]
    zero_rows = metric_all[np.isclose(metric_all["gamma"], 0.0)]
    if zero_rows.empty:
      # Preserve a fully auditable incomplete artifact when every gamma=0 cell failed.
      # Gamma-vs-zero contrasts are then undefined and remain absent.
      pass
    for gamma in gammas:
      if gamma == 0.0:
        continue
      current = metric_all[np.isclose(metric_all["gamma"], gamma)]
      for arm in primary_arms:
        paired = current[current.arm == arm].merge(
          zero_rows[zero_rows.arm == arm], on=pair_keys,
          suffixes=("", "_reference"), validate="one_to_one")
        for _, row in paired.iterrows():
          for metric in metric_names:
            contrast_rows.append({
              "contrast_type": "gamma_vs_zero", "gamma": gamma, "reference_gamma": 0.0,
              "arm": arm, "reference_arm": arm, "replicate": int(row.replicate),
              "validation_size": int(row.validation_size),
              "absolute_bandwidth_sd": float(row.absolute_bandwidth_sd), "metric": metric,
              "value": float(row[metric]), "reference_value": float(row[f"{metric}_reference"]),
              "difference": float(row[metric] - row[f"{metric}_reference"]),
            })
    arm_pairs = (("D_map_plugin_ipw", "A_catalogue_only"),
                 ("D_map_plugin_ipw", "C_known_probability_ipw"),
                 ("C_known_probability_ipw", "B_parent_oracle"),
                 ("E_posterior_mean_density", "D_map_plugin_ipw"))
    for gamma in gammas:
      at_gamma = metric_all[np.isclose(metric_all.gamma, gamma)]
      for arm, reference_arm in arm_pairs:
        paired = at_gamma[at_gamma.arm == arm].merge(
          at_gamma[at_gamma.arm == reference_arm], on=pair_keys,
          suffixes=("", "_reference"), validate="one_to_one")
        for _, row in paired.iterrows():
          for metric in metric_names:
            contrast_rows.append({
              "contrast_type": "arm_pair", "gamma": gamma, "reference_gamma": float("nan"),
              "arm": arm, "reference_arm": reference_arm, "replicate": int(row.replicate),
              "validation_size": int(row.validation_size),
              "absolute_bandwidth_sd": float(row.absolute_bandwidth_sd), "metric": metric,
              "value": float(row[metric]), "reference_value": float(row[f"{metric}_reference"]),
              "difference": float(row[metric] - row[f"{metric}_reference"]),
            })
    contrast_columns = ["contrast_type", "gamma", "reference_gamma", "arm", "reference_arm",
      "replicate", "validation_size", "absolute_bandwidth_sd", "metric", "value",
      "reference_value", "difference"]
    contrasts = pd.DataFrame(contrast_rows, columns=contrast_columns)
    contrasts.to_csv(root/"paired_contrasts_all_cells.csv", index=False)
    group_columns = ["contrast_type", "gamma", "reference_gamma", "arm", "reference_arm",
      "validation_size", "absolute_bandwidth_sd", "metric"]
    if len(contrasts):
      summary = contrasts.groupby(group_columns, dropna=False).difference.agg(
        n_pairs="count", mean_difference="mean", sd_difference="std").reset_index()
      summary["mcse"] = summary.sd_difference / summary.n_pairs.pow(0.5)
    else:
      summary = pd.DataFrame(columns=group_columns+["n_pairs","mean_difference","sd_difference","mcse"])
    summary.to_csv(root/"paired_contrast_summary.csv", index=False)
    cfg["attempted_cells"]=len(ledger)
    cfg["successful_cells"]=sum(x["status"]=="SUCCESS" for x in ledger)
    cfg["failed_cells"]=sum(x["status"]=="FAILED" for x in ledger)
    (root/"config.json").write_text(json.dumps(cfg,indent=2,sort_keys=True)+"\n")
    files=sorted(p for p in root.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS" if cfg["failed_cells"] == 0 else "INCOMPLETE","file_count":len(files),"files":{str(p.relative_to(root)):
      {"bytes":p.stat().st_size,"sha256":_hash(p)} for p in files}}
    (root/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return {"status":"COMPLETED_WITH_FAILURES" if cfg["failed_cells"] else "PASS",
      "planned_cells":len(ledger),"successful_cells":cfg["successful_cells"],
      "failed_cells":cfg["failed_cells"],"metric_rows":sum(len(x) for x in metric_frames),
      "manifest_files":len(files)}

def main():
 p=argparse.ArgumentParser(description=__doc__)
 p.add_argument("--master-seed",type=int,default=20261201)
 p.add_argument("--seed-stride",type=int,default=10)
 p.add_argument("--replicates",type=int,default=20)
 p.add_argument("--gammas",type=lambda s:tuple(float(x) for x in s.split(",")),default=DEFAULT_GAMMAS)
 p.add_argument("--n-parent",type=int,default=4000)
 p.add_argument("--validation-frame",type=int,default=1500)
 p.add_argument("--validation-sizes",type=lambda s:tuple(int(x) for x in s.split(",")),default=(100,300,1000))
 p.add_argument("--beta",type=float,default=-3.9)
 p.add_argument("--quadrature-order",type=int,default=5)
 p.add_argument("--grid-points",type=int,default=401)
 p.add_argument("--scott-multiplier",type=float,default=0.7)
 p.add_argument("--absolute-widths",type=lambda s:tuple(float(x) for x in s.split(",")),default=DEFAULT_WIDTHS)
 p.add_argument("--output-dir",type=Path,required=True)
 print(json.dumps(run(p.parse_args()),indent=2,sort_keys=True))
if __name__=="__main__": main()
