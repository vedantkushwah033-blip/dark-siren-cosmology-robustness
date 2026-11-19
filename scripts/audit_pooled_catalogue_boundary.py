"""Run the existing expanded-prior boundary audit over two independent artifacts.

Each artifact is independently reconstructed on [30,110] using 561/1121 H0
grids. The wrapper pools only after checking identical design and disjoint seed
pairs. It is still a catalogue-only diagnostic, not an all-arm calibration test.
"""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from audit_catalogue_only_prior_boundary import audit

def seed(master,rep,stream):
    return int(np.random.SeedSequence([master,rep,stream]).generate_state(1,dtype=np.uint32)[0])

def run(root_a:Path,root_b:Path,out:Path):
    roots=[root_a.resolve(),root_b.resolve()]
    configs=[json.loads((r/"config.json").read_text()) for r in roots]
    keys=["gammas","validation_sizes","n_parent","n_detected","validation_frame","curvature_beta",
          "h0_true","h0_prior_support","h0_points","scott_multiplier","quadrature_order",
          "sigma_mpc","detection_threshold_mpc"]
    for key in keys:
        if configs[0][key]!=configs[1][key]: raise AssertionError(f"design mismatch: {key}")
    pairs=[]
    for c in configs:
        pairs.extend((seed(c["master_seed"],rep,0),seed(c["master_seed"],rep,1))
                     for rep in range(int(c["replicates"])))
    if len(pairs)!=len(set(pairs)): raise AssertionError("replicate seed pairs overlap")
    frames=[]; reports=[]
    for i,root in enumerate(roots):
        csv=out.with_name(out.stem+f"_stream{i+1}.csv")
        report=audit(root,csv)
        if report["status"]!="PASS": raise AssertionError(f"boundary audit failed for stream {i+1}")
        frame=pd.read_csv(csv); frame["master_seed"]=int(configs[i]["master_seed"])
        frames.append(frame); reports.append(report)
    pooled=pd.concat(frames,ignore_index=True)
    expected=sum(int(c["replicates"]) for c in configs)*len(configs[0]["gammas"])
    if len(pooled)!=expected: raise AssertionError(f"expected {expected} replicate-gamma rows, got {len(pooled)}")
    summary=pooled.groupby("gamma").agg(
        replicates=("replicate","count"),
        mean_median=("median_1121","mean"),
        mean_q025=("q025_1121","mean"),
        mean_mass_30_50=("mass_30_50_1121","mean"),
        mean_mass_30_55=("mass_30_55_1121","mean"),
        mean_mass_70_90=("mass_70_90_1121","mean"),
        max_median_grid_abs_diff=("median_grid_abs_diff","max"),
        max_q025_grid_abs_diff=("q025_grid_abs_diff","max"),
        max_mass_grid_abs_diff=("mass_30_50_grid_abs_diff","max")).reset_index()
    pooled.to_csv(out.with_suffix(".csv"),index=False)
    summary.to_csv(out.with_name(out.stem+"_summary.csv"),index=False)
    result={"status":"PASS","individual_reports":reports,"pooled_replicate_gamma_rows":len(pooled),
        "unique_seed_pairs":len(pairs),"summary_by_gamma":summary.to_dict(orient="records"),
        "max_median_grid_abs_diff":float(pooled.median_grid_abs_diff.max()),
        "max_q025_grid_abs_diff":float(pooled.q025_grid_abs_diff.max()),
        "max_mass_30_50_grid_abs_diff":float(pooled.mass_30_50_grid_abs_diff.max()),
        "interpretation":"Catalogue-only posterior boundary diagnostic; not all-arm calibration or observational inference."}
    out.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__);p.add_argument("artifact_a",type=Path);p.add_argument("artifact_b",type=Path)
    p.add_argument("--output-prefix",type=Path,default=Path("pooled_catalogue_boundary"));a=p.parse_args()
    print(json.dumps(run(a.artifact_a,a.artifact_b,a.output_prefix),indent=2))
if __name__=="__main__":main()
