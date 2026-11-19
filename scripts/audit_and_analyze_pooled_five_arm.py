"""Independent pooled audit/readout for two disjoint five-arm 100-replicate artifacts.

Checks each artifact's hashes, lifecycle denominator, child audit attestations,
deterministic seeds, and archived child/top-level posterior agreement. Requires
matching frozen designs and disjoint replicate seed pairs before pooling.
Writes pooled cell metrics with Wilson coverage intervals and paired t intervals.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import norm, t

ARMS = ["A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw"]
CONTRASTS = [("E_uncertainty_integrated_ipw","D_map_plugin_ipw","E-D"),
             ("E_uncertainty_integrated_ipw","C_known_probability_ipw","E-C"),
             ("D_map_plugin_ipw","C_known_probability_ipw","D-C"),
             ("B_parent_oracle","C_known_probability_ipw","B-C")]
DESIGN_KEYS = ["gammas","validation_sizes","n_parent","n_detected","validation_frame",
               "curvature_beta","h0_true","h0_prior_support","h0_points","scott_multiplier",
               "quadrature_order","sigma_mpc","detection_threshold_mpc"]

def sha(path):
    h=hashlib.sha256()
    with Path(path).open("rb") as f:
        for block in iter(lambda:f.read(1<<20),b""): h.update(block)
    return h.hexdigest()

def seed(master, rep, stream):
    return int(np.random.SeedSequence([master,rep,stream]).generate_state(1,dtype=np.uint32)[0])

def wilson(k,n):
    if n == 0: return float("nan"),float("nan")
    z=float(norm.ppf(.975)); p=k/n; den=1+z*z/n
    c=(p+z*z/(2*n))/den; h=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(c-h),float(c+h)

def audit_one(root):
    root=Path(root).resolve()
    cfg=json.loads((root/"config.json").read_text())
    manifest=json.loads((root/"manifest.json").read_text())
    for rel,meta in manifest["files"].items():
        p=root/rel
        if not p.is_file() or p.stat().st_size != int(meta["bytes"]) or sha(p)!=meta["sha256"]:
            raise AssertionError(f"manifest mismatch: {root.name}/{rel}")
    if len(manifest["files"]) != int(manifest["file_count"]):
        raise AssertionError("manifest file count mismatch")
    planned=pd.read_csv(root/"planned_cells.csv")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    summaries=pd.read_csv(root/"replicate_posterior_summaries.csv")
    aggregates=pd.read_csv(root/"aggregate_metrics.csv")
    contrasts=pd.read_csv(root/"paired_contrasts.csv")
    reps=int(cfg["replicates"]); gammas=list(map(float,cfg["gammas"]))
    sizes=list(map(int,cfg["validation_sizes"]))
    expected=reps*len(gammas)*len(sizes)*len(ARMS)
    if len(planned)!=expected or len(ledger)!=expected or len(summaries)!=expected:
        raise AssertionError("planned/ledger/summary denominator mismatch")
    if planned.cell_id.duplicated().any() or ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate cell ID")
    if set(planned.cell_id)!=set(ledger.cell_id): raise AssertionError("planned/ledger IDs differ")
    if not (ledger.status=="SUCCESS").all() or not (ledger.audit_status=="PASS").all():
        raise AssertionError("non-success lifecycle cells found")
    if len(aggregates)!=len(gammas)*len(sizes)*len(ARMS): raise AssertionError("aggregate row count mismatch")
    if len(contrasts)!=reps*len(gammas)*len(sizes)*8: raise AssertionError("paired contrast row count mismatch")
    for rep in range(reps):
        aseed=seed(cfg["master_seed"],rep,0); vseed=seed(cfg["master_seed"],rep,1)
        block=ledger[ledger.replicate==rep]
        if set(block.analysis_seed)!={aseed} or set(block.validation_seed)!={vseed}:
            raise AssertionError("top-level seed schedule mismatch")
        for gamma in gammas:
            suffix=f"g{gamma:+g}".replace("+","p").replace("-","m")
            child=root/"replicates"/f"rep{rep:03d}_{suffix}"
            report=json.loads((child/"audit_report.json").read_text())
            cc=json.loads((child/"config.json").read_text())
            if report.get("status")!="PASS": raise AssertionError(f"child audit not PASS: {child}")
            if (int(cc["seed"])!=aseed or int(cc["validation_seed"])!=vseed or
                float(cc["curvature_gamma"])!=gamma or int(cc["n_parent"])!=int(cfg["n_parent"]) or
                int(cc["n_detected"])!=int(cfg["n_detected"]) or
                list(map(int,cc["validation_sizes"]))!=sizes):
                raise AssertionError(f"child config mismatch: {child}")
            child_summary=pd.read_csv(child/"results"/"posterior_summaries.csv")
            top=summaries[(summaries.replicate==rep)&(summaries.gamma==gamma)]
            cols=["validation_size","arm","status","h0_true","median_error","posterior_median",
                  "q025","q05","q16","q84","q95","q975"]
            left=child_summary[cols].sort_values(["validation_size","arm"]).reset_index(drop=True)
            right=top[cols].sort_values(["validation_size","arm"]).reset_index(drop=True)
            pd.testing.assert_frame_equal(left,right,check_dtype=False,check_exact=False,rtol=2e-10,atol=2e-10)
    # Recompute aggregate point estimates from replicate summaries.
    for (gamma,nval),_ in ledger.groupby(["gamma","validation_size"]):
        for arm in ARMS:
            g=summaries[(summaries.gamma==gamma)&(summaries.validation_size==nval)&(summaries.arm==arm)]
            row=aggregates[(aggregates.gamma==gamma)&(aggregates.validation_size==nval)&(aggregates.arm==arm)].iloc[0]
            e=g.median_error.to_numpy(float); cov=((g.q16<=g.h0_true)&(g.h0_true<=g.q84)).to_numpy()
            if not np.isclose(row.success_only_bias,e.mean(),rtol=1e-10,atol=1e-10): raise AssertionError("bias mismatch")
            if not np.isclose(row.success_only_rmse,np.sqrt(np.mean(e*e)),rtol=1e-10,atol=1e-10): raise AssertionError("RMSE mismatch")
            if not np.isclose(row.success_only_68pct_coverage,cov.mean(),rtol=1e-10,atol=1e-10): raise AssertionError("coverage mismatch")
    return cfg,summaries,{"status":"PASS","manifest_entries":len(manifest["files"]),
        "planned_cells":len(planned),"successful_cells":int((ledger.status=="SUCCESS").sum()),
        "child_audits":reps*len(gammas),"posterior_rows":len(summaries),
        "aggregate_rows":len(aggregates),"paired_contrast_rows":len(contrasts)}

def pooled_readout(roots, output_prefix):
    audited=[audit_one(r) for r in roots]
    c0,c1=audited[0][0],audited[1][0]
    for key in DESIGN_KEYS:
        if c0[key]!=c1[key]: raise AssertionError(f"frozen design mismatch: {key}")
    if c0["master_seed"]==c1["master_seed"]: raise AssertionError("master seeds must differ")
    seeds=[]
    for cfg,_,_ in audited:
        seeds.extend((seed(cfg["master_seed"],r,0),seed(cfg["master_seed"],r,1)) for r in range(int(cfg["replicates"])))
    if len(seeds)!=len(set(seeds)): raise AssertionError("replicate seed pairs overlap")
    frames=[]
    for cfg,frame,_ in audited:
        f=frame.copy(); f["stream"]=int(cfg["master_seed"]); frames.append(f)
    data=pd.concat(frames,ignore_index=True)
    total=int(c0["replicates"])+int(c1["replicates"])
    truth=float(c0["h0_true"])
    cell_rows=[]
    for (gamma,nval,arm),g in data.groupby(["gamma","validation_size","arm"],sort=True):
        err=g.median_error.to_numpy(float); covered=((g.q16<=truth)&(truth<=g.q84)).to_numpy()
        width=(g.q84-g.q16).to_numpy(float); lo,hi=wilson(int(covered.sum()),len(covered))
        cell_rows.append({"gamma":float(gamma),"validation_size":int(nval),"arm":arm,"n":len(g),
            "bias":float(err.mean()),"rmse":float(np.sqrt(np.mean(err**2))),
            "coverage68":float(covered.mean()),"coverage_wilson95_low":lo,"coverage_wilson95_high":hi,
            "mean_width68":float(width.mean()),"median_error_sd":float(err.std(ddof=1))})
    cells=pd.DataFrame(cell_rows)
    cells.to_csv(output_prefix.with_name(output_prefix.name+"_cells.csv"),index=False)
    paired=[]
    for (gamma,nval),group in data.groupby(["gamma","validation_size"],sort=True):
        byarm={arm:frame.set_index(["stream","replicate"]) for arm,frame in group.groupby("arm")}
        for left,right,label in CONTRASTS:
            z=byarm[left].join(byarm[right],how="inner",lsuffix="_left",rsuffix="_right")
            metrics={
                "median_error":z.median_error_left-z.median_error_right,
                "absolute_error":z.median_error_left.abs()-z.median_error_right.abs(),
                "squared_error":z.median_error_left**2-z.median_error_right**2,
                "coverage68":((z.q16_left<=truth)&(truth<=z.q84_left)).astype(float)-((z.q16_right<=truth)&(truth<=z.q84_right)).astype(float),
                "width68":(z.q84_left-z.q16_left)-(z.q84_right-z.q16_right)}
            for metric,series in metrics.items():
                x=series.to_numpy(float); n=len(x); mean=float(x.mean()); se=float(x.std(ddof=1)/np.sqrt(n))
                q=float(t.ppf(.975,n-1))
                paired.append({"gamma":float(gamma),"validation_size":int(nval),"contrast":label,
                    "metric":metric,"n":n,"mean_difference":mean,"se":se,
                    "t95_low":mean-q*se,"t95_high":mean+q*se})
    pd.DataFrame(paired).to_csv(output_prefix.with_name(output_prefix.name+"_paired.csv"),index=False)
    result={"status":"PASS","individual_audits":[a[2] for a in audited],
        "pooled_replicates":total,"pooled_posterior_rows":len(data),
        "unique_seed_pairs":len(seeds),"pooled_cell_rows":len(cells),"paired_interval_rows":len(paired),
        "interpretation":"Synthetic paired engineering benchmark only; no observational inference or universal estimator ranking."}
    output_prefix.with_suffix(".json").write_text(json.dumps(result,indent=2)+"\n")
    return result

def main():
    p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_a",type=Path); p.add_argument("artifact_b",type=Path)
    p.add_argument("--output-prefix",type=Path,default=Path("pooled_five_arm")); a=p.parse_args()
    print(json.dumps(pooled_readout([a.artifact_a,a.artifact_b],a.output_prefix),indent=2))
if __name__=="__main__": main()
