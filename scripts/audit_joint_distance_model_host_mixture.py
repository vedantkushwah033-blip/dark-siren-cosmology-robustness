"""Independent consistency audit for joint distance/host-mixture CSV artifacts."""
from __future__ import annotations
import argparse, csv, json
from pathlib import Path
import numpy as np
try:\n    from scripts.run_joint_distance_model_host_mixture import stable_seed, TRUTHS, HOSTS, WEIGHTS\nexcept ModuleNotFoundError:\n    from run_joint_distance_model_host_mixture import stable_seed, TRUTHS, HOSTS, WEIGHTS

KEYS=("true_H0","host_z1","host_z2","weight1","weight2")
ARMS=("lcdm_two_host","lowz_two_host","lcdm_collapsed","lowz_collapsed")
PAIRS=(
 ("distance_mismatch_two_host","lowz_two_host","lcdm_two_host"),
 ("distance_mismatch_collapsed","lowz_collapsed","lcdm_collapsed"),
 ("host_collapse_lcdm","lcdm_collapsed","lcdm_two_host"),
 ("host_collapse_lowz","lowz_collapsed","lowz_two_host"),
)

def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))

def audit(root, seed=20261025, atol=1e-10):
    root=Path(root)
    raw=read_csv(root/"per_run.csv")
    summaries=read_csv(root/"summary.csv")
    contrasts=read_csv(root/"paired_contrasts.csv")
    manifest=json.loads((root/"experiment_manifest.json").read_text())
    assert len(raw)==32400, f"raw rows={len(raw)}"
    assert len(summaries)==108, f"summary rows={len(summaries)}"
    assert len(contrasts)==135, f"contrast rows={len(contrasts)}"
    numeric=("true_H0","host_z1","host_z2","weight1","weight2","effective_redshift","replicate",
             "observed_distance_mpc","sigma_mpc","median","lower_68","upper_68","lower_95","upper_95",
             "error_H0","absolute_error_H0","covered_68","covered_95","width_68","width_95")
    for row in raw:
        for k in numeric:
            assert np.isfinite(float(row[k])), (k,row)
        assert 40 <= float(row["median"]) <= 100
        assert float(row["lower_95"]) <= float(row["lower_68"]) <= float(row["median"]) <= float(row["upper_68"]) <= float(row["upper_95"])
    grouped={}
    for row in raw:
        cfg=tuple(float(row[k]) for k in KEYS)
        rep=int(row["replicate"])
        grouped.setdefault(cfg,{}).setdefault(rep,{})[row["inference_arm"]]=row
    assert len(grouped)==27
    for cfg,reps in grouped.items():
        assert len(reps)==300
        for rep,arms in reps.items():
            assert set(arms)==set(ARMS)
            assert len({float(r["observed_distance_mpc"]) for r in arms.values()})==1
            assert len({float(r["true_H0"]) for r in arms.values()})==1
    smap={(tuple(float(r[k]) for k in KEYS),r["arm"]):r for r in summaries}
    c_map={(tuple(float(r[k]) for k in KEYS),r["contrast"]):r for r in contrasts}
    max_summary=max_contrast=max_boot=0.0
    for cfg,reps in grouped.items():
        for arm in ARMS:
            rr=[v[arm] for v in reps.values()]
            e=np.array([float(r["error_H0"]) for r in rr])
            expected={"n":len(e),"bias":float(e.mean()),"mae":float(np.abs(e).mean()),
              "rmse":float(np.sqrt(np.mean(e*e))),"coverage_68":float(np.mean([int(r["covered_68"]) for r in rr])),
              "coverage_95":float(np.mean([int(r["covered_95"]) for r in rr])),
              "mean_width_68":float(np.mean([float(r["width_68"]) for r in rr])),
              "mean_width_95":float(np.mean([float(r["width_95"]) for r in rr]))}
            got=smap[(cfg,arm)]
            for k,v in expected.items():
                diff=abs(float(got[k])-v)
                max_summary=max(max_summary,diff)
                assert diff<=atol,(cfg,arm,k,got[k],v,diff)
        deltas={}
        for name,a,b in PAIRS:
            d=np.array([float(v[a]["absolute_error_H0"])-float(v[b]["absolute_error_H0"]) for v in reps.values()])
            deltas[name]=d
            rng=np.random.default_rng(stable_seed(seed+1,[cfg,name]))
            boot=np.array([d[rng.integers(0,len(d),len(d))].mean() for _ in range(5000)])
            expected={"n":len(d),"mean_delta":float(d.mean()),"bootstrap_95_low":float(np.quantile(boot,.025)),
                      "bootstrap_95_high":float(np.quantile(boot,.975)),"bootstrap_draws":5000}
            got=c_map[(cfg,name)]
            for k,v in expected.items():
                diff=abs(float(got[k])-v)
                if k in ("bootstrap_95_low","bootstrap_95_high"): max_boot=max(max_boot,diff)
                else: max_contrast=max(max_contrast,diff)
                assert diff<=atol,(cfg,name,k,got[k],v,diff)
        interaction=deltas["distance_mismatch_collapsed"]-deltas["distance_mismatch_two_host"]
        rng=np.random.default_rng(stable_seed(seed+1,[cfg,"interaction"]))
        boot=np.array([interaction[rng.integers(0,len(interaction),len(interaction))].mean() for _ in range(5000)])
        expected={"n":len(interaction),"mean_delta":float(interaction.mean()),
                  "bootstrap_95_low":float(np.quantile(boot,.025)),"bootstrap_95_high":float(np.quantile(boot,.975)),
                  "bootstrap_draws":5000}
        got=c_map[(cfg,"distance_mismatch_interaction")]
        for k,v in expected.items():
            diff=abs(float(got[k])-v)
            if k in ("bootstrap_95_low","bootstrap_95_high"): max_boot=max(max_boot,diff)
            else: max_contrast=max(max_contrast,diff)
            assert diff<=atol,(cfg,"interaction",k,got[k],v,diff)
    return {"status":"PASS","raw_rows":len(raw),"summary_rows":len(summaries),"contrast_rows":len(contrasts),
            "cells":len(grouped),"replicates_per_cell":300,"max_summary_difference":max_summary,
            "max_contrast_difference":max_contrast,"max_bootstrap_endpoint_difference":max_boot,
            "manifest_status":manifest.get("status")}

def main():
    p=argparse.ArgumentParser();p.add_argument("artifact_dir");p.add_argument("--seed",type=int,default=20261025)
    a=p.parse_args();print(json.dumps(audit(a.artifact_dir,a.seed),indent=2))
if __name__=="__main__": main()
