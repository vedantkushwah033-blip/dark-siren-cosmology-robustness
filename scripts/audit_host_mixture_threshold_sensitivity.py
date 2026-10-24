"""Independent audit and prespecified grid-resolution check for factorial CSV outputs."""
from __future__ import annotations
import argparse, csv, json, math
from pathlib import Path
import numpy as np
from dark_siren_robustness.hardened_inference import posterior_h0, posterior_summary
from dark_siren_robustness.selection import posterior_h0_selected

METHODS=("two_host_uncorrected","two_host_corrected","collapsed_uncorrected","collapsed_corrected")
KEYS=("true_H0","host_z1","host_z2","weight1","weight2","detection_max_mpc")
TOL=2e-8


def read_csv(path):
    with Path(path).open(newline="",encoding="utf-8") as f:
        return list(csv.DictReader(f))


def audit(input_dir, grid_sensitivity=True):
    root=Path(input_dir)
    raw=read_csv(root/"per_run.csv")
    summaries=read_csv(root/"summary.csv")
    manifest=json.loads((root/"experiment_manifest.json").read_text(encoding="utf-8"))
    if len(raw)!=manifest["raw_rows"] or len(raw)!=manifest["cell_count"]*manifest["replicates_per_cell"]:
        raise ValueError("Raw row count disagrees with manifest.")
    cells={}
    for r in raw:
        key=tuple(float(r[k]) for k in KEYS)
        cells.setdefault(key,[]).append(r)
        d=float(r["observed_distance_mpc"])
        if not math.isfinite(d) or not 0 < d < float(r["detection_max_mpc"]):
            raise ValueError("Non-finite or out-of-range generated distance.")
        for method in METHODS:
            for suffix in ("median","lower_68","upper_68","lower_95","upper_95","error","absolute_error"):
                v=float(r[f"{method}_{suffix}"])
                if not math.isfinite(v): raise ValueError(f"Non-finite raw result: {method}_{suffix}")
            if abs(float(r[f"{method}_error"])-(float(r[f"{method}_median"])-float(r["true_H0"])))>TOL:
                raise ValueError("Stored error does not match median minus truth.")
            if abs(float(r[f"{method}_absolute_error"])-abs(float(r[f"{method}_error"])))>TOL:
                raise ValueError("Stored absolute error mismatch.")
    if len(cells)!=108 or any(len(v)!=manifest["replicates_per_cell"] for v in cells.values()):
        raise ValueError("Expected 108 cells with equal replication.")
    smap={(tuple(float(r[k]) for k in KEYS),r["method"]):r for r in summaries}
    checked_metrics=0
    for key,rows in cells.items():
        truth=key[0]
        for method in METHODS:
            err=np.array([float(r[f"{method}_error"]) for r in rows])
            hit68=np.array([float(r[f"{method}_lower_68"])<=truth<=float(r[f"{method}_upper_68"]) for r in rows])
            hit95=np.array([float(r[f"{method}_lower_95"])<=truth<=float(r[f"{method}_upper_95"]) for r in rows])
            got=smap.get((key,method))
            if got is None: raise ValueError(f"Missing summary: {method}, {key}")
            expected={"bias":err.mean(),"mae":np.abs(err).mean(),"rmse":np.sqrt(np.mean(err**2)),
                      "coverage_68":hit68.mean(),"coverage_95":hit95.mean(),
                      "mean_width_68":np.mean([float(r[f"{method}_upper_68"])-float(r[f"{method}_lower_68"]) for r in rows]),
                      "mean_width_95":np.mean([float(r[f"{method}_upper_95"])-float(r[f"{method}_lower_95"]) for r in rows])}
            for metric,val in expected.items():
                if metric not in got or abs(float(got[metric])-float(val))>TOL:
                    raise ValueError(f"Summary mismatch {metric} for {method} cell {key}")
                checked_metrics+=1
        for correction in ("uncorrected","corrected"):
            name=f"paired_collapsed_minus_two_host_absolute_error_{correction}"
            got=smap.get((key,name))
            if got is None: raise ValueError(f"Missing paired contrast summary {name}")
            delta=np.array([float(r[f"collapsed_{correction}_absolute_error"])-float(r[f"two_host_{correction}_absolute_error"]) for r in rows])
            if abs(float(got["mean_paired_change"])-delta.mean())>TOL:
                raise ValueError(f"Paired mean mismatch for {key} / {correction}")
            # Reproduce the declared cell bootstrap deterministically.
            seed=int(got["bootstrap_seed"]); draws=int(got["bootstrap_draws"])
            rng=np.random.default_rng(seed); boot=np.empty(draws)
            for start in range(0,draws,128):
                stop=min(start+128,draws)
                idx=rng.integers(0,len(delta),size=(stop-start,len(delta)))
                boot[start:stop]=delta[idx].mean(axis=1)
            low,high=np.quantile(boot,[.025,.975])
            if abs(float(got["bootstrap_95_low"])-low)>TOL or abs(float(got["bootstrap_95_high"])-high)>TOL:
                raise ValueError(f"Bootstrap endpoints mismatch for {key} / {correction}")
    report={"status":"PASS","raw_rows":len(raw),"summary_rows":len(summaries),"unique_cells":len(cells),
            "replicates_per_cell":manifest["replicates_per_cell"],"recomputed_method_metrics":checked_metrics,
            "paired_contrasts_recomputed":len(cells)*2,"bootstrap_intervals_recomputed":len(cells)*2,
            "checks":["dimensions","detection bounds","finite values","error arithmetic","bias/MAE/RMSE/coverage/widths","paired means","paired bootstrap endpoints"]}
    if grid_sensitivity:
        # Prespecified subset: truth 70, z=(.02,.04), weights=(.7,.3), thresholds 90/120/150.
        # Compare the 2,001-point grid against the primary 4,001-point grid on identical raw data.
        grid_rows=[]
        selected=[r for r in raw if float(r["true_H0"])==70 and float(r["host_z1"])==.02 and
                  float(r["host_z2"])==.04 and float(r["weight1"])==.7 and float(r["weight2"])==.3]
        for threshold in (90.0,120.0,150.0):
            subset=[r for r in selected if float(r["detection_max_mpc"])==threshold]
            if len(subset)!=manifest["replicates_per_cell"]: raise ValueError("Grid sensitivity subset incomplete.")
            grid2001=np.linspace(25,125,2001); grid4001=np.linspace(25,125,4001)
            for method in METHODS:
                diffs={k:[] for k in ("median","lower_68","upper_68","lower_95","upper_95")}
                for r in subset:
                    z=np.array([float(r["host_z1"]),float(r["host_z2"])])
                    w=np.array([float(r["weight1"]),float(r["weight2"])])
                    if method.startswith("collapsed"):
                        z=np.array([float(r["effective_redshift"])]); w=np.array([1.0])
                    d=float(r["observed_distance_mpc"])
                    for grid,label in ((grid2001,"2001"),(grid4001,"4001")):
                        if method.endswith("corrected"):
                            p=posterior_h0_selected(grid,z,d,20.0,threshold,weights=w,prior_min=50,prior_max=90)
                        else:
                            p=posterior_h0(grid,z,d,20.0,weights=w,prior_min=50,prior_max=90)
                        s=posterior_summary(grid,p)
                        for k in diffs:
                            if label=="2001": lo=s[k]
                            else: hi=s[k]
                            if label=="2001": pass
                        if label=="2001": s2001=s
                        else: s4001=s
                    for k in diffs: diffs[k].append(float(s2001[k])-float(s4001[k]))
                grid_rows.append({"true_H0":70,"host_z1":.02,"host_z2":.04,"weight1":.7,"weight2":.3,
                    "detection_max_mpc":threshold,"method":method,"n":len(subset),
                    **{f"mean_{k}_difference_grid2001_minus_4001":float(np.mean(v)) for k,v in diffs.items()},
                    **{f"max_abs_{k}_difference":float(np.max(np.abs(v))) for k,v in diffs.items()}})
        with (root/"grid_sensitivity.csv").open("w",newline="",encoding="utf-8") as f:
            writer=csv.DictWriter(f,fieldnames=list(grid_rows[0]));writer.writeheader();writer.writerows(grid_rows)
        report["grid_sensitivity_rows"]=len(grid_rows)
        report["grid_sensitivity_subset"]="H0=70, z=(0.02,0.04), weights=(0.7,0.3), thresholds=90/120/150"
        report["grid_sensitivity"]="computed from same observations using 2001- and 4001-point grids"
    (root/"independent_audit_report.json").write_text(json.dumps(report,indent=2)+"\n",encoding="utf-8")
    print(json.dumps(report,indent=2))
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input_dir")
    p.add_argument("--skip-grid-sensitivity",action="store_true")
    a=p.parse_args()
    audit(a.input_dir,not a.skip_grid_sensitivity)

if __name__=="__main__": main()
