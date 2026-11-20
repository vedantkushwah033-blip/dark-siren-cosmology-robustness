"""Exact H0-grid refinement for the 30-event synthetic five-arm design.

Unlike the PCHIP-only audit, this runner reevaluates the selected likelihood at
every H0 point using the original simulator, KDE, and selection-likelihood path.
All grid resolutions use the same synthetic seeds and are independently audited.
"""
from __future__ import annotations
import hashlib, json, subprocess, sys
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.interpolate import PchipInterpolator

SEED, VALIDATION_SEED = 2026112101, 2026112102
GAMMAS, DEGREES, GRIDS = (-1.0, 0.0, 1.0), (1, 2), (281, 561, 1121)
VALIDATION_SIZES, ORDERS = (100, 500), (9,)
ARMS = ("A_catalogue_only", "B_parent_oracle", "C_known_probability_ipw",
        "D_map_plugin_ipw", "E_uncertainty_integrated_ipw")

def sha(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda:f.read(1024*1024),b""): h.update(block)
    return h.hexdigest()

def posterior(grid, loglike):
    shifted=np.asarray(loglike,float)-np.max(loglike)
    raw=np.exp(shifted); area=float(trapezoid(raw,grid))
    if not np.isfinite(area) or area<=0: raise ValueError("posterior normalization failed")
    density=raw/area
    cdf=np.r_[0.0,cumulative_trapezoid(density,grid)]; cdf/=cdf[-1]
    q=np.interp([.025,.16,.5,.84,.975],cdf,grid)
    return density,q

def compare_frames(reference: Path, current: Path):
    checked=0
    for rel in ("observable/catalogue_observed.csv","observable/gw_events_observed.csv",
                "observable/validation_n100.csv","observable/validation_n500.csv",
                "truth/parent_truth.csv","truth/gw_event_truth.csv",
                "truth/validation_frame_truth.csv"):
        a=pd.read_csv(reference/rel); b=pd.read_csv(current/rel)
        pd.testing.assert_frame_equal(a,b,check_dtype=False,check_exact=True)
        checked+=1
    return checked

def run(out: Path):
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise ValueError("output directory must be empty")
    runner=Path(__file__).with_name("run_quadratic_five_arm_smoke.py")
    auditor=Path(__file__).with_name("audit_quadratic_five_arm_smoke.py")
    records=[]; inputs_checked=0; child_dirs={}
    for gamma in GAMMAS:
        reference=None
        for degree in DEGREES:
            for points in GRIDS:
                label=f"gamma_{gamma:+.0f}_degree{degree}_grid{points}"
                child=out/label
                command=[sys.executable,str(runner),
                    "--seed",str(SEED),"--validation-seed",str(VALIDATION_SEED),
                    "--n-parent","1200","--n-detected","30","--validation-frame","500",
                    "--validation-sizes","100,500","--curvature-gamma",str(gamma),
                    "--curvature-beta","-3.9","--inclusion-degree",str(degree),
                    "--h0-true","70","--h0-prior-min","30","--h0-prior-max","110",
                    "--h0-points",str(points),"--scott-multiplier","0.7",
                    "--quadrature-order","9","--sigma-mpc","140",
                    "--detection-threshold-mpc","250","--output-dir",str(child)]
                proc=subprocess.run(command,capture_output=True,text=True)
                (out/f"{label}.stdout.txt").write_text(proc.stdout)
                (out/f"{label}.stderr.txt").write_text(proc.stderr)
                if proc.returncode: raise RuntimeError(f"{label} runner failed: {proc.stderr[-3000:]}")
                audit_path=child/"audit_report.json"
                audited=subprocess.run([sys.executable,str(auditor),str(child),"--output",str(audit_path)],
                    capture_output=True,text=True)
                (out/f"{label}.audit.stdout.txt").write_text(audited.stdout)
                (out/f"{label}.audit.stderr.txt").write_text(audited.stderr)
                if audited.returncode: raise RuntimeError(f"{label} independent audit failed: {audited.stderr[-3000:]}")
                report=json.loads(audit_path.read_text())
                if report.get("status")!="PASS": raise AssertionError(f"{label} audit not PASS")
                cfg=json.loads((child/"config.json").read_text())
                expected={"seed":SEED,"validation_seed":VALIDATION_SEED,"n_parent":1200,
                    "n_detected":30,"validation_frame":500,"validation_sizes":[100,500],
                    "inclusion_model_degree":degree,"curvature_gamma":gamma,
                    "curvature_beta":-3.9,"h0_true":70.0,"h0_prior_support":[30.0,110.0],
                    "h0_grid_points":points,"quadrature_order":9,"scott_multiplier":0.7,
                    "sigma_mpc":140.0,"detection_threshold_mpc":250.0}
                for key,value in expected.items():
                    if cfg.get(key)!=value: raise AssertionError(f"{label} config mismatch {key}")
                if reference is None: reference=child
                else: inputs_checked+=compare_frames(reference,child)
                child_dirs[(gamma,degree,points)]=child
                with np.load(child/"results"/"log_likelihoods.npz") as archive:
                    grid=np.asarray(archive["H0_grid"],float)
                    for n in VALIDATION_SIZES:
                        for arm in ARMS:
                            curve=np.asarray(archive[f"{arm}_n{n}"],float)
                            density,q=posterior(grid,curve)
                            records.append({"gamma":gamma,"degree":degree,"grid_points":points,
                                "validation_size":n,"arm":arm,"q025":q[0],"q16":q[1],
                                "median":q[2],"q84":q[3],"q975":q[4],
                                "independent_audit":"PASS"})
    posterior_df=pd.DataFrame(records)
    posterior_df.to_csv(out/"posterior_summaries.csv",index=False)
    comparisons=[]; common=np.linspace(30,110,2241)
    for gamma in GAMMAS:
        for degree in DEGREES:
            for n in VALIDATION_SIZES:
                for arm in ARMS:
                    curves={}
                    for points in GRIDS:
                        child=child_dirs[(gamma,degree,points)]
                        with np.load(child/"results"/"log_likelihoods.npz") as archive:
                            x=np.asarray(archive["H0_grid"],float)
                            y=np.asarray(archive[f"{arm}_n{n}"],float)
                            d,q=posterior(x,y)
                            interp=PchipInterpolator(x,y-np.max(y),extrapolate=False)
                            raw=np.exp(interp(common)); raw/=trapezoid(raw,common)
                            curves[points]=(raw,q)
                    for lower,higher in ((281,561),(561,1121),(281,1121)):
                        low,ql=curves[lower]; high,qh=curves[higher]
                        comparisons.append({"gamma":gamma,"degree":degree,"validation_size":n,
                            "arm":arm,"lower_grid":lower,"higher_grid":higher,
                            "posterior_density_l1_common_grid":float(trapezoid(np.abs(low-high),common)),
                            "median_abs_shift":abs(float(qh[2]-ql[2])),
                            "q025_abs_shift":abs(float(qh[0]-ql[0])),
                            "q975_abs_shift":abs(float(qh[4]-ql[4]))})
    comparison_df=pd.DataFrame(comparisons)
    comparison_df.to_csv(out/"exact_grid_comparisons.csv",index=False)
    summary=comparison_df.groupby(["arm","lower_grid","higher_grid"]).agg(
        max_density_l1=("posterior_density_l1_common_grid","max"),
        max_median_shift=("median_abs_shift","max"),
        max_q025_shift=("q025_abs_shift","max"),
        max_q975_shift=("q975_abs_shift","max")).reset_index()
    summary.to_csv(out/"exact_grid_comparison_summary.csv",index=False)
    config={"status":"exact_h0_grid_refinement_synthetic_pilot_not_calibration",
        "seeds":[SEED,VALIDATION_SEED],"parent":1200,"detected_events":30,
        "gamma":list(GAMMAS),"degrees":list(DEGREES),"quadrature_orders":list(ORDERS),
        "grid_points":list(GRIDS),"validation_sizes":list(VALIDATION_SIZES),
        "child_runs":len(child_dirs),"independent_child_audits":"all PASS",
        "paired_input_file_checks":inputs_checked,
        "scope":"exact selected-likelihood recomputation at each H0 grid point; one fixed synthetic realization per gamma"}
    (out/"config.json").write_text(json.dumps(config,indent=2)+"\n")
    files=sorted(p for p in out.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS","file_count":len(files),"files":{
        str(p.relative_to(out)):{"bytes":p.stat().st_size,"sha256":sha(p)} for p in files}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    print(summary.to_string(index=False))
    print(json.dumps({"status":"PASS","child_runs":len(child_dirs),
        "posterior_rows":len(posterior_df),"comparison_rows":len(comparison_df),
        "paired_input_checks":inputs_checked,"manifest_files":len(files)},indent=2))
    return config

def main():
    import argparse
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir",type=Path,required=True)
    args=parser.parse_args(); run(args.output_dir.resolve())
if __name__=="__main__": main()
