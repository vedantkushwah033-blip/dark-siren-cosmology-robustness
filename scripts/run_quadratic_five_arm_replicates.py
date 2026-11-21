"""Paired multi-replicate quadratic-logit five-arm 30-event engineering pilot.

Within each replicate, the parent population, GW events, validation-frame
sampling design, and validation IDs are paired across curvature values. Each
child is independently audited before its posterior summaries enter aggregates.
This is an engineering pilot, not confirmatory calibration.
"""
from __future__ import annotations

import argparse
from datetime import datetime, timezone
import hashlib
import json
import subprocess
import sys
from pathlib import Path

import numpy as np
import pandas as pd

from dark_siren_robustness.cell_ledger_store import (
    finalize_interrupted_cells,
    initialize_checkpoint,
    mark_cell_started,
    read_checkpoint,
    record_cell_outcomes,
)

ARMS = [
    "A_catalogue_only",
    "B_parent_oracle",
    "C_known_probability_ipw",
    "D_map_plugin_ipw",
    "E_uncertainty_integrated_ipw",
]


def _sha256(path: Path) -> str:
    h = hashlib.sha256()
    with path.open("rb") as f:
        for block in iter(lambda: f.read(1024 * 1024), b""):
            h.update(block)
    return h.hexdigest()


def _seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([master, replicate, stream]).generate_state(1, dtype=np.uint32)[0])


def _h0_truth(master: int, replicate: int, prior_min: float, prior_max: float) -> tuple[int, float]:
    """Draw one independent H0 truth from a dedicated per-replicate stream."""
    if not np.isfinite(prior_min) or not np.isfinite(prior_max) or prior_min >= prior_max:
        raise ValueError("H0 prior bounds must be finite and increasing")
    seed = _seed(master, replicate, 2)
    truth = float(np.random.default_rng(seed).uniform(prior_min, prior_max))
    return seed, truth


def _paired_contrasts(summaries: pd.DataFrame) -> pd.DataFrame:
    rows = []
    for (rep, gamma, nval), group in summaries.groupby(
        ["replicate", "gamma", "validation_size"], sort=True
    ):
        by_arm = {str(row.arm): row for row in group.itertuples(index=False)}
        if not all(a in by_arm and by_arm[a].status == "SUCCESS" for a in ARMS):
            continue
        for arm in ARMS[1:]:
            for reference in ("A_catalogue_only", "C_known_probability_ipw"):
                x, y = by_arm[arm], by_arm[reference]
                rows.append({
                    "replicate": int(rep), "gamma": float(gamma),
                    "validation_size": int(nval), "arm": arm,
                    "reference_arm": reference,
                    "median_error_difference": float(x.median_error-y.median_error),
                    "absolute_error_difference": float(abs(x.median_error)-abs(y.median_error)),
                    "squared_error_difference": float(x.median_error**2-y.median_error**2),
                    "coverage_68_difference": int(x.q16 <= x.h0_true <= x.q84)
                        - int(y.q16 <= y.h0_true <= y.q84),
                })
    return pd.DataFrame(rows, columns=[
        "replicate","gamma","validation_size","arm","reference_arm",
        "median_error_difference","absolute_error_difference",
        "squared_error_difference","coverage_68_difference",
    ])


def _aggregate(summaries: pd.DataFrame, ledger: pd.DataFrame, replicates: int,
               max_abs_error: float, max_interval_width: float = 40.0) -> pd.DataFrame:
    rows = []
    for (gamma, nval), planned in ledger.groupby(["gamma","validation_size"], sort=True):
        for arm in ARMS:
            group = summaries.loc[
                (summaries.gamma == gamma) & (summaries.validation_size == nval)
                & (summaries.arm == arm) & (summaries.status == "SUCCESS")
            ]
            errors = group.median_error.to_numpy(float)
            widths = (group.q84-group.q16).to_numpy(float)
            covered = ((group.q16.to_numpy(float) <= group.h0_true.to_numpy(float))
                       & (group.h0_true.to_numpy(float) <= group.q84.to_numpy(float)))
            failures = replicates-len(group)
            rows.append({
                "gamma":float(gamma),"validation_size":int(nval),"arm":arm,
                "planned_replicates":replicates,"successful_replicates":len(group),
                "failed_or_non_evaluable":failures,
                "success_only_bias":float(errors.mean()) if len(errors) else np.nan,
                "success_only_rmse":float(np.sqrt(np.mean(errors**2))) if len(errors) else np.nan,
                "success_only_68pct_coverage":float(covered.mean()) if len(covered) else np.nan,
                "success_only_mean_68pct_width":float(widths.mean()) if len(widths) else np.nan,
                "failure_inclusive_rmse_lower":float(np.sqrt(np.sum(errors**2)/replicates)),
                "failure_inclusive_rmse_upper":float(np.sqrt((np.sum(errors**2)+failures*max_abs_error**2)/replicates)),
                "failure_inclusive_coverage_lower":float(np.sum(covered)/replicates),
                "failure_inclusive_coverage_upper":float((np.sum(covered)+failures)/replicates),
                "failure_inclusive_width_lower":float(np.sum(widths)/replicates),
                "failure_inclusive_width_upper":float((np.sum(widths)+failures*max_interval_width)/replicates),
            })
    return pd.DataFrame(rows)


def run(args: argparse.Namespace) -> dict:
    out = args.output_dir.resolve()
    out.mkdir(parents=True, exist_ok=True)
    checkpoint_path = out / "cell_lifecycle_checkpoint.json"
    resume = bool(getattr(args, "resume", False))
    if any(out.iterdir()) and not resume:
        raise ValueError("output directory must be empty unless --resume is specified")
    if resume and not checkpoint_path.is_file():
        raise ValueError("--resume requires an existing cell_lifecycle_checkpoint.json")
    if args.replicates < 2 or args.master_seed < 0:
        raise ValueError("need at least two replicates and a nonnegative master seed")
    inclusion_degree = int(getattr(args, "inclusion_degree", 1))
    if inclusion_degree not in (1, 2):
        raise ValueError("inclusion degree must be 1 or 2")
    if not args.gammas or len(set(args.gammas)) != len(args.gammas):
        raise ValueError("gamma values must be nonempty and unique")
    if sorted(args.gammas) != args.gammas:
        raise ValueError("gamma values must be increasing")
    if not args.validation_sizes or sorted(set(args.validation_sizes)) != args.validation_sizes:
        raise ValueError("validation sizes must be unique and increasing")
    if args.validation_sizes[0] < 5 or args.validation_sizes[-1] > args.validation_frame:
        raise ValueError("validation sizes must be >=5 and no larger than validation frame")
    prior_min=float(getattr(args,"h0_prior_min",50.0))
    prior_max=float(getattr(args,"h0_prior_max",90.0))
    if not np.isfinite(prior_min) or not np.isfinite(prior_max) or prior_min>=prior_max:
        raise ValueError("H0 prior bounds must be finite and increasing")
    prior_predictive=bool(getattr(args,"prior_predictive_h0",False))
    if not prior_predictive and not prior_min<=args.h0_true<=prior_max:
        raise ValueError("H0 truth must lie within the declared prior support")
    h0_schedule={}
    for rep in range(args.replicates):
        if prior_predictive:
            truth_seed, truth = _h0_truth(args.master_seed, rep, prior_min, prior_max)
        else:
            truth_seed, truth = -1, float(args.h0_true)
        h0_schedule[rep]={"h0_truth_seed":truth_seed,"h0_true":truth}
    # Posterior medians and truths both lie in prior support; the full support
    # width is therefore a conservative finite bound for any absolute error.
    max_abs_error=(prior_max-prior_min if prior_predictive else
        max(abs(args.h0_true-prior_min),abs(prior_max-args.h0_true)))
    max_interval_width=prior_max-prior_min
    runner=Path(__file__).with_name("run_quadratic_five_arm_smoke.py")
    auditor=Path(__file__).with_name("audit_quadratic_five_arm_smoke.py")
    planned=[]
    for rep in range(args.replicates):
        analysis_seed=_seed(args.master_seed,rep,0)
        validation_seed=_seed(args.master_seed,rep,1)
        event_seed=_seed(args.master_seed,rep,3)
        validation_sample_seed=_seed(args.master_seed,rep,4)
        h0_truth_seed=int(h0_schedule[rep]["h0_truth_seed"])
        h0_true=float(h0_schedule[rep]["h0_true"])
        for gamma in args.gammas:
            cell=f"rep{rep:03d}_g{gamma:+g}".replace("+","p").replace("-","m")
            for nval in args.validation_sizes:
                for arm in ARMS:
                    planned.append({"cell_id":f"{cell}_n{nval}_{arm}","replicate":rep,
                        "gamma":gamma,"validation_size":nval,"arm":arm,
                        "analysis_seed":analysis_seed,"event_seed":event_seed,
                        "validation_seed":validation_seed,
                        "validation_sample_seed":validation_sample_seed,
                        "h0_truth_seed":h0_truth_seed,"h0_true":h0_true,
                        "n_parent":args.n_parent,"n_detected":args.n_detected,
                        "validation_frame":args.validation_frame,
                        "inclusion_model_degree":inclusion_degree,
                        "h0_prior_min":prior_min,"h0_prior_max":prior_max,
                        "h0_points":args.h0_points,"scott_multiplier":args.scott_multiplier,
                        "quadrature_order":args.quadrature_order,
                        "sigma_mpc":args.sigma_mpc,
                        "detection_threshold_mpc":args.detection_threshold_mpc})
    if resume:
        checkpoint=read_checkpoint(checkpoint_path)
        if checkpoint["planned_cells"] != planned:
            raise ValueError("resume design differs from immutable checkpoint plan")
        finalize_interrupted_cells(checkpoint_path, "Recovered at runner startup before retry.")
    else:
        initialize_checkpoint(checkpoint_path, planned)

    summary_parts=[]
    diagnostic_parts=[]
    audit_rows=[]
    for rep in range(args.replicates):
        analysis_seed=_seed(args.master_seed,rep,0)
        validation_seed=_seed(args.master_seed,rep,1)
        event_seed=_seed(args.master_seed,rep,3)
        validation_sample_seed=_seed(args.master_seed,rep,4)
        h0_truth_seed=int(h0_schedule[rep]["h0_truth_seed"])
        h0_true=float(h0_schedule[rep]["h0_true"])
        for gamma in args.gammas:
            cell=f"rep{rep:03d}_g{gamma:+g}".replace("+","p").replace("-","m")
            child=out/"replicates"/cell
            keys=[(rep,float(gamma),int(nval),str(arm))
                for nval in args.validation_sizes for arm in ARMS]
            current_rows=read_checkpoint(checkpoint_path)["cell_ledger"]
            group_rows=[row for row in current_rows
                if int(row["replicate"])==rep and float(row["gamma"])==float(gamma)]
            terminal={"SUCCESS","NON_EVALUABLE","FAILED"}
            if group_rows and all(row["status"] in terminal for row in group_rows):
                if any(row["status"] in {"SUCCESS","NON_EVALUABLE"} for row in group_rows):
                    if not all(row["status"] in {"SUCCESS","NON_EVALUABLE"}
                        and row["audit_status"]=="PASS" for row in group_rows):
                        raise ValueError("completed child has mixed success and failed cell outcomes")
                    report=json.loads((child/"audit_report.json").read_text())
                    sm=pd.read_csv(child/"results"/"posterior_summaries.csv")
                    sm.insert(0,"replicate",rep); sm.insert(1,"gamma",gamma)
                    sm.insert(2,"analysis_seed",analysis_seed); sm.insert(3,"validation_seed",validation_seed)
                    sm["audit_status"]="PASS"
                    summary_parts.append(sm)
                    wd=pd.read_csv(child/"results"/"weight_diagnostics.csv")
                    wd.insert(0,"replicate",rep); wd.insert(1,"gamma",gamma)
                    diagnostic_parts.append(wd)
                    audit_rows.append({"replicate":rep,"gamma":gamma,"status":"PASS",
                        "manifest_files_verified":report["manifest_files_verified"],
                        "posterior_summary_rows_reconstructed":report["posterior_summary_rows_reconstructed"],
                        "maximum_log_likelihood_absolute_difference":report["maximum_log_likelihood_absolute_difference"]})
                else:
                    audit_rows.append({"replicate":rep,"gamma":gamma,"status":"NOT_PASS",
                        "manifest_files_verified":0,"posterior_summary_rows_reconstructed":0,
                        "maximum_log_likelihood_absolute_difference":np.nan})
                continue
            if any(row["status"] in terminal for row in group_rows):
                raise ValueError("child checkpoint contains a partially terminal cell block")
            now=datetime.now(timezone.utc).isoformat(timespec="microseconds")
            for key in keys:
                mark_cell_started(checkpoint_path,key,now)
            status="SUCCESS"
            audit_status="NOT_RUN"
            error=""
            try:
                cmd=[sys.executable,str(runner),
                    "--seed",str(analysis_seed),"--event-seed",str(event_seed),
                    "--validation-seed",str(validation_seed),
                    "--validation-sample-seed",str(validation_sample_seed),
                    "--n-parent",str(args.n_parent),"--n-detected",str(args.n_detected),
                    "--validation-frame",str(args.validation_frame),
                    "--validation-sizes",",".join(map(str,args.validation_sizes)),
                    "--curvature-gamma",str(gamma),"--curvature-beta",str(args.curvature_beta),
                    "--inclusion-degree",str(inclusion_degree),
                    "--h0-true",str(h0_true),"--h0-prior-min",str(prior_min),
                    "--h0-prior-max",str(prior_max),"--h0-points",str(args.h0_points),
                    "--scott-multiplier",str(args.scott_multiplier),
                    "--quadrature-order",str(args.quadrature_order),
                    "--sigma-mpc",str(args.sigma_mpc),
                    "--detection-threshold-mpc",str(args.detection_threshold_mpc),
                    "--output-dir",str(child)]
                proc=subprocess.run(cmd,capture_output=True,text=True)
                child.parent.mkdir(parents=True,exist_ok=True)
                (child.parent/(cell+".stdout.txt")).write_text(proc.stdout)
                (child.parent/(cell+".stderr.txt")).write_text(proc.stderr)
                if proc.returncode:
                    raise RuntimeError("child runner failed: "+proc.stderr[-2500:])
                audit_path=child/"audit_report.json"
                audit=subprocess.run([sys.executable,str(auditor),str(child),"--output",str(audit_path)],
                    capture_output=True,text=True)
                (child.parent/(cell+".audit.stdout.txt")).write_text(audit.stdout)
                (child.parent/(cell+".audit.stderr.txt")).write_text(audit.stderr)
                if audit.returncode:
                    audit_status="NOT_PASS"
                    raise RuntimeError("independent child audit failed: "+audit.stderr[-2500:])
                report=json.loads(audit_path.read_text())
                if report.get("status")!="PASS":
                    raise RuntimeError("child auditor did not return PASS")
                cfg=json.loads((child/"config.json").read_text())
                if (int(cfg["seed"])!=analysis_seed or int(cfg["validation_seed"])!=validation_seed
                    or float(cfg["curvature_gamma"])!=float(gamma)
                    or int(cfg["n_detected"])!=args.n_detected
                    or int(cfg.get("inclusion_model_degree",1))!=inclusion_degree):
                    raise RuntimeError("child configuration differs from frozen replicate design")
                audit_status="PASS"
                sm=pd.read_csv(child/"results"/"posterior_summaries.csv")
                if len(sm)!=len(args.validation_sizes)*len(ARMS):
                    raise RuntimeError("child posterior-summary row count is incomplete")
                if set(sm["status"].astype(str)) != {"SUCCESS"}:
                    raise RuntimeError("child posterior rows are not all SUCCESS; arm-level NON_EVALUABLE rules are not yet implemented")
                sm.insert(0,"replicate",rep)
                sm.insert(1,"gamma",gamma)
                sm.insert(2,"analysis_seed",analysis_seed)
                sm.insert(3,"validation_seed",validation_seed)
                sm["audit_status"]="PASS"
                wd=pd.read_csv(child/"results"/"weight_diagnostics.csv")
                wd.insert(0,"replicate",rep); wd.insert(1,"gamma",gamma)
                completed=datetime.now(timezone.utc).isoformat(timespec="microseconds")
                record_cell_outcomes(checkpoint_path,[{"key":key,"outcome":{
                    "status":"SUCCESS","audit_status":"PASS","completed_at":completed,
                    "audit_report_path":str(audit_path.relative_to(out)),
                }} for key in keys])
                summary_parts.append(sm)
                diagnostic_parts.append(wd)
                audit_rows.append({"replicate":rep,"gamma":gamma,"status":"PASS",
                    "manifest_files_verified":report["manifest_files_verified"],
                    "posterior_summary_rows_reconstructed":report["posterior_summary_rows_reconstructed"],
                    "maximum_log_likelihood_absolute_difference":report["maximum_log_likelihood_absolute_difference"]})
            except Exception as exc:
                status="FAILED"
                error=f"{type(exc).__name__}: {exc}"[:4000]
                failure_audit_status="NOT_PASS" if audit_status=="PASS" else audit_status
                reason="CHILD_AUDIT_FAILED" if failure_audit_status=="NOT_PASS" else "CHILD_PROCESS_ERROR"
                completed=datetime.now(timezone.utc).isoformat(timespec="microseconds")
                record_cell_outcomes(checkpoint_path,[{"key":key,"outcome":{
                    "status":"FAILED","audit_status":failure_audit_status,
                    "reason_code":reason,"reason_detail":error,
                    "completed_at":completed,
                }} for key in keys])
                audit_rows.append({"replicate":rep,"gamma":gamma,"status":"NOT_PASS",
                    "manifest_files_verified":0,"posterior_summary_rows_reconstructed":0,
                    "maximum_log_likelihood_absolute_difference":np.nan})
    checkpoint=read_checkpoint(checkpoint_path)
    planned_df=pd.DataFrame(checkpoint["planned_cells"])
    ledger_rows=checkpoint["cell_ledger"]
    for row in ledger_rows:
        row["error"]=row.get("reason_detail","")
    ledger_df=pd.DataFrame(ledger_rows)
    summary_columns=["replicate","gamma","analysis_seed","validation_seed","validation_size","arm",
        "status","h0_true","median_error","posterior_median","q025","q05","q16","q84","q95","q975","audit_status"]
    diagnostic_columns=["replicate","gamma","validation_size","arm"]
    summaries=pd.concat(summary_parts,ignore_index=True) if summary_parts else pd.DataFrame(columns=summary_columns)
    diagnostics=pd.concat(diagnostic_parts,ignore_index=True) if diagnostic_parts else pd.DataFrame(columns=diagnostic_columns)
    contrasts=_paired_contrasts(summaries)
    aggregates=_aggregate(summaries,ledger_df,args.replicates,max_abs_error,max_interval_width)
    planned_df.to_csv(out/"planned_cells.csv",index=False)
    ledger_df.to_csv(out/"replicate_ledger.csv",index=False)
    summaries.to_csv(out/"replicate_posterior_summaries.csv",index=False)
    diagnostics.to_csv(out/"replicate_weight_diagnostics.csv",index=False)
    contrasts.to_csv(out/"paired_contrasts.csv",index=False)
    aggregates.to_csv(out/"aggregate_metrics.csv",index=False)
    pd.DataFrame(audit_rows).to_csv(out/"replicate_audits.csv",index=False)
    config={"experiment":"quadratic_logit_matched_five_arm_30event_prior_predictive_stage0" if prior_predictive else "quadratic_logit_matched_five_arm_30event_multireplicate_engineering_pilot",
        "status":"stage0_prior_predictive_engineering_smoke_not_confirmatory" if prior_predictive else "engineering_pilot_not_confirmatory_not_observational",
        "master_seed":args.master_seed,"replicates":args.replicates,
        "gammas":args.gammas,"n_parent":args.n_parent,"n_detected":args.n_detected,
        "validation_frame":args.validation_frame,"validation_sizes":args.validation_sizes,
        "curvature_beta":args.curvature_beta,"inclusion_model_degree":inclusion_degree,
        "h0_truth_mode":"prior_predictive_uniform" if prior_predictive else "fixed",
        "h0_true":None if prior_predictive else float(args.h0_true),
        "h0_prior_support":[prior_min,prior_max],
        "h0_truth_schedule":[{"replicate":r,**h0_schedule[r]} for r in range(args.replicates)],
        "h0_points":args.h0_points,"scott_multiplier":args.scott_multiplier,
        "quadrature_order":args.quadrature_order,"sigma_mpc":args.sigma_mpc,
        "detection_threshold_mpc":args.detection_threshold_mpc,
        "seed_schedule":[{"replicate":r,"population_seed":_seed(args.master_seed,r,0),
            "validation_frame_seed":_seed(args.master_seed,r,1),
            "h0_truth_seed":int(h0_schedule[r]["h0_truth_seed"]),
            "event_seed":_seed(args.master_seed,r,3),
            "validation_sample_seed":_seed(args.master_seed,r,4)} for r in range(args.replicates)],
        "failure_bounds":{"max_abs_median_error":max_abs_error,"max_interval_width":max_interval_width},
        "claim_limit":"Paired engineering pilot only; no estimator ranking, calibration, or observational inference."}
    (out/"config.json").write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
    files=sorted(p for p in out.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS","file_count":len(files),"files":{
        str(p.relative_to(out)):{"bytes":p.stat().st_size,"sha256":_sha256(p)} for p in files}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    counts=ledger_df.status.value_counts().to_dict()
    return {"status":"PASS" if (ledger_df.status=="SUCCESS").all() and (ledger_df.audit_status=="PASS").all() else "INCOMPLETE",
        "planned_cells":len(ledger_df),"status_counts":counts,
        "successful_replicate_gamma_runs":int((pd.DataFrame(audit_rows).status=="PASS").sum()),
        "posterior_summary_rows":len(summaries),"output_dir":str(out)}


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--master-seed",type=int,default=2026111929)
    p.add_argument("--replicates",type=int,default=10)
    p.add_argument("--gammas",type=lambda x:[float(v) for v in x.split(",")],default=[-1.0,0.0,1.0])
    p.add_argument("--n-parent",type=int,default=1200)
    p.add_argument("--n-detected",type=int,default=30)
    p.add_argument("--validation-frame",type=int,default=500)
    p.add_argument("--validation-sizes",type=lambda x:[int(v) for v in x.split(",")],default=[100,300,500])
    p.add_argument("--curvature-beta",type=float,default=-3.9)
    p.add_argument("--inclusion-degree",type=int,choices=(1,2),default=1)
    p.add_argument("--h0-true",type=float,default=70.0)
    p.add_argument("--prior-predictive-h0",action="store_true",
        help="draw one H0 truth per independent replicate from Uniform(prior-min, prior-max)")
    p.add_argument("--h0-prior-min",type=float,default=50.0)
    p.add_argument("--h0-prior-max",type=float,default=90.0)
    p.add_argument("--h0-points",type=int,default=281)
    p.add_argument("--scott-multiplier",type=float,default=0.7)
    p.add_argument("--quadrature-order",type=int,default=5)
    p.add_argument("--sigma-mpc",type=float,default=140.0)
    p.add_argument("--detection-threshold-mpc",type=float,default=250.0)
    p.add_argument("--output-dir",type=Path,required=True)
    p.add_argument("--resume",action="store_true",help="resume unfinished child blocks from the persisted cell checkpoint")
    a=p.parse_args()
    print(json.dumps(run(a),indent=2,sort_keys=True))


if __name__=="__main__":
    main()
