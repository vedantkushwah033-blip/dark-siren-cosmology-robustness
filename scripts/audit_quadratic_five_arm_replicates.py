"""Independent audit of paired quadratic-logit five-arm replicate orchestration.

This file does not import the multi-replicate runner. It checks the aggregate
manifest, planned-cell/ledger reconciliation, deterministic seeds, child audit
reports, cross-curvature shared inputs, nested validation samples, and
reconstructs all aggregate metrics and paired contrasts.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd

ARMS=[
    "A_catalogue_only","B_parent_oracle","C_known_probability_ipw",
    "D_map_plugin_ipw","E_uncertainty_integrated_ipw",
]


def _sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()


def _seed(master: int, replicate: int, stream: int) -> int:
    return int(np.random.SeedSequence([master,replicate,stream]).generate_state(1,dtype=np.uint32)[0])


def _expected_contrasts(summaries: pd.DataFrame) -> pd.DataFrame:
    rows=[]
    for (rep,gamma,nval),group in summaries.groupby(["replicate","gamma","validation_size"],sort=True):
        by_arm={str(r.arm):r for r in group.itertuples(index=False)}
        if not all(a in by_arm and by_arm[a].status=="SUCCESS" for a in ARMS):
            continue
        for arm in ARMS[1:]:
            for ref in ("A_catalogue_only","C_known_probability_ipw"):
                x,y=by_arm[arm],by_arm[ref]
                rows.append({"replicate":int(rep),"gamma":float(gamma),"validation_size":int(nval),
                    "arm":arm,"reference_arm":ref,
                    "median_error_difference":float(x.median_error-y.median_error),
                    "absolute_error_difference":float(abs(x.median_error)-abs(y.median_error)),
                    "squared_error_difference":float(x.median_error**2-y.median_error**2),
                    "coverage_68_difference":int(x.q16<=x.h0_true<=x.q84)-int(y.q16<=y.h0_true<=y.q84)})
    return pd.DataFrame(rows,columns=["replicate","gamma","validation_size","arm","reference_arm",
        "median_error_difference","absolute_error_difference","squared_error_difference","coverage_68_difference"])


def _expected_aggregate(summaries: pd.DataFrame, ledger: pd.DataFrame,
                        replicates: int, max_abs_error: float,
                        max_interval_width: float = 40.0) -> pd.DataFrame:
    rows=[]
    for (gamma,nval),_ in ledger.groupby(["gamma","validation_size"],sort=True):
        for arm in ARMS:
            g=summaries[(summaries.gamma==gamma)&(summaries.validation_size==nval)
                &(summaries.arm==arm)&(summaries.status=="SUCCESS")]
            errors=g.median_error.to_numpy(float)
            widths=(g.q84-g.q16).to_numpy(float)
            covered=((g.q16.to_numpy(float)<=g.h0_true.to_numpy(float))
                &(g.h0_true.to_numpy(float)<=g.q84.to_numpy(float)))
            failures=replicates-len(g)
            rows.append({"gamma":float(gamma),"validation_size":int(nval),"arm":arm,
                "planned_replicates":replicates,"successful_replicates":len(g),
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
                "failure_inclusive_width_upper":float((np.sum(widths)+failures*max_interval_width)/replicates)})
    return pd.DataFrame(rows)


def _assert_frame_equal(left: pd.DataFrame,right: pd.DataFrame,label: str) -> None:
    try:
        pd.testing.assert_frame_equal(left.reset_index(drop=True),right.reset_index(drop=True),
            check_dtype=False,check_exact=True)
    except AssertionError as exc:
        raise AssertionError(f"paired input mismatch: {label}: {exc}") from exc


def audit(root: Path) -> dict:
    root=root.resolve()
    cfg=json.loads((root/"config.json").read_text())
    manifest=json.loads((root/"manifest.json").read_text())
    checked=0
    for rel,meta in manifest["files"].items():
        path=root/rel
        if not path.is_file() or path.stat().st_size!=int(meta["bytes"]) or _sha(path)!=meta["sha256"]:
            raise AssertionError(f"aggregate manifest mismatch: {rel}")
        checked+=1
    if checked!=int(manifest["file_count"]):
        raise AssertionError("aggregate manifest file count mismatch")

    planned=pd.read_csv(root/"planned_cells.csv")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    summaries=pd.read_csv(root/"replicate_posterior_summaries.csv")
    diagnostics=pd.read_csv(root/"replicate_weight_diagnostics.csv")
    aggregates=pd.read_csv(root/"aggregate_metrics.csv")
    contrasts=pd.read_csv(root/"paired_contrasts.csv")
    reps=int(cfg["replicates"]); gammas=[float(x) for x in cfg["gammas"]]
    sizes=[int(x) for x in cfg["validation_sizes"]]
    expected_n=reps*len(gammas)*len(sizes)*len(ARMS)
    if len(planned)!=expected_n or len(ledger)!=expected_n:
        raise AssertionError("planned or ledger row count does not match frozen design")
    if planned.cell_id.duplicated().any() or ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate planned/ledger cell ID")
    if set(planned.cell_id)!=set(ledger.cell_id):
        raise AssertionError("planned cell IDs and ledger cell IDs differ")
    identity=["replicate","gamma","validation_size","arm","analysis_seed","validation_seed"]
    joined=planned[["cell_id",*identity]].merge(ledger[["cell_id",*identity]],on="cell_id",
        how="inner",suffixes=("_planned","_ledger"),validate="one_to_one")
    for name in identity:
        if not np.array_equal(joined[f"{name}_planned"].to_numpy(),joined[f"{name}_ledger"].to_numpy()):
            raise AssertionError(f"planned/ledger mismatch: {name}")
    allowed={"SUCCESS","FAILED"}
    if not set(ledger.status).issubset(allowed):
        raise AssertionError("unknown lifecycle status")
    if not set(ledger.audit_status).issubset({"PASS","NOT_RUN","NOT_PASS"}):
        raise AssertionError("unknown child audit status")

    child_dirs={}
    child_audits=0
    for (rep,gamma),group in ledger.groupby(["replicate","gamma"],sort=True):
        if group.status.nunique()!=1 or group.audit_status.nunique()!=1:
            raise AssertionError("all arms/sizes in a child run must share one lifecycle status")
        status=str(group.status.iloc[0]); audit_status=str(group.audit_status.iloc[0])
        if status=="SUCCESS":
            if audit_status!="PASS":
                raise AssertionError("successful child lacks independent audit PASS")
            rel=f"replicates/rep{int(rep):03d}_g{float(gamma):+g}".replace("+","p").replace("-","m")
            d=root/rel
            report=json.loads((d/"audit_report.json").read_text())
            if report.get("status")!="PASS":
                raise AssertionError(f"child audit report is not PASS: {rel}")
            child_cfg=json.loads((d/"config.json").read_text())
            if (int(child_cfg["seed"])!=_seed(int(cfg["master_seed"]),int(rep),0)
                or int(child_cfg["validation_seed"])!=_seed(int(cfg["master_seed"]),int(rep),1)
                or float(child_cfg["curvature_gamma"])!=float(gamma)
                or int(child_cfg["n_detected"])!=int(cfg["n_detected"])
                or int(child_cfg["n_parent"])!=int(cfg["n_parent"])
                or list(map(float,child_cfg.get("h0_prior_support",[50.0,90.0])))!=list(map(float,cfg.get("h0_prior_support",[50.0,90.0])))):
                raise AssertionError(f"child configuration/seed mismatch: {rel}")
            if list(map(int,child_cfg["validation_sizes"]))!=sizes:
                raise AssertionError(f"child validation sizes mismatch: {rel}")
            child_summary=pd.read_csv(d/"results"/"posterior_summaries.csv")
            expected=summaries[(summaries.replicate==rep)&(summaries.gamma==gamma)]
            cols=["validation_size","arm","status","h0_true","median_error","posterior_median",
                "q025","q05","q16","q84","q95","q975"]
            left=child_summary[cols].sort_values(["validation_size","arm"]).reset_index(drop=True)
            right=expected[cols].sort_values(["validation_size","arm"]).reset_index(drop=True)
            pd.testing.assert_frame_equal(left,right,check_dtype=False,check_exact=False,
                rtol=2e-10,atol=2e-10)
            child_dirs[(int(rep),float(gamma))]=d
            child_audits+=1
        else:
            if audit_status=="PASS" or not str(group.error.iloc[0]):
                raise AssertionError("failed child must retain error and cannot claim audit PASS")

    # Within replicate, parent redshifts, GW events, validation frame and sampled
    # validation IDs/redshifts are paired across curvature; retention differs.
    pairing_blocks=0
    for rep in range(reps):
        available=[(rep,g) for g in gammas if (rep,g) in child_dirs]
        if len(available)<2:
            continue
        base=child_dirs[available[0]]
        parent0=pd.read_csv(base/"truth"/"parent_truth.csv")[["object_id","redshift"]]
        events0=pd.read_csv(base/"observable"/"gw_events_observed.csv")
        event_truth0=pd.read_csv(base/"truth"/"gw_event_truth.csv")
        frame0=pd.read_csv(base/"truth"/"validation_frame_truth.csv")[["frame_id","redshift"]]
        for key in available[1:]:
            d=child_dirs[key]
            _assert_frame_equal(parent0,pd.read_csv(d/"truth"/"parent_truth.csv")[["object_id","redshift"]],"parent")
            _assert_frame_equal(events0,pd.read_csv(d/"observable"/"gw_events_observed.csv"),"GW events")
            _assert_frame_equal(event_truth0,pd.read_csv(d/"truth"/"gw_event_truth.csv"),"GW event truth")
            _assert_frame_equal(frame0,pd.read_csv(d/"truth"/"validation_frame_truth.csv")[["frame_id","redshift"]],"validation frame")
        for n in sizes:
            reference=pd.read_csv(base/"observable"/f"validation_n{n}.csv")[["validation_id","redshift"]]
            for key in available[1:]:
                current=pd.read_csv(child_dirs[key]/"observable"/f"validation_n{n}.csv")[["validation_id","redshift"]]
                _assert_frame_equal(reference,current,f"validation IDs/redshifts n={n}")
            pairing_blocks+=1

    # Independently reconstruct the top-level seed schedule.
    for row in planned[["replicate","analysis_seed","validation_seed"]].drop_duplicates().itertuples(index=False):
        if int(row.analysis_seed)!=_seed(int(cfg["master_seed"]),int(row.replicate),0):
            raise AssertionError("analysis seed schedule mismatch")
        if int(row.validation_seed)!=_seed(int(cfg["master_seed"]),int(row.replicate),1):
            raise AssertionError("validation seed schedule mismatch")

    expected_contrasts=_expected_contrasts(summaries)
    pd.testing.assert_frame_equal(
        contrasts.sort_values(list(contrasts.columns)).reset_index(drop=True),
        expected_contrasts.sort_values(list(expected_contrasts.columns)).reset_index(drop=True),
        check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-12)
    max_error=float(cfg["failure_bounds"]["max_abs_median_error"])
    max_width=float(cfg["failure_bounds"].get("max_interval_width",40.0))
    expected_aggregates=_expected_aggregate(summaries,ledger,reps,max_error,max_width)
    pd.testing.assert_frame_equal(
        aggregates.sort_values(list(aggregates.columns)).reset_index(drop=True),
        expected_aggregates.sort_values(list(expected_aggregates.columns)).reset_index(drop=True),
        check_dtype=False,check_exact=False,rtol=1e-10,atol=1e-12)

    return {"status":"PASS","planned_cells":len(planned),
        "status_counts":ledger.status.value_counts().to_dict(),
        "successful_child_runs":child_audits,
        "aggregate_manifest_entries_verified":checked,
        "cross_curvature_paired_validation_blocks":pairing_blocks,
        "paired_contrast_rows_reconstructed":len(expected_contrasts),
        "aggregate_metric_rows_reconstructed":len(expected_aggregates),
        "interpretation":"Lifecycle, pairing, child-audit attestations and aggregate arithmetic only; not calibration evidence."}


def main():
    p=argparse.ArgumentParser(); p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
    a=p.parse_args(); result=audit(a.artifact_dir); rendered=json.dumps(result,indent=2,sort_keys=True)+"\n"
    if a.output: a.output.write_text(rendered)
    print(rendered)


if __name__=="__main__": main()
