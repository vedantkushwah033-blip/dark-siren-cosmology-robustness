"""Independent orchestration audit for the quadratic-logit H0 calibration pilot."""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd


def sha(path: Path) -> str:
    return hashlib.sha256(path.read_bytes()).hexdigest()



def _assert_summary_matches_ledger(ledger: pd.DataFrame, summary: pd.DataFrame) -> None:
    """Recompute all reported cell counts and metrics from the replicate ledger."""
    keys = ["gamma", "n_validation"]
    if summary.duplicated(keys).any():
        raise AssertionError("duplicate calibration-summary cell")
    expected = ledger[keys].drop_duplicates()
    observed = summary[keys].drop_duplicates()
    paired = expected.merge(observed, on=keys, how="outer", indicator=True)
    if not (paired["_merge"] == "both").all():
        raise AssertionError("calibration summary cell keys do not match ledger")
    count_columns = {"planned": None, "success": "SUCCESS",
                     "non_evaluable": "NON_EVALUABLE", "failed": "FAILED",
                     "not_run": "NOT_RUN"}
    metric_columns = {
        "bias": ("median_bias", "mean"), "rmse": ("median_bias", "rmse"),
        "coverage68": ("covered68", "mean"), "coverage90": ("covered90", "mean"),
        "coverage95": ("covered95", "mean"), "mean_width68": ("width68", "mean"),
        "median_min_kish_ess": ("min_catalogue_kish_ess", "median"),
        "median_max_normalized_weight": ("max_normalized_weight", "median"),
    }
    for row in summary.itertuples(index=False):
        group = ledger[(ledger.gamma.astype(float) == float(row.gamma)) &
                       (ledger.n_validation.astype(int) == int(row.n_validation))]
        for column, status in count_columns.items():
            actual = len(group) if status is None else int((group.status == status).sum())
            if int(getattr(row, column)) != actual:
                raise AssertionError(f"summary {column} mismatch for gamma={row.gamma}, n={row.n_validation}")
        successful = group[group.status == "SUCCESS"]
        for output, (source, operation) in metric_columns.items():
            stored = getattr(row, output)
            if successful.empty:
                if pd.notna(stored):
                    raise AssertionError(f"summary {output} must be null with no successful cells")
                continue
            values = successful[source].to_numpy(dtype=float)
            if not np.isfinite(values).all():
                raise AssertionError(f"non-finite ledger values in {source}")
            if operation == "mean":
                actual = float(np.mean(values))
            elif operation == "median":
                actual = float(np.median(values))
            else:
                actual = float(np.sqrt(np.mean(values ** 2)))
            if not np.isfinite(float(stored)) or not np.isclose(float(stored), actual, rtol=1e-10, atol=1e-12):
                raise AssertionError(f"summary {output} mismatch for gamma={row.gamma}, n={row.n_validation}: stored={stored}, recomputed={actual}")



def _assert_validation_samples_paired(ledger: pd.DataFrame, dirs: dict[str, Path]) -> int:
    """Require identical sampled validation IDs/redshifts across gamma.

    Retention labels are deliberately excluded: gamma changes the inclusion
    mechanism, so those observed outcomes should differ across curvature cells.
    """
    paired_blocks = 0
    for (_, _), block in ledger.groupby(["replicate", "n_validation"]):
        usable = block[block.cell_id.isin(dirs)]
        reference = None
        for row in usable.itertuples(index=False):
            observed = pd.read_csv(dirs[row.cell_id] / "observable/validation_observed.csv")
            sampled_design = observed[["validation_id", "redshift"]]
            if reference is None:
                reference = sampled_design
            else:
                pd.testing.assert_frame_equal(sampled_design, reference)
        if len(usable) > 1:
            paired_blocks += 1
    return paired_blocks


def audit(root: Path) -> dict:
    root=root.resolve()
    planned=pd.read_csv(root/"planned_cells.csv")
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    summary=pd.read_csv(root/"calibration_summary.csv")
    manifest=json.loads((root/"manifest.json").read_text())
    bad=[name for name,digest in manifest.items()
         if not (root/name).is_file() or sha(root/name)!=digest]
    if bad: raise AssertionError(f"aggregate manifest failures: {bad}")
    if planned.cell_id.duplicated().any() or ledger.cell_id.duplicated().any():
        raise AssertionError("duplicate cell identifiers")
    if set(planned.cell_id)!=set(ledger.cell_id):
        raise AssertionError("planned cells and ledger do not match")
    identity=["replicate","gamma","n_validation","analysis_seed","validation_seed","h0_true"]
    joined=planned[["cell_id",*identity]].merge(
        ledger[["cell_id",*identity]],on="cell_id",how="inner",
        suffixes=("_planned","_ledger"),validate="one_to_one")
    for name in identity:
        left=joined[f"{name}_planned"].to_numpy()
        right=joined[f"{name}_ledger"].to_numpy()
        if not np.allclose(left.astype(float),right.astype(float),rtol=0,atol=0):
            raise AssertionError(f"planned/ledger provenance mismatch: {name}")
    allowed={"SUCCESS","NON_EVALUABLE","FAILED","NOT_RUN"}
    if not set(ledger.status).issubset(allowed):
        raise AssertionError("unknown failure-ledger status")
    if (ledger.status=="NOT_RUN").any():
        raise AssertionError("pilot ended with cells not attempted")
    audited=0
    dirs={}
    for row in ledger.itertuples(index=False):
        d=Path(row.artifact_dir)
        if row.status in {"SUCCESS","NON_EVALUABLE"}:
            report=json.loads((d/"audit_report.json").read_text())
            if report.get("status")!="PASS": raise AssertionError(f"child audit failed: {row.cell_id}")
            cfg=json.loads((d/"config.json").read_text())
            if not np.isclose(cfg["inference"]["curvature_gamma"],row.gamma):
                raise AssertionError(f"gamma mismatch: {row.cell_id}")
            if int(cfg["validation_generator_config"]["n_validation"])!=int(row.n_validation):
                raise AssertionError(f"validation-size mismatch: {row.cell_id}")
            if int(cfg["analysis_generator_config"]["seed"])!=int(row.analysis_seed):
                raise AssertionError(f"analysis-seed mismatch: {row.cell_id}")
            if int(cfg["validation_generator_config"]["seed"])!=int(row.validation_seed):
                raise AssertionError(f"validation-seed mismatch: {row.cell_id}")
            if cfg["analysis_generator_config"]["n_catalogue_observed"]!=len(pd.read_csv(d/"observable/analysis_catalogue.csv")):
                raise AssertionError(f"catalogue count mismatch: {row.cell_id}")
            dirs[row.cell_id]=d; audited+=1
        elif row.status=="FAILED":
            # Failed attempts are retained in the ledger; partial raw files remain available.
            if not row.error: raise AssertionError(f"failed cell has no error: {row.cell_id}")
    paired=0
    for (rep,gamma),block in ledger.groupby(["replicate","gamma"]):
        usable=block[block.cell_id.isin(dirs)]
        if len(usable)<2: continue
        usable=usable.sort_values("n_validation")
        first=dirs[usable.iloc[0].cell_id]
        ref_frame=pd.read_csv(first/"truth/validation_frame_truth.csv")
        ref_catalogue=pd.read_csv(first/"observable/analysis_catalogue.csv")
        ref_events=pd.read_csv(first/"observable/gw_events_observed.csv")
        prior_ids=None
        for row in usable.itertuples(index=False):
            d=dirs[row.cell_id]
            pd.testing.assert_frame_equal(pd.read_csv(d/"observable/analysis_catalogue.csv"),ref_catalogue)
            pd.testing.assert_frame_equal(pd.read_csv(d/"observable/gw_events_observed.csv"),ref_events)
            frame=pd.read_csv(d/"truth/validation_frame_truth.csv")
            cols=["frame_id","redshift","selected_for_validation"]
            pd.testing.assert_frame_equal(frame[cols],ref_frame[cols])
            ids=pd.read_csv(d/"observable/validation_observed.csv").validation_id.to_numpy()
            if prior_ids is not None and not np.array_equal(prior_ids,ids[:len(prior_ids)]):
                raise AssertionError("validation SRS samples are not nested prefixes")
            prior_ids=ids
        paired+=1
    # Parent, GW events, and validation-frame sampling must be invariant across gamma.
    # Retention labels/probabilities are intentionally excluded because gamma changes them.
    for rep,block in ledger.groupby("replicate"):
        usable=block[block.cell_id.isin(dirs)]
        if usable.empty: continue
        reference=None
        for row in usable.itertuples(index=False):
            d=dirs[row.cell_id]
            parent=pd.read_csv(d/"truth/analysis_parent_truth.csv")[["object_id","redshift"]]
            events=pd.read_csv(d/"observable/gw_events_observed.csv")
            frame=pd.read_csv(d/"truth/validation_frame_truth.csv")
            frame=frame[["frame_id","redshift","selected_for_validation"]]
            current=(parent,events,frame)
            if reference is None: reference=current
            else:
                pd.testing.assert_frame_equal(parent,reference[0])
                pd.testing.assert_frame_equal(events,reference[1])
                pd.testing.assert_frame_equal(frame,reference[2])
    paired_validation_blocks = _assert_validation_samples_paired(ledger, dirs)
    counts=ledger.status.value_counts().to_dict()
    _assert_summary_matches_ledger(ledger, summary)
    return {"status":"PASS","planned_cells":len(planned),"status_counts":counts,
            "child_audits_verified":audited,"paired_nested_blocks_checked":paired,
            "aggregate_manifest_entries_verified":len(manifest),
            "summary_metrics_recomputed":True,
            "interpretation":"Arithmetic, pairing, nested design and artifact integrity only; not calibration evidence."}


def main():
    p=argparse.ArgumentParser(); p.add_argument("artifact_dir",type=Path); p.add_argument("--output",type=Path)
    a=p.parse_args(); result=audit(a.artifact_dir)
    payload=json.dumps(result,indent=2)
    if a.output: a.output.write_text(payload+"\n")
    print(payload)

if __name__=="__main__": main()
