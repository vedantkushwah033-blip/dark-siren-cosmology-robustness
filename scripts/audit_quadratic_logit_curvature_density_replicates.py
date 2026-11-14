"""Independent outer auditor for paired quadratic-logit curvature density pilots."""
from __future__ import annotations
import argparse, hashlib, importlib.util, json
from pathlib import Path
import numpy as np
import pandas as pd

def _hash(path):
 h=hashlib.sha256()
 with Path(path).open("rb") as f:
  for block in iter(lambda:f.read(1<<20),b""): h.update(block)
 return h.hexdigest()

def _load_child_auditor():
 path=Path(__file__).with_name("audit_quadratic_logit_curvature_density_smoke.py")
 spec=importlib.util.spec_from_file_location("quadratic_logit_child_raw_auditor",path)
 module=importlib.util.module_from_spec(spec)
 assert spec.loader is not None
 spec.loader.exec_module(module)
 return module.audit

def _compare_aggregate(root,filename,frames,keys,values):
 expected=pd.concat(frames,ignore_index=True) if frames else pd.DataFrame(columns=keys+values)
 actual=pd.read_csv(root/filename)
 if len(actual)!=len(expected): raise AssertionError(f"{filename}: aggregate row count mismatch")
 if len(expected)==0: return 0
 if expected.duplicated(keys).any() or actual.duplicated(keys).any():
  raise AssertionError(f"{filename}: duplicate aggregate key")
 expected=expected.sort_values(keys).reset_index(drop=True)
 actual=actual.sort_values(keys).reset_index(drop=True)
 for key in keys:
  if actual[key].tolist()!=expected[key].tolist(): raise AssertionError(f"{filename}: aggregate key mismatch {key}")
 for col in values:
  if pd.api.types.is_numeric_dtype(expected[col]):
   if not np.allclose(actual[col].to_numpy(float),expected[col].to_numpy(float),rtol=2e-10,atol=2e-10,equal_nan=True):
    raise AssertionError(f"{filename}: aggregate values mismatch {col}")
  elif actual[col].tolist()!=expected[col].tolist():
   raise AssertionError(f"{filename}: aggregate values mismatch {col}")
 return len(expected)

def _assert_frame_equivalent(root, filename, expected, key_columns, numeric_columns, categorical_columns=()):
 actual=pd.read_csv(root/filename)
 if len(actual)!=len(expected): raise AssertionError(f"{filename}: row count mismatch")
 if len(expected)==0: return 0
 sort_keys=[c for c in key_columns if c!="reference_gamma"]
 left=expected.copy(); right=actual.copy()
 left["_sort_reference_gamma"]=left.get("reference_gamma",pd.Series(index=left.index,dtype=float)).fillna(-999.0)
 right["_sort_reference_gamma"]=right.get("reference_gamma",pd.Series(index=right.index,dtype=float)).fillna(-999.0)
 sort_keys=sort_keys+["_sort_reference_gamma"]
 left=left.sort_values(sort_keys).reset_index(drop=True)
 right=right.sort_values(sort_keys).reset_index(drop=True)
 for col in key_columns:
  if col=="reference_gamma":
   if not np.allclose(left[col].to_numpy(float),right[col].to_numpy(float),rtol=0,atol=0,equal_nan=True):
    raise AssertionError(f"{filename}: key mismatch {col}")
  elif pd.api.types.is_numeric_dtype(left[col]):
   if not np.array_equal(left[col].to_numpy(),right[col].to_numpy()):
    raise AssertionError(f"{filename}: key mismatch {col}")
  elif left[col].fillna("").astype(str).tolist()!=right[col].fillna("").astype(str).tolist():
   raise AssertionError(f"{filename}: key mismatch {col}")
 for col in categorical_columns:
  if left[col].fillna("").astype(str).tolist()!=right[col].fillna("").astype(str).tolist():
   raise AssertionError(f"{filename}: categorical mismatch {col}")
 for col in numeric_columns:
  if not np.allclose(left[col].to_numpy(float),right[col].to_numpy(float),rtol=2e-10,atol=2e-10,equal_nan=True):
   raise AssertionError(f"{filename}: numeric mismatch {col}")
 return len(expected)

def audit(root:Path):
 root=Path(root).resolve()
 manifest=json.loads((root/"manifest.json").read_text())
 checked=0
 for rel,meta in manifest["files"].items():
  path=root/rel
  if not path.is_file() or path.stat().st_size!=int(meta["bytes"]) or _hash(path)!=meta["sha256"]:
   raise AssertionError(f"outer manifest verification failed: {rel}")
  checked+=1
 actual_files={str(path.relative_to(root)) for path in root.rglob("*") if path.is_file() and path.relative_to(root).as_posix()!="manifest.json"}
 listed_files=set(manifest["files"])
 if actual_files!=listed_files:
  missing=sorted(listed_files-actual_files)
  unlisted=sorted(actual_files-listed_files)
  raise AssertionError(f"outer manifest inventory mismatch: missing={missing[:5]}, unmanifested={unlisted[:5]}")
 if checked!=int(manifest["file_count"]): raise AssertionError("outer manifest file count mismatch")
 cfg=json.loads((root/"config.json").read_text())
 if manifest.get("status") not in {"PASS", "INCOMPLETE"}: raise AssertionError("invalid outer manifest status")
 if cfg.get("experiment")!="paired_quadratic_logit_curvature_density_pilot": raise AssertionError("wrong experiment identifier")
 if cfg.get("status")!="DEVELOPMENT_PILOT_NOT_CONFIRMATORY_NOT_BANDWIDTH_SELECTION": raise AssertionError("claim status mismatch")
 if cfg.get("analysis_model")!="degree_one_logistic": raise AssertionError("analysis model mismatch")
 if float(cfg.get("coordinate_center",float("nan")))!=0.205 or float(cfg.get("coordinate_scale",float("nan")))!=0.195: raise AssertionError("coordinate specification mismatch")
 if cfg.get("coefficient_prior")!={"family":"independent_gaussian","mean":[0.0,0.0],"sd":[5.0,2.5]}: raise AssertionError("coefficient prior mismatch")
 gammas=[float(x) for x in cfg["gammas"]]
 reps=int(cfg["replicates"]); stride=int(cfg["seed_stride"]); master=int(cfg["master_seed"])
 if reps<2 or len(gammas)<2 or len(set(gammas))!=len(gammas) or stride<2: raise AssertionError("invalid planned design")
 if master<0 or not np.isfinite(gammas).all() or not any(g==0.0 for g in gammas):
  raise AssertionError("invalid seed or gamma reference design")
 widths=[float(x) for x in cfg["absolute_widths"]]
 sizes=[int(x) for x in cfg["validation_sizes"]]
 if len(widths)<2 or widths!=sorted(set(widths)) or not np.isfinite(widths).all() or any(x<=0 for x in widths):
  raise AssertionError("invalid absolute bandwidth design")
 if not sizes or sizes!=sorted(set(sizes)) or sizes[0]<5 or sizes[-1]>int(cfg["validation_frame"]):
  raise AssertionError("invalid nested validation-size design")
 if int(cfg["n_parent"])<100 or int(cfg["validation_frame"])<sizes[-1]:
  raise AssertionError("invalid parent or validation-frame size")
 if int(cfg["quadrature_order"])<3 or int(cfg["grid_points"])<51:
  raise AssertionError("invalid numerical resolution")
 if not np.isfinite(float(cfg["beta"])) or not np.isfinite(float(cfg["scott_multiplier"])) or float(cfg["scott_multiplier"])<=0:
  raise AssertionError("invalid model or bandwidth parameter")
 ledger=pd.read_csv(root/"cell_ledger.csv").fillna("")
 expected=[(r,g) for r in range(reps) for g in gammas]
 got=list(zip(ledger.replicate.astype(int),ledger.gamma.astype(float)))
 if got!=expected: raise AssertionError("cell ledger does not enumerate exact replicate x gamma design")
 if set(ledger.status)-{"SUCCESS","FAILED"}: raise AssertionError("invalid lifecycle status")
 if set(ledger.audit_status)-{"PENDING_OUTER_AUDIT","NOT_RUN"}: raise AssertionError("invalid child audit lifecycle status")
 if not (ledger.loc[ledger.status=="SUCCESS","audit_status"]=="PENDING_OUTER_AUDIT").all(): raise AssertionError("successful cells must await the independent outer child audit")
 if not (ledger.loc[ledger.status=="FAILED","audit_status"]=="NOT_RUN").all(): raise AssertionError("failed cells cannot claim a child audit")
 if int(cfg["attempted_cells"])!=len(expected): raise AssertionError("attempted cell count mismatch")
 if int(cfg["successful_cells"])!=int((ledger.status=="SUCCESS").sum()) or int(cfg["failed_cells"])!=int((ledger.status=="FAILED").sum()):
  raise AssertionError("config totals disagree with cell ledger")
 child_audit=_load_child_auditor()
 metric_frames=[]; diagnostic_frames=[]; scott_frames=[]; child_reports=[]
 paired_frames_by_replicate={}
 for row in ledger.to_dict("records"):
  rid=int(row["replicate"]); gamma=float(row["gamma"])
  seed=master+rid*stride
  if int(row["population_seed"])!=seed or int(row["validation_seed"])!=seed+1: raise AssertionError("paired seed schedule mismatch")
  label=f"gamma_{gamma:+.2f}".replace("+","p").replace("-","m").replace(".","_")
  child=root/"cells"/label/f"replicate_{rid:04d}"
  if row["status"]=="FAILED":
   if int(row["metric_rows"])!=0 or not row["error"] or not child.is_dir(): raise AssertionError("failed cell lacks preserved lifecycle record")
   continue
  child_cfg=json.loads((child/"config.json").read_text())
  if int(child_cfg["population_seed"])!=seed or int(child_cfg["validation_seed"])!=seed+1: raise AssertionError("child seed mismatch")
  if float(child_cfg["gamma"])!=gamma or int(child_cfg["n_parent"])!=int(cfg["n_parent"]): raise AssertionError("child DGP/config mismatch")
  if list(child_cfg["validation_sizes"])!=list(cfg["validation_sizes"]): raise AssertionError("child validation-size mismatch")
  if list(child_cfg["absolute_bandwidth_candidates"])!=list(cfg["absolute_widths"]): raise AssertionError("child width-grid mismatch")
  if int(child_cfg["validation_frame"])!=int(cfg["validation_frame"]) or float(child_cfg["beta"])!=float(cfg["beta"]): raise AssertionError("child validation-frame or beta mismatch")
  if int(child_cfg["quadrature_order"])!=int(cfg["quadrature_order"]) or int(child_cfg["grid_points"])!=int(cfg["grid_points"]): raise AssertionError("child numerical-resolution mismatch")
  if float(child_cfg["scott_multiplier"])!=float(cfg["scott_multiplier"]): raise AssertionError("child Scott multiplier mismatch")
  if float(child_cfg["coordinate_center"])!=float(cfg["coordinate_center"]) or float(child_cfg["coordinate_scale"])!=float(cfg["coordinate_scale"]): raise AssertionError("child coordinate specification mismatch")
  if child_cfg["coefficient_prior"]!=cfg["coefficient_prior"]: raise AssertionError("child coefficient prior mismatch")
  if child_cfg.get("analysis_model")!=cfg["analysis_model"] or child_cfg.get("arm_C")!=cfg["arm_C"] or child_cfg.get("arm_D_E")!=cfg["arm_D_E"]: raise AssertionError("child analysis or arm-label mismatch")
  # Verify the common-random-number design from persisted raw inputs, not merely
  # from repeated seed labels: the parent/frame draws and validation subsample
  # must be identical across gamma cells within each replicate.
  parent_raw=pd.read_csv(child/"truth"/"parent_truth.csv")
  validation_frame_raw=pd.read_csv(child/"truth"/"validation_frame_truth.csv")
  validation_raw=pd.read_csv(child/"observable"/"validation_max.csv")
  current_common={
   "parent_object_id":parent_raw.object_id.to_numpy(),
   "parent_redshift":parent_raw.redshift.to_numpy(float),
   "validation_frame_id":validation_frame_raw.frame_id.to_numpy(),
   "validation_frame_redshift":validation_frame_raw.redshift.to_numpy(float),
   "validation_selected":validation_frame_raw.selected_for_validation.to_numpy(),
   "validation_id":validation_raw.validation_id.to_numpy(),
   "validation_redshift":validation_raw.redshift.to_numpy(float),
  }
  if rid in paired_frames_by_replicate:
   reference_common=paired_frames_by_replicate[rid]
   for key,values in current_common.items():
    reference_values=reference_common[key]
    if np.issubdtype(np.asarray(values).dtype,np.number) and np.issubdtype(np.asarray(reference_values).dtype,np.number):
     equal=np.array_equal(values,reference_values) if not np.issubdtype(np.asarray(values).dtype,np.floating) else np.allclose(values,reference_values,rtol=0,atol=1e-14)
    else:
     equal=np.array_equal(values,reference_values)
    if not equal: raise AssertionError(f"common-random-number raw input mismatch within replicate {rid}: {key}")
  else:
   paired_frames_by_replicate[rid]=current_common
  report=child_audit(child)
  child_reports.append({"replicate":rid,"gamma":gamma,**report})
  for fn,frames in [("density_metrics.csv",metric_frames),("weight_diagnostics.csv",diagnostic_frames),("scott_widths.csv",scott_frames)]:
   frame=pd.read_csv(child/"results"/fn)
   if fn=="density_metrics.csv" and len(frame)!=int(row["metric_rows"]): raise AssertionError("cell metric count differs from ledger")
   frame.insert(0,"replicate",rid); frame.insert(1,"gamma",gamma); frames.append(frame)
 metric_rows=_compare_aggregate(root,"density_metrics_all_cells.csv",metric_frames,
  ["replicate","gamma","validation_size","absolute_bandwidth_sd","arm","node"],
  ["iae","ise","boundary_iae","status","kish_ess"])
 diag_rows=_compare_aggregate(root,"weight_diagnostics_all_cells.csv",diagnostic_frames,
  ["replicate","gamma","validation_size","arm","node"],
  ["kish_ess","kish_ess_fraction","maximum_normalized_weight","top_1_percent_mass","top_5_percent_mass"])
 scott_rows=_compare_aggregate(root,"scott_widths_all_cells.csv",scott_frames,
  ["replicate","gamma","validation_size","arm","node"],["realized_bandwidth_redshift"])

 # Reconstruct paired gamma-vs-zero and within-gamma estimator contrasts from
 # the independently audited child metric records, not the saved contrast CSV.
 metric_columns=["gamma","replicate","validation_size","absolute_bandwidth_sd","arm","node","iae","ise","boundary_iae","status","kish_ess"]
 metric_all=pd.concat(metric_frames,ignore_index=True) if metric_frames else pd.DataFrame(columns=metric_columns)
 primary_arms=("A_catalogue_only","B_parent_oracle","C_known_probability_ipw",
               "D_map_plugin_ipw","E_posterior_mean_density")
 metric_all=metric_all[metric_all.arm.isin(primary_arms) & (metric_all.node==-1)].copy()
 metric_names=("iae","ise","boundary_iae")
 zero_rows=metric_all[np.isclose(metric_all.gamma,0.0)]
 # Any failed planned replicate x gamma cell breaks the preregistered
 # common-random-number grid. The runner therefore suppresses every paired
 # contrast, including within-gamma arm pairs, for an incomplete design.
 pair_keys=["replicate","validation_size","absolute_bandwidth_sd"]
 expected_contrasts=[]
 for gamma in gammas:
  if gamma==0.0: continue
  current=metric_all[np.isclose(metric_all.gamma,gamma)]
  for arm in primary_arms:
   paired=current[current.arm==arm].merge(zero_rows[zero_rows.arm==arm],
    on=pair_keys,suffixes=("","_reference"),validate="one_to_one")
   for _,row in paired.iterrows():
    for metric in metric_names:
     expected_contrasts.append({"contrast_type":"gamma_vs_zero","gamma":gamma,
      "reference_gamma":0.0,"arm":arm,"reference_arm":arm,"replicate":int(row.replicate),
      "validation_size":int(row.validation_size),"absolute_bandwidth_sd":float(row.absolute_bandwidth_sd),
      "metric":metric,"value":float(row[metric]),"reference_value":float(row[f"{metric}_reference"]),
      "difference":float(row[metric]-row[f"{metric}_reference"])})
 arm_pairs=(("D_map_plugin_ipw","A_catalogue_only"),
  ("D_map_plugin_ipw","C_known_probability_ipw"),
  ("C_known_probability_ipw","B_parent_oracle"),
  ("E_posterior_mean_density","D_map_plugin_ipw"))
 for gamma in gammas:
  at_gamma=metric_all[np.isclose(metric_all.gamma,gamma)]
  for arm,reference_arm in arm_pairs:
   paired=at_gamma[at_gamma.arm==arm].merge(at_gamma[at_gamma.arm==reference_arm],
    on=pair_keys,suffixes=("","_reference"),validate="one_to_one")
   for _,row in paired.iterrows():
    for metric in metric_names:
     expected_contrasts.append({"contrast_type":"arm_pair","gamma":gamma,
      "reference_gamma":float("nan"),"arm":arm,"reference_arm":reference_arm,
      "replicate":int(row.replicate),"validation_size":int(row.validation_size),
      "absolute_bandwidth_sd":float(row.absolute_bandwidth_sd),"metric":metric,
      "value":float(row[metric]),"reference_value":float(row[f"{metric}_reference"]),
      "difference":float(row[metric]-row[f"{metric}_reference"])})
 contrast_columns=["contrast_type","gamma","reference_gamma","arm","reference_arm",
  "replicate","validation_size","absolute_bandwidth_sd","metric","value","reference_value","difference"]
 expected_contrasts_frame=pd.DataFrame(expected_contrasts,columns=contrast_columns)
 # A failed planned gamma/replicate cell breaks the common-random-number
 # grid. Partial intersections must not be summarized as paired evidence.
 if not (ledger.status == "SUCCESS").all():
  expected_contrasts_frame=expected_contrasts_frame.iloc[0:0].copy()
 contrast_rows=_assert_frame_equivalent(root,"paired_contrasts_all_cells.csv",
  expected_contrasts_frame,
  ["contrast_type","gamma","reference_gamma","arm","reference_arm","replicate",
   "validation_size","absolute_bandwidth_sd","metric"],
  ["value","reference_value","difference"])
 group_columns=["contrast_type","gamma","reference_gamma","arm","reference_arm",
  "validation_size","absolute_bandwidth_sd","metric"]
 if len(expected_contrasts_frame):
  expected_summary=expected_contrasts_frame.groupby(group_columns,dropna=False).difference.agg(
   n_pairs="count",mean_difference="mean",sd_difference="std").reset_index()
  expected_summary["mcse"]=expected_summary.sd_difference/expected_summary.n_pairs.pow(0.5)
 else:
  expected_summary=pd.DataFrame(columns=group_columns+["n_pairs","mean_difference","sd_difference","mcse"])
 summary_rows=_assert_frame_equivalent(root,"paired_contrast_summary.csv",expected_summary,
  group_columns,["n_pairs","mean_difference","sd_difference","mcse"])
 failed_count=int((ledger.status=="FAILED").sum())
 expected_manifest_status="PASS" if failed_count==0 else "INCOMPLETE"
 if manifest.get("status")!=expected_manifest_status: raise AssertionError("outer manifest status disagrees with lifecycle ledger")
 return {"status":"PASS" if failed_count==0 else "INCOMPLETE","manifest_files_verified":checked,"planned_cells":len(expected),
  "successful_cells":len(child_reports),"failed_cells":failed_count,
  "child_raw_audits_passed":len(child_reports),"aggregate_metric_rows_reconstructed":metric_rows,
  "aggregate_weight_diagnostic_rows_reconstructed":diag_rows,
  "aggregate_scott_width_rows_reconstructed":scott_rows,
  "paired_contrast_rows_reconstructed":contrast_rows,
  "paired_contrast_summary_rows_reconstructed":summary_rows,
  "scope":"Independent child raw-input audits plus exact outer lifecycle and aggregation reconstruction; synthetic density stage only."}

def main():
 p=argparse.ArgumentParser(description=__doc__); p.add_argument("artifact_dir",type=Path)
 p.add_argument("--output",type=Path); a=p.parse_args()
 report=audit(a.artifact_dir)
 if a.output: a.output.write_text(json.dumps(report,indent=2,sort_keys=True)+"\n")
 print(json.dumps(report,indent=2,sort_keys=True))
 if report["status"] != "PASS": raise SystemExit(2)
if __name__=="__main__": main()
