"""Engineering smoke for common absolute-bandwidth five-arm density recovery.

This stage evaluates density recovery only. It deliberately does not evaluate
GW likelihoods, H0 posteriors, calibration, or estimator superiority.
"""
from __future__ import annotations
import argparse, hashlib, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.integrate import trapezoid

from dark_siren_robustness.generative_selection import (
    generate_synthetic_population, generate_independent_validation_sample,
    latent_galaxy_intensity,
)
from dark_siren_robustness.five_arm_comparison import build_five_arm_weights, weight_concentration_diagnostics
from dark_siren_robustness.catalogue_density import estimate_host_density, reflected_weighted_kde_bandwidth

Z_MIN, Z_MAX = 0.01, 0.40
DEFAULT_WIDTHS = (0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.040)
BOUNDARY_LOWER, BOUNDARY_UPPER = 0.049, 0.361

def _hash(path: Path) -> str:
    h=hashlib.sha256()
    with path.open("rb") as f:
        for b in iter(lambda:f.read(1<<20), b""): h.update(b)
    return h.hexdigest()

def _metric(density, truth, grid):
    delta=np.asarray(density)-truth
    boundary=(grid <= BOUNDARY_LOWER) | (grid >= BOUNDARY_UPPER)
    left=grid <= BOUNDARY_LOWER
    right=grid >= BOUNDARY_UPPER
    return {
      "iae":float(trapezoid(np.abs(delta),grid)),
      "ise":float(trapezoid(delta**2,grid)),
      "boundary_iae":float(trapezoid(np.abs(delta[left]),grid[left])+trapezoid(np.abs(delta[right]),grid[right])),
    }

def run(args):
    out=args.output_dir.resolve()
    out.mkdir(parents=True,exist_ok=True)
    if any(out.iterdir()): raise ValueError("output directory must be empty")
    widths=tuple(sorted(set(float(x) for x in args.absolute_widths)))
    if len(widths)<2 or any(not np.isfinite(x) or x<=0 for x in widths): raise ValueError("need >=2 positive finite widths")
    sizes=list(args.validation_sizes)
    if not sizes or sorted(set(sizes))!=sizes or sizes[0]<5 or sizes[-1]>args.validation_frame: raise ValueError("invalid validation sizes")
    if args.n_parent<100 or args.quadrature_order<3: raise ValueError("invalid parent size or quadrature order")
    (out/"observable").mkdir(); (out/"truth").mkdir(); (out/"results").mkdir()
    pop=generate_synthetic_population(seed=args.seed,n_parent=args.n_parent,n_detected=1,h0_true=70.0,
        sigma_mpc=140.0,detection_threshold_mpc=250.0,z_min=Z_MIN,z_max=Z_MAX,
        catalogue_floor=args.catalogue_floor,catalogue_midpoint=0.22,catalogue_width=0.05)
    val=generate_independent_validation_sample(seed=args.validation_seed,n_validation_frame=args.validation_frame,
        n_validation=sizes[-1],z_min=Z_MIN,z_max=Z_MAX,catalogue_floor=args.catalogue_floor,
        catalogue_midpoint=0.22,catalogue_width=0.05)
    parent=pop["parent_truth"]; catalogue=pop["catalogue_observed"]
    parent.to_csv(out/"truth"/"parent_truth.csv",index=False)
    catalogue.to_csv(out/"observable"/"catalogue_observed.csv",index=False)
    val["validation_frame_truth"].to_csv(out/"truth"/"validation_frame_truth.csv",index=False)
    val["validation_observed"].to_csv(out/"observable"/"validation_max.csv",index=False)
    pz=parent.redshift.to_numpy(float); cz=catalogue.redshift.to_numpy(float)
    base_p=1.0+1.5*(pz-Z_MIN)/(Z_MAX-Z_MIN)
    base_c=1.0+1.5*(cz-Z_MIN)/(Z_MAX-Z_MIN)
    grid=np.linspace(Z_MIN,Z_MAX,args.grid_points)
    merger=1.0+1.5*(grid-Z_MIN)/(Z_MAX-Z_MIN)
    truth_raw=latent_galaxy_intensity(grid)*merger
    truth=truth_raw/float(trapezoid(truth_raw,grid))
    archive={"redshift_grid":grid,"merger_weight_grid":merger,"truth_density":truth,
             "parent_redshift":pz,"catalogue_redshift":cz,"parent_base_weights":base_p,
             "catalogue_base_weights":base_c,"absolute_widths":np.asarray(widths)}
    rows=[]; diag=[]; scott_rows=[]
    scott_multiplier=float(getattr(args, "scott_multiplier", 0.7))
    if not np.isfinite(scott_multiplier) or scott_multiplier <= 0:
        raise ValueError("scott_multiplier must be finite and positive")
    for n in sizes:
        v=val["validation_observed"].iloc[:n]
        v.to_csv(out/"observable"/f"validation_n{n}.csv",index=False)
        arms=build_five_arm_weights(pz,cz,v.redshift.to_numpy(float),v.retained.to_numpy(int),
          catalogue_floor=args.catalogue_floor,catalogue_midpoint=0.22,catalogue_width=0.05,
          degree=1,coordinate_center=0.205,coordinate_scale=0.195,
          quadrature_order=args.quadrature_order,parent_base_weights=base_p,catalogue_base_weights=base_c)
        archive[f"posterior_nodes_n{n}"]=arms.inclusion_posterior.nodes
        archive[f"posterior_weights_n{n}"]=arms.inclusion_posterior.weights
        archive[f"posterior_mode_n{n}"]=arms.inclusion_posterior.mode
        raw_inputs={
          "A_catalogue_only":(cz,arms.arm_a_catalogue_only),
          "B_parent_oracle":(pz,arms.arm_b_parent_oracle),
          "C_known_probability_ipw":(cz,arms.arm_c_known_probability),
          "D_map_plugin_ipw":(cz,arms.arm_d_map_plugin),
        }
        # Record the physical widths induced by the shared-Scott rule on
        # the exact same raw arm weights; these inform the absolute-width grid.
        for name,(sample,weights) in raw_inputs.items():
            scott_rows.append({
              "validation_size":n,"arm":name,"node":-1,
              "realized_bandwidth_redshift":reflected_weighted_kde_bandwidth(
                sample,lower_bound=Z_MIN,upper_bound=Z_MAX,
                sample_weights=weights,scott_multiplier=scott_multiplier)
            })
        for node,weights in enumerate(arms.arm_e_node_weights):
            scott_rows.append({
              "validation_size":n,"arm":"E_uncertainty_integrated_ipw","node":node,
              "realized_bandwidth_redshift":reflected_weighted_kde_bandwidth(
                cz,lower_bound=Z_MIN,upper_bound=Z_MAX,
                sample_weights=weights,scott_multiplier=scott_multiplier)
            })
        for name,(sample,weights) in raw_inputs.items():
            archive[f"weights_{name}_n{n}"]=weights
            diag.append({"validation_size":n,"arm":name,"node":-1,**weight_concentration_diagnostics(weights)})
        archive[f"weights_E_nodes_n{n}"]=arms.arm_e_node_weights
        archive[f"weights_E_quadrature_n{n}"]=arms.arm_e_quadrature_weights
        e_diag=[weight_concentration_diagnostics(w) for w in arms.arm_e_node_weights]
        diag.append({"validation_size":n,"arm":"E_uncertainty_integrated_ipw","node":-1,
          "kish_ess":min(x["kish_ess"] for x in e_diag),
          "kish_ess_fraction":min(x["kish_ess_fraction"] for x in e_diag),
          "maximum_normalized_weight":max(x["maximum_normalized_weight"] for x in e_diag),
          "top_1_percent_mass":max(x["top_1_percent_mass"] for x in e_diag),
          "top_5_percent_mass":max(x["top_5_percent_mass"] for x in e_diag)})
        for width in widths:
            per_arm=[]
            for name,(sample,weights) in raw_inputs.items():
                d=estimate_host_density(grid,sample,sample_weights=weights,merger_weight=merger,absolute_bandwidth_sd=width)
                archive[f"density_{name}_n{n}_w{width:.6f}"]=d
                per_arm.append((name,d,None))
            node_dens=[]
            for j,weights in enumerate(arms.arm_e_node_weights):
                d=estimate_host_density(grid,cz,sample_weights=weights,merger_weight=merger,absolute_bandwidth_sd=width)
                node_dens.append(d)
                rows.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":width,
                  "arm":"E_uncertainty_integrated_ipw","node":j,**_metric(d,truth,grid),
                  "status":"SUCCESS","kish_ess":e_diag[j]["kish_ess"]})
            node_dens=np.asarray(node_dens)
            archive[f"density_E_nodes_n{n}_w{width:.6f}"]=node_dens
            q=arms.arm_e_quadrature_weights
            mean_density=np.tensordot(q,node_dens,axes=(0,0))
            archive[f"density_E_posterior_mean_n{n}_w{width:.6f}"]=mean_density
            mean_metric=_metric(mean_density,truth,grid)
            rows.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":width,
              "arm":"E_posterior_mean_density","node":-1,**mean_metric,"status":"SUCCESS",
              "kish_ess":float("nan")})
            weighted_metric={k:float(np.dot(q,[ _metric(d,truth,grid)[k] for d in node_dens])) for k in ("iae","ise","boundary_iae")}
            rows.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":width,
              "arm":"E_posterior_weighted_node_metric","node":-1,**weighted_metric,"status":"SUCCESS",
              "kish_ess":float("nan")})
            for name,d,_ in per_arm:
                weights=raw_inputs[name][1]
                rows.append({"replicate":0,"validation_size":n,"absolute_bandwidth_sd":width,
                  "arm":name,"node":-1,**_metric(d,truth,grid),"status":"SUCCESS",
                  "kish_ess":float(weight_concentration_diagnostics(weights)["kish_ess"])})
    np.savez_compressed(out/"results"/"density_reconstruction_inputs.npz",**archive)
    pd.DataFrame(rows).to_csv(out/"results"/"density_metrics.csv",index=False)
    pd.DataFrame(diag).to_csv(out/"results"/"weight_diagnostics.csv",index=False)
    pd.DataFrame(scott_rows).to_csv(out/"results"/"scott_widths.csv",index=False)
    config={"schema_version":1,"experiment":"absolute_physical_bandwidth_five_arm_density_smoke",
      "status":"ENGINEERING_SMOKE_NOT_SELECTION_NOT_CONFIRMATORY",
      "seed":args.seed,"validation_seed":args.validation_seed,"n_parent":args.n_parent,
      "n_catalogue":len(catalogue),"validation_frame":args.validation_frame,"validation_sizes":sizes,
      "catalogue_floor":args.catalogue_floor,"catalogue_midpoint":0.22,"catalogue_width":0.05,
      "scott_multiplier":scott_multiplier,
      "coordinate_center":0.205,"coordinate_scale":0.195,"quadrature_order":args.quadrature_order,
      "grid_points":args.grid_points,"redshift_bounds":[Z_MIN,Z_MAX],
      "absolute_bandwidth_candidates":list(widths),
      "boundary_regions":[[Z_MIN,BOUNDARY_LOWER],[BOUNDARY_UPPER,Z_MAX]],
      "analytic_target":"normalized latent_galaxy_intensity(z) * (1+1.5*(z-z_min)/(z_max-z_min))",
      "arm_E_metric":"posterior-weighted expectation of node-specific IAE/ISE/boundary IAE; node rows retained",
      "limitations":["one engineering replicate only","posterior nodes are consumed from the existing inclusion-posterior fitter; this smoke does not independently validate that posterior",
        "no training/holdout bandwidth selection","no H0 likelihood or calibration","synthetic low-z Euclidean methodology only"]}
    (out/"config.json").write_text(json.dumps(config,indent=2,sort_keys=True)+"\n")
    files=sorted(p for p in out.rglob("*") if p.is_file() and p.name!="manifest.json")
    manifest={"status":"PASS","file_count":len(files),"files":{str(p.relative_to(out)):{
      "bytes":p.stat().st_size,"sha256":_hash(p)} for p in files}}
    (out/"manifest.json").write_text(json.dumps(manifest,indent=2,sort_keys=True)+"\n")
    return {"status":"PASS","replicate_count":1,"validation_sizes":sizes,"width_count":len(widths),
      "metric_rows":len(rows),"scott_width_rows":len(scott_rows),"manifest_files":len(files)}
def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--seed",type=int,default=20261114); p.add_argument("--validation-seed",type=int,default=20261115)
    p.add_argument("--n-parent",type=int,default=1200); p.add_argument("--validation-frame",type=int,default=1500)
    p.add_argument("--validation-sizes",type=lambda x:tuple(int(v) for v in x.split(",")),default=(100,300,1000))
    p.add_argument("--catalogue-floor",type=float,default=0.10); p.add_argument("--quadrature-order",type=int,default=5)
    p.add_argument("--grid-points",type=int,default=501)
    p.add_argument("--scott-multiplier",type=float,default=0.7)
    p.add_argument("--absolute-widths",type=lambda x:tuple(float(v) for v in x.split(",")),default=DEFAULT_WIDTHS)
    p.add_argument("--output-dir",type=Path,required=True)
    a=p.parse_args(); print(json.dumps(run(a),indent=2,sort_keys=True))
if __name__=="__main__": main()
