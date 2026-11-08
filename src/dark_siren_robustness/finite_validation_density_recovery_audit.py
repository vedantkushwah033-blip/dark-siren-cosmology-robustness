"""Independent auditor for finite-validation density-recovery artifacts.

This module intentionally does not import the benchmark runner, inclusion
fitter, validation splitter, or production KDE estimator.
"""
from __future__ import annotations

import argparse
import hashlib
import json
from pathlib import Path

import numpy as np
import pandas as pd
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.optimize import brentq
from scipy.special import expit
from scipy.stats import gaussian_kde


Z_MIN, Z_MAX = 0.01, 0.40
TOL = 2e-8


def _lambda(z):
    z = np.asarray(z, dtype=float)
    return 1.0 + 2.2*np.exp(-0.5*((z-0.105)/0.025)**2) + 1.6*np.exp(-0.5*((z-0.275)/0.040)**2)


def _center():
    z = np.linspace(Z_MIN, Z_MAX, 30001)
    x = (z-Z_MIN)/(Z_MAX-Z_MIN)
    return float(trapezoid(_lambda(z)*x*x,z)/trapezoid(_lambda(z),z))


def _alpha(gamma, retention, center):
    z=np.linspace(Z_MIN,Z_MAX,30001); x=(z-Z_MIN)/(Z_MAX-Z_MIN); lam=_lambda(z)
    den=trapezoid(lam,z)
    def mean(a):
        return float(trapezoid(lam*expit(a-x+gamma*(x*x-center)),z)/den)
    return float(brentq(lambda a: mean(a)-retention,-30,30))


def _truth(grid):
    merger=1.0+1.5*(grid-Z_MIN)/(Z_MAX-Z_MIN)
    raw=_lambda(grid)*merger
    return raw/trapezoid(raw,grid)


def _independent_kde(grid, sample, probability, multiplier, merger):
    sample=np.asarray(sample,dtype=float)
    if probability is None:
        weights=np.ones(sample.size)
    else:
        weights=1.0/np.asarray(probability,dtype=float)
    scaled=weights/np.max(weights)
    augmented=np.concatenate((sample,2*grid[0]-sample,2*grid[-1]-sample))
    augmented_weights=np.tile(scaled/3.0,3)
    neff=float(np.sum(scaled)**2/np.sum(scaled**2))
    factor=float(multiplier)*neff**(-1.0/5.0)
    kde=gaussian_kde(augmented,weights=augmented_weights,bw_method=lambda _kde:factor)
    raw=np.asarray(kde(grid))*merger
    return raw/trapezoid(raw,grid)


def audit_artifact(artifact_dir: Path) -> dict:
    attempts=pd.read_csv(artifact_dir/"attempt_metrics.csv")
    units=pd.read_csv(artifact_dir/"parent_units.csv")
    summary=json.loads((artifact_dir/"summary.json").read_text(encoding="utf-8"))
    manifest=json.loads((artifact_dir/"manifest.json").read_text(encoding="utf-8"))
    for name, record in manifest["files"].items():
        artifact_path=artifact_dir/record["path"]
        if not artifact_path.is_file():
            raise AssertionError(f"Manifest file missing: {record['path']}")
        payload=artifact_path.read_bytes()
        if len(payload)!=int(record["bytes"]):
            raise AssertionError(f"Manifest byte count mismatch: {record['path']}")
        if hashlib.sha256(payload).hexdigest()!=record["sha256"]:
            raise AssertionError(f"Manifest checksum mismatch: {record['path']}")
    grid_spec=summary["grid"]
    grid=np.linspace(float(grid_spec[0]),float(grid_spec[1]),int(grid_spec[2]))
    merger=1.0+1.5*(grid-Z_MIN)/(Z_MAX-Z_MIN)
    target=_truth(grid)
    if not np.allclose(target,np.asarray(summary["analytic_target_density"]),rtol=1e-10,atol=1e-12):
        raise AssertionError("Stored analytic target does not match independent DGP reconstruction")
    target_json=json.dumps(target.tolist(),separators=(",",":"))
    target_hash=hashlib.sha256(target_json.encode("utf-8")).hexdigest()
    if target_hash != summary["analytic_target_sha256"]:
        raise AssertionError("Analytic target checksum mismatch")

    required={"replicate","gamma","target_retention","n_validation","object_id","redshift",
              "true_inclusion_probability","retention_uniform","validation_rank","retained","role"}
    if not required.issubset(units.columns):
        raise AssertionError(f"Missing raw-unit columns: {sorted(required-units.columns)}")
    if units.duplicated(["replicate","gamma","target_retention","n_validation","object_id"]).any():
        raise AssertionError("Duplicate parent object key")
    center=_center()
    checked_partitions=0
    for key, group in units.groupby(["replicate","gamma","target_retention","n_validation"],sort=False):
        if len(group)!=int(summary["n_parent"]) or group.object_id.nunique()!=len(group):
            raise AssertionError(f"Incomplete parent frame: {key}")
        group=group.sort_values("object_id")
        if not np.array_equal(group.object_id.to_numpy(),np.arange(len(group))):
            raise AssertionError(f"Parent IDs are not complete and contiguous: {key}")
        gamma=float(key[1]); retention=float(key[2]); nval=int(key[3])
        alpha=_alpha(gamma,retention,center)
        x=(group.redshift.to_numpy()-Z_MIN)/(Z_MAX-Z_MIN)
        p=expit(alpha-x+gamma*(x*x-center))
        if not np.allclose(p,group.true_inclusion_probability.to_numpy(),rtol=2e-9,atol=2e-10):
            raise AssertionError(f"Truth probability reconstruction failed: {key}")
        y=(group.retention_uniform.to_numpy()<p).astype(int)
        if not np.array_equal(y,group.retained.to_numpy(dtype=int)):
            raise AssertionError(f"Retention labels do not match stored uniforms: {key}")
        rank=group.validation_rank.to_numpy(dtype=int)
        expected_role=np.where(rank<nval,"VALIDATION",
            np.where(y==1,"ANALYSIS_CATALOGUE","NON_RETAINED_HOLDOUT"))
        if not np.array_equal(expected_role,group.role.to_numpy()):
            raise AssertionError(f"Disjoint partition reconstruction failed: {key}")
        if int(np.sum(expected_role=="VALIDATION"))!=nval:
            raise AssertionError(f"Validation sample size mismatch: {key}")
        checked_partitions+=1

    # Nested validation sets must agree across n values for each parent scenario.
    for key, block in units.groupby(["replicate","gamma","target_retention"],sort=False):
        validation_sets={}
        for nval,g in block.groupby("n_validation"):
            validation_sets[int(nval)]=set(g.loc[g.role=="VALIDATION","object_id"].astype(int))
        ordered=sorted(validation_sets)
        for small,large in zip(ordered,ordered[1:]):
            if not validation_sets[small].issubset(validation_sets[large]):
                raise AssertionError(f"Validation samples are not nested: {key}")

    attempt_keys=["replicate","gamma","target_retention","n_validation"]
    checked_densities=0
    for row in attempts.itertuples(index=False):
        key=(row.replicate,row.gamma,row.target_retention,row.n_validation)
        mask=np.ones(len(units),dtype=bool)
        for col,val in zip(attempt_keys,key):
            mask &= units[col].to_numpy()==val
        frame=units.loc[mask].sort_values("object_id")
        if frame.empty:
            raise AssertionError(f"Attempt has no parent-unit records: {key}")
        z=frame.redshift.to_numpy(dtype=float)
        y=frame.retained.to_numpy(dtype=int)
        role=frame.role.to_numpy()
        valmask=role=="VALIDATION"
        analysis=role=="ANALYSIS_CATALOGUE"
        gamma=float(row.gamma); nval=int(row.n_validation)
        sample=None; prob=None
        if row.status!="SUCCESS":
            if row.density_json and isinstance(row.density_json,str):
                raise AssertionError(f"Failed attempt unexpectedly contains density: {key}")
            if row.arm.startswith("estimated_degree_") and row.fit_status not in {"SUCCESS","FAILED","NON_EVALUABLE"}:
                raise AssertionError(f"Failed estimated arm has inconsistent fit status: {key}")
            continue
        if row.arm=="catalogue_only":
            sample=z[analysis]
        elif row.arm=="known_inclusion_ipw":
            sample=z[analysis]
            prob=frame.true_inclusion_probability.to_numpy()[analysis]
        elif row.arm=="parent_sample_kde":
            sample=z
        elif row.arm.startswith("estimated_degree_"):
            sample=z[analysis]
            beta=np.asarray(json.loads(row.coefficients_json),dtype=float)
            center_fit=float(row.fit_center); scale_fit=float(row.fit_scale)
            xv=(z[valmask]-center_fit)/scale_fit
            design=np.column_stack([xv**j for j in range(len(beta))])
            pval=expit(design@beta)
            if not np.all(np.isfinite(pval)):
                raise AssertionError(f"Non-finite fitted validation probabilities: {key}")
            xa=(sample-center_fit)/scale_fit
            design_a=np.column_stack([xa**j for j in range(len(beta))])
            prob=expit(design_a@beta)
            # Independently reconstruct the observed-information covariance.
            info=design.T@((pval*(1-pval))[:,None]*design)
            cov=np.linalg.inv(info)
            stored_cov=np.asarray(json.loads(row.fit_covariance_json),dtype=float)
            if not np.allclose(cov,stored_cov,rtol=2e-6,atol=2e-8):
                raise AssertionError(f"Fitted covariance reconstruction failed: {key}")
        else:
            raise AssertionError(f"Unknown arm: {row.arm}")

        density=np.asarray(json.loads(row.density_json),dtype=float)
        reconstructed=_independent_kde(
            grid,sample,prob,float(row.bandwidth_multiplier),merger
        )
        if density.shape!=grid.shape or not np.allclose(density,reconstructed,rtol=TOL,atol=TOL):
            raise AssertionError(f"Independent KDE reconstruction failed: {key}/{row.arm}")
        if not np.isclose(trapezoid(density,grid),1.0,rtol=1e-9,atol=1e-9):
            raise AssertionError(f"Density normalization failed: {key}/{row.arm}")
        iae=float(trapezoid(np.abs(density-target),grid))
        if not np.isclose(iae,float(row.iae),rtol=TOL,atol=TOL):
            raise AssertionError(f"IAE reconstruction failed: {key}/{row.arm}")
        checked_densities+=1

    if len(attempts)!=int(summary["attempt_count"]):
        raise AssertionError("Attempt count differs from summary")
    observed_counts={str(k):int(v) for k,v in attempts.status.value_counts().items()}
    if observed_counts!=summary["status_counts"]:
        raise AssertionError("Attempt status ledger differs from summary")
    return {
        "status":"PASS",
        "audit_scope":"raw partition, DGP probabilities, retention labels, nested validation, fitted covariance, KDE, IAE, summary ledger",
        "checked_parent_partitions":checked_partitions,
        "checked_density_estimates":checked_densities,
        "attempt_count":int(len(attempts)),
        "status_counts":observed_counts,
        "manifest_files_verified":int(len(manifest["files"])),
        "limitations":[
            "Numerical artifact reconstruction is not scientific validation of the DGP.",
            "No H0 likelihood, posterior calibration, or coverage is audited.",
            "The auditor independently implements the declared KDE protocol but shares SciPy's gaussian_kde primitive.",
        ],
    }


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--artifact-dir",type=Path,required=True)
    parser.add_argument("--output",type=Path,default=None)
    args=parser.parse_args()
    report=audit_artifact(args.artifact_dir)
    encoded=json.dumps(report,indent=2,sort_keys=True)+"\n"
    if args.output:
        args.output.write_text(encoded,encoding="utf-8")
    print(encoded)


if __name__=="__main__":
    main()
