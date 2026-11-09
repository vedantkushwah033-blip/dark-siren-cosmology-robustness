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
from scipy.stats import gaussian_kde, t as student_t


Z_MIN, Z_MAX = 0.01, 0.40
TOL = 2e-8


def _lambda(z):
    z = np.asarray(z, dtype=float)
    return 1.0 + 2.2*np.exp(-0.5*((z-0.105)/0.025)**2) + 1.6*np.exp(-0.5*((z-0.275)/0.040)**2)


def _center():
    z = np.linspace(Z_MIN, Z_MAX, 20001)
    x = (z-Z_MIN)/(Z_MAX-Z_MIN)
    return float(trapezoid(_lambda(z)*x*x,z)/trapezoid(_lambda(z),z))


def _alpha(gamma, retention, center):
    z=np.linspace(Z_MIN,Z_MAX,20001); x=(z-Z_MIN)/(Z_MAX-Z_MIN); lam=_lambda(z)
    den=trapezoid(lam,z)
    def mean(a):
        return float(trapezoid(lam*expit(a-x+gamma*(x*x-center)),z)/den)
    return float(brentq(lambda a: mean(a)-retention,-30,30))


def _truth(grid):
    merger=1.0+1.5*(grid-Z_MIN)/(Z_MAX-Z_MIN)
    raw=_lambda(grid)*merger
    return raw/trapezoid(raw,grid)


def _expected_rng_frame(seed, replicate, n_parent):
    """Independently reconstruct the three declared, disjoint RNG streams."""
    parent_grid=np.linspace(Z_MIN,Z_MAX,20001)
    intensity=_lambda(parent_grid)
    cdf=np.r_[0.0,cumulative_trapezoid(intensity,parent_grid)]
    cdf/=cdf[-1]
    parent_rng=np.random.default_rng(
        np.random.SeedSequence([int(seed),int(replicate),1])
    )
    expected_redshift=np.interp(parent_rng.random(n_parent),cdf,parent_grid)
    split_rng=np.random.default_rng(
        np.random.SeedSequence([int(seed),int(replicate),2])
    )
    validation_order=split_rng.permutation(n_parent)
    expected_rank=np.empty(n_parent,dtype=np.int64)
    expected_rank[validation_order]=np.arange(n_parent,dtype=np.int64)
    retention_rng=np.random.default_rng(
        np.random.SeedSequence([int(seed),int(replicate),3])
    )
    expected_uniform=retention_rng.random(n_parent)
    return expected_redshift,expected_uniform,expected_rank


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
    original_probability=scaled/np.sum(scaled)
    original_mean=float(original_probability@sample)
    original_denom=1.0-float(original_probability@original_probability)
    original_variance=float(
        np.sum(original_probability*(sample-original_mean)**2)/original_denom
    )
    augmented_probability=augmented_weights/np.sum(augmented_weights)
    augmented_mean=float(augmented_probability@augmented)
    augmented_denom=1.0-float(augmented_probability@augmented_probability)
    augmented_variance=float(
        np.sum(augmented_probability*(augmented-augmented_mean)**2)/augmented_denom
    )
    if min(original_variance,augmented_variance)<=0:
        raise AssertionError("Independent KDE covariance is undefined")
    factor=float(multiplier)*neff**(-1.0/5.0)*np.sqrt(
        original_variance/augmented_variance
    )
    kde=gaussian_kde(augmented,weights=augmented_weights,bw_method=lambda _kde:factor)
    raw=np.asarray(kde(grid))*merger
    return raw/trapezoid(raw,grid)


def _manual_gaussian_kde(grid, sample, probability, multiplier, merger):
    """Direct weighted Gaussian-kernel sum, independent of scipy.gaussian_kde."""
    sample=np.asarray(sample,dtype=float)
    weights=np.ones(sample.size,dtype=float) if probability is None else 1.0/np.asarray(probability,dtype=float)
    if sample.ndim!=1 or sample.size<2 or weights.shape!=sample.shape:
        raise AssertionError("Manual KDE received invalid sample or weights")
    if not np.all(np.isfinite(weights)) or np.any(weights<=0):
        raise AssertionError("Manual KDE received nonpositive/nonfinite weights")
    scaled=weights/np.max(weights)
    neff=float(np.sum(scaled)**2/np.sum(scaled**2))
    augmented=np.concatenate((sample,2.0*grid[0]-sample,2.0*grid[-1]-sample))
    augmented_weights=np.tile(scaled/3.0,3)
    probability_aug=augmented_weights/np.sum(augmented_weights)
    denominator=1.0-float(probability_aug@probability_aug)
    if not np.isfinite(denominator) or denominator<=0:
        raise AssertionError("Manual KDE weighted covariance is undefined")
    original_probability=scaled/np.sum(scaled)
    original_mean=float(original_probability@sample)
    original_denominator=1.0-float(original_probability@original_probability)
    if not np.isfinite(original_denominator) or original_denominator<=0:
        raise AssertionError("Original-sample weighted covariance is undefined")
    original_variance=float(
        np.sum(original_probability*(sample-original_mean)**2)/original_denominator
    )
    bandwidth=float(multiplier)*neff**(-0.2)*np.sqrt(original_variance)
    if not np.isfinite(bandwidth) or bandwidth<=0:
        raise AssertionError("Manual KDE bandwidth is invalid")
    raw=np.empty(grid.size,dtype=float)
    normalizer=np.sqrt(2.0*np.pi)*bandwidth
    for start in range(0,grid.size,128):
        stop=min(start+128,grid.size)
        standardized=(grid[start:stop,None]-augmented[None,:])/bandwidth
        raw[start:stop]=np.exp(-0.5*standardized**2)@probability_aug/normalizer
    weighted=raw*np.asarray(merger,dtype=float)
    area=float(trapezoid(weighted,grid))
    if not np.isfinite(area) or area<=0:
        raise AssertionError("Manual KDE normalization failed")
    return weighted/area


def _reconstruct_weight_diagnostics(probability, n):
    weights=np.ones(n,dtype=float) if probability is None else 1.0/np.asarray(probability,dtype=float)
    scaled=weights/np.max(weights)
    normalized=scaled/np.sum(scaled)
    entropy=-float(np.sum(normalized*np.log(normalized)))
    ordered=np.sort(normalized)[::-1]
    return {
        "kish_ess":float(np.sum(scaled)**2/np.sum(scaled**2)),
        "entropy_ess":float(np.exp(entropy)),
        "max_normalized_weight":float(np.max(normalized)),
        "top_1pct_weight_mass":float(np.sum(ordered[:max(1,int(np.ceil(0.01*n)))])),
        "top_5pct_weight_mass":float(np.sum(ordered[:max(1,int(np.ceil(0.05*n)))])),
    }


def _reconstruct_paired_contrasts(attempts, *, n_replicates, seed, bootstrap_resamples=2000):
    results=[]
    rng=np.random.default_rng(seed)
    keys=["gamma","target_retention","n_validation","bandwidth_multiplier"]
    for key,group in attempts.groupby(keys,sort=True):
        selected=group.loc[group.arm.isin(
            ["catalogue_only","estimated_degree_1_ipw","estimated_degree_2_ipw"]
        )]
        wide=selected.pivot(index="replicate",columns="arm",values="iae")
        if "catalogue_only" not in wide:
            continue
        for degree in (1,2):
            arm=f"estimated_degree_{degree}_ipw"
            if arm not in wide:
                continue
            paired=wide[["catalogue_only",arm]].dropna()
            diff=(paired[arm]-paired["catalogue_only"]).to_numpy(dtype=float)
            n=int(diff.size)

            # Independently reconstruct attempt-status overlap as well as
            # conditional paired accuracy; failed attempts remain visible.
            status_wide=selected.pivot(index="replicate",columns="arm",values="status")
            cat_status=status_wide["catalogue_only"] if "catalogue_only" in status_wide else pd.Series(index=status_wide.index,dtype=object)
            est_status=status_wide[arm] if arm in status_wide else pd.Series(index=status_wide.index,dtype=object)
            both_present=cat_status.notna() & est_status.notna()
            if "n_analysis_catalogue" in selected:
                catalogue_sizes=(
                    selected.loc[selected.arm=="catalogue_only"]
                    .set_index("replicate")["n_analysis_catalogue"]
                    .reindex(status_wide.index)
                    .dropna()
                    .astype(float)
                )
            else:
                catalogue_sizes=pd.Series(dtype=float)
            cat_success=cat_status.eq("SUCCESS")
            est_success=est_status.eq("SUCCESS")
            item={
                "gamma":float(key[0]),"target_retention":float(key[1]),
                "n_validation":int(key[2]),"bandwidth_multiplier":float(key[3]),
                "degree":degree,"contrast":f"{arm}_minus_catalogue_only",
                "n_pairs":n,"n_replicates_planned":int(n_replicates),
                "n_replicates_with_both_attempts":int(both_present.sum()),
                "mean_n_analysis_catalogue":float(catalogue_sizes.mean()) if len(catalogue_sizes) else None,
                "sd_n_analysis_catalogue":float(catalogue_sizes.std(ddof=1)) if len(catalogue_sizes)>1 else None,
                "min_n_analysis_catalogue":int(catalogue_sizes.min()) if len(catalogue_sizes) else None,
                "max_n_analysis_catalogue":int(catalogue_sizes.max()) if len(catalogue_sizes) else None,
                "n_catalogue_size_records":int(len(catalogue_sizes)),
                "n_catalogue_only_success":int(cat_success.sum()),
                "n_catalogue_only_failed":int(cat_status.eq("FAILED").sum()),
                "n_catalogue_only_non_evaluable":int(cat_status.eq("NON_EVALUABLE").sum()),
                "n_catalogue_only_not_run":int(cat_status.eq("NOT_RUN").sum()),
                "n_estimated_success":int(est_success.sum()),
                "n_estimated_failed":int(est_status.eq("FAILED").sum()),
                "n_estimated_non_evaluable":int(est_status.eq("NON_EVALUABLE").sum()),
                "n_estimated_not_run":int(est_status.eq("NOT_RUN").sum()),
                "n_both_success":int((cat_success & est_success).sum()),
                "n_estimated_not_success_catalogue_success":int((cat_success & ~est_success).sum()),
                "n_catalogue_not_success_estimated_success":int((~cat_success & est_success).sum()),
                "n_both_not_success":int((~cat_success & ~est_success).sum()),
                "conditional_pair_fraction_of_planned":float(n/n_replicates) if n_replicates else None,
                "mean_difference":float(np.mean(diff)) if n else None,
                "monte_carlo_se":float(np.std(diff,ddof=1)/np.sqrt(n)) if n>1 else None,
                "bootstrap_95_interval":None,"paired_t_95_interval":None,
            }
            if n>=2:
                draws=rng.integers(0,n,size=(bootstrap_resamples,n))
                boot=diff[draws].mean(axis=1)
                item["bootstrap_95_interval"]=[
                    float(np.quantile(boot,0.025)),float(np.quantile(boot,0.975))
                ]
                se=float(np.std(diff,ddof=1)/np.sqrt(n))
                crit=float(student_t.ppf(0.975,n-1))
                mean=float(np.mean(diff))
                item["paired_t_95_interval"]=[mean-crit*se,mean+crit*se]
            results.append(item)
    return results


def audit_artifact(artifact_dir: Path) -> dict:
    attempts=pd.read_csv(artifact_dir/"attempt_metrics.csv")
    units=pd.read_csv(artifact_dir/"parent_units.csv")
    summary=json.loads((artifact_dir/"summary.json").read_text(encoding="utf-8"))
    manifest=json.loads((artifact_dir/"manifest.json").read_text(encoding="utf-8"))
    config=json.loads((artifact_dir/"config.json").read_text(encoding="utf-8"))
    for key,value in config.items():
        if key not in summary or summary[key]!=value:
            raise AssertionError(f"Immutable config differs from summary: {key}")
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
    attempt_key_columns=["replicate","gamma","target_retention","n_validation",
                         "degree","bandwidth_multiplier","arm"]
    if attempts.duplicated(attempt_key_columns).any():
        raise AssertionError("Duplicate attempt key")

    # Reconstruct the complete planned factorial attempt set. A CSV that has
    # internally consistent summaries but silently omits difficult attempts
    # must fail the independent audit.
    expected_attempt_keys=set()
    for rep in range(int(summary["n_replicates"])):
        for gamma in summary["gamma_values"]:
            for retention in summary["retention_values"]:
                for nval in summary["validation_sizes"]:
                    for multiplier in summary["bandwidth_multipliers"]:
                        for arm in ("catalogue_only","known_inclusion_ipw","parent_sample_kde"):
                            expected_attempt_keys.add((
                                rep,float(gamma),float(retention),int(nval),0,
                                float(multiplier),arm,
                            ))
                        for degree in summary["degrees"]:
                            expected_attempt_keys.add((
                                rep,float(gamma),float(retention),int(nval),int(degree),
                                float(multiplier),f"estimated_degree_{int(degree)}_ipw",
                            ))
    actual_attempt_keys={
        (int(row.replicate),float(row.gamma),float(row.target_retention),
         int(row.n_validation),int(row.degree),float(row.bandwidth_multiplier),
         str(row.arm))
        for row in attempts.itertuples(index=False)
    }
    missing_attempts=expected_attempt_keys-actual_attempt_keys
    unexpected_attempts=actual_attempt_keys-expected_attempt_keys
    if missing_attempts or unexpected_attempts:
        raise AssertionError(
            "Attempt design mismatch: "
            f"{len(missing_attempts)} missing, {len(unexpected_attempts)} unexpected"
        )
    allowed_statuses={"SUCCESS","NON_EVALUABLE","FAILED","NOT_RUN"}
    if not set(attempts.status).issubset(allowed_statuses):
        raise AssertionError("Unknown attempt lifecycle status")
    expected_rng_protocol=(
        "numpy SeedSequence([seed, replicate, stream_id]); "
        "stream_id 1=parent redshifts, 2=validation permutation, 3=retention uniforms"
    )
    if summary.get("rng_stream_protocol")!=expected_rng_protocol:
        raise AssertionError("Unknown or missing independent RNG stream protocol")
    center=_center()
    checked_partitions=0
    expected_frames={
        rep:_expected_rng_frame(
            int(summary["seed_start"]),int(rep),int(summary["n_parent"])
        )
        for rep in range(int(summary["n_replicates"]))
    }
    for key, group in units.groupby(["replicate","gamma","target_retention","n_validation"],sort=False):
        if len(group)!=int(summary["n_parent"]) or group.object_id.nunique()!=len(group):
            raise AssertionError(f"Incomplete parent frame: {key}")
        group=group.sort_values("object_id")
        if not np.array_equal(group.object_id.to_numpy(),np.arange(len(group))):
            raise AssertionError(f"Parent IDs are not complete and contiguous: {key}")
        expected_z,expected_uniform,expected_rank=expected_frames[int(key[0])]
        if not np.allclose(
            group.redshift.to_numpy(dtype=float),expected_z,rtol=0.0,atol=5e-15
        ):
            raise AssertionError(f"Parent redshifts do not match their declared RNG stream: {key}")
        if not np.allclose(
            group.retention_uniform.to_numpy(dtype=float),
            expected_uniform,rtol=0.0,atol=5e-15,
        ):
            raise AssertionError(f"Retention uniforms do not match their independent RNG stream: {key}")
        if not np.array_equal(
            group.validation_rank.to_numpy(dtype=int),expected_rank
        ):
            raise AssertionError(f"Validation permutation does not match its independent RNG stream: {key}")
        rank=group.validation_rank.to_numpy(dtype=int)
        if not np.array_equal(np.sort(rank),np.arange(len(group))):
            raise AssertionError(f"Validation ranks are not a permutation: {key}")
        if np.any((group.redshift.to_numpy()<Z_MIN)|(group.redshift.to_numpy()>Z_MAX)):
            raise AssertionError(f"Parent redshift outside declared support: {key}")
        if np.any((group.retention_uniform.to_numpy()<0)|(group.retention_uniform.to_numpy()>=1)):
            raise AssertionError(f"Retention uniform outside [0,1): {key}")
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

    # Common random numbers require parent frames, retention uniforms and
    # validation ranks to be identical across gamma/retention scenarios.
    for rep, block in units.groupby("replicate",sort=False):
        reference=None
        for _, scenario in block.groupby(["gamma","target_retention","n_validation"],sort=False):
            ordered=scenario.sort_values("object_id")
            values=ordered[["redshift","retention_uniform","validation_rank"]].to_numpy()
            if reference is None:
                reference=values
            elif not np.array_equal(reference,values):
                raise AssertionError(f"Common-random-number frame mismatch in replicate {rep}")

    # Nested validation sets must agree across n values for each parent scenario.
    for key, block in units.groupby(["replicate","gamma","target_retention"],sort=False):
        validation_sets={}
        for nval,g in block.groupby("n_validation"):
            validation_sets[int(nval)]=set(g.loc[g.role=="VALIDATION","object_id"].astype(int))
        ordered=sorted(validation_sets)
        for small,large in zip(ordered,ordered[1:]):
            if not validation_sets[small].issubset(validation_sets[large]):
                raise AssertionError(f"Validation samples are not nested: {key}")

    # Index parent partitions once. Re-scanning the entire raw-unit table for
    # every attempt would make the full confirmatory audit quadratic in the
    # number of scenarios and prohibitively expensive.
    unit_group_lookup={
        (int(key[0]),float(key[1]),float(key[2]),int(key[3])):
            group.sort_values("object_id")
        for key,group in units.groupby(
            ["replicate","gamma","target_retention","n_validation"],sort=False
        )
    }
    checked_densities=0
    manual_checks={}
    manual_validation=min(int(value) for value in summary["validation_sizes"])
    manual_bandwidth=min(float(value) for value in summary["bandwidth_multipliers"])
    for row in attempts.itertuples(index=False):
        if int(row.seed)!=int(summary["seed_start"]):
            raise AssertionError("Attempt root-seed label differs from immutable design seed")
        key=(int(row.replicate),float(row.gamma),float(row.target_retention),int(row.n_validation))
        frame=unit_group_lookup.get(key)
        if frame is None or frame.empty:
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
            if row.arm.startswith("estimated_degree_") and row.fit_status not in {"SUCCESS","FAILED","NON_EVALUABLE","NOT_RUN"}:
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
            degree=int(row.degree)
            if row.fit_status!="SUCCESS" or beta.shape!=(degree+1,):
                raise AssertionError(f"Successful estimated arm has inconsistent fit record: {key}")
            validation_z=z[valmask]
            validation_y=y[valmask].astype(float)
            center_fit=float(row.fit_center); scale_fit=float(row.fit_scale)
            expected_center=float(np.mean(validation_z))
            expected_scale=float(np.std(validation_z,ddof=0))
            if expected_scale<=0:
                expected_scale=1.0
            if not np.isclose(center_fit,expected_center,rtol=1e-12,atol=1e-12):
                raise AssertionError(f"Fitted redshift center mismatch: {key}")
            if not np.isclose(scale_fit,expected_scale,rtol=1e-12,atol=1e-12):
                raise AssertionError(f"Fitted redshift scale mismatch: {key}")
            xv=(validation_z-center_fit)/scale_fit
            design=np.column_stack([xv**j for j in range(degree+1)])
            if np.linalg.matrix_rank(design)!=degree+1:
                raise AssertionError(f"Successful fit has rank-deficient validation design: {key}")
            pval=expit(design@beta)
            if not np.all(np.isfinite(pval)) or np.any((pval<0)|(pval>1)):
                raise AssertionError(f"Invalid fitted validation probabilities: {key}")
            score=design.T@(pval-validation_y)
            score_norm=float(np.linalg.norm(score,ord=np.inf))
            if not np.isclose(score_norm,float(row.fit_gradient_inf_norm),rtol=2e-5,atol=2e-8):
                raise AssertionError(f"Fitted logistic score mismatch: {key}")
            if score_norm>1e-6:
                raise AssertionError(f"Successful fit fails independent score tolerance: {key}")
            xa=(sample-center_fit)/scale_fit
            design_a=np.column_stack([xa**j for j in range(degree+1)])
            prob=expit(design_a@beta)
            if not np.all(np.isfinite(prob)) or np.any((prob<0)|(prob>1)):
                raise AssertionError(f"Invalid analysis-catalogue fitted probabilities: {key}")
            # Independently reconstruct observed information and its diagnostics.
            info=design.T@((pval*(1-pval))[:,None]*design)
            eigenvalues=np.linalg.eigvalsh(info)
            if np.min(eigenvalues)<=0:
                raise AssertionError(f"Successful fit has non-positive information: {key}")
            cov=np.linalg.inv(info)
            stored_cov=np.asarray(json.loads(row.fit_covariance_json),dtype=float)
            if not np.allclose(cov,stored_cov,rtol=2e-6,atol=2e-8):
                raise AssertionError(f"Fitted covariance reconstruction failed: {key}")
            hessian_condition=float(np.max(eigenvalues)/np.min(eigenvalues))
            if not np.isclose(hessian_condition,float(row.fit_hessian_condition_number),rtol=2e-5,atol=2e-8):
                raise AssertionError(f"Fitted information condition number mismatch: {key}")
            if not np.isclose(float(np.max(np.abs(beta))),float(row.fit_max_abs_coefficient),rtol=2e-5,atol=2e-8):
                raise AssertionError(f"Maximum coefficient diagnostic mismatch: {key}")
            if int(np.linalg.matrix_rank(design))!=int(row.fit_design_rank):
                raise AssertionError(f"Fitted design-rank diagnostic mismatch: {key}")
            design_condition=float(np.linalg.cond(design))
            if not np.isclose(design_condition,float(row.fit_design_condition_number),rtol=2e-5,atol=2e-8):
                raise AssertionError(f"Fitted design-condition diagnostic mismatch: {key}")
            objective=float(np.sum(np.logaddexp(0.0,design@beta)-validation_y*(design@beta)))
            if not np.isclose(objective,float(row.fit_objective),rtol=2e-8,atol=2e-8):
                raise AssertionError(f"Fitted logistic objective mismatch: {key}")
            if int(row.n_analysis_catalogue)!=int(np.sum(analysis)):
                raise AssertionError(f"Analysis-catalogue size mismatch: {key}")
            if int(row.n_validation_actual)!=int(np.sum(valmask)):
                raise AssertionError(f"Validation sample size mismatch: {key}")
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
        ise=float(trapezoid((density-target)**2,grid))
        max_error=float(np.max(np.abs(density-target)))
        if not np.isclose(iae,float(row.iae),rtol=TOL,atol=TOL):
            raise AssertionError(f"IAE reconstruction failed: {key}/{row.arm}")
        if (int(row.n_validation)==manual_validation
            and np.isclose(float(row.bandwidth_multiplier),manual_bandwidth)):
            check_key=(float(row.gamma),str(row.arm))
            if check_key not in manual_checks:
                direct=_manual_gaussian_kde(
                    grid,sample,prob,float(row.bandwidth_multiplier),merger
                )
                direct_error=float(np.max(np.abs(density-direct)))
                if not np.allclose(density,direct,rtol=2e-10,atol=2e-10):
                    raise AssertionError(f"Independent direct Gaussian-kernel spot check failed: {key}/{row.arm}")
                manual_checks[check_key]=direct_error
        if not np.isclose(ise,float(row.ise),rtol=TOL,atol=TOL):
            raise AssertionError(f"ISE reconstruction failed: {key}/{row.arm}")
        if not np.isclose(max_error,float(row.max_abs_density_error),rtol=TOL,atol=TOL):
            raise AssertionError(f"Maximum density error reconstruction failed: {key}/{row.arm}")
        weight_metrics=_reconstruct_weight_diagnostics(prob,len(sample))
        for metric,value in weight_metrics.items():
            if not np.isclose(value,float(getattr(row,metric)),rtol=TOL,atol=TOL):
                raise AssertionError(f"Weight diagnostic reconstruction failed: {key}/{row.arm}/{metric}")
        checked_densities+=1

    expected_attempt_count=(
        int(summary["n_replicates"])*len(summary["gamma_values"])*
        len(summary["retention_values"])*len(summary["validation_sizes"])*
        len(summary["bandwidth_multipliers"])*(3+len(summary["degrees"]))
    )
    if len(attempts)!=expected_attempt_count or len(attempts)!=int(summary["attempt_count"]):
        raise AssertionError("Attempt count differs from complete factorial design or summary")
    observed_counts={str(k):int(v) for k,v in attempts.status.value_counts().items()}
    if observed_counts!=summary["status_counts"]:
        raise AssertionError("Attempt status ledger differs from summary")
    success=attempts.loc[attempts.status=="SUCCESS"]
    reconstructed_arm_means={str(k):float(v.iae.mean())
                              for k,v in success.groupby("arm",sort=True)}
    if reconstructed_arm_means.keys()!=summary["conditional_success_mean_iae"].keys():
        raise AssertionError("Conditional-success arm summary keys differ")
    for arm,value in reconstructed_arm_means.items():
        if not np.isclose(value,summary["conditional_success_mean_iae"][arm],
                          rtol=1e-10,atol=1e-12):
            raise AssertionError(f"Conditional-success IAE summary mismatch: {arm}")
    reconstructed_contrasts=_reconstruct_paired_contrasts(
        attempts,n_replicates=int(summary["n_replicates"]),
        seed=int(summary["seed_start"])+9109,
    )
    stored_contrasts=summary["paired_ipw_minus_catalogue_only_contrasts"]
    if len(reconstructed_contrasts)!=len(stored_contrasts):
        raise AssertionError("Paired contrast count differs from independent reconstruction")
    for expected,stored in zip(reconstructed_contrasts,stored_contrasts):
        for field in ("gamma","target_retention","n_validation","bandwidth_multiplier",
                      "degree","n_pairs","n_replicates_planned",
                      "n_replicates_with_both_attempts","min_n_analysis_catalogue",
                      "max_n_analysis_catalogue","n_catalogue_size_records","n_catalogue_only_success",
                      "n_catalogue_only_failed","n_catalogue_only_non_evaluable","n_catalogue_only_not_run",
                      "n_estimated_success","n_estimated_failed","n_estimated_non_evaluable","n_estimated_not_run",
                      "n_both_success","n_estimated_not_success_catalogue_success",
                      "n_catalogue_not_success_estimated_success","n_both_not_success"):
            if expected[field]!=stored[field]:
                raise AssertionError(f"Paired contrast key mismatch: {field}")
        for field in ("mean_n_analysis_catalogue","sd_n_analysis_catalogue",
                      "mean_difference","monte_carlo_se","conditional_pair_fraction_of_planned"):
            a,b=expected[field],stored[field]
            if a is None or b is None:
                if a is not None or b is not None:
                    raise AssertionError(f"Paired contrast mismatch: {field}")
            elif not np.isclose(a,b,rtol=1e-10,atol=1e-12):
                raise AssertionError(f"Paired contrast mismatch: {field}")
        for field in ("bootstrap_95_interval","paired_t_95_interval"):
            a,b=expected[field],stored[field]
            if a is not None and not np.allclose(a,b,rtol=1e-10,atol=1e-12):
                raise AssertionError(f"Paired contrast interval mismatch: {field}")
    return {
        "status":"PASS",
        "audit_scope":"independent RNG streams, raw partition, DGP probabilities, retention labels, nested validation, fitted covariance, KDE, IAE, summary ledger",
        "checked_rng_streams":int(len(expected_frames)),
        "checked_parent_partitions":checked_partitions,
        "checked_density_estimates":checked_densities,
        "attempt_count":int(len(attempts)),
        "status_counts":observed_counts,
        "manifest_files_verified":int(len(manifest["files"])),
        "manual_gaussian_kde_crosschecks":int(len(manual_checks)),
        "manual_max_density_abs_difference":max(manual_checks.values()) if manual_checks else None,
        "manual_crosscheck_cells":[
            {"gamma":key[0],"arm":key[1],"max_abs_difference":value}
            for key,value in sorted(manual_checks.items())
        ],
        "limitations":[
            "Numerical artifact reconstruction is not scientific validation of the DGP.",
            "No H0 likelihood, posterior calibration, or coverage is audited.",
            "Bulk KDE reconstruction uses SciPy gaussian_kde; one deterministic direct Gaussian-kernel cross-check is run per curvature-by-arm cell.",
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
