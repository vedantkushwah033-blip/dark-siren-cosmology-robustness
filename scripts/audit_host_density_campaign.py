"""Independent replay and posterior reconstruction of all host-density screen fits.

This audit does not import the package's generator, KDE estimator, or likelihood.
It independently regenerates the parent/catalogue/event data, estimates reflected
weighted KDEs, reconstructs selected H0 posteriors, checks artifact hashes, and
re-aggregates all saved arm-level summaries. All scientific claims remain within
the declared synthetic low-redshift toy model.
"""
from __future__ import annotations
import argparse, csv, hashlib, json, math
from pathlib import Path
import numpy as np
from scipy.integrate import cumulative_trapezoid, simpson
from scipy.special import ndtr

C=299_792.458
SIGMA=140.0
THRESHOLD=250.0
H0_TRUE=70.0
N_PARENT=6000
N_DETECTED=40
Z_LO,Z_HI=0.01,0.40
SAMPLER_N=20001
ARMS=("oracle","catalogue_only","known_completeness")
CONFIGS=("bw_scott","bw_0.05","bw_0.1","bw_0.2","bw_0.5","bw_0.8","bw_1.2",
         "grid_coarse","grid_fine","misspec_mid_0.18","misspec_mid_0.26",
         "misspec_width_0.10","misspec_width_0.025")
TOL={"posterior_mean":0.05,"posterior_median":0.10,"ci68_low":0.10,"ci68_high":0.10,
     "ci68_width":0.10,"l1_host_density_error":0.025,"min_detection_efficiency":0.003,
     "posterior_integral":0.002}


def latent_intensity(z):
    z=np.asarray(z,dtype=float)
    return 1.0+2.2*np.exp(-0.5*((z-0.105)/0.025)**2)+1.6*np.exp(-0.5*((z-0.275)/0.040)**2)


def inclusion(z,floor=0.10,midpoint=0.22,width=0.05):
    return floor+(1.0-floor)/(1.0+np.exp((np.asarray(z)-midpoint)/width))


def draw_parent_redshifts(n,rng,zmin=Z_LO,zmax=Z_HI):
    grid=np.linspace(zmin,zmax,SAMPLER_N)
    intensity=latent_intensity(grid)
    increments=0.5*(intensity[:-1]+intensity[1:])*np.diff(grid)
    cdf=np.r_[0.0,np.cumsum(increments)]
    cdf/=cdf[-1]
    return np.interp(rng.random(n),cdf,grid)


def simulate(seed):
    """Independent implementation of the campaign's fully specified RNG sequence."""
    rng=np.random.default_rng(int(seed))
    parent=draw_parent_redshifts(N_PARENT,rng)
    retained=rng.random(N_PARENT)<inclusion(parent)
    observed=parent[retained]
    if observed.size<2:
        raise RuntimeError("independent replay retained too few catalogue objects")
    merger=1.0+1.5*(parent-Z_LO)/(Z_HI-Z_LO)
    host_prob=merger/merger.sum()
    detected=[]
    attempts=0
    max_attempts=max(100_000,N_DETECTED*2_000)
    while len(detected)<N_DETECTED and attempts<max_attempts:
        batch=min(max(1000,2*(N_DETECTED-len(detected))),max_attempts-attempts)
        ids=rng.choice(N_PARENT,size=batch,replace=True,p=host_prob)
        distance=rng.normal(C*parent[ids]/H0_TRUE,SIGMA)
        accepted=distance[distance>=THRESHOLD]
        if accepted.size:
            detected.extend(accepted[:N_DETECTED-len(detected)].tolist())
        attempts+=batch
    if len(detected)!=N_DETECTED:
        raise RuntimeError("independent replay failed to produce the requested detected-event count")
    return parent,observed,np.asarray(detected,dtype=float)


def reflected_weighted_kde(zg,sample,probabilities,factor):
    """Independent reflected Gaussian-kernel sum; input probabilities are inclusion probs."""
    weights=np.ones(sample.size) if probabilities is None else 1.0/np.asarray(probabilities,dtype=float)
    scaled=weights/np.max(weights)
    total=scaled.sum()
    mean=float(np.sum(scaled*sample)/total)
    variance=float(np.sum(scaled*(sample-mean)**2)/total)
    if not np.isfinite(variance) or variance<=0:
        raise ValueError("independent KDE has invalid weighted variance")
    neff=float(total*total/np.sum(scaled*scaled))
    bw=(neff**(-0.2) if factor is None else float(factor))*math.sqrt(variance)
    if not np.isfinite(bw) or bw<=0:
        raise ValueError("independent KDE bandwidth is not positive and finite")
    centers=np.concatenate((sample,2.0*zg[0]-sample,2.0*zg[-1]-sample))
    center_weights=np.tile(scaled/total/3.0,3)
    density=np.empty(zg.size,dtype=float)
    norm=1.0/(math.sqrt(2.0*math.pi)*bw)
    for start in range(0,zg.size,32):
        stop=min(zg.size,start+32)
        u=(zg[start:stop,None]-centers[None,:])/bw
        density[start:stop]=(np.exp(-0.5*u*u)@center_weights)*norm
    return density,weights


def selected_posterior(zg,density,distances,h0g):
    """Independent selected likelihood/posterior using Simpson quadrature."""
    density=density/simpson(density,x=zg)
    means=C*zg[None,:]/h0g[:,None]
    loglike=np.zeros(h0g.size,dtype=float)
    normalizer=SIGMA*math.sqrt(2.0*math.pi)
    for distance in distances:
        integrand=density[None,:]*np.exp(-0.5*((distance-means)/SIGMA)**2)/normalizer
        event_likelihood=simpson(integrand,x=zg,axis=1)
        if np.any(event_likelihood<=0) or not np.all(np.isfinite(event_likelihood)):
            raise ArithmeticError("independent event likelihood is invalid")
        loglike+=np.log(event_likelihood)
    beta=simpson(density[None,:]*ndtr((means-THRESHOLD)/SIGMA),x=zg,axis=1)
    if np.any(beta<=0) or not np.all(np.isfinite(beta)):
        raise ArithmeticError("independent detection efficiency is invalid")
    loglike-=distances.size*np.log(beta)
    p=np.exp(loglike-loglike.max())
    p/=simpson(p,x=h0g)
    cdf=np.r_[0.0,cumulative_trapezoid(p,x=h0g)]
    cdf/=cdf[-1]
    q16,q50,q84=np.interp([0.16,0.50,0.84],cdf,h0g)
    return {"posterior_mean":float(simpson(h0g*p,x=h0g)),
            "posterior_median":float(q50),"ci68_low":float(q16),"ci68_high":float(q84),
            "ci68_width":float(q84-q16),"posterior_integral":float(simpson(p,x=h0g)),
            "min_detection_efficiency":float(np.min(beta))}


def wilson(k,n,z=1.959964):
    p=k/n
    den=1.0+z*z/n
    center=(p+z*z/(2*n))/den
    half=z*math.sqrt(p*(1.0-p)/n+z*z/(4*n*n))/den
    return [float(max(0.0,center-half)),float(min(1.0,center+half))]


def verify_close(saved,actual,field,tolerance,context,failures,maximum):
    delta=abs(float(saved)-float(actual))
    maximum[field]=max(maximum[field],delta)
    if not np.isfinite(delta) or delta>tolerance:
        failures.append({"context":context,"field":field,"saved":float(saved),"independent":float(actual),
                         "absolute_difference":float(delta),"tolerance":tolerance})


def _read_csv(path):
    with path.open(newline="",encoding="utf-8") as handle:
        return list(csv.DictReader(handle))


def audit(campaign_dir,*,audit_z_points=401,audit_h0_points=81,only_configs=None,replicate_limit=None):
    root=Path(campaign_dir)
    manifest=json.loads((root/"campaign_manifest.json").read_text(encoding="utf-8"))
    summary=json.loads((root/"campaign_summary.json").read_text(encoding="utf-8"))
    names=tuple(summary.get("configuration_names",()))
    if names!=CONFIGS:
        raise ValueError("campaign configuration list differs from locked contract")
    nrep=int(summary.get("replicates_per_configuration",-1))
    seed_start=int(summary.get("seed_start",-1))
    seed_end=int(summary.get("seed_end",-1))
    if nrep!=40 or seed_start!=20261029 or seed_end!=20261068:
        raise ValueError("campaign replicate/seed contract mismatch")
    if summary.get("status")!="synthetic_exploratory_host_density_campaign_not_real_H0_measurement":
        raise ValueError("unexpected scientific scope label")
    expected_files={f"{name}_replicates.csv" for name in CONFIGS}|{"campaign_summary.json"}
    if set(manifest.get("files",{}))!=expected_files:
        raise ValueError("manifest must cover campaign_summary.json and all 13 replicate CSVs")
    for name,record in manifest["files"].items():
        path=root/name
        data=path.read_bytes()
        if len(data)!=int(record.get("bytes",-1)) or hashlib.sha256(data).hexdigest()!=record.get("sha256"):
            raise ValueError(f"SHA-256/byte-count mismatch: {name}")
    if manifest.get("replicates_per_configuration")!=nrep or manifest.get("seed_start")!=seed_start or manifest.get("seed_end")!=seed_end:
        raise ValueError("manifest and summary seed metadata disagree")

    selected=list(CONFIGS if only_configs is None else only_configs)
    if not selected or len(set(selected))!=len(selected) or not set(selected).issubset(CONFIGS):
        raise ValueError("only_configs must be a non-empty unique subset of the locked configuration list")
    n_audit=nrep if replicate_limit is None else int(replicate_limit)
    if n_audit<1 or n_audit>nrep:
        raise ValueError("replicate_limit must be in [1, 40]")
    if audit_z_points<201 or audit_h0_points<41:
        raise ValueError("independent audit grid must be at least 201 redshift × 41 H0 points")

    zg=np.linspace(Z_LO,Z_HI,audit_z_points)
    h0g=np.linspace(50.0,90.0,audit_h0_points)
    merger=1.0+1.5*(zg-zg[0])/(zg[-1]-zg[0])
    truth=latent_intensity(zg)*merger
    truth/=simpson(truth,x=zg)
    simulations={seed_start+r:simulate(seed_start+r) for r in range(n_audit)}
    failures=[]
    borderline=[]
    maximum={key:0.0 for key in TOL}
    raw_by_config={}
    perconfig={}
    checked=0

    for name in selected:
        rows=_read_csv(root/f"{name}_replicates.csv")
        if len(rows)!=nrep*len(ARMS):
            raise ValueError(f"{name}: expected {nrep*len(ARMS)} rows, got {len(rows)}")
        lookup={}
        for row in rows:
            key=(int(row["replicate"]),row["arm"])
            if key in lookup:
                raise ValueError(f"{name}: duplicate replicate/arm row {key}")
            lookup[key]=row
        expected_keys={(rep,arm) for rep in range(nrep) for arm in ARMS}
        if set(lookup)!=expected_keys:
            raise ValueError(f"{name}: replicate/arm keys are incomplete")

        cfg=summary["results"][name]["config"]
        factor=cfg.get("bandwidth_factor")
        midpoint=float(cfg.get("supplied_midpoint",0.22))
        width=float(cfg.get("supplied_width",0.05))
        reconstructed=[]
        for rep in range(n_audit):
            seed=seed_start+rep
            parent,observed,distances=simulations[seed]
            probs=inclusion(observed,midpoint=midpoint,width=width)
            samples={"oracle":(parent,None),"catalogue_only":(observed,None),
                     "known_completeness":(observed,probs)}
            for arm in ARMS:
                sample,probability=samples[arm]
                galaxy_density,weights=reflected_weighted_kde(zg,sample,probability,factor)
                density=galaxy_density*merger
                density/=simpson(density,x=zg)
                metrics=selected_posterior(zg,density,distances,h0g)
                l1=float(simpson(np.abs(density-truth),x=zg))
                weight_ratio=float(np.max(weights)/np.mean(weights))
                ess=float(weights.sum()**2/np.sum(weights*weights))
                ess_fraction=ess/weights.size
                saved=lookup[(rep,arm)]
                context={"configuration":name,"replicate":rep,"arm":arm}
                metrics["l1_host_density_error"]=l1
                for field,tolerance in TOL.items():
                    verify_close(saved[field],metrics[field],field,tolerance,context,failures,maximum)
                for field,value in (("n_parent",N_PARENT),("n_catalogue_observed",observed.size),("n_detected",N_DETECTED)):
                    if int(saved[field])!=int(value):
                        failures.append({"context":context,"field":field,"saved":int(saved[field]),"independent":int(value)})
                if int(saved["seed"])!=seed:
                    failures.append({"context":context,"field":"seed","saved":int(saved["seed"]),"independent":seed})
                saved_covered=str(saved["covered_68"]).strip().lower() in ("true","1")
                independent_covered=metrics["ci68_low"]<=H0_TRUE<=metrics["ci68_high"]
                if saved_covered!=independent_covered:
                    edge=min(abs(float(saved["ci68_low"])-H0_TRUE),abs(float(saved["ci68_high"])-H0_TRUE),
                             abs(metrics["ci68_low"]-H0_TRUE),abs(metrics["ci68_high"]-H0_TRUE))
                    if edge<=max(TOL["ci68_low"],TOL["ci68_high"]):
                        borderline.append({"context":context,"saved":saved_covered,"independent":bool(independent_covered),
                                           "nearest_endpoint_distance":float(edge)})
                    else:
                        failures.append({"context":context,"field":"covered_68","saved":saved_covered,
                                         "independent":bool(independent_covered),"nearest_endpoint_distance":float(edge)})
                reconstructed.append({"replicate":rep,"arm":arm,"posterior_mean":metrics["posterior_mean"],
                    "covered_68":bool(metrics["ci68_low"]<=H0_TRUE<=metrics["ci68_high"]),
                    "ci68_width":metrics["ci68_width"],"l1_host_density_error":l1,
                    "max_weight_ratio":weight_ratio,"ess_fraction":ess_fraction})
                checked+=1
        raw_by_config[name]=rows
        perconfig[name]={"rows_independently_reconstructed":len(reconstructed),
                         "failures_after_configuration":len(failures)}
        print(f"independently reconstructed {name}: {len(reconstructed)} rows; failures={len(failures)}",flush=True)

    # Independent re-aggregation of all saved headline statistics and fixed-seed
    # bootstrap intervals. This checks raw CSV -> summary rather than trusting JSON.
    aggregate_checks=0
    max_aggregate_diff=0.0
    for name in selected:
        saved_arms={item["arm"]:item for item in summary["results"][name]["arms"]}
        rows=raw_by_config[name]
        byarm={arm:[row for row in rows if row["arm"]==arm] for arm in ARMS}
        errors={arm:np.array([float(row["posterior_mean"])-H0_TRUE for row in byarm[arm]]) for arm in ARMS}
        bootstrap_rng=np.random.default_rng(1)
        for arm in ARMS:
            grp=byarm[arm]
            e=errors[arm]
            n=len(e)
            covered=sum(str(row["covered_68"]).strip().lower() in ("true","1") for row in grp)
            expected={"n":n,"bias":float(e.mean()),"bias_se":float(e.std(ddof=1)/np.sqrt(n)),
                "rmse":float(np.sqrt(np.mean(e*e))),"coverage68":covered/n,"cov_wilson95":wilson(covered,n),
                "width":float(np.mean([float(row["ci68_width"]) for row in grp])),
                "l1_density":float(np.mean([float(row["l1_host_density_error"]) for row in grp])),
                "max_weight_ratio_max":float(np.max([float(row["max_weight_ratio"]) for row in grp])),
                "ess_fraction_min":float(np.min([float(row["ess_fraction"]) for row in grp]))}
            saved=saved_arms[arm]
            for field in ("bias","bias_se","rmse","coverage68","width","l1_density","max_weight_ratio_max","ess_fraction_min"):
                delta=abs(float(saved[field])-expected[field])
                max_aggregate_diff=max(max_aggregate_diff,delta); aggregate_checks+=1
                if not np.isfinite(delta) or delta>1e-9:
                    failures.append({"context":{"configuration":name,"arm":arm},"field":"aggregate_"+field,
                                    "saved":saved[field],"independent":expected[field],"absolute_difference":delta})
            delta=max(abs(float(a)-float(b)) for a,b in zip(saved["cov_wilson95"],expected["cov_wilson95"]))
            max_aggregate_diff=max(max_aggregate_diff,delta); aggregate_checks+=1
            if delta>1e-9:
                failures.append({"context":{"configuration":name,"arm":arm},"field":"aggregate_cov_wilson95",
                                "saved":saved["cov_wilson95"],"independent":expected["cov_wilson95"],"absolute_difference":delta})
            if int(saved["n"])!=n:
                failures.append({"context":{"configuration":name,"arm":arm},"field":"aggregate_n"})
            if arm!="oracle":
                err_by_rep={int(row["replicate"]):float(row["posterior_mean"])-H0_TRUE for row in grp}
                oracle_by_rep={int(row["replicate"]):float(row["posterior_mean"])-H0_TRUE for row in byarm["oracle"]}
                d=np.array([abs(err_by_rep[i])-abs(oracle_by_rep[i]) for i in range(n)])
                contrast=float(d.mean())
                boot=np.empty(4000)
                for start in range(0,4000,500):
                    stop=min(start+500,4000)
                    ids=bootstrap_rng.integers(0,n,size=(stop-start,n))
                    boot[start:stop]=d[ids].mean(axis=1)
                ci=[float(np.quantile(boot,0.025)),float(np.quantile(boot,0.975))]
                delta=abs(float(saved["paired_abs_err_minus_oracle"])-contrast)
                max_aggregate_diff=max(max_aggregate_diff,delta); aggregate_checks+=1
                if delta>1e-9:
                    failures.append({"context":{"configuration":name,"arm":arm},"field":"aggregate_paired_abs_err",
                                    "saved":saved["paired_abs_err_minus_oracle"],"independent":contrast,"absolute_difference":delta})
                delta=max(abs(float(a)-float(b)) for a,b in zip(saved["paired_ci95"],ci))
                max_aggregate_diff=max(max_aggregate_diff,delta); aggregate_checks+=1
                if delta>1e-9:
                    failures.append({"context":{"configuration":name,"arm":arm},"field":"aggregate_paired_ci95",
                                    "saved":saved["paired_ci95"],"independent":ci,"absolute_difference":delta})

    result={"status":"PASS" if not failures else "FAIL","selected_configurations":selected,
        "configurations_audited":len(selected),"campaign_replicates_per_configuration":nrep,
        "replicates_independently_reconstructed_per_configuration":n_audit,
        "raw_rows_independently_reconstructed":checked,"expected_raw_rows":len(selected)*n_audit*len(ARMS),
        "manifest_files_checked":len(manifest["files"]),"saved_aggregate_fields_checked":aggregate_checks,
        "maximum_aggregate_absolute_difference":max_aggregate_diff,"maximum_row_reconstruction_differences":maximum,
        "coverage_flag_borderline_disagreement_count":len(borderline),"coverage_flag_borderline_disagreements":borderline[:100],
        "independent_grid":{"z_points":audit_z_points,"h0_points":audit_h0_points},"tolerances":TOL,
        "per_configuration":perconfig,"failure_count":len(failures),"failures":failures[:100],
        "claim_scope":"independent_replay_and_posterior_reconstruction_for_declared_synthetic_model_only"}
    return result


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("campaign_dir",type=Path)
    p.add_argument("--audit-z-points",type=int,default=401)
    p.add_argument("--audit-h0-points",type=int,default=81)
    p.add_argument("--only-configs",nargs="+",default=None)
    p.add_argument("--replicate-limit",type=int,default=None)
    p.add_argument("--report",type=Path,default=None)
    a=p.parse_args()
    result=audit(a.campaign_dir,audit_z_points=a.audit_z_points,audit_h0_points=a.audit_h0_points,
                 only_configs=a.only_configs,replicate_limit=a.replicate_limit)
    text=json.dumps(result,indent=2,sort_keys=True)
    print(text,flush=True)
    if a.report:
        a.report.parent.mkdir(parents=True,exist_ok=True)
        a.report.write_text(text+"\n",encoding="utf-8")
    raise SystemExit(0 if result["status"]=="PASS" else 1)


if __name__=="__main__":
    main()
