"""Post-hoc uncertainty summaries for the paired quadratic-logit H0 pilot.

Consumes the aggregate replicate_ledger.csv. Coverage intervals are Wilson
intervals; continuous-parameter uncertainty uses t intervals and paired
bootstrap intervals. All outputs remain descriptive at pilot scale.
"""
from __future__ import annotations
import argparse, json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t


def mean_t_interval(values: np.ndarray) -> tuple[float, float, float]:
    x=np.asarray(values,dtype=float); x=x[np.isfinite(x)]
    if len(x)==0: return (float("nan"),float("nan"),float("nan"))
    mean=float(x.mean())
    if len(x)<2: return (mean,float("nan"),float("nan"))
    half=float(t.ppf(.975,len(x)-1)*x.std(ddof=1)/np.sqrt(len(x)))
    return mean,mean-half,mean+half


def wilson(k: int,n: int,z: float=1.959963984540054) -> tuple[float,float]:
    if n==0: return (float("nan"),float("nan"))
    p=k/n; den=1+z*z/n
    centre=(p+z*z/(2*n))/den
    half=z*np.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return float(max(0,centre-half)),float(min(1,centre+half))


def bootstrap_mean_interval(values: np.ndarray, rng: np.random.Generator,
                            draws: int=10000) -> tuple[float,float]:
    x=np.asarray(values,dtype=float); x=x[np.isfinite(x)]
    if len(x)<2: return (float("nan"),float("nan"))
    idx=rng.integers(0,len(x),size=(draws,len(x)))
    means=x[idx].mean(axis=1)
    return tuple(map(float,np.quantile(means,[.025,.975])))



def _assert_pairing_keys(paired: pd.DataFrame, suffixes: tuple[str, str]) -> None:
    """Reject replicate-ID matches that do not share the generating design."""
    left, right = suffixes
    for key in ("analysis_seed", "validation_seed", "h0_true"):
        a = paired[f"{key}_{left}"].to_numpy()
        b = paired[f"{key}_{right}"].to_numpy()
        if not np.array_equal(a, b):
            raise ValueError(f"paired replicate mismatch in {key}")


def analyze(root: Path, output: Path, seed: int=20261119) -> dict:
    root=root.resolve(); output.mkdir(parents=True,exist_ok=True)
    ledger=pd.read_csv(root/"replicate_ledger.csv")
    required={"replicate","gamma","n_validation","status","h0_true","h0_median",
              "h0_q16","h0_q84","covered68","covered90","covered95",
              "min_catalogue_kish_ess","max_normalized_weight","analysis_seed","validation_seed"}
    missing=required-set(ledger.columns)
    if missing: raise ValueError(f"ledger missing columns: {sorted(missing)}")
    if ledger.duplicated(["replicate","gamma","n_validation"]).any():
        raise ValueError("duplicate replicate/gamma/validation-size cells")
    rng=np.random.default_rng(seed)
    rows=[]
    for (gamma,nval),g in ledger.groupby(["gamma","n_validation"],sort=True):
        ok=g[g.status=="SUCCESS"].copy(); n=len(ok)
        bias=(ok.h0_median-ok.h0_true).to_numpy(float)
        mean,lo,hi=mean_t_interval(bias)
        rmse=float(np.sqrt(np.mean(bias*bias))) if n else float("nan")
        rlo,rhi=bootstrap_mean_interval(bias*bias,rng)
        if n >= 2 and np.isfinite([rlo, rhi]).all():
            rlo,rhi=(float(np.sqrt(max(0,rlo))),float(np.sqrt(max(0,rhi))))
        else:
            # Undefined bootstrap bounds must remain undefined for n < 2.
            rlo=rhi=float("nan")
        planned=len(g)
        unknown=planned-n
        row={"gamma":gamma,"n_validation":nval,"planned":planned,"success":n,
             "non_success":unknown,
             "non_evaluable":int((g.status=="NON_EVALUABLE").sum()),
             "failed":int((g.status=="FAILED").sum()),
             "not_run":int((g.status=="NOT_RUN").sum()),
             "mean_median_bias":mean,"bias_ci95_low":lo,"bias_ci95_high":hi,
             "median_bias":float(np.median(bias)) if n else float("nan"),
             "rmse_posterior_median":rmse,"rmse_bootstrap95_low":rlo,
             "rmse_bootstrap95_high":rhi,
             "mean_width68":float((ok.h0_q84-ok.h0_q16).mean()) if n else float("nan"),
             "median_min_kish_ess":float(ok.min_catalogue_kish_ess.median()) if n else float("nan"),
             "median_max_normalized_weight":float(ok.max_normalized_weight.median()) if n else float("nan")}
        for level in (68,90,95):
            col=f"covered{level}"; k=int(ok[col].sum()) if n else 0
            low,high=wilson(k,n)
            # Wilson intervals describe SUCCESS cells; all-attempt bounds
            # additionally treat each non-success as uncovered/covered.
            row.update({f"coverage{level}":k/n if n else float("nan"),
                        f"coverage{level}_wilson95_low":low,
                        f"coverage{level}_wilson95_high":high,
                        f"coverage{level}_all_attempt_lower":k/planned if planned else float("nan"),
                        f"coverage{level}_all_attempt_upper":(k+unknown)/planned if planned else float("nan")})
        rows.append(row)
    cell=pd.DataFrame(rows); cell.to_csv(output/"cell_uncertainty_summary.csv",index=False)

    paired=[]
    success=ledger[ledger.status=="SUCCESS"]
    gammas=sorted(ledger.gamma.unique())
    if 0.0 in gammas:
        for gamma in gammas:
            if gamma==0: continue
            for nval in sorted(ledger.n_validation.unique()):
                left=success[(success.gamma==gamma)&(success.n_validation==nval)]
                right=success[(success.gamma==0)&(success.n_validation==nval)]
                m=left.merge(right,on="replicate",suffixes=("_gamma","_zero"),validate="one_to_one")
                if m.empty: continue
                _assert_pairing_keys(m, ("gamma", "zero"))
                for metric,values in {
                    "signed_median_bias":m.h0_median_gamma-m.h0_median_zero,
                    "absolute_median_error":abs(m.h0_median_gamma-m.h0_true_gamma)-abs(m.h0_median_zero-m.h0_true_zero),
                    "width68":(m.h0_q84_gamma-m.h0_q16_gamma)-(m.h0_q84_zero-m.h0_q16_zero),
                }.items():
                    x=np.asarray(values,dtype=float); avg,lo,hi=mean_t_interval(x)
                    blo,bhi=bootstrap_mean_interval(x,rng)
                    paired.append({"contrast":f"gamma_{gamma:g}_minus_gamma_0",
                        "gamma":gamma,"n_validation":nval,"metric":metric,
                        "paired_replicates":len(x),"mean_difference":avg,
                        "paired_t95_low":lo,"paired_t95_high":hi,
                        "paired_bootstrap95_low":blo,"paired_bootstrap95_high":bhi})
    pd.DataFrame(paired).to_csv(output/"paired_gamma_contrasts.csv",index=False)

    size_effects=[]
    for gamma,g in success.groupby("gamma"):
        sizes=sorted(g.n_validation.unique())
        if len(sizes)<2: continue
        low=g[g.n_validation==sizes[0]]
        high=g[g.n_validation==sizes[-1]]
        m=high.merge(low,on="replicate",suffixes=("_high","_low"),validate="one_to_one")
        if m.empty: continue
        _assert_pairing_keys(m, ("high", "low"))
        for metric,values in {
            "signed_median_bias":m.h0_median_high-m.h0_median_low,
            "absolute_median_error":abs(m.h0_median_high-m.h0_true_high)-abs(m.h0_median_low-m.h0_true_low),
            "width68":(m.h0_q84_high-m.h0_q16_high)-(m.h0_q84_low-m.h0_q16_low),
        }.items():
            x=np.asarray(values,dtype=float); avg,lo,hi=mean_t_interval(x)
            blo,bhi=bootstrap_mean_interval(x,rng)
            size_effects.append({"gamma":gamma,"comparison":f"{sizes[-1]}_minus_{sizes[0]}",
                "metric":metric,"paired_replicates":len(x),"mean_difference":avg,
                "paired_t95_low":lo,"paired_t95_high":hi,
                "paired_bootstrap95_low":blo,"paired_bootstrap95_high":bhi})
    pd.DataFrame(size_effects).to_csv(output/"paired_validation_size_contrasts.csv",index=False)
    report={"status":"PASS","ledger_rows":len(ledger),
       "status_counts":{str(k):int(v) for k,v in ledger.status.value_counts().items()},
       "cell_summaries":len(cell),"paired_gamma_contrasts":len(paired),
       "paired_validation_size_contrasts":len(size_effects),
       "uncertainty_note":"Coverage uses Wilson binomial intervals. Bias uses t intervals; RMSE and paired contrasts use percentile bootstrap intervals. With 10 replicates these are descriptive and unstable, not confirmatory.",
       "coverage_note":"Wilson intervals describe SUCCESS cells. All-attempt lower/upper bounds count every non-success as uncovered/covered, respectively; failures, non-evaluable and not-run cells remain in the planned denominator.",
       "sbc_note":"This is fixed-truth repeated-sampling coverage, not joint simulation-based calibration."}
    (output/"posthoc_report.json").write_text(json.dumps(report,indent=2)+"\n")
    return report


def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("artifact_dir",type=Path); p.add_argument("--output-dir",type=Path,default=Path("h0_pilot_posthoc"))
    p.add_argument("--seed",type=int,default=20261119)
    a=p.parse_args(); print(json.dumps(analyze(a.artifact_dir,a.output_dir,a.seed),indent=2))


if __name__=="__main__": main()
