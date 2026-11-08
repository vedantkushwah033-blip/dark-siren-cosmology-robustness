"""Predeclared replicate-level analysis for the fixed-catalogue confirmatory design."""
from __future__ import annotations
import argparse,json
from pathlib import Path
import numpy as np
import pandas as pd
from scipy.stats import t as student_t

PRIMARY_GAMMA_NEG=-1.2
PRIMARY_GAMMA_POS=1.2
PRIMARY_N_VALIDATION=300
PRIMARY_DEGREE=1
PRIMARY_BANDWIDTH=0.15

def _paired_difference(group,arm):
 wide=group.pivot(index="replicate",columns="arm",values="iae")
 if "catalogue_only" not in wide or arm not in wide:return pd.Series(dtype=float)
 return (wide[arm]-wide["catalogue_only"]).dropna()

def _holm_adjust(pvalues):
 """Holm adjustment over the full predeclared family.

 Unavailable tests are assigned p=1 *for adjustment only*; their raw p-value
 remains null and their contrast remains explicitly non-evaluable. This keeps
 missing tests from silently shrinking the multiplicity family.
 """
 values=np.asarray(pvalues,float)
 if values.ndim!=1:raise ValueError("pvalues must be one-dimensional")
 if np.any(np.isfinite(values)&((values<0)|(values>1))):
  raise ValueError("finite p-values must lie in [0,1]")
 safe=np.where(np.isfinite(values),values,1.0)
 out=np.empty(safe.shape,float)
 order=np.argsort(safe);running=0.0;m=len(order)
 for rank,idx in enumerate(order):
  adjusted=min(1.0,(m-rank)*safe[idx])
  running=max(running,adjusted);out[idx]=running
 return out

def analyze_fixed_catalogue_attempts(attempts:pd.DataFrame,*,planned_replicates:int,
        bootstrap_resamples:int=10000,bootstrap_seed:int=20261110)->dict:
 required={"replicate","gamma","target_retention","n_validation","bandwidth_multiplier",
           "degree","arm","status","iae"}
 if not required.issubset(attempts.columns):
  raise ValueError(f"Missing attempt columns: {sorted(required-set(attempts.columns))}")
 if planned_replicates<1 or bootstrap_resamples<1000:
  raise ValueError("planned_replicates>=1 and bootstrap_resamples>=1000 required")
 primary=attempts[(attempts.n_validation==PRIMARY_N_VALIDATION)&
  (np.isclose(attempts.bandwidth_multiplier,PRIMARY_BANDWIDTH))&
  (attempts.degree.isin([0,PRIMARY_DEGREE]))&
  (attempts.arm.isin(["catalogue_only",f"estimated_degree_{PRIMARY_DEGREE}_ipw"]))&
  (attempts.gamma.isin([PRIMARY_GAMMA_NEG,PRIMARY_GAMMA_POS]))]
 deltas={}
 for gamma in (PRIMARY_GAMMA_NEG,PRIMARY_GAMMA_POS):
  g=primary[primary.gamma==gamma]
  d=_paired_difference(g,f"estimated_degree_{PRIMARY_DEGREE}_ipw")
  deltas[gamma]=d
 if all(len(deltas[g]) for g in deltas):
  paired=pd.concat([deltas[PRIMARY_GAMMA_NEG].rename("negative"),
                    deltas[PRIMARY_GAMMA_POS].rename("positive")],axis=1).dropna()
 else:
  paired=pd.DataFrame(columns=["negative","positive"])
 interaction=(paired.positive-paired.negative).to_numpy(float)
 primary_result={"estimand":"Delta(+1.2)-Delta(-1.2), degree-1, n_validation=300, Scott multiplier=0.15",
  "n_replicates_planned":int(planned_replicates),"n_complete_pairs":int(len(interaction)),
  "n_incomplete_replicates":int(planned_replicates-len(interaction)),
  "mean_interaction":float(interaction.mean()) if len(interaction) else None,
  "monte_carlo_se":float(interaction.std(ddof=1)/np.sqrt(len(interaction))) if len(interaction)>1 else None,
  "bootstrap_95_interval":None,"paired_t_95_interval":None,
  "failure_aware_cell_counts":{}}
 for gamma in (PRIMARY_GAMMA_NEG,PRIMARY_GAMMA_POS):
  for arm in ("catalogue_only","estimated_degree_1_ipw"):
   cell=attempts[(attempts.gamma==gamma)&(attempts.n_validation==PRIMARY_N_VALIDATION)&
    (np.isclose(attempts.bandwidth_multiplier,PRIMARY_BANDWIDTH))&(attempts.arm==arm)]
   primary_result["failure_aware_cell_counts"][f"gamma={gamma}:{arm}"]={
    "n_planned":int(len(cell)),"n_success":int((cell.status=="SUCCESS").sum()),
    "n_failed":int((cell.status=="FAILED").sum()),
    "n_non_evaluable":int((cell.status=="NON_EVALUABLE").sum()),
    "n_not_run":int((cell.status=="NOT_RUN").sum())}
 if len(interaction)>=2:
  rng=np.random.default_rng(bootstrap_seed)
  draws=rng.integers(0,len(interaction),size=(bootstrap_resamples,len(interaction)))
  boot=interaction[draws].mean(axis=1)
  primary_result["bootstrap_95_interval"]=[float(np.quantile(boot,.025)),float(np.quantile(boot,.975))]
  se=float(interaction.std(ddof=1)/np.sqrt(len(interaction)))
  crit=float(student_t.ppf(.975,len(interaction)-1));mean=float(interaction.mean())
  primary_result["paired_t_95_interval"]=[mean-crit*se,mean+crit*se]
 secondary=[]
 for gamma in sorted(attempts.gamma.unique()):
  for nval in (100,300,1000):
   for degree in (1,2):
    arm=f"estimated_degree_{degree}_ipw"
    g=attempts[(attempts.gamma==gamma)&(attempts.n_validation==nval)&
     (np.isclose(attempts.bandwidth_multiplier,PRIMARY_BANDWIDTH))&
     (attempts.arm.isin(["catalogue_only",arm]))]
    diff=_paired_difference(g,arm).to_numpy(float)
    n=len(diff);mean=float(diff.mean()) if n else None
    se=float(diff.std(ddof=1)/np.sqrt(n)) if n>1 else None
    if n>1:
     crit=float(student_t.ppf(.975,n-1));lo=mean-crit*se;hi=mean+crit*se
     pval=float(2*student_t.sf(abs(mean/se),n-1)) if se>0 else (0.0 if mean else 1.0)
    else:lo=hi=pval=None
    secondary.append({"gamma":float(gamma),"n_validation":int(nval),"degree":int(degree),
     "n_pairs":int(n),"n_replicates_planned":int(planned_replicates),
     "mean_iae_difference":mean,"monte_carlo_se":se,
     "paired_t_95_interval":[lo,hi] if lo is not None else None,"p_value":pval})
 adjusted=_holm_adjust([x["p_value"] if x["p_value"] is not None else np.nan for x in secondary])
 for item,p in zip(secondary,adjusted):
  item["holm_adjusted_p_value"]=float(p) if np.isfinite(p) else None
 return {"analysis":"predeclared_fixed_analysis_catalogue_confirmatory_analysis",
  "primary_interaction":primary_result,"secondary_contrasts":secondary,
  "secondary_multiplicity":"Holm adjustment across the full predeclared 18-test degree x curvature x validation-size family; unavailable tests retain null raw p-values and are assigned p=1 only for multiplicity adjustment, so missing contrasts cannot shrink the family. Intervals are pointwise paired-t intervals.",
  "n_secondary_tests_predeclared":18,
  "n_secondary_tests_evaluable":int(sum(x["p_value"] is not None for x in secondary)),
  "n_secondary_tests_non_evaluable":int(sum(x["p_value"] is None for x in secondary)),
  "limitations":["All estimates concern synthetic density recovery only.",
   "The primary interaction is pilot-informed and must not be described as independent hypothesis generation.",
   "Bootstrap resamples replicate IDs, preserving the within-replicate curvature pairing.",
   "Failures are reported separately; no numeric IAE penalty is imputed.",
   "The Holm family always contains all 18 predeclared contrasts; unavailable tests have null raw p-values and use p=1 only as a conservative adjustment placeholder."]}

def main():
 p=argparse.ArgumentParser(description=__doc__);p.add_argument("--attempts-csv",type=Path,required=True)
 p.add_argument("--planned-replicates",type=int,default=100);p.add_argument("--output",type=Path,required=True)
 a=p.parse_args();frame=pd.read_csv(a.attempts_csv)
 result=analyze_fixed_catalogue_attempts(frame,planned_replicates=a.planned_replicates)
 a.output.parent.mkdir(parents=True,exist_ok=True)
 a.output.write_text(json.dumps(result,indent=2,sort_keys=True)+"\n",encoding="utf-8")
 print(json.dumps(result,indent=2,sort_keys=True))
if __name__=="__main__":main()
