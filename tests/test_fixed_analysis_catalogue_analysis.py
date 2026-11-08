import pandas as pd
from dark_siren_robustness.fixed_analysis_catalogue_analysis import analyze_fixed_catalogue_attempts

def test_primary_curvature_interaction_and_holm_secondary_family():
 rows=[]
 for rep in range(5):
  for gamma in (-1.2,0.0,1.2):
   for nval in (100,300,1000):
    base=1.0+0.01*rep+0.02*gamma
    for arm,degree,delta in [
     ("catalogue_only",0,0.0),
     ("estimated_degree_1_ipw",1,(-0.1+0.01*rep if gamma==-1.2 and nval==300 else
       0.2+0.02*rep if gamma==1.2 and nval==300 else -0.03)),
     ("estimated_degree_2_ipw",2,0.01),
    ]:
     rows.append({"replicate":rep,"gamma":gamma,"target_retention":0.25,
      "n_validation":nval,"bandwidth_multiplier":0.15,"degree":degree,"arm":arm,
      "status":"SUCCESS","iae":base+delta})
 result=analyze_fixed_catalogue_attempts(pd.DataFrame(rows),planned_replicates=5,
  bootstrap_resamples=1000,bootstrap_seed=3)
 primary=result["primary_interaction"]
 assert primary["n_complete_pairs"]==5
 assert abs(primary["mean_interaction"]-0.32)<1e-12
 assert len(result["secondary_contrasts"])==18
 assert all(x["n_pairs"]==5 for x in result["secondary_contrasts"])
 adjusted=[x["holm_adjusted_p_value"] for x in result["secondary_contrasts"]]
 assert all(x is not None for x in adjusted)
 assert all(0<=x<=1 for x in adjusted)


def test_unavailable_secondary_contrast_does_not_shrink_holm_family():
 rows=[]
 for rep in range(5):
  for gamma in (-1.2,0.0,1.2):
   for nval in (100,300,1000):
    base=1.0+0.01*rep+0.02*gamma
    rows.append({"replicate":rep,"gamma":gamma,"target_retention":0.25,
     "n_validation":nval,"bandwidth_multiplier":0.15,"degree":0,
     "arm":"catalogue_only","status":"SUCCESS","iae":base})
    for degree in (1,2):
     # Deliberately make one predeclared contrast non-evaluable.
     if gamma==-1.2 and nval==100 and degree==1:
      continue
     rows.append({"replicate":rep,"gamma":gamma,"target_retention":0.25,
      "n_validation":nval,"bandwidth_multiplier":0.15,"degree":degree,
      "arm":f"estimated_degree_{degree}_ipw","status":"SUCCESS",
      "iae":base+0.02*degree})
 result=analyze_fixed_catalogue_attempts(pd.DataFrame(rows),planned_replicates=5,
  bootstrap_resamples=1000,bootstrap_seed=9)
 assert len(result["secondary_contrasts"])==18
 unavailable=[x for x in result["secondary_contrasts"]
  if x["gamma"]==-1.2 and x["n_validation"]==100 and x["degree"]==1]
 assert len(unavailable)==1
 assert unavailable[0]["p_value"] is None
 assert unavailable[0]["holm_adjusted_p_value"]==1.0
 assert result["n_secondary_tests_predeclared"]==18
 assert result["n_secondary_tests_evaluable"]==17
 assert result["n_secondary_tests_non_evaluable"]==1
 assert all(x["holm_adjusted_p_value"] is not None
  for x in result["secondary_contrasts"])
