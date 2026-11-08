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
 assert primary["mean_interaction"]==0.32
 assert len(result["secondary_contrasts"])==18
 assert all(x["n_pairs"]==5 for x in result["secondary_contrasts"])
 adjusted=[x["holm_adjusted_p_value"] for x in result["secondary_contrasts"]]
 assert all(x is not None for x in adjusted)
 assert all(0<=x<=1 for x in adjusted)
