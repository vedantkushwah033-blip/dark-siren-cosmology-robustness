"""Paired bootstrap comparison: unknown-host minus known-host at fixed H0."""
from __future__ import annotations
import argparse, csv
from pathlib import Path
import numpy as np

METRICS={
 "median_error":lambda r:float(r["median_error"]),
 "absolute_error":lambda r:abs(float(r["median_error"])),
 "covers_68":lambda r:float(r["covers_68"]),
 "covers_95":lambda r:float(r["covers_95"]),
 "width_68":lambda r:float(r["width_68"]),
 "width_95":lambda r:float(r["width_95"]),
}
def main():
 p=argparse.ArgumentParser()
 p.add_argument("--input",default="results/matched-fixed-h0-coverage/per_run.csv")
 p.add_argument("--output",default="results/matched-fixed-h0-coverage/paired_comparison.csv")
 p.add_argument("--bootstrap",type=int,default=10000)
 p.add_argument("--seed",type=int,default=20261023)
 a=p.parse_args()
 if a.bootstrap<1000:p.error("--bootstrap must be at least 1000")
 rows=list(csv.DictReader(open(a.input,newline="",encoding="utf-8")))
 by={}
 for r in rows:by.setdefault((float(r["true_H0"]),int(r["replicate"])),{})[r["scenario"]]=r
 rng=np.random.default_rng(a.seed);out=[]
 for h0 in sorted({k[0] for k in by}):
  pairs=[v for (h,rep),v in sorted(by.items()) if h==h0]
  if not pairs or any("known_host" not in x or "unknown_host" not in x for x in pairs):
   raise ValueError(f"unpaired scenario rows at H0={h0}")
  for metric,fn in METRICS.items():
   delta=np.array([fn(x["unknown_host"])-fn(x["known_host"]) for x in pairs])
   draws=np.empty(a.bootstrap)
   for start in range(0,a.bootstrap,250):
    stop=min(start+250,a.bootstrap)
    idx=rng.integers(0,len(delta),size=(stop-start,len(delta)))
    draws[start:stop]=delta[idx].mean(axis=1)
   lo,hi=np.quantile(draws,[.025,.975])
   out.append({"true_H0":h0,"metric":metric,"n_paired":len(pairs),
    "delta_unknown_minus_known":float(delta.mean()),"bootstrap_95_low":float(lo),
    "bootstrap_95_high":float(hi),"bootstrap_draws":a.bootstrap,"seed":a.seed})
  e0=np.array([float(x["known_host"]["median_error"]) for x in pairs])
  e1=np.array([float(x["unknown_host"]["median_error"]) for x in pairs])
  d=np.sqrt(np.mean(e1**2))-np.sqrt(np.mean(e0**2))
  draws=np.empty(a.bootstrap)
  for start in range(0,a.bootstrap,250):
   stop=min(start+250,a.bootstrap)
   idx=rng.integers(0,len(e0),size=(stop-start,len(e0)))
   draws[start:stop]=np.sqrt(np.mean(e1[idx]**2,axis=1))-np.sqrt(np.mean(e0[idx]**2,axis=1))
  lo,hi=np.quantile(draws,[.025,.975])
  out.append({"true_H0":h0,"metric":"rmse_median","n_paired":len(pairs),
    "delta_unknown_minus_known":float(d),"bootstrap_95_low":float(lo),
    "bootstrap_95_high":float(hi),"bootstrap_draws":a.bootstrap,"seed":a.seed})
 path=Path(a.output);path.parent.mkdir(parents=True,exist_ok=True)
 with path.open("w",newline="",encoding="utf-8") as f:
  w=csv.DictWriter(f,fieldnames=list(out[0]));w.writeheader();w.writerows(out)
 print(f"Wrote {len(out)} paired comparisons to {path}")
 for r in out:
  if r["metric"] in ("rmse_median","covers_68","median_error"):print(r)
if __name__=="__main__":main()
