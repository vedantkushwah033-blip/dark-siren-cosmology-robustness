"""Joint scorecard for coverage, interval width, bias, and RMSE.

Runs the prior-support synthetic experiment once and joins all performance
metrics for each truth × prior × method. This is a diagnostic scorecard, not
an automatic prior-selection rule.
"""
from __future__ import annotations
import argparse, csv, importlib.util
from pathlib import Path

def load_module(name, path):
    spec=importlib.util.spec_from_file_location(name,path)
    if spec is None or spec.loader is None: raise RuntimeError(f"Cannot load {path}")
    module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)
    return module

HERE=Path(__file__).resolve().parent
prior=load_module("prior_support_source",HERE/"run_prior_support_sensitivity.py")
width=load_module("prior_width_source",HERE/"run_prior_interval_width_sensitivity.py")

def build_scorecard(summaries,width_summaries):
    lookup={(r["true_H0"],r["prior_min"],r["prior_max"],r["method"]):r for r in width_summaries}
    output=[]
    for row in summaries:
        key=(row["true_H0"],row["prior_min"],row["prior_max"],row["method"])
        w=lookup[key]
        out=dict(row)
        out["abs_coverage_68_error_from_nominal"]=abs(row["coverage_68"]-0.68)
        out["abs_coverage_95_error_from_nominal"]=abs(row["coverage_95"]-0.95)
        for field in ("mean_width_68","median_width_68","p90_width_68",
                      "mean_width_95","median_width_95","p90_width_95"):
            out[field]=w[field]
        output.append(out)
    expected=len(prior.TRUTHS)*len(prior.PRIORS)*len(prior.METHODS)
    if len(output)!=expected: raise RuntimeError(f"Expected {expected} scorecard rows, got {len(output)}")
    if len({(r["true_H0"],r["prior_min"],r["prior_max"],r["method"]) for r in output})!=expected:
        raise RuntimeError("Duplicate or missing scorecard keys")
    return output

def write_csv(path,rows):
    if not rows: raise RuntimeError("Refusing to write empty scorecard")
    fields=list(dict.fromkeys(k for r in rows for k in r))
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f); reread=list(reader)
        if reader.fieldnames!=fields or len(reread)!=len(rows):
            raise RuntimeError("Scorecard CSV round-trip validation failed")
        for row in reread:
            for name in ("bias","rmse","coverage_68","coverage_95","mean_width_68","mean_width_95"):
                if not row[name]: raise RuntimeError(f"Missing {name}")
    return fields

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261024)
    p.add_argument("--grid-size",type=int,default=4001)
    p.add_argument("--output",default="results/joint-prior-scorecard.csv")
    a=p.parse_args()
    _,summaries=prior.run_audit(a.replicates,a.seed,a.grid_size)
    rows,width_summaries=width.summarize_widths(prior.run_audit(a.replicates,a.seed,a.grid_size)[0])
    score=build_scorecard(summaries,width_summaries)
    out=Path(a.output);out.parent.mkdir(parents=True,exist_ok=True)
    fields=write_csv(out,score)
    print(f"PASS: {len(score)} rows; {len(fields)} columns; output={out}")
    for r in score:
        print({k:r[k] for k in ("true_H0","prior_min","prior_max","method","bias","rmse","coverage_68","abs_coverage_68_error_from_nominal","mean_width_68","coverage_95","mean_width_95")})
if __name__=="__main__": main()
