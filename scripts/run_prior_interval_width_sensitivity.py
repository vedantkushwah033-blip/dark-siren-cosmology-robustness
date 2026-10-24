"""Paired posterior interval-width sensitivity audit.

Reuses the same synthetic observations and posterior summaries as the prior-support
audit. Reports whether changes in coverage accompany changes in interval width.
"""
from __future__ import annotations
import argparse, csv, importlib.util
from pathlib import Path
import numpy as np

BASE = Path(__file__).resolve().parent / "run_prior_support_sensitivity.py"
SPEC = importlib.util.spec_from_file_location("prior_support_source", BASE)
if SPEC is None or SPEC.loader is None:
    raise RuntimeError("Unable to load prior-support audit.")
source = importlib.util.module_from_spec(SPEC)
SPEC.loader.exec_module(source)

def summarize_widths(rows):
    baseline=(50.0,90.0)
    lookup={(r["true_H0"],r["replicate"],r["prior_min"],r["prior_max"]):r for r in rows}
    summary=[]; paired=[]
    for truth in source.TRUTHS:
        for prior_min,prior_max in source.PRIORS:
            cell=[r for r in rows if r["true_H0"]==truth and r["prior_min"]==prior_min and r["prior_max"]==prior_max]
            for method in source.METHODS:
                item={"true_H0":truth,"prior_min":prior_min,"prior_max":prior_max,"method":method,"n":len(cell)}
                for level in ("68","95"):
                    widths=np.array([r[f"{method}_upper_{level}"]-r[f"{method}_lower_{level}"] for r in cell],dtype=float)
                    item[f"mean_width_{level}"]=float(widths.mean())
                    item[f"median_width_{level}"]=float(np.median(widths))
                    item[f"p90_width_{level}"]=float(np.quantile(widths,.90))
                summary.append(item)
            if (prior_min,prior_max)==baseline: continue
            for row in cell:
                base=lookup[(truth,row["replicate"],*baseline)]
                for method in source.METHODS:
                    rec={"true_H0":truth,"replicate":row["replicate"],"method":method,
                         "baseline_prior_min":baseline[0],"baseline_prior_max":baseline[1],
                         "comparison_prior_min":prior_min,"comparison_prior_max":prior_max}
                    for level in ("68","95"):
                        bw=base[f"{method}_upper_{level}"]-base[f"{method}_lower_{level}"]
                        cw=row[f"{method}_upper_{level}"]-row[f"{method}_lower_{level}"]
                        rec[f"baseline_width_{level}"]=bw
                        rec[f"comparison_width_{level}"]=cw
                        rec[f"width_difference_{level}"]=cw-bw
                        rec[f"absolute_width_difference_{level}"]=abs(cw-bw)
                    paired.append(rec)
    return summary,paired

def write_csv(path,rows):
    if not rows: raise RuntimeError(f"Refusing to write empty artifact: {path.name}")
    fields=list(dict.fromkeys(k for row in rows for k in row))
    with path.open("w",newline="",encoding="utf-8") as f:
        writer=csv.DictWriter(f,fieldnames=fields);writer.writeheader();writer.writerows(rows)
    with path.open(newline="",encoding="utf-8") as f:
        reader=csv.DictReader(f); loaded=list(reader)
        if reader.fieldnames!=fields or len(loaded)!=len(rows):
            raise RuntimeError(f"CSV validation failed: {path.name}")
    return len(rows)

def main():
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates",type=int,default=300)
    p.add_argument("--seed",type=int,default=20261024)
    p.add_argument("--grid-size",type=int,default=4001)
    p.add_argument("--output-dir",default="results/prior-interval-width")
    a=p.parse_args()
    rows,_=source.run_audit(a.replicates,a.seed,a.grid_size)
    summary,paired=summarize_widths(rows)
    out=Path(a.output_dir);out.mkdir(parents=True,exist_ok=True)
    n1=write_csv(out/"interval_width_summary.csv",summary)
    n2=write_csv(out/"paired_interval_width_differences.csv",paired)
    expected_summary=len(source.TRUTHS)*len(source.PRIORS)*len(source.METHODS)
    expected_paired=len(source.TRUTHS)*(len(source.PRIORS)-1)*a.replicates*len(source.METHODS)
    if n1!=expected_summary or n2!=expected_paired:
        raise RuntimeError(f"Unexpected artifact sizes: summary={n1}, paired={n2}")
    print(f"PASS: {n1} summary rows and {n2} paired interval-width rows")
    for row in summary: print(row)
if __name__=="__main__": main()
