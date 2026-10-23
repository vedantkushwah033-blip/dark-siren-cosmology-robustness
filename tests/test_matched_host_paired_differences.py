import csv, subprocess, sys
from pathlib import Path

def test_paired_comparison_detects_controlled_difference(tmp_path):
 root=Path(__file__).resolve().parents[1]
 src=tmp_path/"per_run.csv"; rows=[]
 for rep in range(30):
  for scenario,err,c68,w68 in [("known_host",2.0,1,10.0),("unknown_host",1.0,1,8.0)]:
   rows.append({"true_H0":70.0,"replicate":rep,"scenario":scenario,
    "median_error":err,"covers_68":c68,"covers_95":c68,"width_68":w68,"width_95":w68*2})
 with src.open("w",newline="") as f:
  w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
 out=tmp_path/"paired.csv"
 subprocess.run([sys.executable,str(root/"scripts/analyze_matched_host_paired_differences.py"),
  "--input",str(src),"--output",str(out),"--bootstrap","1000"],check=True,capture_output=True,text=True)
 with out.open(newline="") as stream:\n  got=list(csv.DictReader(stream))
 assert len(got)==7
 rmse=next(r for r in got if r["metric"]=="rmse_median")
 assert float(rmse["delta_unknown_minus_known"])<0
 assert float(rmse["bootstrap_95_high"])<0
