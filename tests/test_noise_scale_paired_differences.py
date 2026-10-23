import csv
import subprocess
import sys
from pathlib import Path

def test_paired_bootstrap_writes_expected_metrics(tmp_path):
    root=Path(__file__).resolve().parents[1]
    src=tmp_path/"per_run.csv"
    rows=[]
    for h0 in (70.0,):
        for rep in range(20):
            for mode,err,cov,width in [
                ("fixed_reported_sigma",2.0,1,10.0),
                ("sigma_depends_on_h0",1.0,1,8.0),
            ]:
                rows.append({"true_H0":h0,"replicate":rep,"mode":mode,
                    "median_error":err,"covers_68":cov,"covers_95":cov,
                    "width_68":width,"width_95":width*2})
    with src.open("w",newline="") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0]));w.writeheader();w.writerows(rows)
    out=tmp_path/"paired.csv"
    subprocess.run([sys.executable,str(root/"scripts/analyze_noise_scale_paired_differences.py"),
        "--input",str(src),"--output",str(out),"--bootstrap","1000"],check=True,capture_output=True,text=True)
    got=list(csv.DictReader(out.open()))
    assert len(got)==6
    rmse=next(r for r in got if r["metric"]=="rmse_median")
    assert float(rmse["delta_dependent_minus_fixed"]) < 0
    assert float(rmse["bootstrap_95_high"]) < 0
