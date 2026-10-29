import csv, importlib.util, json
from pathlib import Path
from dark_siren_robustness.host_density_benchmark import run_benchmark

spec=importlib.util.spec_from_file_location("independent_host_density_audit",Path(__file__).resolve().parents[1]/"scripts"/"audit_host_density_independent.py")
module=importlib.util.module_from_spec(spec); spec.loader.exec_module(module)


def test_reconstruction_passes_and_detects_self_consistent_tampering(tmp_path):
    frame,summary=run_benchmark(n_replicates=1,n_parent=1500,n_detected=15,bandwidth_factor=.2)
    frame.to_csv(tmp_path/"replicate_metrics.csv",index=False)
    (tmp_path/"summary.json").write_text(json.dumps(summary),encoding="utf-8")
    result=module.audit(tmp_path,n_reps=1,z_points=801,h0_points=321)
    assert result["status"]=="PASS"
    rows=list(csv.DictReader((tmp_path/"replicate_metrics.csv").open(newline="",encoding="utf-8")))
    for row in rows:
        if row["arm"]=="known_completeness":
            row["posterior_mean"]=str(float(row["posterior_mean"])+.5)
    with (tmp_path/"replicate_metrics.csv").open("w",newline="",encoding="utf-8") as handle:
        writer=csv.DictWriter(handle,fieldnames=list(rows[0]));writer.writeheader();writer.writerows(rows)
    result=module.audit(tmp_path,n_reps=1,z_points=801,h0_points=321)
    assert result["status"]=="FAIL" and any(x["field"]=="posterior_mean" for x in result["failures"])
