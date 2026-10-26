import csv, hashlib, importlib.util, json
from pathlib import Path
import pytest
from scipy.integrate import cumulative_trapezoid, trapezoid
from scipy.stats import norm

ROOT=Path(__file__).resolve().parents[1]
def load(name,path):
    spec=importlib.util.spec_from_file_location(name,ROOT/"scripts"/path)
    mod=importlib.util.module_from_spec(spec); spec.loader.exec_module(mod); return mod
RUN=load("posterior_runner","run_posterior_calibration.py")
AUDIT=load("posterior_auditor","audit_posterior_calibration.py")

def write_csv(path,rows):
    with path.open("w",newline="",encoding="utf-8") as f:
        w=csv.DictWriter(f,fieldnames=list(rows[0])); w.writeheader(); w.writerows(rows)

def fixture(root):
    raw,summaries=RUN.run(replicates=100,seed=1771)
    write_csv(root/"per_run.csv",raw); write_csv(root/"summary.csv",summaries)
    meta={"audit_status":"PENDING","base_seed":1771,"replicates_per_truth_noise_cell":100,
          "raw_records":len(raw),"summary_rows":len(summaries),
          "sha256":{n:hashlib.sha256((root/n).read_bytes()).hexdigest() for n in ("per_run.csv","summary.csv")}}
    (root/"run_metadata.json").write_text(json.dumps(meta,indent=2)+"\n")


def test_posterior_uses_likelihood_truncated_to_positive_distances():
    # Independent direct calculation for one arm. The generator conditions on
    # positive distances, so an untruncated Gaussian mixture is misspecified.
    import numpy as np
    distance, sigma, truth = 25.0, 40.0, 65.0
    support = (20.0, 220.0)
    h = RUN.GRID[(RUN.GRID >= support[0]) & (RUN.GRID <= support[1])]
    mu1, mu2 = RUN.C * RUN.Z1 / h, RUN.C * RUN.Z2 / h
    likelihood = (
        0.7 * norm.pdf(distance, loc=mu1, scale=sigma) / norm.cdf(mu1 / sigma)
        + 0.3 * norm.pdf(distance, loc=mu2, scale=sigma) / norm.cdf(mu2 / sigma)
    )
    density = likelihood / trapezoid(likelihood, h)
    cdf = np.concatenate(([0.0], cumulative_trapezoid(density, h)))
    cdf /= cdf[-1]
    expected_median = float(np.interp(0.5, cdf, h))
    actual = RUN.posterior_summary(
        distance, sigma, (0.7, 0.3), support, "uniform_H0", truth
    )
    assert actual["median"] == pytest.approx(expected_median, abs=1e-10)

def test_posterior_calibration_design_and_summary_metrics():
    raw,summaries=RUN.run(replicates=100,seed=1771)
    assert len(raw)==4*3*100
    assert len(summaries)==4*3*2*2*2
    assert {r["true_H0"] for r in raw}==set(RUN.TRUTHS)
    assert all(0 <= int(value) <= 1 for r in raw for key, value in r.items() if key.endswith(("_covered68", "_covered95")))
    assert all(0<=r["coverage68"]<=1 and 0<=r["coverage95"]<=1 for r in summaries)

def test_independent_audit_reconstructs_posterior_outputs(tmp_path):
    fixture(tmp_path)
    report=AUDIT.audit(tmp_path)
    assert report["audit"]=="PASS"
    assert report["raw_records"]==1200 and report["summary_rows"]==96

def test_audit_rejects_tampered_summary_after_hash_refresh(tmp_path):
    fixture(tmp_path)
    path=tmp_path/"summary.csv"
    with path.open(newline="",encoding="utf-8") as f: rows=list(csv.DictReader(f))
    rows[0]["rmse_median"]=str(float(rows[0]["rmse_median"])+0.5)
    write_csv(path,rows)
    mp=tmp_path/"run_metadata.json"; meta=json.loads(mp.read_text())
    meta["sha256"]["summary.csv"]=hashlib.sha256(path.read_bytes()).hexdigest()
    mp.write_text(json.dumps(meta,indent=2)+"\n")
    with pytest.raises(ValueError,match="rmse_median mismatch"): AUDIT.audit(tmp_path)


@pytest.mark.parametrize(
    ("field", "value", "message"),
    [
        ("replicates_per_truth_noise_cell", 100.5, "must be an integer >= 100"),
        ("replicates_per_truth_noise_cell", True, "must be an integer >= 100"),
        ("replicates_per_truth_noise_cell", 99, "must be an integer >= 100"),
        ("base_seed", -1, "must be a non-negative integer"),
        ("base_seed", 1.5, "must be a non-negative integer"),
        ("base_seed", False, "must be a non-negative integer"),
    ],
)
def test_audit_rejects_malformed_numeric_metadata(tmp_path, field, value, message):
    fixture(tmp_path)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata[field] = value
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match=message):
        AUDIT.audit(tmp_path)


@pytest.mark.parametrize("replicates", [99, 100.5, True, "100"])
def test_runner_rejects_invalid_replicate_types_and_counts(replicates):
    with pytest.raises(ValueError, match="replicates must be an integer of at least 100"):
        RUN.run(replicates=replicates, seed=1771)


@pytest.mark.parametrize("seed", [-1, 1.5, True, "1771"])
def test_runner_rejects_invalid_seed_types_and_values(seed):
    with pytest.raises(ValueError, match="seed must be a non-negative integer"):
        RUN.run(replicates=100, seed=seed)


@pytest.mark.parametrize(
    ("field", "value"),
    [
        ("raw_records", 1199),
        ("raw_records", 1200.0),
        ("raw_records", True),
        ("summary_rows", 95),
        ("summary_rows", 96.0),
        ("summary_rows", False),
    ],
)
def test_audit_rejects_incorrect_or_noninteger_declared_row_counts(tmp_path, field, value):
    fixture(tmp_path)
    metadata_path = tmp_path / "run_metadata.json"
    metadata = json.loads(metadata_path.read_text(encoding="utf-8"))
    metadata[field] = value
    metadata_path.write_text(json.dumps(metadata, indent=2) + "\n", encoding="utf-8")
    with pytest.raises(ValueError, match=f"{field} metadata mismatch"):
        AUDIT.audit(tmp_path)
