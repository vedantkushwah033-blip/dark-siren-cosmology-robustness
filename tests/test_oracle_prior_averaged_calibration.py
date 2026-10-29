import importlib.util
from pathlib import Path

spec=importlib.util.spec_from_file_location("oracle_bias_diagnostic",Path(__file__).resolve().parents[1]/"scripts"/"diagnose_oracle_bias.py")
module=importlib.util.module_from_spec(spec);spec.loader.exec_module(module)


def test_exact_model_is_prior_averaged_calibrated():
    r=module.run("prior",300,11,401,81)
    assert abs(r["mean_bias"])<4*r["mean_bias_se"]
    assert .58<r["cov68"]<.78 and .84<r["cov90"]<.96


def test_negative_control_dropping_selection_normalisation_is_detected():
    r=module.run("prior",300,11,401,81,key="unconditional_log_likelihood")
    assert abs(r["mean_bias"])>4*r["mean_bias_se"] or r["cov68"]<.58
