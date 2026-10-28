from __future__ import annotations

import json
import os
from pathlib import Path
import subprocess
import sys

import pandas as pd

ROOT = Path(__file__).resolve().parents[1]
RUNNER = ROOT / "scripts" / "run_catalogue_thinning_audit.py"
AUDITOR = ROOT / "scripts" / "audit_catalogue_thinning_audit.py"


def run_script(script: Path, *args: str, env: dict[str, str]) -> subprocess.CompletedProcess[str]:
    return subprocess.run(
        [sys.executable, str(script), *args],
        cwd=ROOT,
        env=env,
        text=True,
        capture_output=True,
        check=False,
    )


def test_catalogue_thinning_runner_and_independent_audit_round_trip(tmp_path: Path) -> None:
    output_dir = tmp_path / "audit"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    generated = run_script(
        RUNNER,
        "--n-parent", "1200",
        "--seed", "731",
        "--output-dir", str(output_dir),
        env=env,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr

    checked = run_script(AUDITOR, str(output_dir), env=env)
    assert checked.returncode == 0, checked.stdout + checked.stderr
    report = json.loads((output_dir / "catalogue_thinning_audit.json").read_text())
    assert report["status"] == "PASS"
    assert report["n_parent"] == 1200
    assert report["strata_reconstructed"] == 4


def test_independent_audit_rejects_modified_raw_records(tmp_path: Path) -> None:
    output_dir = tmp_path / "tampered-records"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    generated = run_script(
        RUNNER,
        "--n-parent", "1200",
        "--seed", "732",
        "--output-dir", str(output_dir),
        env=env,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr

    records_path = output_dir / "catalogue_thinning_records.csv"
    records = pd.read_csv(records_path)
    records.loc[0, "retained"] = 1 - int(records.loc[0, "retained"])
    records.to_csv(records_path, index=False, float_format="%.12g")

    checked = run_script(AUDITOR, str(output_dir), env=env)
    assert checked.returncode != 0
    assert "SHA-256" in checked.stdout + checked.stderr


def test_independent_audit_rejects_misdeclared_completeness_model(tmp_path: Path) -> None:
    output_dir = tmp_path / "tampered-summary"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    generated = run_script(
        RUNNER,
        "--n-parent", "1200",
        "--seed", "733",
        "--output-dir", str(output_dir),
        env=env,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr

    summary_path = output_dir / "catalogue_thinning_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["completeness_model"]["scale_mpc"] = 999.0
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    checked = run_script(AUDITOR, str(output_dir), env=env)
    assert checked.returncode != 0
    assert "parameters differ" in checked.stdout + checked.stderr

def test_independent_audit_rejects_tampered_binomial_standard_deviation(tmp_path: Path) -> None:
    output_dir = tmp_path / "tampered-binomial-sd"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    generated = run_script(
        RUNNER,
        "--n-parent", "1200",
        "--seed", "734",
        "--output-dir", str(output_dir),
        env=env,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr

    summary_path = output_dir / "catalogue_thinning_summary.json"
    summary = json.loads(summary_path.read_text())
    summary["strata"][0]["binomial_sd"] += 0.5
    summary_path.write_text(json.dumps(summary, indent=2) + "\n")

    checked = run_script(AUDITOR, str(output_dir), env=env)
    assert checked.returncode != 0
    assert "binomial_sd" in checked.stdout + checked.stderr

def test_independent_audit_rejects_nonbinary_raw_retained_flag(tmp_path: Path) -> None:
    output_dir = tmp_path / "nonbinary-retained-flag"
    env = os.environ.copy()
    env["PYTHONPATH"] = str(ROOT / "src") + os.pathsep + env.get("PYTHONPATH", "")
    generated = run_script(
        RUNNER,
        "--n-parent", "1200",
        "--seed", "735",
        "--output-dir", str(output_dir),
        env=env,
    )
    assert generated.returncode == 0, generated.stdout + generated.stderr

    records_path = output_dir / "catalogue_thinning_records.csv"
    records = pd.read_csv(records_path)
    # Pandas 3.0 rejects lossy assignment into an integer column; promote
    # before injecting the deliberately malformed fractional flag.
    records["retained"] = records["retained"].astype(float)
    records.loc[0, "retained"] = 0.5
    records.to_csv(records_path, index=False, float_format="%.12g")

    checked = run_script(AUDITOR, str(output_dir), env=env)
    assert checked.returncode != 0
    assert "binary values" in checked.stdout + checked.stderr

