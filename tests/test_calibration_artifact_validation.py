import csv
import json
import math
from pathlib import Path
import subprocess
import sys

import pytest

ROOT = Path(__file__).resolve().parents[1]
VALIDATOR = ROOT / "scripts" / "validate_calibration_artifacts.py"
LEVELS = (0.50, 0.68, 0.80, 0.90, 0.95)


def write_artifacts(root: Path, *, n: int = 10, bad_mcse: bool = False,
                    decreasing_width: bool = False, wrong_replicates: bool = False) -> None:
    root.mkdir(parents=True, exist_ok=True)
    summary = {"requested_replicates": n}
    for scenario in ("unknown_host", "known_host"):
        summary[scenario] = {
            "n": n,
            "coverage_68": 0.7,
            "coverage_95": 0.95,
            "coverage_68_mc_se": math.sqrt(0.7 * 0.3 / n),
            "coverage_95_mc_se": math.sqrt(0.95 * 0.05 / n),
        }
    (root / "summary.json").write_text(json.dumps(summary), encoding="utf-8")

    with (root / "coverage_curve.csv").open("w", newline="", encoding="utf-8") as f:
        fields = ["scenario", "nominal_level", "replicates", "coverage",
                  "coverage_mc_se", "mean_interval_width"]
        writer = csv.DictWriter(f, fieldnames=fields)
        writer.writeheader()
        for scenario in ("unknown_host", "known_host"):
            for i, level in enumerate(LEVELS):
                coverage = level
                width = float(i + 1)
                if decreasing_width and scenario == "unknown_host" and level == 0.90:
                    width = 2.5
                rep_count = n - 1 if wrong_replicates and scenario == "known_host" and level == 0.68 else n
                mcse = math.sqrt(coverage * (1.0 - coverage) / n)
                if bad_mcse and scenario == "unknown_host" and level == 0.68:
                    mcse += 0.01
                writer.writerow({
                    "scenario": scenario, "nominal_level": level, "replicates": rep_count,
                    "coverage": coverage, "coverage_mc_se": mcse,
                    "mean_interval_width": width,
                })
    for name in ("per_run_unknown_host.csv", "per_run_known_host.csv"):
        with (root / name).open("w", newline="", encoding="utf-8") as f:
            writer = csv.writer(f)
            writer.writerow(["replicate"])
            writer.writerows([[i] for i in range(n)])


def run_validator(root: Path) -> subprocess.CompletedProcess:
    return subprocess.run(
        [sys.executable, str(VALIDATOR), str(root)],
        capture_output=True, text=True, check=False,
    )


def test_validator_accepts_consistent_artifacts(tmp_path):
    write_artifacts(tmp_path)
    result = run_validator(tmp_path)
    assert result.returncode == 0, result.stderr


@pytest.mark.parametrize(
    ("kwargs", "message"),
    [
        ({"bad_mcse": True}, "MC standard error inconsistent"),
        ({"decreasing_width": True}, "interval width decreases"),
        ({"wrong_replicates": True}, "expected 10 replicates"),
    ],
)
def test_validator_rejects_inconsistent_coverage_curve(tmp_path, kwargs, message):
    write_artifacts(tmp_path, **kwargs)
    result = run_validator(tmp_path)
    assert result.returncode != 0
    assert message in result.stderr
