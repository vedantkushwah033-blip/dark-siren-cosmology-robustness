"""Command-line entry point for the exploratory host-density H0 benchmark."""
from __future__ import annotations
import argparse, json
from pathlib import Path
from dark_siren_robustness.host_density_benchmark import run_benchmark

def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--replicates", type=int, default=10)
    p.add_argument("--n-parent", type=int, default=6000)
    p.add_argument("--n-detected", type=int, default=40)
    p.add_argument("--seed", type=int, default=20261028)
    p.add_argument("--h0-points", type=int, default=161)
    p.add_argument("--z-points", type=int, default=801)
    p.add_argument("--output-dir", type=Path, default=Path("artifacts/host_density_benchmark"))
    a = p.parse_args()
    rows, summary = run_benchmark(n_replicates=a.replicates, n_parent=a.n_parent,
        n_detected=a.n_detected, seed=a.seed, h0_points=a.h0_points, z_points=a.z_points)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(a.output_dir/"replicate_metrics.csv", index=False)
    (a.output_dir/"summary.json").write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
if __name__ == "__main__": main()
