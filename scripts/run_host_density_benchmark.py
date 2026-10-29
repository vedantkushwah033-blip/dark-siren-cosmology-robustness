"""Command-line entry point for the exploratory host-density H0 benchmark."""
from __future__ import annotations
import argparse, hashlib, json, os
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
    p.add_argument("--bandwidth-factor", type=float, default=None,
                   help="Physical factor times weighted SD of the original sample; default is weighted Scott.")
    p.add_argument("--supplied-midpoint", type=float, default=0.22,
                   help="Midpoint assumed in the inference correction; generator truth remains 0.22.")
    p.add_argument("--supplied-width", type=float, default=0.05,
                   help="Width assumed in the inference correction; generator truth remains 0.05.")
    p.add_argument("--output-dir", type=Path, default=Path("artifacts/host_density_benchmark"))
    a = p.parse_args()
    rows, summary = run_benchmark(n_replicates=a.replicates, n_parent=a.n_parent,
        n_detected=a.n_detected, seed=a.seed, h0_points=a.h0_points, z_points=a.z_points,
        bandwidth_factor=a.bandwidth_factor, supplied_midpoint=a.supplied_midpoint, supplied_width=a.supplied_width)
    a.output_dir.mkdir(parents=True, exist_ok=True)
    rows.to_csv(a.output_dir/"replicate_metrics.csv", index=False)
    summary_path = a.output_dir/"summary.json"
    summary_path.write_text(json.dumps(summary, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    file_hashes = {}
    for path in (a.output_dir/"replicate_metrics.csv", summary_path):
        file_hashes[path.name] = {"sha256": hashlib.sha256(path.read_bytes()).hexdigest(),
                                  "bytes": path.stat().st_size}
    manifest = {
        "schema_version": 1,
        "git_commit": os.environ.get("GITHUB_SHA"),
        "replicate_count": int(summary["n_replicates"]),
        "seed_start": int(summary["seed_start"]),
        "seed_end": int(summary["seed_start"]) + int(summary["n_replicates"]) - 1,
        "configuration": {key: summary[key] for key in (
            "n_parent", "n_detected", "h0_true_km_s_mpc", "sigma_mpc",
            "detection_threshold_mpc", "h0_grid", "redshift_grid", "bandwidth_factor",
            "supplied_inclusion_midpoint", "supplied_inclusion_width")},
        "files": file_hashes,
    }
    (a.output_dir/"manifest.json").write_text(json.dumps(manifest, indent=2, sort_keys=True)+"\n", encoding="utf-8")
    print(json.dumps(summary, indent=2, sort_keys=True))
if __name__ == "__main__": main()
