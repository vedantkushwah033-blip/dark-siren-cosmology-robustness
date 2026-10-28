# Independent audit for the host-density-to-H0 benchmark

The benchmark runner writes `replicate_metrics.csv` and `summary.json`. The audit script validates the scientific-scope label, design metadata, paired replicate/arm structure, seed mapping, posterior normalization, interval ordering and width, coverage flags, and independently recomputes the arm-level bias, RMSE, empirical 68% coverage, interval width, and mean posterior mean from the raw CSV.

Run:

```bash
python scripts/audit_host_density_benchmark.py artifacts/host_density_benchmark
```

A PASS verifies internal consistency and reproducibility of the saved summaries; it does **not** validate the toy likelihood as a realistic gravitational-wave detector model, prove completeness correction is robust to a misspecified selection law, or establish a real measurement of H0. The audit also reports SHA-256 hashes for the two artifacts. Those hashes identify the audited files; they are not a digital signature.
