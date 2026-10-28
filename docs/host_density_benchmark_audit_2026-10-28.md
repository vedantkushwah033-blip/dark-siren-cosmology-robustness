# Audit scope for the host-density-to-H0 benchmark

The benchmark runner writes `replicate_metrics.csv` and `summary.json`. The audit script validates the scientific-scope label, design metadata, paired replicate/arm structure, seed mapping, consistency of recorded interval summaries, coverage flags, and recomputes arm-level bias, RMSE, empirical 68% coverage, interval width, and mean posterior mean from the raw CSV.

Run:

```bash
python scripts/audit_host_density_benchmark.py artifacts/host_density_benchmark
```

## What PASS currently means

A PASS means that the saved CSV and JSON artifacts are internally consistent under the checks implemented by the auditor. The auditor reports SHA-256 hashes for the two files. These hashes identify the audited files; they are not digital signatures.

## Important limitation: posterior grid is not independently audited

The current raw schema records posterior summaries and a scalar `posterior_integral`, but does not preserve the posterior density vector on the H0 grid. Therefore the auditor can check that the *recorded scalar* is approximately one, but cannot independently integrate the actual posterior density or reconstruct its mean, median, and credible interval from grid values.

The output explicitly reports `posterior_grid_independently_reconstructed: false`. Do not describe this as independent verification of posterior normalization.

A stronger schema should save either (a) the posterior density values and corresponding H0 grid per replicate and arm, or (b) a deterministic, versioned representation sufficient to reconstruct them exactly from preserved inputs. The independent audit should then recompute normalization, mean, quantiles, interval width, coverage, and aggregate metrics. Add adversarial tests that alter the posterior vector while leaving the summary columns unchanged.

## Scientific scope

A PASS does **not** validate the toy likelihood as a realistic gravitational-wave detector model, prove completeness correction is robust to a misspecified selection law, or establish a real measurement of H0. The experiment is a synthetic methodological study. Coverage estimates are noisy when replicate counts are small; report the numerator, denominator, and Monte Carlo uncertainty.

