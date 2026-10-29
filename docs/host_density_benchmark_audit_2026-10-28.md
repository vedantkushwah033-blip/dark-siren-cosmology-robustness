# Audit scope for the host-density-to-H0 benchmark

The benchmark runner writes `replicate_metrics.csv` and `summary.json`. Each replicate/arm row now includes the normalized posterior density values as JSON on the declared uniform H0 grid. The audit script independently reconstructs that grid, integrates it using a pure-Python trapezoid implementation, recomputes the posterior mean and equal-tail quantiles, checks interval width and coverage flags, and recomputes arm-level bias, RMSE, empirical 68% coverage, interval width, and mean posterior mean from raw records.

Run:

```bash
python scripts/audit_host_density_benchmark.py artifacts/host_density_benchmark
```

## What PASS currently means

A PASS means that the saved CSV and JSON artifacts are internally consistent under the checks implemented by the auditor. The auditor reports SHA-256 hashes for the two files. These hashes identify the audited files; they are not digital signatures.

## Posterior-grid audit and its boundary

The current raw schema preserves the normalized posterior density vector as JSON for every replicate and arm. The auditor reconstructs the uniform H0 grid from the declared metadata and independently recalculates the trapezoid integral, posterior mean, median, equal-tail 68% interval, interval width, and coverage consistency. Adversarial tests alter the saved density vector while leaving summary columns unchanged and require rejection.

The artifact audit is still not an independent scientific reimplementation of the likelihood or the KDE. It verifies internal consistency of saved posterior vectors and summaries, not whether the likelihood is a realistic gravitational-wave model or whether the posterior density was generated from the correct scientific assumptions. The file hashes identify the audited files; they are not digital signatures.

## Scientific scope

A PASS does **not** validate the toy likelihood as a realistic gravitational-wave detector model, prove completeness correction is robust to a misspecified selection law, or establish a real measurement of H0. The experiment is a synthetic methodological study. Coverage estimates are noisy when replicate counts are small; report the numerator, denominator, and Monte Carlo uncertainty.



## Bandwidth convention correction (2026-10-29)

The benchmark runner now uses the weighted Scott rule by default (`bandwidth_factor=None`). An explicit factor is defined by the current estimator as `h = factor × weighted_SD`, not as a multiplier of Scott's bandwidth. The former hard-coded `bandwidth_factor=0.8` therefore represented a broad physical bandwidth under the corrected convention and is no longer the benchmark default. Previously generated benchmark artifacts and numerical summaries must be regenerated before comparison or interpretation; the source change alone does not validate new scientific conclusions.
