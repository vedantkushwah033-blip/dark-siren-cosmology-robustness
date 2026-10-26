# Independent audit of distance-uncertainty sensitivity artifacts

The sensitivity runner writes `replicates.csv` and `summary.json`. Before interpreting a result, run the independent reconstruction audit:

```bash
python scripts/run_distance_uncertainty_sensitivity.py --replicates 300 --seed 20261026 --bootstrap-resamples 2000 --output-dir artifacts/distance_uncertainty_sensitivity
python scripts/audit_distance_uncertainty_sensitivity.py artifacts/distance_uncertainty_sensitivity
```

The auditor does not import the runner or the likelihood module. It independently:
- checks the expected truth × fractional-error × replicate × arm design and uniqueness of every record;
- reconstructs true distances, generated uncertainty scales, and observations from the recorded standardized noise;
- verifies the exact paired observations for Arms A and B;
- recomputes the three Gaussian log-likelihood forms and posterior grid summaries;
- reconstructs coverage, Wilson intervals, pull summaries, error, RMSE, and interval-width summaries from the raw rows;
- independently reconstructs the paired replicate-level A-minus-B coverage, absolute-error, squared-error, and interval-width contrasts and their bootstrap intervals;
- records SHA-256 hashes of the raw CSV and summary JSON in `audit_report.json`.

The audit is a consistency check, not proof that the model is physically complete or that a likelihood is universally correct. It does not turn the low-redshift one-known-host toy model into an LVK analysis. The default 300 replicates per cell are exploratory; use the prespecified larger run for confirmatory comparisons. Report prior-boundary truths (H0=50 and 90) separately from interior truths.

The CI run uses only two replicates per cell to check code execution and the audit's ability to reconstruct the artifacts. These smoke-test outputs are not scientific results.
