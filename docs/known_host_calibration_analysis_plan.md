# Known-host H₀ calibration surface: analysis plan

**Status:** prespecified analysis plan; no calibration-surface results are reported in this document.  
**Companion protocol:** [Known-host calibration protocol](known_host_calibration_protocol.md).

## 1. Objective and scope

Evaluate the repeated-sampling behaviour of the repository's known-host, low-redshift, matched-model simulator and inference pipeline. This is a software/statistical calibration diagnostic under the stated synthetic model, not an observational estimate of H₀ and not evidence by itself that realistic dark-siren analyses are calibrated.

The planned design has nine truth values (50, 55, 60, 65, 70, 75, 80, 85, 90 km s⁻¹ Mpc⁻¹), four fractional distance uncertainties (0.05, 0.10, 0.15, 0.20), and 300 replicates in each cell by default. That is 36 cells and 10,800 raw records.

## 2. Validation gates before interpreting results

No scientific interpretation should begin until all of these conditions hold:

1. The designated workflow run is completed successfully.
2. Its artifact contains raw per-run observations, cell summaries, metadata, independent audit output, paired-comparison outputs, and rendered report.
3. The independent audit reports `PASS`, verifies the exact design and row counts, and independently reconstructs posterior summaries from raw observations.
4. The paired analysis runs only after that audit and verifies that its input hashes still match the audited files.
5. The report renderer succeeds and its input hashes match the audited files.
6. The final commit and workflow run ID are recorded alongside the artifact.
7. Any deviations from the default replicate count, seed, simulator revision, prior, or posterior grid are explicitly recorded.

A CI pass or a unit-test pass alone does not satisfy these gates.

## 3. Primary cell-level summaries

Report all 36 truth × uncertainty cells, not just a grand average.

### Coverage

For central 68% and 95% posterior intervals, report the number of intervals containing the generating truth, empirical coverage, and 95% Wilson binomial interval. With 300 replicates per cell, the approximate binomial standard error near nominal coverage is √[p(1-p)/300]; quote the observed standard error and Wilson interval rather than treating empirical coverage as exact.

Compare empirical coverage with the corresponding nominal level (0.68 or 0.95). Describe a discrepancy as a finite-simulation diagnostic, not proof of universal under- or over-coverage. The two coverage levels are nested and must not be treated as independent outcomes.

### Point-estimation error

Report median bias, mean posterior-median error, mean posterior-mean error, and RMSE separately. Keep signed mean error distinct from absolute error and RMSE; do not substitute one for another. Report the scale in km s⁻¹ Mpc⁻¹.

### Interval and posterior spread

Report mean 68% and 95% interval widths, mean posterior standard deviation, and median-pull diagnostics as available in the independently audited summary. Interpret interval width jointly with coverage and error; narrower intervals are not automatically better if coverage deteriorates.

### Prior boundaries

Report the two boundary truths (50 and 90) separately from the seven interior truths. A uniform prior on [50, 90] truncates posterior support at the boundary, so pooling boundary and interior cells can conceal materially different behaviour. Include the fraction of posterior medians near either prior edge, but do not use that metric alone to declare a cell invalid.

## 4. Comparisons across uncertainty levels

The intended common-random-number pairing is within a fixed true-H₀ value and replicate across the four uncertainty levels. Preserve that pairing in any comparisons. Do not describe different true-H₀ values as fully paired host-level realizations solely because the same seed integer was reused.

For each true H₀, use the 0.05 uncertainty cell as the reference and report paired replicate-level differences for:
- absolute posterior-median error;
- 68% and 95% posterior interval width;
- the 68% and 95% truth-inclusion indicators.

The implementation is `scripts/analyze_known_host_paired_comparisons.py`. It consumes the raw per-run CSV only after verifying the audit status and SHA-256 hashes, checks for complete replicate/noise blocks, and writes `paired_comparisons.csv`, `paired_comparisons.json`, and `paired_comparisons_report.md` to the experiment output directory.

For every truth, comparison noise level, and metric, the analysis reports the mean and median paired difference. Difference direction is **comparison noise level minus the 5% reference**. Its percentile bootstrap interval is for the mean paired difference; replicate IDs are resampled within each fixed-truth stratum, which preserves the pairing for each contrast. The default workflow uses 4,000 bootstrap resamples and records the deterministic bootstrap seed. These comparisons are descriptive; no multiplicity-adjusted confirmatory claims are made. The analysis does not use the four rows as if they were independent.

The main expected qualitative relationship is that larger distance uncertainty will generally widen posteriors and reduce precision. Do not assume coverage must move monotonically with uncertainty; measure it. Do not turn an expected trend into a pass/fail threshold after seeing results.

## 5. Monte Carlo uncertainty and reporting language

- Display per-cell sample size and uncertainty intervals alongside coverage.
- Avoid ranking cells by small differences that are compatible with Monte Carlo noise.
- Treat any across-cell trend as descriptive unless supported by a prespecified model or paired procedure.
- Do not claim nominal coverage merely because a nominal value lies inside a broad interval; state the estimate and interval.
- If an audit fails, investigate the failure and preserve the failed artifact/log. Do not selectively remove cells or raw records to obtain a pass.

## 6. Missingness, failures, and deviations

The expected number of raw records is exactly 10,800 for the default design. A missing cell, duplicated replicate key, unexpected seed, altered design grid, non-finite observation, invalid uncertainty, inconsistent standardized noise, summary mismatch, or checksum mismatch is an audit failure requiring investigation.

If the default run fails, correct the cause and rerun the same prespecified configuration when appropriate. If the configuration changes, identify the new run as a deviation or sensitivity analysis rather than silently replacing the default result. Report failed runs where they materially affect reproducibility.

## 7. Limits on conclusions

Even perfect empirical coverage in this experiment applies only to the specified matched low-redshift, known-host, Gaussian-distance model, prior support, grid, and simulation implementation. It does not test uncertain host identity, realistic cosmological distance-redshift relations, detector selection, catalogue incompleteness, peculiar velocities, calibration systematics, or a mismatch between the true and fitted likelihood.

The next stage should isolate these complications in separately prespecified experiments, with each experiment's estimand and audit conditions stated before examining its outcomes. Any paper should distinguish software verification, simulation calibration, and astrophysical validity as three different claims.

## 8. Expected Monte Carlo resolution at the default sample size

With 300 independent replicate blocks per truth, the binomial standard error of an empirical coverage estimate is approximately 0.0269 at coverage 0.68 and 0.0126 at coverage 0.95. These correspond to roughly 2.7 and 1.3 percentage points, respectively, for one standard error. Near the nominal values, a 95% Wilson interval will therefore typically span about ten percentage points for 68% coverage and about five percentage points for 95% coverage; the exact interval depends on the observed success count.

This is a planning calculation, not a result from the unrun calibration experiment. It means that a cell estimate such as 0.65 versus 0.68 should not be described as a meaningful difference on its own with n=300. Report the actual success count and Wilson interval. If the first run suggests a discrepancy worth resolving, increase the replicate count in a separately identified follow-up run (within the workflow's supported limit) rather than quietly changing the original protocol or stopping early based on the observed result.

For paired comparisons across uncertainty levels, the independent units are replicate IDs within a fixed truth stratum, not the four rows as if they were independent. Any bootstrap confidence interval for a paired metric must resample those replicate IDs as blocks, retaining all four uncertainty-level records for each sampled ID. Coverage indicators across uncertainty levels are correlated by design; do not calculate uncertainty for their differences using an independent-proportions formula.
