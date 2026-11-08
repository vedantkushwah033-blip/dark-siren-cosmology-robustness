# Host-density H0-grid sensitivity protocol

**Scope:** synthetic, low-redshift Euclidean toy benchmark only. This is a numerical sensitivity check, not an H0 measurement, detector validation, or posterior-calibration result.

## Design

The paired experiment runs the same benchmark seeds and the same parent-galaxy/event realizations at H0 grid resolutions 161, 321, 641, 1,281, and 3,201 points. The redshift grid, population generator, noise model, detection threshold, KDE bandwidth, and three benchmark arms remain fixed. Each resolution writes its own schema-v2 raw artifact and summary.

The three arms are:
- `oracle`: parent-truth redshifts, a diagnostic reference only;
- `catalogue_only`: observed catalogue without completeness correction;
- `known_completeness`: observed catalogue with the declared inclusion probabilities treated as known.

## Independent audit

For every resolution, the existing schema-v2 artifact auditor recomputes posterior normalization, mean, median, 16th/84th percentiles, interval width, and aggregate metrics from the saved posterior density vectors. The grid-sensitivity auditor then verifies checksums, record counts, paired replicate/arm keys, identical seeds, and every value in the comparison table.

The default numerical tolerance is 0.01 km s^-1 Mpc^-1 for the maximum absolute difference in posterior mean, median, and 68% interval endpoints/width relative to the finest grid. The two finest resolutions must both meet this tolerance relative to the finest grid before the experiment reports `PASS`. The report also records the coarsest grid meeting the tolerance and raw CSV storage bytes at each resolution.

## Observed exploratory result (workflow run 42437607983)

The configured CI experiment used 5 paired replicates, 3,000 parent galaxies, 20 detected events per replicate, and a fixed 301-point redshift grid. The independent artifact audit reported `PASS` at the predeclared 0.01 km s^-1 Mpc^-1 tolerance.

| H0 grid points | Maximum absolute difference vs. 3,201 across mean, median, 68% endpoints/width (km s^-1 Mpc^-1) | Raw replicate CSV bytes |
|---:|---:|---:|
| 161 | 0.006320 | 68,168 |
| 321 | 0.001592 | 135,536 |
| 641 | 0.000404 | 275,451 |
| 1,281 | 0.000083 | 564,557 |
| 3,201 | 0 | 1,380,926 |

The 161-point grid was the coarsest tested resolution meeting the tolerance. Storage increased by about 20.3 times from 161 to 3,201 points because the schema-v2 artifact preserves full posterior density vectors. The workflow's grid-sensitivity tests passed on Python 3.10, 3.11, and 3.12; the full configured experiment and independent audit ran on Python 3.12.

This is a small exploratory numerical check, not a reliable estimate of repeated-sampling behavior: it has only five replicates. The 161-point result is not a general recommendation for other likelihood widths, priors, or model configurations. Re-run convergence checks when those change.

## Interpretation limits

- Passing means only that the saved numerical summaries are stable under the tested H0-grid refinements and tolerance.
- It does not establish calibration, coverage, correctness of the event likelihood, robustness to KDE bandwidth, or scientific validity of the toy selection law.
- A failure is reported as `UNDERRESOLVED`; do not silently relax the tolerance after seeing results.
- The grid experiment is a computational diagnostic and must not be presented as evidence of estimator superiority.
