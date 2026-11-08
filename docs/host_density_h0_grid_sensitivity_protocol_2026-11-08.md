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

## Interpretation limits

- Passing means only that the saved numerical summaries are stable under the tested H0-grid refinements and tolerance.
- It does not establish calibration, coverage, correctness of the event likelihood, robustness to KDE bandwidth, or scientific validity of the toy selection law.
- A failure is reported as `UNDERRESOLVED`; do not silently relax the tolerance after seeing results.
- The grid experiment is a computational diagnostic and must not be presented as evidence of estimator superiority.
