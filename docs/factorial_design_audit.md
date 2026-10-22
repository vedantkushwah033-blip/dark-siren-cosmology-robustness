# Factorial experiment: design audit and interpretation

## Source-verified selection behavior

The pinned simulator's random, redshift-dependent, and sky-dependent selection functions all select exactly `round(N * completeness)` galaxies without replacement. The faint-selection function also keeps exactly that number, choosing the highest luminosity proxies. Thus, for the current `N=500` and target completeness `0.7`, each mechanism retains 350 of 500 galaxies (70%) by construction. This is not merely an expected or approximate completeness.

The mechanisms differ in *which* galaxies are retained, not in the number retained. Random selection is a size-matched baseline; redshift-dependent selection favors smaller redshifts through a score proportional to `1/(z + 1e-9)`; faint selection retains the most luminous objects; and sky selection applies lower sampling scores in a defined RA patch while maintaining the same sample size.

## Current factorial matrix

The nine scenarios vary completeness mechanism and host weighting, including:
- complete catalogues with matched host/inference weighting;
- 70% random or redshift-dependent catalogues;
- luminosity-weighted host generation with matched luminosity inference weights;
- luminosity-weighted host generation with uniform inference weights.

The contrast script compares eight prespecified pairs over seven metrics. Its bootstrap resamples each scenario independently because the experiment currently uses different seed blocks per scenario. The resulting contrasts are independent-scenario comparisons, not paired-run effects.

## Statistical interpretation

For each reported coverage proportion, include the number of simulations and its Monte Carlo standard error, approximately `sqrt(p * (1-p) / n)`. For example, with 500 runs and observed coverage near 0.68, the Monte Carlo standard error is about 2.1 percentage points. A difference of only one or two percentage points between scenarios would therefore be weak evidence on its own at this run count.

For bias, RMSE, and interval widths, report the scenario estimate plus a bootstrap interval. Avoid claiming a scenario is better merely because its point estimate is numerically smaller. State whether comparisons were prespecified and acknowledge that multiple exploratory contrasts can produce chance findings.

## Scope and limitations

These experiments are synthetic and use a simplified low-redshift relation `d_L ≈ cz/H0`, a Gaussian distance likelihood, a finite H0 grid, and a generated galaxy population. They do not constitute an analysis of real LVK events or a measurement of the real Hubble constant. The missing-host term uses hidden simulated galaxies and is an oracle validation device, not a directly observable correction.

Before interpreting outputs, confirm that the manual workflow completed successfully and archive its artifact bundle. The automated unit-test workflow validates code behavior, not the scientific hypotheses.
