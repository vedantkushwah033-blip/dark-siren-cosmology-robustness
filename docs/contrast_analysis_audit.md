# Contrast-analysis implementation audit

## Scope

This is a source-code audit of the current analysis script, not a report of simulation findings. The factorial experiment must still be run and its artifacts inspected before scientific conclusions are drawn.

## Confirmed behavior

- The contrast script computes eight named scenario comparisons across seven metrics, for 56 contrast rows.
- Each contrast is scenario B minus scenario A.
- Scenario seed blocks differ, so independent scenario bootstrap resampling is consistent with the current simulation design; the contrasts are not paired-run effects.
- Coverage contrasts are differences in empirical coverage proportions.
- The RMSE summary is calculated as the square root of mean squared posterior-median error within each scenario. The contrast table compares mean squared error, not the difference in RMSE itself.
- The bootstrap intervals are percentile intervals and are not adjusted for the 56 comparisons.

## Small implementation issue to fix

In `scripts/analyze_factorial_contrasts.py`, the `pull_median` calculation uses a truthiness fallback: `row.get("pull_median") or calculated_value`.

If the CSV contains a valid pull value of exactly zero, Python treats it as false and the script replaces it with the calculated value. This is a narrow edge case, but zero is a valid statistical value. The robust implementation should distinguish a missing/empty value from a numeric zero (for example, parse the field and use it whenever it is not empty).

## Reporting recommendations

1. Label the contrast metric as `squared_median_error` or mean squared error; do not describe it as an RMSE difference.
2. If an RMSE contrast is desired, compute the scenario-level RMSE in each bootstrap resample and subtract those two RMSEs. This is a different estimand from the difference in mean squared error.
3. State that the 56 intervals are exploratory and unadjusted for multiplicity.
4. Do not interpret statistical intervals as uncertainty in the real-world value of H0; they describe Monte Carlo variability under the specified synthetic model.
5. Preserve scenario labels, seed blocks, simulator commit, run ID, and artifact hashes alongside the output.

## Bottom line

The contrast script has a transparent, reproducible structure, but the meaning of each metric must remain explicit. Fix the zero-value fallback before relying on `pull_median` when that field is supplied by input data.