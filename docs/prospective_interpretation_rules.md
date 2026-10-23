# Prospective interpretation rules for the synthetic dark-siren experiments

## Purpose and timing

Use this checklist before interpreting the first completed factorial run. It is intended to reduce selective interpretation; it is not a formal preregistration and does not turn the simulation into observational evidence. If the analysis plan changes after results are inspected, record the change and label the affected analysis exploratory.

## 1. Validate the run before interpreting metrics

Do not interpret results unless all of the following hold:

- The GitHub Actions run completed successfully, and its run ID, main-repository commit SHA, and pinned original-simulator commit SHA are recorded.
- The factorial output contains nine scenario summaries and exactly nine times the requested number of per-run rows.
- Every scenario has the requested replicate count, expected scenario label, and finite posterior summaries.
- For 70% scenarios, each simulated catalogue retains exactly 350 of 500 galaxies under the current design.
- Credible-interval endpoints are ordered, widths are non-negative, and reported coverage counts agree with the per-run records.
- Missing-host and full-population-oracle artifacts pass their own validations; their results are not silently substituted for the main factorial results.

A passing unit-test workflow is evidence about code behavior, not evidence that this experiment has run or that its scientific conclusions are valid.

## 2. Interpret the nine scenarios conditionally

The scenario matrix is targeted, not a complete factorial design. Describe contrasts as conditional comparisons between the two named scenarios. Do not claim that a contrast estimates a universal main effect of completeness, missingness mechanism, or host weighting.

For every contrast, state:
- scenario A and scenario B;
- the metric and direction (B minus A);
- the point estimate and bootstrap interval;
- the replicate count and relevant Monte Carlo uncertainty;
- the assumptions that differ between A and B.

Do not generalize a comparison beyond the tested completeness level, generated population, weighting rule, distance model, and likelihood.

## 3. Coverage and uncertainty

Report 68% and 95% coverage as proportions with replicate counts and Monte Carlo standard errors. At 500 independent replicates, coverage near 0.68 has a binomial Monte Carlo standard error of about 0.021; coverage near 0.95 has a standard error of about 0.0097.

Treat a small coverage difference as inconclusive when it is comparable to Monte Carlo noise. Do not use nominal coverage alone to rank methods: also report bias, RMSE, and interval width. A method with wider intervals may attain higher coverage simply by being less informative.

## 4. Bias, RMSE, and interval width

Report the exact estimand used by the code. In particular, distinguish:
- median bias: the median of replicate-level posterior-median errors;
- RMSE of posterior medians: the square root of the mean squared replicate-level posterior-median error;
- mean credible-interval width: the average width across replicates.

Do not label a statistic as posterior-mean bias if it is calculated from posterior medians. Interpret uncertainty intervals alongside point estimates, and retain null or unfavorable findings in the report.

## 5. Contrasts and multiplicity

The current analysis contains eight prespecified scenario pairs across seven metrics (56 reported comparisons). The bootstrap intervals are exploratory and are not multiplicity-adjusted. Do not treat an interval excluding zero as a confirmatory discovery in isolation. Discuss the pattern across related metrics and contrasts, identify the comparisons planned before seeing results, and state explicitly that multiplicity was not controlled.

If formal confirmatory claims are later desired, define a primary estimand and multiplicity procedure before rerunning or examining the confirmatory results.

## 6. Oracle controls

The missing-host oracle and full-population oracle are simulation-only benchmarks that use hidden population information. They can diagnose the gap between a specified idealized correction and an incomplete-catalogue analysis; they cannot establish that the required hidden information is available in real observations.

Report these controls separately from the main nine-scenario matrix. Do not call an oracle result a deployable observational method.

## 7. Allowed and disallowed conclusions

Allowed:
- “Under this synthetic setup, scenario B had a higher/lower estimated coverage or error metric than scenario A, with the stated uncertainty.”
- “The contrast was inconclusive at the current simulation precision.”
- “This result motivates a follow-up simulation using a more realistic selection function.”

Not supported by this experiment alone:
- claims that the method recovers the real Universe's Hubble constant;
- claims about actual LVK detections or a real galaxy survey;
- claims that one weighting method is universally optimal;
- claims that a simulation-only oracle correction is available for observational data;
- causal or factorial main-effect claims that the targeted scenario matrix cannot identify.

## 8. Required reporting sentence

Every public summary of this experiment should state that it is a synthetic simulation study using simplified modelling assumptions, not a measurement from real gravitational-wave data. Include the run ID and code revisions when reporting numerical findings.
