# Results audit and reporting checklist

## Status

This checklist is for use after the manual `factorial-sensitivity` workflow completes. It is not evidence that the experiment has run. Do not fill in numerical findings until the raw artifacts are available and independently checked.

## 1. Record the run identity

Record all of the following before analysis:

- Git commit SHA of this repository.
- Git commit SHA of the original simulator (expected pin: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`).
- GitHub Actions run ID and final conclusion.
- Requested replicate count and base seed.
- Python and dependency versions from the run environment.
- Artifact download date and a checksum or immutable copy of the original artifact bundle.

Do not cite a mutable `main` branch as the sole reproducibility identifier.

## 2. Verify artifact completeness

For the main factorial experiment, verify:

- `scenario_summary.csv`, `scenario_summary.json`, and `per_run_results.csv` exist and are non-empty.
- There are exactly nine main scenarios.
- Each scenario has the requested number of valid replicate records.
- The per-run CSV has exactly nine times the requested replicate count.
- All scenario labels agree between summary and per-run files.
- Each 70% target-completeness record has realized completeness 0.70, given the current fixed-size mask design.
- Posterior summaries and interval endpoints are finite.
- Every credible interval has lower endpoint less than or equal to upper endpoint.
- Failed, excluded, or invalid runs are counted and explained rather than silently dropped.

For the missing-host control, verify six scenarios and the expected number of records. For the full-population oracle benchmark, verify two scenarios and expected counts. Check that the paired oracle summaries are identical where the design says the posterior is computed from the same full simulated population.

## 3. Interpret coverage with Monte Carlo uncertainty

For observed coverage proportion \(\hat p = k/n\), report the numerator, denominator, proportion, and a binomial uncertainty interval (the workflow uses Wilson intervals in its relevant analysis).

As a rough diagnostic, the Monte Carlo standard error is

\[
\operatorname{SE}(\hat p) \approx \sqrt{\hat p(1-\hat p)/n}.
\]

At \(n=500\) and coverage near 0.68, this is approximately 0.021, or 2.1 percentage points. Differences of only one or two percentage points are therefore not persuasive by themselves. Do not interpret an interval overlapping or excluding a reference value without also reporting the design, uncertainty, and effect size.

## 4. Report the full metric set

For every scenario, report:

- valid replicate count and any failed-run count;
- empirical 68% and 95% interval coverage, with uncertainty;
- bias and RMSE of the posterior median relative to the known simulation truth;
- mean and/or median 68% and 95% interval width;
- relevant posterior diagnostics and pull summaries, where defined;
- realized completeness and the selection mechanism;
- host-generation weights and inference weights.

For each contrast, state the direction and units of the effect, the interval estimate, the exact contrast definition, and whether it was pre-specified. The current contrast analysis has multiple exploratory comparisons and does not automatically adjust all intervals for multiplicity. Avoid presenting a list of intervals as if each were a confirmatory test.

## 5. Keep the scientific claims proportional

Permissible claims must describe the tested simulation design. For example:

- “Under the specified synthetic selection mechanism, scenario A had a different estimated coverage from scenario B.”
- “The estimated contrast was small relative to its Monte Carlo uncertainty.”
- “The result was inconclusive at the tested replicate count.”

Do not claim that a correction works for real galaxy surveys, that the simulation establishes the real value of \(H_0\), or that an oracle control is operationally available for an incomplete survey. The inverse-redshift missingness rule is a synthetic stress test, not a calibrated survey selection function. The fixed-size masks condition on the retained count.

## 6. Separate findings from follow-up work

Do not begin the flat-\(\Lambda\)CDM comparative experiment until the baseline factorial artifacts have been validated and archived. That follow-up should compare matched model generation/inference against a deliberately mismatched distance model and must state that the initial design fixes \(\Omega_m=0.3\), rather than marginalizing over it.

## 7. Final sign-off

- [ ] Run completed successfully.
- [ ] Artifacts archived with run ID and commit SHAs.
- [ ] Counts and scenario labels validated.
- [ ] Invalid and failed runs accounted for.
- [ ] Coverage uncertainty reported.
- [ ] Bias, RMSE, interval widths, and contrasts checked.
- [ ] Multiple-comparison limitations stated.
- [ ] Scientific language restricted to the synthetic design.
- [ ] Raw data, scripts, configuration, and analysis outputs retained.
- [ ] Independent review of code and interpretation requested before external scientific claims.
