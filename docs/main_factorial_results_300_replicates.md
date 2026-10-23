# Main paired factorial experiment: 300-replicate results

## Execution and validation

The corrected paired factorial workflow completed successfully.

- Workflow: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39959387009
- Artifact: `paired-factorial-main-300-replicates` (artifact ID `12503494036`)
- Pinned simulator commit: `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`
- Design: 9 scenarios × 300 replicates = 2,700 per-run records.
- Paired contrasts: 49 rows, each based on 300 paired universes.
- Diagnostic figures were generated and included in the workflow artifact.

The earlier run `39958960480` failed only at a validation check that incorrectly expected 490 contrast rows instead of 49. The simulation had produced 9 scenario summaries, 2,700 per-run records, and 49 contrast rows. The corrected run completed validation and artifact upload successfully.

## Scenario summaries

All scenarios use synthetic data with fixed generating value $H_0=70$ km s$^{-1}$ Mpc$^{-1}$. The reported values describe this simulator and this design only.

| Scenario | Median bias summary | RMSE, posterior median | RMSE, posterior mean | 68% coverage | 95% coverage | Mean 68% width | Host observed |
|---|---:|---:|---:|---:|---:|---:|---:|
| Complete, uniform generation/inference | 4.150 | 6.198 | 4.823 | 93.7% | 99.7% | 22.458 | 100.0% |
| Random 70%, uniform generation/inference | 4.076 | 6.251 | 4.865 | 94.0% | 99.7% | 22.422 | 70.7% |
| Redshift-dependent 70%, uniform generation/inference | 3.228 | 5.861 | 4.587 | 93.7% | 99.7% | 22.849 | 71.0% |
| Complete, luminosity generation/matched inference | 3.493 | 6.004 | 4.693 | 95.3% | 100.0% | 22.357 | 100.0% |
| Complete, luminosity generation/uniform inference | 3.456 | 6.055 | 4.697 | 95.0% | 100.0% | 22.409 | 100.0% |
| Random 70%, luminosity generation/matched inference | 3.495 | 6.038 | 4.729 | 95.0% | 100.0% | 22.297 | 68.3% |
| Random 70%, luminosity generation/uniform inference | 3.510 | 6.059 | 4.707 | 95.0% | 100.0% | 22.380 | 68.3% |
| Redshift-dependent 70%, luminosity generation/matched inference | 2.546 | 5.666 | 4.467 | 95.3% | 100.0% | 22.723 | 73.0% |
| Redshift-dependent 70%, luminosity generation/uniform inference | 2.737 | 5.746 | 4.487 | 95.0% | 100.0% | 22.765 | 73.0% |

Bias, RMSE, and interval widths are in km s$^{-1}$ Mpc$^{-1}$. The `median_bias` column is the simulator's reported median-bias summary; it is not the mean signed error. Coverage is the empirical proportion of intervals containing the fixed simulated truth.

## Interpretation

1. **Random catalogue incompleteness:** In the uniform-generation/uniform-inference block, random 70% catalogue completeness produced RMSE of 6.251 versus 6.198 for the complete catalogue. Empirical 68% coverage was 94.0% versus 93.7%, while mean 68% interval width was very similar (22.422 versus 22.458). These marginal summaries suggest a small change in this run, but paired contrasts and uncertainty should be used before claiming an effect.
2. **Redshift-dependent selection:** In the same block, the redshift-dependent selection scenario had lower posterior-median RMSE (5.861) and a lower median-bias summary (3.228) than the complete scenario. This does not imply that incompleteness improves real dark-siren inference: it is a result of the particular synthetic redshift selection, host distribution, and inference model.
3. **Host-weight mismatch:** In complete catalogues with luminosity-generated hosts, matched inference had RMSE 6.004 compared with 6.055 for uniform inference. The difference is small in marginal summaries. The paired contrast analysis, not comparison of rounded scenario summaries alone, is the intended way to quantify this.
4. **Coverage is not nominal here:** The empirical 68% interval coverage is around 93.7%–95.3%, and 95% interval coverage is 99.7%–100%. This shows that these Bayesian intervals are not achieving their nominal frequentist coverage under this fixed-truth simulation setup. Calibration needs to be interpreted alongside the separate matched-model simulation-based calibration; the two diagnostics answer different questions.
5. **Host visibility:** The realized host-observed fractions are stochastic and do not always equal the nominal 70% catalogue target. The true host is a single draw and its visibility depends on the generated population and mask. In luminosity-host scenarios, observed-host fractions are about 68.3% under random masking and 73.0% under the redshift-dependent mask.

## Precision and limitations

With 300 replicates, an observed coverage near 94% has a binomial Monte Carlo standard error of roughly 1.4 percentage points. Therefore small differences in coverage should not be overinterpreted. The experiment is still a synthetic single-event study with simplified distance inference, synthetic luminosity proxies, and a toy redshift-dependent selection rule. It does not establish effects in real LVK events or calibrated galaxy surveys.

The 49 paired contrast rows and their bootstrap intervals are available in the workflow artifact and should be examined before making contrast-specific claims. The scenario table alone is descriptive.

## Reproducibility artifacts

Download the artifact from the successful workflow page above. It contains scenario summaries, per-run results, paired contrasts, metadata, and generated diagnostic figures. The plotting script is `scripts/plot_paired_factorial_results.py`.

## Next steps

1. Inspect the paired-contrast bootstrap intervals and ensure the contrast labels are scientifically correct.
2. Audit the unexpectedly high 68% and 95% coverage against the exact interval construction and fixed-truth simulation design.
3. Add repeated-truth calibration or a clearly separated coverage grid over multiple fixed $H_0$ values.
4. Freeze the analysis and update manuscript figures and discussion only after these checks.
