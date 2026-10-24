# Coverage-transition map: validated results

## Provenance and validation

The coverage-transition workflow completed successfully for Python 3.10, 3.11, and 3.12, with 5,400 per-run rows and 36 summary rows in each artifact. The full test suite also passed on all three Python versions for commit `10ad5e5a3127c1fc31b9ee9888800ad0f300953b`.

- Workflow: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40075566569
- Full suite: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40075566548
- Configuration: truths 52, 55, 60, 65, 70, 75, 80, 85, 88; priors [50,90] and [40,100]; 300 observations per truth; seed 20261025; 4,001-point grid.

The three runtime artifacts agree to floating-point precision; the largest absolute difference among summary values checked was about 1.8×10⁻¹⁵. CSVs were inspected directly. The summary table below is from the Python 3.10 artifact; tiny runtime-level rounding differences do not affect interpretation.

## Fixed-truth coverage results

Coverage is the fraction of the 300 repeated simulations in which the equal-tailed posterior interval contained the fixed true value. Wilson intervals are included in the CSV artifacts.

| Truth | Prior | Method | Bias | RMSE | 68% coverage | 95% coverage |
|---:|---|---|---:|---:|---:|---:|
|52|[50,90]|Uncorrected|+13.320|14.343|0.0%|89.3%|
|52|[50,90]|Corrected|+11.336|12.496|0.0%|91.7%|
|52|[40,100]|Uncorrected|+10.221|11.952|73.0%|95.7%|
|52|[40,100]|Corrected|+7.331|9.810|79.7%|96.7%|
|55|[50,90]|Uncorrected|+11.695|13.064|56.0%|94.3%|
|55|[50,90]|Corrected|+9.882|11.464|62.7%|96.7%|
|55|[40,100]|Uncorrected|+9.890|12.440|71.7%|94.7%|
|55|[40,100]|Corrected|+7.136|10.772|74.3%|96.3%|
|60|[50,90]|Uncorrected|+7.965|9.652|72.3%|99.3%|
|60|[50,90]|Corrected|+6.219|8.399|77.7%|99.7%|
|60|[40,100]|Uncorrected|+7.199|10.656|78.0%|98.0%|
|60|[40,100]|Corrected|+4.812|9.718|84.3%|99.0%|
|65|[50,90]|Uncorrected|+4.917|7.507|84.0%|100.0%|
|65|[50,90]|Corrected|+3.413|6.884|91.7%|100.0%|
|65|[40,100]|Uncorrected|+5.449|9.871|82.7%|99.0%|
|65|[40,100]|Corrected|+3.238|9.520|85.3%|99.3%|
|70|[50,90]|Uncorrected|+2.220|5.800|95.7%|100.0%|
|70|[50,90]|Corrected|+0.895|5.734|98.0%|100.0%|
|70|[40,100]|Uncorrected|+3.053|8.584|87.7%|100.0%|
|70|[40,100]|Corrected|+0.733|9.212|91.0%|100.0%|
|75|[50,90]|Uncorrected|−1.887|6.118|98.3%|100.0%|
|75|[50,90]|Corrected|−3.229|7.077|96.7%|100.0%|
|75|[40,100]|Uncorrected|−0.892|9.268|92.3%|100.0%|
|75|[40,100]|Corrected|−3.243|11.139|94.0%|100.0%|
|80|[50,90]|Uncorrected|−6.178|8.465|85.3%|100.0%|
|80|[50,90]|Corrected|−7.542|9.834|81.3%|100.0%|
|80|[40,100]|Uncorrected|−5.519|11.104|95.0%|100.0%|
|80|[40,100]|Corrected|−8.187|14.038|94.3%|100.0%|
|85|[50,90]|Uncorrected|−11.206|13.087|53.7%|100.0%|
|85|[50,90]|Corrected|−12.605|14.621|42.3%|99.7%|
|85|[40,100]|Uncorrected|−10.088|15.051|88.7%|100.0%|
|85|[40,100]|Corrected|−12.316|17.677|76.0%|100.0%|
|88|[50,90]|Uncorrected|−13.625|15.233|7.7%|96.7%|
|88|[50,90]|Corrected|−14.974|16.738|5.7%|88.0%|
|88|[40,100]|Uncorrected|−12.323|16.799|76.3%|100.0%|
|88|[40,100]|Corrected|−14.577|19.619|57.3%|99.7%|

Bias is the mean posterior-median estimate minus the truth. RMSE is computed over the posterior-median estimation errors. (H_0) units are km s⁻¹ Mpc⁻¹.

## Interval-width tradeoff

Mean credible-interval widths below are calculated from the per-run artifact, not inferred from the coverage summary. Each cell uses 300 observations.

| Truth | Prior | Method | Mean 68% width | Mean 95% width |
|---:|---|---|---:|---:|
|52|[50,90]|Uncorrected|26.29|37.65|
|52|[50,90]|Corrected|25.13|37.51|
|52|[40,100]|Uncorrected|37.91|54.39|
|52|[40,100]|Corrected|35.75|54.31|
|55|[50,90]|Uncorrected|25.12|37.24|
|55|[50,90]|Corrected|24.20|37.18|
|55|[40,100]|Uncorrected|35.97|52.93|
|55|[40,100]|Corrected|34.17|52.99|
|60|[50,90]|Uncorrected|23.96|36.43|
|60|[50,90]|Corrected|23.39|36.48|
|60|[40,100]|Uncorrected|33.67|51.38|
|60|[40,100]|Corrected|32.60|51.52|
|65|[50,90]|Uncorrected|22.46|35.47|
|65|[50,90]|Corrected|22.30|35.68|
|65|[40,100]|Uncorrected|31.54|49.42|
|65|[40,100]|Corrected|31.06|49.68|
|70|[50,90]|Uncorrected|21.40|34.80|
|70|[50,90]|Corrected|21.69|35.17|
|70|[40,100]|Uncorrected|30.71|48.16|
|70|[40,100]|Corrected|30.85|48.50|
|75|[50,90]|Uncorrected|21.34|34.35|
|75|[50,90]|Corrected|21.68|34.77|
|75|[40,100]|Uncorrected|30.46|47.24|
|75|[40,100]|Corrected|30.47|47.65|
|80|[50,90]|Uncorrected|21.41|34.36|
|80|[50,90]|Corrected|21.94|34.80|
|80|[40,100]|Uncorrected|30.99|47.18|
|80|[40,100]|Corrected|30.93|47.61|
|85|[50,90]|Uncorrected|21.45|33.78|
|85|[50,90]|Corrected|21.77|34.23|
|85|[40,100]|Uncorrected|30.21|46.11|
|85|[40,100]|Corrected|29.96|46.60|
|88|[50,90]|Uncorrected|21.25|33.38|
|88|[50,90]|Corrected|21.60|33.86|
|88|[40,100]|Uncorrected|29.66|45.42|
|88|[40,100]|Corrected|29.53|45.94|

Widening the prior from [50,90] to [40,100] increases mean 68% interval widths by roughly 9–13 units in this setup. Coverage often improves near the edges, but the larger uncertainty is an important cost.

## Interpretation

1. **Boundary behavior is severe.** With the [50,90] prior, 68% coverage is 0/300 at truth 52 and 23/300 at truth 88 for uncorrected inference; corrected coverage is 0/300 and 17/300 respectively. The Wilson intervals in the CSV quantify Monte Carlo uncertainty around these fractions.
2. **Widening the prior often improves edge coverage.** At truth 52, 68% coverage under the wider prior is 73.0% uncorrected and 79.7% corrected. At truth 88 it is 76.3% and 57.3%, respectively. The improvement is not identical across methods.
3. **Correction is not uniformly better.** It lowers RMSE at truths 52–65 in both prior settings in this run, but at truths 75–88 it generally increases RMSE. At truth 70, the effect depends on the prior.
4. **95% coverage can be misleading without width.** Several settings achieve 100% empirical 95% coverage, but their credible intervals are broad. Coverage alone does not imply high precision.
5. **Paired observations matter.** The comparison uses the same detected observations across prior supports and inference methods, reducing noise due to unrelated Monte Carlo samples when comparing conditions.

## Scope and caveats

This is a simplified synthetic low-redshift host-mixture model, not an analysis of real LVK events or a real galaxy catalogue. It tests nine truth values and two prior intervals, with 300 repetitions each. Results are conditional on the simulation assumptions, detection interval, host mixture, noise model, and estimator. They do not establish universal superiority of either likelihood, nor do they justify selecting a prior after looking at these outcomes.

Fixed-truth coverage asks whether intervals contain a specified truth over repeated data sets. Prior-predictive simulation-based calibration (SBC) asks a different question about posterior ranks when truths are drawn from the prior. Both are useful and should not be conflated.

## Reproduction

```bash
python scripts/run_coverage_transition_map.py --replicates 300 --seed 20261025 --grid-size 4001 --output-dir results/coverage-transition
python -m pytest tests/test_coverage_transition_map.py -q
```

The workflow uploads `per_run.csv` and `summary.csv` for each supported Python version.
