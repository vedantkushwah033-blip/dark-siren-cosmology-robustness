# Detection-threshold misspecification: results

## Status and reproducibility

The experiment completed successfully on Python 3.10, 3.11, and 3.12 in GitHub Actions run [40080814841](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40080814841), commit `b835a650d2`. Each matrix job uploaded a 6,000-row `per_run.csv`. The full test suite also passed at that commit in [run 40080814829](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40080814829).

The CSV artifact from the successful run contained the expected 6,000 rows (5 truths × 300 replicates × 4 assumed thresholds). I independently recomputed the summary metrics and paired bootstrap contrasts from the per-run observations. A follow-up workflow-file correction fixes the multiline artifact path so future runs also upload `summary.csv` and `paired_threshold_contrasts.csv`; that later commit must be validated separately.

## Design

Observations are generated under the true detection interval 0 < distance < 120 Mpc, then analyzed using assumed upper thresholds 120, 150, 180, and 210 Mpc. The synthetic host mixture has redshifts 0.02 and 0.04 with weights 0.7 and 0.3; Gaussian distance noise has sigma 20 Mpc. The prior is [50,90], with 300 paired observations per true H0 value (52, 60, 70, 80, 88), seed 20261027, and a 4,001-point posterior grid.

## Recomputed performance summaries

Coverage is the fraction of simulation intervals containing the generating truth. Widths are in the same distance-parameter units as H0. Values are Monte Carlo estimates from 300 replicates per cell.

| True H0 | Assumed max (Mpc) | Bias | MAE | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
|52|120|+11.174|11.174|12.285|0.047|0.897|22.593|36.968|
|52|150|+12.217|12.217|13.043|0.000|0.847|22.204|36.712|
|52|180|+11.697|11.697|12.624|0.000|0.857|22.328|36.844|
|52|210|+12.095|12.095|13.099|0.000|0.870|23.386|37.076|
|60|120|+6.622|7.087|9.158|0.783|0.987|22.575|36.409|
|60|150|+7.306|7.358|9.256|0.717|0.980|21.965|35.988|
|60|180|+6.924|7.090|9.155|0.763|0.987|22.293|36.172|
|60|210|+7.488|7.645|9.725|0.740|0.980|22.969|36.347|
|70|120|−0.336|5.788|6.721|0.950|1.000|22.025|35.583|
|70|150|+0.057|5.204|6.076|0.950|1.000|21.379|35.053|
|70|180|−0.168|5.566|6.457|0.950|1.000|21.801|35.281|
|70|210|+0.462|5.695|6.579|0.943|1.000|22.137|35.382|
|80|120|−7.451|7.960|10.211|0.763|1.000|21.032|34.285|
|80|150|−7.289|7.723|9.709|0.760|1.000|20.532|33.735|
|80|180|−7.347|7.887|10.005|0.767|1.000|20.941|33.986|
|80|210|−6.733|7.384|9.567|0.807|1.000|21.004|33.971|
|88|120|−13.764|13.764|15.558|0.087|0.840|19.976|33.020|
|88|150|−13.686|13.686|15.248|0.073|0.840|19.621|32.526|
|88|180|−13.658|13.658|15.390|0.113|0.853|19.950|32.775|
|88|210|−13.106|13.106|14.880|0.117|0.910|19.890|32.670|

## Paired contrasts against the correctly specified 120 Mpc threshold

For each wrong threshold, the contrast is wrong-threshold minus 120-Mpc result, using the same 300 observations. Intervals below are percentile-bootstrap 95% intervals for the mean paired difference (5,000 resamples; seed 20261028). These were recomputed from the downloaded per-run artifact.

Selected absolute-error and 68%-coverage contrasts:

| True H0 | Assumed max | Δ MAE (95% bootstrap interval) | Δ 68% coverage (95% bootstrap interval) |
|---:|---:|---:|---:|
|52|150|+1.044 [+0.983, +1.104]|−0.047 [−0.073, −0.023]|
|52|180|+0.523 [+0.479, +0.567]|−0.047 [−0.070, −0.023]|
|52|210|+0.922 [+0.912, +0.931]|−0.047 [−0.073, −0.023]|
|60|150|+0.271 [+0.177, +0.364]|−0.067 [−0.097, −0.040]|
|60|180|+0.003 [−0.049, +0.055]|−0.020 [−0.037, −0.007]|
|60|210|+0.558 [+0.483, +0.629]|−0.043 [−0.067, −0.020]|
|70|150|−0.584 [−0.642, −0.528]|+0.000 [+0.000, +0.000]|
|70|180|−0.222 [−0.257, −0.189]|+0.000 [+0.000, +0.000]|
|70|210|−0.093 [−0.181, −0.004]|−0.007 [−0.017, +0.000]|
|80|150|−0.236 [−0.303, −0.173]|−0.003 [−0.010, +0.000]|
|80|180|−0.072 [−0.106, −0.040]|+0.003 [+0.000, +0.010]|
|80|210|−0.576 [−0.629, −0.522]|+0.043 [+0.023, +0.067]|
|88|150|−0.078 [−0.145, −0.013]|−0.013 [−0.027, −0.003]|
|88|180|−0.106 [−0.136, −0.076]|+0.027 [+0.010, +0.047]|
|88|210|−0.658 [−0.682, −0.633]|+0.030 [+0.013, +0.050]|

For the coverage contrast, the indicator differences are discrete and bootstrap intervals can be degenerate at zero for cells where all paired coverage indicators agree. These should be interpreted as Monte Carlo summaries, not as proof of equivalence.

## Interpretation

1. **The correctly specified threshold is not uniformly best on finite-sample point-estimation metrics.** For true H0 = 70, some deliberately incorrect thresholds reduce MAE/RMSE in this run. This can happen by chance or because a misspecified likelihood shifts estimates toward the truth for a particular parameter configuration. It does not justify using an incorrect selection function.
2. **Threshold misspecification can reduce coverage.** At true H0 = 52, all three wrong thresholds reduce 68% coverage to zero in this particular sample; at H0 = 60, coverage also declines.
3. **Prior-boundary effects are prominent.** At H0 = 52 and 88, coverage is poor even with the correct threshold (4.7% and 8.7% for 68% intervals). The truth is close to a boundary of the [50,90] prior, and the broader experiment has already shown that prior support strongly affects coverage. This experiment does not isolate threshold error from those boundary effects.
4. **The 95% coverage values of 1.0 at some truths are not proof of perfect calibration.** They are empirical outcomes from 300 repetitions and require uncertainty intervals to interpret rigorously.
5. The effect is not monotone in assumed threshold. For example, the H0 = 60 RMSE is 9.158 at 120 Mpc, 9.256 at 150, 9.155 at 180, and 9.725 at 210. A simple “larger mismatch always causes worse performance” claim is not supported.

## Limitations and next steps

This remains a simplified Gaussian distance-threshold model, not a realistic gravitational-wave detector-selection pipeline. It excludes SNR-dependent efficiency, detector duty cycles, waveform dependence, multiple-event population inference, redshift uncertainty, and realistic galaxy-catalogue incompleteness.

The most important next experiment should address prior-boundary sensitivity and threshold misspecification jointly using a predeclared factorial design. Report coverage Wilson intervals as well as bootstrap contrasts, and validate the selection normalization against an independent numerical integral. Avoid selecting the assumed threshold based on which result looks best.

## Reproduction

```bash
python scripts/run_threshold_misspecification.py --replicates 300 --seed 20261027 --grid-size 4001 --output-dir results/threshold-misspecification
python -m pytest tests/test_threshold_misspecification.py -q
```
