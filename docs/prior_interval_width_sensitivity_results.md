# Prior interval-width sensitivity audit: validated results

**Workflow:** [Prior interval-width sensitivity audit](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40074452995)  
**Commit:** `a947c1db7f179caf55b2dff60dc23dd8a5e10b9a`  
**Status:** experiment, focused tests, and artifact uploads passed on Python 3.10, 3.11, and 3.12.  
**Design:** 300 paired synthetic observations per fixed truth; seed 20261024; 4,001-point posterior grid.

This is a controlled synthetic model, not a real gravitational-wave or galaxy-catalogue analysis.

## Design and outputs

Truths are (H_0=52,70,88); prior supports are [50,90], [40,100], [30,110]. Both uncorrected and selection-corrected posteriors are analyzed for the same observations under all priors. The model uses two host redshifts (0.02, 0.04), weights (0.7, 0.3), Gaussian distance noise of 20 Mpc, and detection interval (0<d_{obs}<180) Mpc.

The CI artifacts contain 18 summary rows (truth × prior × method) and 3,600 paired contrast rows (truth × replicate × wider prior × method). Python-version artifacts were inspected and have identical summary values.

## Mean credible-interval widths

Widths are in the units of (H_0), conventionally km s(^{-1}) Mpc(^{-1}).

| True (H_0) | Prior | Method | Mean 68% width | Median 68% width | Mean 95% width | Median 95% width |
|---:|---|---|---:|---:|---:|---:|
| 52 | [50,90] | Uncorrected | 26.154 | 23.927 | 37.569 | 37.865 |
| 52 | [50,90] | Corrected | 25.142 | 23.328 | 37.456 | 37.505 |
| 52 | [40,100] | Uncorrected | 37.583 | 38.370 | 54.082 | 55.669 |
| 52 | [40,100] | Corrected | 35.556 | 34.488 | 54.022 | 55.500 |
| 52 | [30,110] | Uncorrected | 45.429 | 47.782 | 64.811 | 66.267 |
| 52 | [30,110] | Corrected | 43.744 | 45.611 | 65.162 | 66.512 |
| 70 | [50,90] | Uncorrected | 21.711 | 22.835 | 35.013 | 36.214 |
| 70 | [50,90] | Corrected | 21.894 | 22.150 | 35.338 | 36.506 |
| 70 | [40,100] | Uncorrected | 31.468 | 30.034 | 48.727 | 49.533 |
| 70 | [40,100] | Corrected | 31.260 | 29.790 | 49.024 | 49.570 |
| 70 | [30,110] | Uncorrected | 38.532 | 38.061 | 58.160 | 58.391 |
| 70 | [30,110] | Corrected | 37.799 | 36.857 | 58.372 | 58.522 |
| 88 | [50,90] | Uncorrected | 21.816 | 21.205 | 33.896 | 35.031 |
| 88 | [50,90] | Corrected | 22.138 | 21.748 | 34.364 | 35.611 |
| 88 | [40,100] | Uncorrected | 30.383 | 28.183 | 46.080 | 45.465 |
| 88 | [40,100] | Corrected | 30.233 | 29.198 | 46.620 | 46.336 |
| 88 | [30,110] | Uncorrected | 36.558 | 34.078 | 55.481 | 54.408 |
| 88 | [30,110] | Corrected | 36.360 | 35.146 | 56.205 | 55.259 |

## Paired change in mean interval width versus [50,90]

| Truth | Method | Wider prior | Increase in mean 68% width | Increase in mean 95% width |
|---:|---|---|---:|---:|
| 52 | Uncorrected | [40,100] | +11.429 | +16.513 |
| 52 | Corrected | [40,100] | +10.413 | +16.565 |
| 52 | Uncorrected | [30,110] | +19.275 | +27.242 |
| 52 | Corrected | [30,110] | +18.601 | +27.706 |
| 70 | Uncorrected | [40,100] | +9.757 | +13.713 |
| 70 | Corrected | [40,100] | +9.366 | +13.686 |
| 70 | Uncorrected | [30,110] | +16.821 | +23.147 |
| 70 | Corrected | [30,110] | +15.905 | +23.034 |
| 88 | Uncorrected | [40,100] | +8.568 | +12.184 |
| 88 | Corrected | [40,100] | +8.095 | +12.256 |
| 88 | Uncorrected | [30,110] | +14.743 | +21.585 |
| 88 | Corrected | [30,110] | +14.222 | +21.841 |

The paired CSV stores per-observation differences and absolute differences. Mean width increases above are differences of means for the same paired observations.

## Interpretation

1. **Wider support consistently produces wider average credible intervals in this experiment.** Mean 68% widths increase by about 8–11 units under [40,100] and 14–19 units under [30,110]. Mean 95% widths increase by about 12–17 and 22–28 units, respectively.
2. **The prior can improve edge-truth coverage at the cost of precision.** The earlier coverage audit found severe undercoverage near the [50,90] boundaries, and wider support markedly improved coverage in several edge cases. This width audit shows that the improvement accompanies substantially wider intervals.
3. **This is not automatically a flaw.** If the baseline prior excludes plausible values, a wider prior can reduce artificial truncation. However, wider intervals are less precise, and a wider prior is not automatically more scientifically appropriate.
4. **Selection correction has smaller effects on width than prior support in these runs.** Corrected and uncorrected mean widths are often close relative to the changes caused by widening the prior. This statement is specific to this synthetic model.
5. **Coverage, bias, RMSE, and interval width must be reported together.** Neither nominal coverage nor narrow intervals alone establish good inference.

## Limitations

- The model uses only two fixed host redshifts and simplified Gaussian distance errors.
- The observations are synthetic and conditioned on a finite detection interval.
- Only three truth values and three prior supports were studied.
- These results do not establish a preferred prior or validate an observational (H_0) measurement.
- Posterior grid stability is a separate numerical question from model adequacy and frequentist calibration.

## Reproduction

```bash
python scripts/run_prior_interval_width_sensitivity.py --replicates 300 --seed 20261024 --grid-size 4001 --output-dir results/prior-interval-width
python -m pytest tests/test_prior_interval_width_sensitivity.py -q
```
