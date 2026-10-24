# Single-host fixed-truth selection control: results

**Status:** workflow passed; artifact inspected.  
**Workflow:** [Toy single-host fixed-truth selection control](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40099867069)  
**Commit:** `15d16e3b06cc9edc204df8d2f5690704a47b7f20`  
**Artifact:** `single-host-fixed-truth-control-results`; SHA-256 `88c39ea9f176c92fa601f50aaaef8823d67a2959f7c56bfb4bd5887c955c0ffa`.

## Design and integrity

The run completed on Python 3.11 with 300 replicates at each of nine fixed truths (52, 55, 60, 65, 70, 75, 80, 85, 88), using one host at redshift 0.03, Gaussian distance noise of 20 Mpc, detection interval (0,180) Mpc, prior [50,90], and a 4,001-point H0 grid. Each detected observation was analysed by both the uncorrected and selection-normalized likelihoods.

- Per-run records: 2,700 (9 truths × 300).
- Summary records: 27 (two method summaries plus one paired absolute-error contrast for each truth).
- Runtime manifest counts matched the CSVs.
- Focused tests passed; the run validated finite values and observation thresholds and uploaded the artifact.
- Artifact digest: `sha256:88c39ea9f176c92fa601f50aaaef8823d67a2959f7c56bfb4bd5887c955c0ffa`.

This is a toy-model control, not a real gravitational-wave measurement or a realistic detector-selection model.

## Fixed-truth results

Bias and RMSE are in the model's H0 units. Coverage is the empirical fraction of equal-tailed posterior intervals containing the known truth. The target nominal coverages are 68% and 95%; these are pointwise repeated-sampling diagnostics, not guaranteed properties of Bayesian credible intervals.

| True H0 | Method | Bias | RMSE | 68% coverage | 95% coverage |
|---:|---|---:|---:|---:|---:|
| 52 | Uncorrected | +7.531 | 8.893 | 27.0% | 89.7% |
| 52 | Corrected | +6.493 | 8.142 | 50.0% | 93.3% |
| 55 | Uncorrected | +6.084 | 8.016 | 65.0% | 94.7% |
| 55 | Corrected | +5.108 | 7.477 | 71.7% | 97.7% |
| 60 | Uncorrected | +3.996 | 7.379 | 77.3% | 97.3% |
| 60 | Corrected | +3.153 | 7.257 | 79.0% | 98.0% |
| 65 | Uncorrected | +2.686 | 7.727 | 71.0% | 97.7% |
| 65 | Corrected | +2.042 | 7.915 | 69.0% | 98.0% |
| 70 | Uncorrected | +1.872 | 7.255 | 74.7% | 99.7% |
| 70 | Corrected | +1.450 | 7.517 | 74.0% | 99.7% |
| 75 | Uncorrected | −1.485 | 7.018 | 80.7% | 98.7% |
| 75 | Corrected | −1.826 | 7.425 | 80.0% | 98.0% |
| 80 | Uncorrected | −3.301 | 6.468 | 87.3% | 99.3% |
| 80 | Corrected | −3.491 | 6.765 | 86.7% | 99.0% |
| 85 | Uncorrected | −6.955 | 8.917 | 71.3% | 96.0% |
| 85 | Corrected | −7.107 | 9.174 | 70.3% | 96.0% |
| 88 | Uncorrected | −9.390 | 10.857 | 20.0% | 91.3% |
| 88 | Corrected | −9.524 | 11.078 | 20.0% | 91.3% |

Each method's 68% and 95% coverage also has a Wilson 95% interval in `summary.csv`. With 300 replicates, uncertainty in the coverage fraction is material, but the large deficits at truths 52 and 88 are not plausibly explained by ordinary binomial Monte Carlo noise alone.

## Paired change in absolute error

The paired metric is corrected minus uncorrected absolute error. Negative values favor the correction; positive values favor the uncorrected estimator for this loss. The 95% intervals are percentile-bootstrap intervals over the 300 matched observations (5,000 draws; seed 20261031).

| True H0 | Mean paired absolute-error change | Bootstrap 95% interval |
|---:|---:|---:|
| 52 | −1.038 | [−1.062, −1.011] |
| 55 | −0.768 | [−0.830, −0.702] |
| 60 | −0.025 | [−0.119, +0.073] |
| 65 | +0.232 | [+0.150, +0.312] |
| 70 | +0.226 | [+0.169, +0.285] |
| 75 | +0.274 | [+0.229, +0.319] |
| 80 | +0.179 | [+0.151, +0.208] |
| 85 | +0.152 | [+0.126, +0.179] |
| 88 | +0.134 | [+0.110, +0.161] |

The correction reduces absolute error at truths 52 and 55 in this experiment, is inconclusive at 60, and increases it at truths 65–88. This is another direct reason not to claim that selection normalization universally improves fixed-truth point estimation.

## What this control clarifies

1. **The two-host mixture is not necessary for edge-truth problems to appear.** In this one-host control, 68% coverage is 27% (uncorrected) and 50% (corrected) at truth 52, and only 20% for both methods at truth 88. The earlier two-host maps also showed severe edge undercoverage. Thus mixture complexity alone cannot explain the boundary behavior.
2. **The correction has truth-dependent effects.** It improves bias and RMSE at some lower truths, but increases RMSE and paired absolute error at most interior-to-upper truths in this run.
3. **Good 95% coverage does not imply precise inference.** Several cells have 95% coverage near or above nominal while median estimates remain biased and 68% coverage can be poor.
4. **This does not isolate a single cause.** Prior boundaries, the selected-data distribution, the broad noise model, posterior geometry, and the finite design can all contribute. Further controls are needed before attributing the pattern to any one mechanism.

## Scope limits and next steps

This model has one fixed host, Gaussian distance noise, a hard finite detection interval, and a low-redshift distance approximation. It omits a realistic detector statistic, injection-based selection, a source-population model, galaxy-catalogue completeness, and cosmological distance-redshift effects. It cannot establish how a real LVK analysis would behave.

The next useful analysis is a matched comparison of single-host and two-host results under the same truths, prior, detection interval, noise, replicate count, and posterior grid. That would improve attribution while preserving the same simulated data design. Any later expansion should be prespecified and should report all truths, intervals, and negative findings.
