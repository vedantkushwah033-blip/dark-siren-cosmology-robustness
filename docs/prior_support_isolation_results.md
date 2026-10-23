# Matched-model prior-support isolation audit

**Status:** completed; exploratory numerical/statistical diagnostic.  
**Workflow:** [Prior-support isolation audit, run 40005022824](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40005022824)  
**Artifact:** `prior-support-isolation-results` (contains `summary.csv` and `per_run.csv`).  
**Design:** 7 fixed truth values × 2 prior supports × 300 replicates = 4,200 posterior fits.

## Question

How much do the apparent bias and interval coverage change when the same simulated distance observations are analysed with a narrower prior, H0 ∈ [50, 90], versus a broader prior, H0 ∈ [30, 120] km s⁻¹ Mpc⁻¹?

This is a deliberately simple matched-model experiment. For each truth value, a single redshift z=0.03 is fixed and the observed distance is drawn from

$$
d_{\rm obs}\sim\mathcal N\!\left(\frac{c z}{H_{0,\rm true}},\; (0.10\,d_{\rm true})^2\right).
$$

The same observed distance is fitted under both priors for each replicate. The likelihood and data-generating distribution both use the same Gaussian distance model. The prior is uniform in H0 within its specified support. The grid has 3,001 points spanning 25–125; posterior density is zero outside each prior interval.

## Results

Coverage is the empirical fraction of central posterior intervals containing the fixed injected truth. Error and width units are km s⁻¹ Mpc⁻¹. Coverage MCSE is approximately the binomial standard error at n=300.

| True H0 | Prior | Mean median error | RMSE | 68% coverage (MCSE) | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 55 | 50–90 | +3.28 | 5.65 | 72.3% (2.6 pp) | 96.3% | 10.79 | 21.00 |
| 55 | 30–120 | +2.03 | 6.23 | 67.7% (2.7 pp) | 95.3% | 12.17 | 24.82 |
| 60 | 50–90 | +2.27 | 6.41 | 70.7% (2.6 pp) | 94.3% | 12.25 | 23.65 |
| 60 | 30–120 | +2.13 | 7.32 | 67.0% (2.7 pp) | 93.3% | 13.23 | 26.89 |
| 65 | 50–90 | +1.72 | 6.67 | 64.3% (2.8 pp) | 95.0% | 13.10 | 24.95 |
| 65 | 30–120 | +2.05 | 7.61 | 62.7% (2.8 pp) | 93.3% | 14.21 | 28.86 |
| 70 | 50–90 | +1.32 | 6.25 | 68.0% (2.7 pp) | 97.3% | 13.33 | 24.85 |
| 70 | 30–120 | +2.24 | 7.76 | 66.0% (2.7 pp) | 94.3% | 15.26 | 30.82 |
| 75 | 50–90 | +0.94 | 5.92 | 67.7% (2.7 pp) | 98.0% | 12.66 | 22.98 |
| 75 | 30–120 | +3.30 | 9.36 | 60.7% (2.8 pp) | 91.3% | 16.50 | 32.77 |
| 80 | 50–90 | −1.46 | 5.24 | 79.3% (2.3 pp) | 97.3% | 11.91 | 21.26 |
| 80 | 30–120 | +1.64 | 8.40 | 69.7% (2.7 pp) | 95.0% | 16.74 | 33.07 |
| 85 | 50–90 | −3.17 | 4.91 | 80.3% (2.3 pp) | 98.7% | 10.52 | 18.55 |
| 85 | 30–120 | +1.99 | 8.20 | 68.7% (2.7 pp) | 97.0% | 17.54 | 34.09 |

“pp” means percentage points. Values come from the workflow output; they are not observational H0 measurements.

## Interpretation

1. **The broader prior gives wider intervals in every truth cell**, as expected when more posterior support is available. That alone does not guarantee better fixed-truth coverage.
2. **The narrower prior shifts the posterior median downward at the high-truth end.** At H0=85, mean signed median error is −3.17 under [50,90], versus +1.99 under [30,120]. This is evidence that prior support can materially affect posterior summaries in this setup.
3. **Coverage varies by truth and prior.** Under the broad prior, empirical 68% coverage ranges from 60.7% to 69.7%; with n=300, the per-cell MCSE is around 2.6–2.8 percentage points. Some cells are close to nominal 68%, while others are lower by more than one MCSE.
4. **The effect is not just a simple boundary-clipping story.** Some 68% intervals do not touch the narrow prior boundary, yet the posterior median and coverage still differ between prior choices. Nonlinearity of d=c z/H0, a prior uniform in H0, finite noise, and interval construction all affect repeated-sampling behavior.
5. A uniform prior in H0 is not equivalent to a uniform prior in distance or in another parameterization. Bayesian credible intervals are not guaranteed to have exactly nominal frequentist coverage at each fixed truth.
6. This experiment uses one fixed redshift and a simple Gaussian distance model. It is not a realistic dark-siren analysis, includes no galaxy catalogue uncertainty, and does not use public gravitational-wave observations.

## Conclusion

Prior support and parameterization can materially influence posterior summaries in this simplified model. However, this experiment does **not** explain every coverage anomaly in the previous host-mixture simulations, and it does not establish that the implementation is calibrated for real gravitational-wave cosmology.

## Next step: simulation-based calibration (SBC)

Run a prior-predictive rank diagnostic, drawing H0 from the same prior used for inference, generating distances from the matched likelihood, and checking the ranks of the injected values against their posterior distributions. This tests the complete Bayesian computation under the assumed prior-predictive model, complementing fixed-truth frequentist coverage. It should be conducted separately for each prior support and should include a test of posterior grid resolution.
