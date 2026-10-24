# Fixed-truth coverage map: results and interpretation

**Status:** completed synthetic experiment; workflow passed.  
**Workflow:** [Toy fixed-truth coverage map](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40072545551)  
**Commit tested:** `f96395906528d47c16d5d1c4d2e254e928b5d292`

## Design

The experiment evaluates pointwise repeated-sampling behavior at seven fixed true values of (H_0): 52, 55, 60, 70, 80, 85, and 88. It uses 500 detected synthetic observations per truth (3,500 total), a uniform prior on [50, 90], a 4,001-point posterior grid, host redshifts 0.02 and 0.04 with weights 0.7 and 0.3, Gaussian distance noise of 20 Mpc, and the toy detection condition (0<d_{obs}<180) Mpc. Corrected inference divides the host-mixture likelihood by the model's detection efficiency; uncorrected inference deliberately omits that normalization.

## Results

Each cell is the empirical fraction of the 500 simulations whose equal-tailed credible interval contains the fixed truth. Values are Monte Carlo estimates, not exact coverage probabilities.

| True (H_0) | Bias, uncorrected | Bias, corrected | RMSE, uncorrected | RMSE, corrected | 68% coverage, uncorrected | 68% coverage, corrected | 95% coverage, uncorrected | 95% coverage, corrected |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 52 | 14.648 | 12.530 | 15.794 | 13.843 | 0.000 | 0.000 | 0.872 | 0.902 |
| 55 | 11.133 | 9.159 | 12.353 | 10.606 | 0.580 | 0.654 | 0.950 | 0.964 |
| 60 | 8.035 | 6.404 | 9.835 | 8.678 | 0.708 | 0.770 | 0.986 | 0.990 |
| 70 | 1.820 | 0.509 | 6.172 | 6.269 | 0.952 | 0.962 | 1.000 | 1.000 |
| 80 | -6.436 | -7.851 | 8.727 | 10.156 | 0.854 | 0.812 | 1.000 | 1.000 |
| 85 | -10.954 | -12.375 | 12.604 | 14.152 | 0.550 | 0.454 | 1.000 | 0.996 |
| 88 | -13.592 | -14.983 | 15.191 | 16.745 | 0.074 | 0.068 | 0.978 | 0.904 |

Bias and RMSE are in the model's (H_0) units. The workflow output also reports Wilson intervals for empirical coverage and paired changes in absolute error.

## Interpretation

1. **Performance is strongly truth-dependent.** The posterior is biased upward at low true values and downward at high true values. The fixed prior [50, 90] and the toy host-distance mixture constrain the posterior; this pattern should not be generalized to real dark-siren analyses.
2. **Selection normalization is not uniformly superior under every fixed truth and metric.** It reduces RMSE at 52, 55, and 60, but increases RMSE at 70, 80, 85, and 88 in this particular simulation. At (H_0=70), the RMSE difference is small relative to the reported simulation variability. The correct conditional likelihood is not guaranteed to minimize fixed-truth RMSE for every parameter value.
3. **Pointwise 68% coverage is far from nominal near the prior edges.** For (H_0=52), neither method's 68% interval covered the truth in any of the 500 replicates. For (H_0=88), the rates were 7.4% and 6.8%. The truth is within the prior, but close to its boundary, and the posterior is systematically displaced. This is a substantive property of the toy setup to investigate, not a reason to tune the generator.
4. **95% coverage is also imperfect near the edges.** In particular, corrected coverage at (H_0=88) was 0.904. Its Wilson interval from the run was approximately [0.875, 0.927], below 0.95, so Monte Carlo noise alone is unlikely to explain the entire discrepancy for this experiment.
5. **Bayesian credible interval mass and pointwise frequentist coverage are different calibration concepts.** Prior-predictive SBC can be well calibrated even while fixed-truth repeated-sampling coverage varies across (H_0). Both diagnostics are useful and answer different questions.

## What this experiment does not establish

It does not establish that selection normalization is generally harmful or helpful, that a realistic gravitational-wave analysis would show the same bias, or that the model is suitable for observational inference. The model has a finite host mixture, a simplified Gaussian distance likelihood, a finite prior, and a toy distance threshold. Follow-up work should isolate prior-boundary effects, compare with analytic/single-host controls, report Monte Carlo uncertainty for all comparisons, and preserve the same generative model when evaluating a correctly specified likelihood.
