# Simulation-based calibration (SBC) rank audit

**Status:** rerun completed after correcting the prior-support CDF calculation.  
**Latest successful workflow:** [SBC rank audit, run 40068188858](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40068188858)  
**Successful workflow:** [SBC rank audit, run 40005653560](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40005653560)  
**Cross-version test suite:** [run 40005653507](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40005653507)  
**Cross-version tests:** [run 40005653555](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40005653555)  
**Artifact:** `sbc-rank-audit-results` (contains per-run ranks and summary diagnostics).

## Purpose and design

Simulation-based calibration checks the Bayesian computation under a prior-predictive model. For each of two prior supports, H0 is sampled uniformly from the same prior used for inference; a distance observation is generated from the likelihood; the posterior is calculated; and the posterior CDF evaluated at the injected H0 is recorded.

If the prior-predictive generator and inference calculation match, these CDF ranks should be approximately Uniform(0, 1), apart from Monte Carlo and finite-grid effects.

The design used 1,000 replicates per prior (2,000 posterior fits total), fixed redshift z=0.03, and a fixed Gaussian distance uncertainty
$$
\sigma = 0.10\,\frac{c(0.03)}{70}.
$$
The fixed uncertainty is important: the inference likelihood treats `gw_sigma` as fixed, so the simulation must not secretly vary it as a function of the sampled H0 unless that dependence is also included in the likelihood.

## Results

| Prior support (km s⁻¹ Mpc⁻¹) | Replicates | Mean CDF rank | Rank variance | KS statistic | KS p-value | 10-bin chi-square p-value |
|---|---:|---:|---:|---:|---:|---:|
| 50–90 | 1,000 | 0.4989 | 0.0793 | 0.0214 | 0.7421 | 0.3361 |
| 30–120 | 1,000 | 0.4934 | 0.0820 | 0.0239 | 0.6065 | 0.9179 |

At the conventional 5% level, neither KS test nor the decile chi-square test rejects uniformity for these runs. Mean ranks are close to 0.5 and the empirical rank variances are close to the Uniform(0,1) variance of 1/12 ≈ 0.0833.

These are diagnostics, not proof of exact calibration. The tests have finite power; a non-significant p-value does not establish that the implementation is correct for all conditions.

## Important correction to the first attempt

The first SBC run used a measurement uncertainty proportional to the true distance, while the inference likelihood treated that uncertainty as fixed and independent of candidate H0. This was a mismatch between the data-generating process and the likelihood, and the first run showed low KS p-values. That result should **not** be interpreted as evidence of a bug in the posterior routine.

The generator was corrected to use a fixed, H0-independent uncertainty, matching the likelihood being fitted. The corrected run produced the results reported above. The original output remains part of the workflow history and is superseded for the SBC question.


## Follow-up correction to the SBC rank calculation

A later code review found that the SBC helper's CDF calculation integrated over the full computational grid, including zero-density points outside the prior support. This could introduce artificial trapezoidal boundary mass, even though the posterior routine itself was normalized on the prior-supported grid. The helper now restricts its CDF integration to the nonzero posterior support, and a regression test checks a constant-density example with analytically known CDF value. The table above has been refreshed from the rerun linked at the top of this report. This correction reinforces why numerical calibration should be audited independently of whether the original tests passed.

## What this validates—and what it does not

- It provides evidence that the posterior CDF computation behaves consistently with the specified one-host, Gaussian-distance, uniform-H0 prior-predictive model in these tested settings.
- It does not show that fixed-truth frequentist coverage must equal nominal coverage at every H0 value.
- It does not validate the unknown-host mixture model, galaxy-catalogue incompleteness, selection effects, a full cosmological distance relation, or inference from real gravitational-wave events.
- The fixed-truth prior-support experiment and prior-predictive SBC answer different questions; both should remain reported.

## Next step

Apply SBC to a deliberately small, analytically inspectable two-host mixture. Then validate the same posterior ranks under multiple grid resolutions. Only after these tests pass should the study extend the calibration checks to catalogue incompleteness and host-weight misspecification.
