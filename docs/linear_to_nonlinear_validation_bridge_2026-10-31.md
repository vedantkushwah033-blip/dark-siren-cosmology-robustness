# Validation bridge: finite-population means to dark-siren posterior functionals

**Status:** design note. The v1 linear pilot completed and its outer artifact audit passed on the recorded source commit. The audit boundary and required v2 inner-attempt accounting are documented separately. This note does not claim confirmatory evidence.

## Why a validation bridge is needed

The project studies uncertainty under incomplete galaxy catalogues. Its downstream quantity is not an ordinary sample mean: it is a nonlinear posterior functional assembled from distance likelihoods, candidate-host contributions, catalogue weights, smoothing, and a selection-normalization term. A method that behaves well for a linear finite-population mean is not thereby calibrated for this posterior.

The validation programme therefore separates increasingly demanding claims rather than treating one successful simulation as universal validation.

## Layer 1 — linear design-based benchmark

Target: \(\bar Y_N=N^{-1}\sum_{i=1}^{N}Y_i\), the exact mean of a realized finite frame sampled by SRSWOR.

The benchmark compares a studentized pseudo-population SRSWOR bootstrap with a design-based Wald interval and an ordinary with-replacement percentile bootstrap. It is a diagnostic test of the finite-population resampling implementation and the audit pipeline. Its literature basis is Wang, Peng & Kim (2022), DOI 10.1111/rssb.12506.

Required evidence before interpreting it:
- exact-head CI passes the dedicated tests;
- the pilot completes with all attempts retained;
- the independent auditor reconstructs frames, samples, seeds, intervals, failures, and summaries;
- coverage is reported with its denominator, Monte Carlo uncertainty, interval width, and failure rate;
- multiple independently generated frames are used before making regime-level generalizations.

## Layer 2 — finite-validation bootstrap for catalogue quantities

The next target is uncertainty induced by finite validation samples and estimated catalogue-inclusion probabilities. This is a different estimand and a different resampling problem. It requires:
- an explicit validation-sampling design and known inclusion probabilities;
- separation of training/validation information from posterior evaluation to prevent leakage;
- checks for positivity, extreme weights, effective sample size, and weight concentration;
- paired repeated-sampling comparisons on identical simulated catalogues;
- an independent replay of the entire resampling and posterior pipeline.

The linear benchmark does not establish any of these properties.

## Layer 3 — nonlinear dark-siren posterior calibration

The final target is a declared posterior functional, for example the posterior mean of H0 under a specified synthetic data-generating process. Calibration must be evaluated against known injected truths over repeated independently generated datasets.

At minimum, report:
- bias and RMSE of the posterior mean;
- empirical coverage and width of posterior credible intervals;
- simulation-based calibration diagnostics where their assumptions apply;
- numerical integration error and posterior-grid resolution sensitivity;
- selection-normalization stability;
- sensitivity to catalogue incompleteness and host-population misspecification;
- failure rates, effective sample size, maximum normalized weight, and reproducibility metadata.

A credible interval's nominal posterior probability is not the same quantity as its repeated-sampling coverage. Both may be informative, but they must not be conflated.

## Claim ladder for the manuscript

| Evidence achieved | Defensible claim | Claim not justified |
|---|---|---|
| Unit tests pass | The implementation satisfies tested software invariants | The interval procedure is calibrated |
| Linear pilot independently audited | Recorded pilot artifacts match the prespecified linear design and calculations | General finite-population coverage or nonlinear posterior validity |
| Multi-frame linear study with precision targets | Operating characteristics for the declared linear regimes and design | Validity for galaxy-catalogue weights or dark-siren posteriors |
| Finite-validation repeated-sampling study | Operating characteristics for the declared validation design and estimands | Real-observatory performance |
| Nonlinear synthetic SBC and repeated-sampling study | Calibration under the specified synthetic generative model | Calibration on real LVK events |
| External data benchmark with versioned inputs and independent review | Performance on that exact external benchmark | Universal detector realism or a new cosmological measurement |

## Reporting rules

1. Keep pilot and confirmatory datasets separate.
2. Freeze confirmatory estimands, regimes, sample sizes, stopping rules, and primary contrasts before execution.
3. Treat frames as a higher-level sampling unit; do not use the total number of within-frame samples as if they were independent frame draws.
4. Preserve failed attempts and report their frequency.
5. Distinguish numerical reproducibility from scientific validity.
6. Do not describe a synthetic H0 estimate as a measurement of the Universe.
7. Do not describe a method as novel when the contribution is an implementation, audit, or application of an established method.

## Current stop/go gate

The v1 PR-triggered workflow passed on Python 3.10, 3.11, and 3.12, and its Python 3.12 artifact audit reported PASS for three frames, 1,800 outer method attempts, and 18 summary cells. The v1 artifact does not serialize an independently checkable inner-bootstrap failure histogram; see the v2 remediation protocol. A later workflow on commit `9b26b23c29e26ce5bcdd1f729749072fd5ac5329` passed Python 3.10 and 3.11 while Python 3.12 continued through the separate catalogue-thinning calibration. The nonlinear dark-siren claims remain bounded by their own synthetic-model assumptions and the separate external-data benchmark gate.
