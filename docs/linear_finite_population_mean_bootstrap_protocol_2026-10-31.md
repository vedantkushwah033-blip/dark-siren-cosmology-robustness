# Linear finite-population mean benchmark: protocol and claim boundary

**Status:** implementation and unit-test layer added; no repeated-sampling coverage result is claimed. This is a calibration benchmark for the resampling primitive, not a result about the nonlinear dark-siren posterior.

## 1. Why this benchmark is separate

The dark-siren experiment's target is a nonlinear posterior-mean functional involving a fitted inclusion model, inverse-probability weights, a weighted KDE, and a selection-aware likelihood. A finite-population-mean bootstrap theorem does not automatically apply to that functional.

Before comparing bootstrap methods on the nonlinear target, the project now includes a transparent linear benchmark whose estimand is the ordinary mean of a fixed finite population. This makes it possible to test the pseudo-population construction, reapplication of SRSWOR, studentization, finite-population correction, and interval inversion in a setting with a known finite-frame target.

## 2. Procedure implemented

For a fixed finite population of size N, let a simple random sample without replacement of size n be observed, with n < N. Let y-bar be the sample mean and s-squared the sample variance (denominator n-1). The estimated design variance of the sample mean is

V-hat = (1 - n/N) s-squared / n.

For each bootstrap replicate:

1. Generate a pseudo-population of size N by drawing multinomial counts with equal cell probabilities 1/n over the n observed records.
2. Expand those counts into a pseudo-population.
3. Draw n positions without replacement from that pseudo-population.
4. Calculate the bootstrap sample mean and the pseudo-population mean.
5. Calculate the bootstrap SRSWOR standard error using the bootstrap sample variance and the same finite-population factor.
6. Store the studentized statistic (bootstrap sample mean - pseudo-population mean) / bootstrap standard error.

The confidence interval inverts the empirical bootstrap-t quantiles around the observed sample mean, using the observed design-based standard error. Replicates with zero or non-finite estimated variance are explicitly counted as failed; they are not silently substituted.

Implementation: `src/dark_siren_robustness/finite_population_mean_bootstrap.py`.

Unit tests check deterministic replay, attempt accounting, studentized-statistic reconstruction, finite-population variance arithmetic, interval inversion, and input validation. Passing these tests establishes software-contract behavior only, not nominal coverage.

## 3. Theoretical reference and limits

Wang, Peng and Kim (2022), *Bootstrap Inference for the Finite Population Mean under Complex Sampling Designs*, Journal of the Royal Statistical Society: Series B 84(4), 1150–1174, DOI: 10.1111/rssb.12506, develop a pseudo-population bootstrap that regenerates a finite population, reapplies the sampling design, and studentizes. Their SRS theory is for a finite-population mean and relies on asymptotic regularity conditions, including a sample fraction bounded away from one, non-degenerate finite-population variance, a moment condition, and a non-lattice limiting distribution.

This repository's implementation is an explicit computational benchmark inspired by that SRS construction. It is not represented as a new theorem or a complete reproduction of every theoretical condition. The empirical coverage study remains to be run and independently checked.

Most importantly, this finite-mean result does **not** prove validity for the nonlinear dark-siren posterior mean. The project will report the linear benchmark and nonlinear posterior-functional experiment as distinct analyses. A successful linear benchmark is a prerequisite sanity check, not a bridge that licenses theoretical claims about the nonlinear estimator.

## 4. Required next experiment

A separate fixed-frame simulation runner should:

- generate several deterministic finite populations with varied skewness, tail behavior, and sample fractions;
- draw repeated exact SRSWOR samples for each frame;
- compare the studentized pseudo-population interval with a conventional Wald interval and a clearly labelled ordinary with-replacement bootstrap diagnostic;
- evaluate coverage of the exact finite-frame mean, interval width, failure rate, and Monte Carlo uncertainty;
- report scenario-level results and avoid treating repeated samples from one frame as independent finite populations;
- preserve seeds, sampled indices, pseudo-population draws or sufficient replay data, and an independent reconstruction audit.

The pilot should determine replication counts and Monte Carlo precision before any confirmatory claim. If coverage is poor in any prespecified regime, that result must be reported rather than tuned away.

## 5. Citation

Wang, Z., Peng, L. & Kim, J. K. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *Journal of the Royal Statistical Society: Series B (Statistical Methodology)*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506
