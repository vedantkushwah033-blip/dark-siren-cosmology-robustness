# Fixed-truth Bayesian credible coverage versus an inverted Gaussian confidence interval

**Status:** independent secondary analysis of an archived synthetic artifact; descriptive, not an astrophysical validation.

## Why this check was added

The 32,400-record distance-uncertainty experiment reports fixed-truth coverage of central Bayesian credible intervals. That quantity is useful, but a Bayesian 68% credible interval is not generally required to contain a fixed frequentist truth in 68% of repeated samples. Treating a deviation from 68% as a model or implementation failure without a matched reference interval would therefore be too strong.

This note adds a reference calculation for Arm A, where the realized Gaussian distance uncertainty is held fixed in the likelihood. It asks whether the same simulated observations support a conventional Gaussian confidence interval with the expected repeated-sampling coverage.

## Input and provenance

- Source: archived **distance_uncertainty_sensitivity_full/replicates.csv** from the numerical-diagnostics Python 3.12 artifact, GitHub Actions run 40630753593.
- Raw CSV SHA-256: **62af121350d62b141f0f93917258e1967cd51b611234e6c4c9839f3088e2a4da**.
- Records: 32,400 total; Arm A contributes 10,800 records.
- Arm A model: \(d_{\mathrm{obs}}=cz/H_{0,\mathrm{true}}+\sigma_A\epsilon\), with \(\epsilon\sim N(0,1)\), \(\sigma_A=fcz/H_{0,\mathrm{true}}\).
- Design: \(H_0\in\{50,55,\ldots,90\}\), four fractional-noise settings \(f\in\{0.05,0.10,0.15,0.20\}\), 300 replicates per cell.
- This is a synthetic one-known-host low-redshift model, not LVK data. Arm A's sigma is generated using the true \(H_0\), so it is an oracle-style conditional benchmark rather than an operational recipe available to an analyst who does not know \(H_0\).

The source CSV and archived posterior intervals were read directly. Coverage was reconstructed independently from each row's observation, realized sigma, redshift, truth, and stored posterior interval. No production inference code or historical artifact was changed for this analysis.

## Reference interval

For a central Gaussian distance interval with nominal level \(1-\alpha\), let \(q=\Phi^{-1}(1-\alpha/2)\). Under the Arm A data-generating model,

\[
[d_{\mathrm{obs}}-q\sigma_A,\ d_{\mathrm{obs}}+q\sigma_A]
\]

has the nominal repeated-sampling coverage for the Gaussian distance mean. Because \(d=cz/H_0\) is monotone decreasing for positive \(H_0\), invert the interval to obtain

\[
H_0\in\left[\frac{cz}{d_{\mathrm{obs}}+q\sigma_A},\ 
\frac{cz}{d_{\mathrm{obs}}-q\sigma_A}\right].
\]

When the lower distance endpoint is non-positive, the upper \(H_0\) endpoint is treated as \(+\infty\), rather than discarding the replicate. For this design, this occurred in four of the 10,800 Arm A records for the 95% interval and none for the 68% interval. This interval is a calibration reference conditional on the known realized sigma; it is not a claim that this oracle interval is operationally available.

## Results

Coverage below is from Arm A only. Wilson intervals are 95% binomial intervals for Monte Carlo uncertainty.

| Truth set / interval | Bayesian equal-tail posterior interval | Inverted Gaussian distance interval |
|---|---:|---:|
| Interior truths \(55\le H_0\le85\), nominal 68% | 6,156 / 8,400 = **73.29%** (Wilson 95%: 72.33–74.22%) | 5,680 / 8,400 = **67.62%** (Wilson 95%: 66.61–68.61%) |
| Interior truths \(55\le H_0\le85\), nominal 95% | 8,158 / 8,400 = **97.12%** (Wilson 95%: 96.74–97.46%) | 7,990 / 8,400 = **95.12%** (Wilson 95%: 94.64–95.56%) |
| Central truth \(H_0=70\), nominal 68% | 858 / 1,200 = **71.50%** (Wilson 95%: 68.88–73.98%) | 796 / 1,200 = **66.33%** |
| Central truth \(H_0=70\), nominal 95% | 1,178 / 1,200 = **98.17%** (Wilson 95%: 97.24–98.79%) | 1,146 / 1,200 = **95.50%** |

The inverted-interval coverage matches the Gaussian noise-generation rule within Monte Carlo uncertainty, including when pooling the seven interior truths. The Bayesian equal-tail intervals are more conservative in this particular setup: their interior-truth coverage is above their nominal levels, including at the central truth \(H_0=70\). This is consistent with the fact that posterior credible coverage at fixed parameter values is not generally identical to frequentist confidence coverage; the nonlinear \(d\leftrightarrow H_0\) transformation and uniform \(H_0\) prior are part of this setup.

The exact prior-endpoint truths \(H_0=50\) and \(H_0=90\) remain a separate issue: the stored central equal-tail posterior intervals do not include these exact endpoints, so they have zero coverage there by construction. Do not pool those boundary cells into a nominal-coverage claim.

## Interpretation and limits

1. **Do not label 73.3% Bayesian 68% coverage as an implementation bug on its own.** The matched Gaussian reference interval behaves as expected, which makes prior/parameterization effects a credible explanation for the difference.
2. **Do not conclude that the posterior is universally calibrated.** This reference check is conditional on Arm A's known realized sigma and on the stated synthetic likelihood. It does not test catalogue incompleteness, host-population misspecification, detector selection, or realistic GW likelihoods.
3. **Do not use the oracle interval as a practical estimator.** Its sigma depends on the true \(H_0\). It is a controlled reference for the repeated-sampling experiment only.
4. The next useful study is a prespecified prior-sensitivity comparison (for example, uniform \(H_0\) versus a clearly defined alternative prior, with identical generated observations and exact same interval estimands), reported by truth and noise fraction. It should be reconciled with the existing unintegrated prior-support experiment in the nine-patch research bundle before implementation, to avoid duplicating work.
5. For a publishable claim, this result needs a scripted, tested reproduction step in CI and a systematic literature comparison. It is not evidence of novelty or astrophysical validity by itself.

## Reproduction sketch

The core calculation uses the stored columns **observed_distance_mpc**, **generated_sigma_mpc**, **host_redshift**, **true_h0**, and **standard_normal**. For each interval level, compute \(q\), form the distance endpoints, invert the positive endpoints to \(H_0\), allow an infinite upper endpoint if needed, and count coverage. Compare against the archived **lower_68**, **upper_68**, **lower_95**, and **upper_95** columns. The raw artifact hash above pins the input for independent replication.
