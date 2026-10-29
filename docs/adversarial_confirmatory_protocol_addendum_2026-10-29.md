# Adversarial protocol addendum: selection conditioning, boundary strata, and confirmatory design

**Status:** preregistration-quality design correction; this document is not a report of completed runs.  
**Applies to:** synthetic dark-siren/cosmology robustness experiments on `research/audit-posterior-grid-integrity`.  
**Scope limit:** the low-redshift linear relation and toy selection model do not validate an LVK analysis or yield a real-data estimate of H0.

## 1. Correct conditioning statement

Let the *pre-detection* event-data density be
\[
q(d\mid H_0,f)=(1-f)q_{\rm cat}(d\mid H_0)+f q_{\rm miss}(d\mid H_0),
\]
where each component integrates to one over the same pre-selection data space and (f) is explicitly the pre-selection missing-host fraction. If detection probability conditional on data (d) is \(\epsilon_{\rm true}(d)\), define
\[
\alpha_k(H_0)=\int_{\mathcal D}\epsilon_{\rm true}(d)q_k(d\mid H_0)\,dd,\qquad
\beta(H_0,f)=(1-f)\alpha_{\rm cat}(H_0)+f\alpha_{\rm miss}(H_0).
\]
For one detected datum, under a model with assumed efficiency \(\epsilon_{\rm model}\), the conditional density is
\[
p(d\mid D=1,H_0,f)=
\frac{\epsilon_{\rm model}(d)q(d\mid H_0,f)}
{\int_{\mathcal D}\epsilon_{\rm model}(u)q(u\mid H_0,f)\,du},
\]
for data in the detection domain. The numerator and denominator must use the same assumed selection law for a correctly specified likelihood. A deliberately mismatched likelihood must use its assumed law consistently in both places and be labelled misspecified. For a hard cut, \(\epsilon(d)=1\) for accepted observations and zero otherwise; do not multiply by an extra nontrivial efficiency factor inside the accepted region.

For (N) conditionally independent detected events and fixed observed count (N), multiply the *per-event conditional density* (N) times; the population selection normalization therefore appears as \(\beta(H_0,f)^{-N}\), not merely once. This statement assumes conditioning on (N); a rate-based Poisson point-process likelihood has a different count/rate term and must be derived separately. Do not mix these formulations.

The detected-event mixture fraction is generally not (f). Its missing-host component is
\[
f_{\rm miss,det}(H_0)=\frac{f\alpha_{\rm miss}(H_0)}{(1-f)\alpha_{\rm cat}(H_0)+f\alpha_{\rm miss}(H_0)}.
\]
Tests should check this identity numerically and verify that the mixture weights sum to one after conditioning.

## 2. Boundary effects: what is and is not established

A posterior truncated to prior support \([50,90]\) can be pulled inward when the likelihood would otherwise place substantial mass beyond a boundary. This is a predictable support effect, but it does **not** by itself prove a particular coverage figure such as 57%, nor imply that all endpoint truths must have the same coverage. Those numbers require a reproducible artifact and a fully specified simulation design.

Report truth strata separately: exact endpoints (50 and 90), near-boundary values (55 and 85), and interior (e.g. 60–80). Include signed bias, RMSE, 68% and 95% coverage with Wilson intervals, mean interval width, posterior mass within one and two grid cells of each prior boundary, and the fraction of posterior mass that would lie outside the analysis support under a wider-support diagnostic where available. Do not pool endpoint and interior results into one headline calibration number. Distinguish prior-support sensitivity from grid truncation: widening a prior is not the same as extending a numerical grid while keeping the prior fixed.

## 3. Design audit for the proposed factorial study

The proposed levels (H_0^{\rm true}\in\{50,55,70,85,90\}), (f_{\rm true}\in\{0.1,0.3,0.5\}), (f_{\rm model}\in\{0.1,0.3,0.5\}), and three cutoff offsets imply 5 × 3 × 3 × 3 = **135 cells**. At 500 replicates per cell, this is 67,500 replicate datasets and 20.25 million event-level observations at 300 events per replicate, before failed-fit reruns and diagnostic arms. Do a small runtime pilot and estimate Monte Carlo precision before launching the full grid. If computationally infeasible, use a predeclared staged design: first validate all correctly specified cells and endpoint/interior contrasts; then run the full crossed misspecification grid. Do not select scenarios based on favourable outcomes.

The protocol must state the true cutoff, the exact selection function, whether selection acts on latent distance or noisy observed distance, whether the 300 events are generated conditional on detection, and how events are sampled from the catalogue/missing-host mixture. A numeric offset alone is not reproducible if the baseline cutoff and selection law are unspecified.

For nominal coverage (p=0.68) and (R=500), the binomial standard error is
\[
\sqrt{p(1-p)/R}=0.02086,
\]
or 2.09 percentage points. A rough 95% Monte Carlo half-width is about 4.1 percentage points near the nominal rate. Report the covered count, denominator, Wilson interval, and exact scenario; do not describe 500 replicates as guaranteeing ±2 percentage-point 95% precision. For a paired bias contrast, calculate its Monte Carlo standard error from the replicate-level paired differences, not from the coverage formula.

## 4. Predeclared estimands and reporting

Primary outcomes per cell:
1. Mean signed posterior-median error (R^{-1}\sum_r(\hat H_{0,r}-H_{0,\rm true})), with Monte Carlo standard error and 95% interval.
2. RMSE of posterior medians.
3. 68% and 95% interval coverage, covered counts, and Wilson intervals.
4. Mean interval widths and fit-failure/non-finite rates.
5. Boundary-stratified results and posterior boundary-mass diagnostics.

Use paired replicate seeds across inference arms when they share the same generated dataset. Preserve the replicate-level raw table and full configuration. If the same seed is reused across different conditions, label those conditions as potentially dependent; do not call them independent replications. Define the bias acceptance rule numerically using the replicate-level standard error and report the interval regardless of pass/fail. The coverage acceptance band is a decision threshold, not a proof of calibration.

## 5. Adversarial controls required before a scientific claim

- **Selection negative control:** omit the H0-dependent normalization on purpose; the audit must detect the expected calibration/bias degradation under a condition where selection matters.
- **Cutoff mismatch:** generate under a declared true cutoff and infer under predeclared lower/equal/higher assumed cutoffs. Confirm that the equal-cutoff arm is the correctly specified control.
- **Mixture control:** test (f=0), (f=1), and intermediate values; reject NaN, negative, or non-normalized component weights.
- **Boundary control:** compare exact endpoints, near-boundary, and interior truths under fixed and widened prior supports.
- **Numerical control:** independently recompute posterior normalization, medians, quantiles, and interval coverage from saved grid likelihoods; repeat at finer grids and report maximum endpoint/summary changes.
- **Integrity control:** fail on missing/duplicate replicate-arm rows, mismatched manifests, altered posterior arrays, or disagreement between saved summaries and independent reconstruction.

## 6. Claim gates

- **Software gate:** tests pass on the exact commit and supported Python versions.
- **Numerical gate:** grid convergence and independent posterior reconstruction pass on the exact artifacts.
- **Synthetic statistical gate:** confirmatory outcomes and Monte Carlo uncertainty are reported, including all preregistered strata and failed fits.
- **Observational gate:** remains **NO-GO** until a documented, version-pinned real-data pipeline, compatible inputs, provenance, and independent end-to-end reproduction exist.
- **Publication gate:** do not call the work publication-ready solely because CI passes or a synthetic experiment yields favourable results.

No value in this document should be reported as an observed project result unless the corresponding immutable run artifact exists and has been independently checked.
