# Independent adjudication of the proposed 48-hour research sprint

**Date:** 5 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Scope:** Mathematical and source-level review of the user-supplied execution package. This document records an adjudication, not an experiment result. No confirmatory simulations, test runs, or manuscript-completion claims are made here.

## Executive decision

The proposed sprint is directionally useful, but its final verdict — **“confirmatory-experiment ready / manuscript draft ready” — is not supported** by the evidence inspected. Several proposed changes would alter the statistical estimator or make assumptions that have not been established. They must not be merged into the primary analysis simply because they appear to stabilize computation.

**Decision: NO-GO for the R=2,000 confirmatory campaign and NO-GO for journal submission at this stage.** A bounded, non-confirmatory paired pilot may proceed after exact-head CI and artifact checks. The primary scientific objective remains to establish when finite-validation inclusion estimation improves or degrades the reconstructed host population and downstream synthetic H0 calibration.

## 1. Ridge regularization is not a numerical patch to the existing V2 estimator

The supplied package proposes replacing the current objective by

\[
\ell_\lambda(\beta)=\ell(\beta)+\frac{\lambda}{2}\|\beta\|_2^2.
\]

This is a legitimate alternative estimator for a declared \(\lambda>0\), but it is **not** the same estimator as the current V2. The current `estimated_inclusion_v2.py` implements ordinary logistic negative log-likelihood under a finite coefficient box, with analytic gradient/Hessian, projected-KKT stopping and damped/backtracked steps. It does not add an L2 penalty. The damping term used to stabilize a Newton linear solve is algorithmic; it is not a ridge penalty in the statistical objective.

The distinction matters:

- A finite box produces a constrained MLE (when the optimization is solved), not an unconstrained MLE under separation.
- Ridge changes the objective and shrinks coefficients. If the intercept is penalized, it also shrinks the baseline inclusion level. A common alternative is a penalty matrix that leaves the intercept unpenalized, but that too is a distinct estimator and must be stated.
- A ridge penalty can produce finite estimates under separation, but it does not establish positivity, correct functional form, calibration, or acceptable inverse-weight variance.
- Selecting \(\lambda\) using the same validation data without accounting for tuning adds another layer of estimation uncertainty.

**Recommended design:** retain V1 and the existing box-constrained V2 as the primary solver comparison. If ridge is scientifically useful, add a separately named sensitivity arm with a pre-specified penalty matrix and fixed/tuned-by-nested-validation rule. Persist \(\lambda\), penalty matrix, convergence diagnostics, and solver identity in every artifact. Do not relabel V2 as ridge.

## 2. Armijo descent does not resolve statistical non-existence

For Bernoulli logistic regression,

\[
\ell(\beta)=\sum_i\{\log(1+e^{x_i^T\beta})-y_i x_i^T\beta\},\quad
g=X^T(\pi-y),\quad H=X^TWX,
\]
where \(W_{ii}=\pi_i(1-\pi_i)\).

The Hessian is positive semidefinite and is positive definite for finite coefficients when the design matrix has full column rank. Under complete or quasi-complete separation, the unpenalized likelihood can approach its supremum along a diverging coefficient sequence; the negative log-likelihood approaches a finite infimum. A monotone line search can control individual steps but cannot create a finite unconstrained MLE where none exists.

The existing V2 therefore needs to be described as a **box-constrained ordinary logistic MLE solver**. Its diagnostics should distinguish at least:
1. interior KKT convergence;
2. boundary KKT convergence;
3. rank-deficient/non-identifiable coefficient directions;
4. detected complete or quasi-complete separation;
5. numerical non-convergence / line-search failure;
6. invalid or non-evaluable input.

A boundary solution is not automatically an optimizer failure: it may be the correct constrained optimum. Conversely, a solver success flag is not evidence that the unconstrained MLE exists. Active-set gradient signs and projected-gradient residuals should both be checked. A formal separation detector (for example, a documented linear-programming test) is needed before assigning a separation label; large coefficients and extreme fitted probabilities alone are only warning signals.

## 3. The proposed FPC covariance is not generally justified

The package proposes

\[
\widehat\Sigma_\beta=(X_{val}^TWX_{val})^{-1}
\left(1-\frac{n_{val}}{N_{parent}}\right).
\]

The scalar finite-population correction (FPC) is not a universal multiplier for the covariance of a logistic-regression coefficient estimator. It may arise in particular design-based linearization calculations under simple random sampling without replacement, but the formula cannot be adopted without specifying:

- the finite population and which variables/outcomes are fixed in it;
- the validation sampling design (SRSWOR, stratified, clustered, unequal probability, or otherwise);
- whether the inclusion indicators are observed for the entire sampled frame;
- any validation sampling weights and nonresponse adjustments;
- the target of inference (finite-population coefficient, superpopulation parameter, or predictive inclusion function);
- the appropriate design-based variance estimator and its finite-sample behavior.

For a superpopulation Bernoulli logistic model, the inverse observed information is a model-based covariance approximation and does not automatically acquire this FPC. For a finite-population design, a design-aware sandwich/linearization or replicate-weight procedure must be derived for the actual sampling design. Multiplying the inverse Hessian by an FPC by analogy is not sufficient.

**Decision:** do not implement the proposed covariance formula until the validation sampling design and estimand are explicit and a derivation/test is available.

## 4. “Method C” requires a declared distribution over the inclusion parameters

The proposed mixture

\[
p(H_0\mid d,\mathcal V)\approx
\frac{1}{B}\sum_{b=1}^B p(H_0\mid d,\beta^{(b)})
\]

is mathematically a mixture of normalized conditional posteriors. It becomes a coherent Bayesian marginal posterior only if the draws come from a justified posterior \(p(\beta\mid\mathcal V)\) under a specified joint model and compatible prior/likelihood factorization. Draws from a frequentist Wald approximation or a bootstrap distribution instead define a different uncertainty-propagation procedure; they must not be labelled Bayesian marginalization.

Further issues with the proposed implementation:

- A normal approximation around a boundary-limited or separated logistic fit can place substantial mass outside the coefficient box and can be a poor approximation to the estimator's sampling distribution.
- The inverse Hessian may be ill-conditioned or singular under rank deficiency.
- A simple FPC multiplier does not repair those issues.
- \(B=50\) is not automatically adequate. The Monte Carlo error of posterior summaries and mixture quantiles must be assessed by increasing \(B\), using independent seeds, and reporting convergence.
- Averaging posterior densities on a common H0 grid is different from averaging posterior means or interval endpoints. Quantiles must be computed from the mixture density, not averaged across component quantiles.

**Safer near-term alternative:** retain the existing validation-resampling sensitivity analysis as a clearly labelled frequentist two-stage uncertainty diagnostic, conditional on the fixed catalogue and GW data. If a Bayesian nuisance-parameter marginalization is desired, specify a joint generative model for validation outcomes, parent inclusion, host population and GW data, then derive the posterior factorization before coding it. Do not call a bootstrap mixture a Bayesian posterior.

## 5. The proposed hard assertion f=0 is conditional, not universal

If the weighted catalogue density has already been normalized to estimate the **entire parent host population**, and \(f\) in a separate mixture denotes the fraction of hosts missing from the catalogue, then adding a second missing-host component can double-correct the same missing population. In that specific construction, a separate missing-host mixture is redundant and the model should not add it.

However, the statement “IPW always implies \(f=0\)” is too broad without defining exactly what the catalogue component and missing component represent. A mixture may represent a different physical or statistical decomposition, or the IPW estimator may target only a subpopulation. The correct invariant is semantic and generative:

> Do not combine a full-parent IPW reconstruction with a second term that reintroduces the same uncatalogued parent mass.

Source inspection of the estimated-inclusion pilot path found the reconstructed host density passed to the event integrals and selection normalization, and no separate hidden-host mixture fraction was evident on that inspected path. This is not a global clearance of every likelihood implementation in the repository. Before adding an assertion, trace every likelihood path and document the estimand represented by each component. An unconditional `assert f == 0` could conceal a model-definition error rather than fix it.

## 6. Confirmatory design: correct the estimands before choosing R

The proposed paired contrast

\[
D_r=|\widehat H_{0,cat}^{(r)}-H_{0,true}|
-|\widehat H_{0,IPW}^{(r)}-H_{0,true}|
\]

is a useful paired performance estimand: positive values favor the IPW arm on absolute posterior-mean error. It is not sufficient alone. At minimum report paired differences in:
- host-density error against the known synthetic parent density;
- posterior-mean error and bias;
- interval width and empirical coverage;
- weight concentration and positivity diagnostics;
- optimizer/separation/boundary outcomes;
- operational success and conditional-on-success performance.

### Truth stratum

If the generating truth is fixed at \(H_0=70\), it is not meaningful to call \(55\leq H_0\leq85\) a generating-truth stratum: every replicate has the same truth. This range may instead be a pre-specified **posterior-grid interior region** or a reporting restriction on posterior mass. These are different estimands. If truth-stratified performance is wanted, vary the generating H0 across predeclared values and report the strata separately.

### Monte Carlo precision

For a coverage estimate near 0.68 with 2,000 independent successful replicates, the binomial Monte Carlo standard error is approximately
\[
\sqrt{0.68(1-0.68)/2000}\approx0.0104,
\]
or roughly ±0.020 at 95% normal scale. This calculation assumes independent replicates and no loss to failure. At 0.95 coverage, the standard error is about 0.0049 (roughly ±0.0096). Failure rates reduce the effective number of evaluable replicates. The required R must also be justified for the paired continuous estimand using a pilot estimate of the standard deviation of \(D_r\) and a predeclared target MCSE/interval width.

Thus R=2,000 is a plausible candidate, not a proof of adequate precision. A 100-replicate pilot is a planning input, not confirmatory evidence; re-estimate precision separately across misspecification and validation-size cells.

### Grid resolution

A fixed \(N_{H_0}=3,201\) is not a convergence certificate. The existing project already found that grid refinement and prior-domain sensitivity are distinct: local grid convergence can be excellent while widening the prior changes posterior means and intervals materially. Use a frozen prior/domain policy and per-scenario refinement checks for posterior means, interval endpoints, normalization, and selection-integral error. Do not use a large grid to mask prior sensitivity.

### Failure accounting

Keep execution state and scientific evaluability distinct. `NOT_RUN` is a design/scheduling state, not an observed replicate outcome; `FAILED` is an attempted computational failure; `NON_EVALUABLE` is a completed run for which a defined estimand cannot be calculated; `SUCCESS` is an evaluable run. Publish denominators for all states and report both unconditional operational performance and conditional-on-success statistical performance. Never silently remove failed replicates.

## 7. Manuscript and artifact status

The proposed package names `manuscript/finite_validation_dark_sirens.md` as a completed manuscript. That path was not found in the branch during this inspection (GitHub returned 404); therefore manuscript completion is **not verified**. Even a present outline would not establish a submission-ready paper. Results, methods, limitations, verified primary references, figure provenance, and artifact hashes must be populated from audited outputs.

The branch workflow file contains a V2 smoke run and a paired V1/V2 comparison script. The smoke has only two replicates and a coarse H0/redshift grid, so it is an integration check, not a performance study. The current workflow definition also contains repeated V2 smoke-step/artifact-path entries; inspect the exact workflow execution and remove accidental duplication only after confirming which step is intended. Do not infer a successful latest run from a workflow definition or from an earlier commit's CI.

## 8. Revised execution order

1. **Exact-head CI and artifact verification:** capture the exact branch SHA, matrix job conclusions, test logs and generated V2/paired artifacts. A previous successful matrix does not validate subsequent commits.
2. **D-01 forensics:** search retained historical CSV/JSON/logs and reproduce the failing rows. If source records are absent, record the cause as unrecoverable; do not manufacture a root cause.
3. **Solver adjudication:** run V1/V2 on identical arrays for linear and quadratic designs; compare objective, fitted probability curve, projected-gradient/KKT residual, boundary status, conditioning, and failure classification. Add formal separation and active-set tests.
4. **Weight/positivity analysis:** evaluate fitted probability quantiles and inverse-weight concentration over validation support and the full target domain; distinguish theoretical coefficient-box extrema from observed pilot extrema. Any clipping/stabilization is a separately named estimator arm with explicit estimand consequences.
5. **Likelihood semantic audit:** trace every host-density and mixture path, ensuring numerator and selection normalization use the same reconstructed population measure.
6. **Uncertainty target:** choose between design-aware frequentist propagation and a coherent Bayesian joint model. Do not conflate them.
7. **Design lock:** predeclare estimands, H0 truths, priors/domain, validation sampling design, failure denominators, solver arms, scenario grid and MC precision targets.
8. **Pilot, then confirmatory run:** proceed only when the pilot supports the computational budget and precision target.
9. **Manuscript:** write results only after independent artifact audit. Keep synthetic scope explicit.

## 9. Current status

| Gate | Status | Reason |
|---|---|---|
| V2 software foundation | Substantially implemented; exact latest CI pending verification | V2 is box-constrained unpenalized logistic NLL, not ridge |
| D-01 historical cause | OPEN / UNKNOWN | Original failure records and traces not recovered |
| V1/V2 paired pilot | Integration path present; artifacts not yet independently reviewed here | Two-replicate smoke is not performance evidence |
| Positivity / weight stability | OPEN | Coefficient box permits extremely small probabilities; no calibrated acceptability rule |
| Double-correction | PARTIAL | No hidden-host mixture seen in inspected pilot path; not global audit |
| Uncertainty propagation | OPEN | FPC/Wald Method C not justified for specified design |
| Confirmatory design | NO-GO | Estimands and precision plan not frozen |
| Manuscript | NOT VERIFIED / NOT SUBMISSION-READY | Named file not found; results and references require audit |

**Bottom line:** The scientifically responsible next achievement is not a large R=2,000 run or a claim of “Newton + ridge closure.” It is a reproducible, paired, artifact-audited demonstration of how the inclusion estimator changes host-density reconstruction and H0 calibration, with the failure and weight-instability mechanisms measured rather than assumed.
