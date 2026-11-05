# Independent adjudication of external research dossiers

**Date:** 2026-11-05  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Purpose:** Reconcile the supplied literature-review and consortium reports against the inspected repository source. This is a source-level adjudication, not a confirmatory statistical result or release approval.

## Executive decision

The defensible project contribution remains a **synthetic, failure-aware benchmark** of finite-validation estimation of catalogue-inclusion probabilities and its downstream effects on host-density reconstruction and toy dark-siren (H_0) inference. The supplied literature review supports a conditional research-gap hypothesis, not a proof that no prior study exists. Novelty wording must remain qualified by the search scope and date.

The project is **not confirmatory-experiment ready**. In particular, the supplied consortium dossier contains material statements that conflict with the inspected implementation. Do not use it as a factual description of V2 or as evidence that D-01 was historically reproduced.

## 1. Repository facts checked on this branch

The following source files were fetched directly from the named branch:

- `src/dark_siren_robustness/estimated_inclusion_v2.py`
- `src/dark_siren_robustness/estimated_inclusion.py`
- `src/dark_siren_robustness/estimated_inclusion_study.py`
- `docs/logistic_inclusion_v2_design_and_validation_2026-11-05.md`
- `docs/critical_statistical_validity_gate_2026-11-05.md`

The branch exists. The V2 source identifies itself as a separately versioned, opt-in **box-constrained logistic negative-log-likelihood optimizer**. The branch documentation states V2 is not wired into the host-density estimator, pilot runner, bootstrap, or artifact schemas.

## 2. Corrections to the supplied consortium dossier

### 2.1 V2 is not ridge-regularized logistic regression

The dossier describes V2 as minimizing a ridge-penalized objective with gradient containing (-\lambda\beta) and Hessian containing (-\lambda I). That is not the implementation fetched from the branch.

The actual V2 objective is the ordinary logistic negative log-likelihood:

[
\ell(\theta)=-\sum_i\{y_i\log p_i+(1-y_i)\log(1-p_i)\}.
]

Its gradient and Hessian are:

[
g(\theta)=X^T(p-y),\qquad H(\theta)=X^T W X.
]

The algorithm uses temporary diagonal damping when solving for a Newton direction, plus projected Armijo backtracking and projected-gradient fallback. The damping is an **optimization safeguard**, not an L2 penalty in the statistical objective. The declared coefficient box changes the estimator to a constrained MLE. A boundary-limited solution is not an unconstrained finite MLE and does not prove positivity.

### 2.2 Historical D-01 cause is not established

The dossier states that D-01 failures occurred primarily because of complete/quasi-complete separation and ill-conditioned Hessians. Those are plausible failure mechanisms for logistic likelihoods, but the source and artifacts inspected here do not establish that they caused the historical D-01 rows.

The V2 design document explicitly leaves D-01 open pending replay of the recorded failures, including the two failures in the earlier 840-row pilot. Therefore:

- separation is a **candidate mechanism**, not a verified historical diagnosis;
- the two historical failed rows have not been shown to be reproduced by V2;
- passing unit tests and CI does not close D-01.

### 2.3 “IPW requires mixture fraction f=0” is too broad as a universal theorem

In the inspected estimated-inclusion pilot path, one IPW-weighted host density is passed into the selection-aware likelihood; the inspected path does not add a separate hidden-host mixture fraction. That path does not demonstrate double correction.

If the IPW reconstruction genuinely targets the full parent host distribution, adding a separate missing-host component for the same omitted population would generally double count it. But the statement (f=0) is valid only for that precisely defined construction. A different model could use weights that reconstruct only a subpopulation and represent a distinct residual population explicitly. The complete generative factorization, target population, and selection normalization must decide the issue—not a blanket rule.

### 2.4 The quadratic centering constant is not 1/3 in the current pilot

The inspected pilot computes its quadratic DGP centering constant as the second moment of standardized redshift under the declared parent intensity, using numerical integration on a 20,001-point grid. It does not simply assume a uniform standardized-redshift distribution and set the constant to (1/3).

The DGP centering convention is separate from the fitted quadratic basis. Production V1 and V2 fit the uncentered basis (a+bx+qx^2), preserving the same declared coefficient box. Do not describe this as orthogonalization: centering by a second moment alone does not generally make the quadratic term orthogonal to the linear term.

### 2.5 “R=2,000 provides high power” is not a demonstrated design result

A replicate count alone does not establish power or Monte Carlo precision. For a paired primary contrast (D_r), the estimated standard error is (s_D/\sqrt R), but (s_D) must be estimated from an appropriate pilot under each primary scenario. The target confidence-interval half-width or minimum scientifically meaningful effect must be specified before selecting (R). Failed and non-evaluable replicates must remain visible in the denominator accounting.

The proposed sign convention
[
D_r=|\widehat H_{0,cat}^{(r)}-H_{0,true}|-|\widehat H_{0,IPW}^{(r)}-H_{0,true}|
]
is coherent: positive values favor IPW on absolute posterior-mean error. It should be frozen only after the primary scenario, estimator version, prior, H0 domain, and treatment of failures are declared.

## 3. IPW first-order expansion: accepted with conditions

For a target parent density (f), inclusion probability (pi), and fitted probability (widehat\pi=\pi(1+\delta)), the normalized population reconstruction is

[
\widehat f(z)=
\frac{f(z)/(1+\delta(z))}
{\mathbb E_f[1/(1+\delta(Z))]}.
]

If (delta) is uniformly small enough for a first-order expansion and the required expectations exist,

[
\frac{\widehat f(z)-f(z)}{f(z)}
\approx -\{\delta(z)-\mathbb E_f[\delta(Z)]\}.
]

This is a useful population-level identity. It does not by itself establish finite-sample KDE behavior, positivity, bootstrap validity, or downstream posterior calibration. When fitted probabilities approach zero, higher-order terms and extreme weights can dominate. The analysis must report probability lower tails and weight concentration, not just average propensity prediction error.

## 4. Bootstrap and uncertainty language

The existing source-level validity gate describes paired validation-record resampling and a separate design-aware SRSWOR exploratory construction. It does not establish nominal coverage for the complete nonlinear estimator. Percentiles of bootstrap posterior means are not automatically (H_0) credible intervals or frequentist confidence intervals.

Keep these outputs explicitly exploratory until a target uncertainty estimand and its repeated-sampling properties are validated. Do not average bootstrap likelihoods or posterior arrays and call the result a Bayesian marginal posterior without a separate derivation.

## 5. Confirmatory gate: required before any large campaign

1. **D-01:** recover original raw failure records and replay them on V1 and V2; if irrecoverable, record the limitation and conduct a new prospective adversarial suite under a new identifier.
2. **V1/V2 parity:** compare linear and quadratic fits on identical validation records, reporting coefficients, objective, predictions, KKT residuals, boundary flags, and failure states. Do not treat solver agreement as downstream validity.
3. **Weight stability:** quantify minimum fitted inclusion probability, maximum raw and normalized weight, Kish and entropy ESS, top-1%/5% weight mass, and density-reconstruction error. Treat clipping/trimming as separate estimators and conduct sensitivity analyses if used.
4. **Generative coherence:** trace the exact host-density object used by the event numerator and (eta(H_0)); verify whether any missing-population term exists in each estimator arm. Keep the fixed-count conditioning explicit.
5. **Uncertainty target:** distinguish validation-sample variability, finite-parent/catalogue variability, event variability, inclusion-model uncertainty, numerical error, and Monte Carlo error.
6. **Design freeze:** predeclare primary estimand, primary scenario, estimator arms, prior/domain, numerical tolerances, failure denominator rules, and Monte Carlo precision target.
7. **Independent audit:** reconstruct summaries from raw replicate-level records, not producer-generated summaries alone.

## 6. Status and permitted claims

| Item | Status |
|---|---|
| Conditional novelty hypothesis | Plausible; literature search is not exhaustive proof of absence |
| IPW first-order population expansion | Derived under small relative error and integrability conditions |
| V2 optimizer | Implemented as damped, box-constrained Newton for ordinary logistic NLL |
| V2 statistical regularization | **Not implemented**; diagonal damping is algorithmic |
| Historical D-01 diagnosis | **Unverified / OPEN** |
| V2 downstream H0 integration | Not implemented on this branch per design document |
| Bootstrap coverage for full pipeline | Not established |
| IPW/mixture double correction | Not present in the inspected pilot path; audit any other path separately |
| Confirmatory simulation readiness | **NO-GO pending gates above** |
| Submission readiness | **NOT READY** |

## 7. Defensible contribution statement

> We conduct a reproducible synthetic benchmark of finite-validation estimation of galaxy-catalogue inclusion probabilities in a simplified dark-siren inference pipeline. The study evaluates how finite validation, nonlinear inclusion-model misspecification, and inverse-probability weight instability affect reconstructed host distributions and downstream (H_0) posterior behavior, while retaining optimizer and replicate failures in an explicit audit trail. The work evaluates established estimators and does not introduce a new estimator or claim observational validation.

This wording should be revisited after the confirmatory experiment and independent artifact audit; until then it describes the intended benchmark, not completed findings.
