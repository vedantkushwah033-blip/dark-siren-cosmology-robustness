# Critical statistical validity gate

**Date:** 5 November 2026  
**Status:** source-level triage; not a scientific release approval  
**Branch:** `research/quadratic-parameterization-audit-2026-11-05`  
**Scope:** bootstrap interpretation, selection normalization, quadratic centering, and possible IPW/mixture double counting.

## Executive decision

The four concerns are important, but the blanket statement that all four are already demonstrated defects is not supported by the source inspected on this branch. The correct disposition is:

| Concern | Source-level finding | Gate disposition |
|---|---|---|
| Validation bootstrap | Existing routines resample paired validation records and recompute each replicate's density and H0 posterior. They do not average bootstrap likelihoods into one posterior. The code explicitly labels the output exploratory and makes no validity claim for the nonlinear functional. | **Not proven invalid as a descriptive resampling diagnostic; not validated as a confidence procedure or full uncertainty propagation.** Do not call its percentile summaries calibrated intervals. |
| \(\beta(H_0)\) normalization | In the fixed-threshold Gaussian toy, \(\beta(H_0)=\int p_{host}(z)P(d\ge d_{th}\mid z,H_0)dz\), and the conditional likelihood for a fixed number \(n\) of detected events is the product of event-mixture numerators divided by \(\beta(H_0)^n\). The implementation uses this formula and the same normalized host density in numerator and efficiency. | **Algebra is internally coherent under the stated toy assumptions; independent numerical/known-answer tests and grid-convergence evidence are still required.** This is not an event-rate likelihood. |
| Quadratic centering | The DGP computes \(c=E_{parent}[x^2]\) from the declared latent parent intensity using trapezoidal quadrature on the parent sampler's grid. It does not substitute 1/3. This is a declared reference-measure centering convention, not a general orthogonalization guarantee. | **Choice is documented and defensible; interpretation and finite-frame sensitivity need testing.** Fixed coefficient boxes make centered/uncentered constrained fits non-equivalent unless bounds are transformed. |
| IPW plus hidden-host mixture | The inspected estimated-inclusion pilot constructs one IPW-weighted catalogue host density and passes it to the selection-aware toy likelihood. No hidden-host mixture fraction is present in that path. | **No double counting demonstrated in this pilot path.** Any separate mixture likelihood or future joint model must be audited from its complete factorization before reuse. |

This gate therefore blocks confirmatory claims and expensive confirmatory simulation, but it does not justify deleting or declaring invalid every existing software component. The current code provides useful synthetic infrastructure and exploratory diagnostics; it has not yet established a validated end-to-end statistical procedure.

## 1. Bootstrap: identify the estimand before judging validity

### What the current code actually does

`posterior_validation_bootstrap.py` implements two distinct resampling constructions:

- `paired_bootstrap`: resamples paired (validation redshift, retention indicator) records with replacement, refits the inclusion model, rebuilds the weighted KDE, and recomputes a normalized H0 posterior for each replicate.
- `design_aware_srswor`: creates a multinomial pseudo-population of size N from the observed validation records, draws an SRSWOR sample, and computes a reference pseudo-population posterior and a resampled posterior for each attempt.

The routine retains posterior-mean draws and per-attempt posterior density arrays. Its percentile interval is a quantile of the replicate posterior means. It does **not** form a single posterior by averaging bootstrap likelihoods. The module explicitly calls the design-aware construction exploratory, says it is not the studentized finite-population-mean interval procedure, and labels the FPC-adjusted standard error diagnostic only.

### What remains unresolved

- Ordinary pairs bootstrap is not automatically design-consistent for a fixed finite parent frame sampled by SRSWOR.
- The multinomial pseudo-population construction is not automatically valid for the nonlinear logistic → inverse probability → KDE → selection-conditioned posterior functional.
- The scalar FPC multiplier is not a derived variance estimator for that nonlinear functional.
- The resampling distribution is conditional on the observed catalogue and GW distances; it does not propagate all outer population, catalogue, event, or model uncertainty.
- Percentiles of posterior means are not H0 posterior credible intervals and are not automatically confidence intervals with nominal repeated-sampling coverage.

### Required gate

Keep these routines as exploratory diagnostics only. A defensible replacement must specify the target (finite-frame conditional variability, repeated SRSWOR variability, or Bayesian coefficient uncertainty), derive or simulate its sampling properties, and demonstrate coverage/calibration against a known synthetic reference. Do not add a likelihood-averaging algorithm merely to satisfy a verbal requirement; the current code does not do that.

## 2. Selection normalization \(\beta(H_0)\)

The toy model declares:

\[
d\mid z,H_0\sim\mathcal N(cz/H_0,\sigma^2),\qquad
D=\{d\ge d_{th}\}.
\]

For a normalized pre-detection host density \(p_{host}(z)\),

\[
P(D\mid H_0)=\beta(H_0)
=\int p_{host}(z)
\Phi\!\left(\frac{cz/H_0-d_{th}}{\sigma}\right)dz.
\]

For n conditionally independent detected events with observed distances \(d_1,\ldots,d_n\), conditioned on the fixed detected count, the likelihood is

\[
L(H_0)=\prod_{j=1}^{n}
\frac{\int p_{host}(z)\,p(d_j\mid z,H_0)\,dz}
{\beta(H_0)}.
\]

Thus the log likelihood subtracts \(n\log\beta(H_0)\). The current `selection_toy.py` computes the event integrals and beta with trapezoidal quadrature on the same redshift grid and same supplied normalized host density. It checks the density integral and beta range.

This derivation is conditional on the declared toy assumptions: fixed n, independent events, a hard distance threshold, Gaussian distance noise, no rate evolution beyond what is already folded into the supplied host density, and no explicit event-rate/count likelihood. It is not automatically the selection normalization for a realistic dark-siren hierarchical model.

### Required gate

1. Add analytic/known-answer tests in limiting regimes (for example, very low and very high threshold relative to the distance distribution).
2. Independently compute beta by a second numerical route on a refined grid.
3. Verify \(\beta\in(0,1]\), convergence of beta and posterior summaries with redshift-grid refinement, and consistency between numerator and selection model.
4. State explicitly that the current likelihood conditions on a fixed number of detected events and omits rate information.
5. Audit any future mixture or population model by deriving its numerator and denominator together; do not patch beta in isolation.

## 3. Quadratic centering

The inclusion DGP uses

\[
\operatorname{logit}\pi(x)=\alpha+\beta x+\gamma(x^2-c),
\qquad c=E_{parent}[x^2].
\]

The implementation calculates c from the latent parent-redshift intensity on a 20,001-point trapezoidal grid matching the parent sampler's numerical reference. For a uniform symmetric x on [-1,1], c=1/3; the actual latent intensity is not uniform, so 1/3 is not the correct declared parent-law moment here.

This centering makes the parent-law expectation of \(x^2-c\) zero. It does **not**, by itself, guarantee that the quadratic column is orthogonal to the linear column; that additionally depends on \(E[x^3]\). Nor is the theoretical intensity moment identical to the realized empirical moment of each finite parent draw.

A separate issue is parameterization constraints. With centered coefficients (alpha,beta,gamma) and uncentered coefficients (a,b,q),

\[
a=\alpha-\gamma c,\quad b=\beta,\quad q=\gamma.
\]

The existing uncentered intercept box \([-25,25]\) maps to the q-dependent centered interval \([-25+qc,25+qc]\). Applying the same fixed numeric intercept box in both bases changes the feasible set.

### Required gate

- Keep the DGP's documented c unless a predeclared alternative is scientifically motivated.
- Test the actual numerical value of c and its convergence under quadrature refinement.
- Report the finite-parent empirical moment as a diagnostic, not as a replacement chosen after outcomes.
- Compare centered and uncentered optimizers with transformed bounds; separately study wider-box sensitivity.
- Do not claim orthogonality merely from centering by \(E[x^2]\).

The algebraic regression tests have been added in `tests/test_quadratic_parameterization_audit.py`; their exact-head CI was pending at the time of this document. They do not yet establish constrained optimizer equivalence.

## 4. IPW and hidden-host mixture

The estimated-inclusion pilot currently estimates one host density from the retained catalogue using estimated inclusion probabilities, then supplies that density to the fixed-count selection-aware toy likelihood. In the inspected path there is no explicit hidden-host mixture term \((1-f)L_{catalogue}+fL_{missing}\). Consequently, the alleged double counting is not established for this specific pipeline.

That does not settle the issue for other likelihood modules, earlier branches, or a proposed extension. The correct audit object is the full generative factorization, including:
- parent host population;
- catalogue-retention mechanism;
- conditioning on the observed catalogue and/or its size;
- whether the IPW reconstruction estimates a pre-retention parent density;
- GW detection and fixed-count conditioning;
- any missing-host mixture component;
- selection efficiency and rate terms.

### Required gate

Do not combine IPW reconstruction with a hidden-host mixture fraction until the estimand and conditional likelihood are derived from one joint model. If the IPW density already estimates the parent host distribution, adding a missing-host component may duplicate the incompleteness correction; if it estimates only a catalogue-conditioned component, a mixture may have a different interpretation. The code must make that distinction explicit and test limiting cases such as complete catalogue retention and zero missing-host fraction.

## 5. Stop/go decision

### Allowed now
- Source-level mathematical review.
- Deterministic unit tests for algebraic identities and numerical integration.
- Reduced smoke tests used only for software verification.
- Exploratory simulations clearly labelled non-confirmatory.
- Development of a separate uncertainty-propagation method after its target is defined.

### Not allowed as a scientific claim yet
- A validated nonlinear finite-population bootstrap interval.
- Calibrated uncertainty propagation for estimated inclusion coefficients through the H0 posterior.
- Confirmatory evidence that estimated inclusion weighting improves H0 inference.
- A realistic LVK or survey result.
- A claim that all four concerns have already been empirically shown to invalidate the current pilot.

### Required next sequence

1. Complete the selection-normalization known-answer and grid-convergence tests.
2. Finish the centered/uncentered constrained-optimizer equivalence test with transformed bounds.
3. Write a joint factorization for any model that includes both IPW and a hidden-host mixture; leave the current no-mixture pilot unchanged until then.
4. Recast bootstrap outputs strictly as exploratory resampling diagnostics and design a known-reference repeated-SRSWOR study for the nonlinear functional.
5. Only after these gates pass, freeze a confirmatory protocol and allocate the expensive simulation budget.

## 6. Evidence boundary

This document records source inspection and mathematical derivation only. It does not claim that the required tests have been implemented or passed, that the full curvature screen has run, or that an independent artifact reconstruction has been completed. No production estimator, likelihood, DGP, configuration, or result was changed by this document.
