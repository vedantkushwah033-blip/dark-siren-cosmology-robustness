# Independent adjudication of the November 2026 V2/IPW audit dossier

**Date:** 5 November 2026  
**Branch:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Purpose:** Correct mathematical overstatements and convert the supplied consortium dossier into testable research gates.  
**Scope:** Source-level and mathematical adjudication; no new simulation or test execution is claimed.

## Executive decision

The dossier is useful as a risk register, but its final classification **“Ready for Confirmatory Simulation” is not supported** by the inspected repository state. The current V2 solver is an opt-in box-constrained ordinary logistic MLE solver; it is not L2-regularized, not Firth estimation, and is not integrated into the pilot. The historical D-01 cause remains unknown. The IPW pilot has substantial observable weight diagnostics, but no empirically calibrated acceptability rule. The proposed R=2,000 design is a candidate, not a justified sample size.

**Decision: NO-GO for the confirmatory campaign until the estimand, truth strata, solver arm, failure policy, uncertainty target, and Monte Carlo precision are frozen and a planning study supports the computational design.** A bounded exploratory planning campaign may proceed after exact-head tests.

## 1. Corrections to the logistic mathematical claims

For binary outcomes and design matrix X, the stable negative log-likelihood is

\[
\ell(\beta)=\sum_i\{\log(1+e^{x_i^T\beta})-y_i x_i^T\beta\},
\quad g=X^T(\pi-y),\quad H=X^TWX.
\]

The Hessian is positive semidefinite. For finite coefficients, each logistic variance is strictly positive, so H is positive definite exactly when X has full column rank. If X is rank deficient, the objective is convex but not strictly convex in coefficient space; multiple coefficient vectors can yield the same linear predictors and fitted probabilities.

### Critical correction: separation does not make the NLL tend to minus infinity

Under complete or quasi-complete separation, there is a direction along which the likelihood increases toward its supremum, while the negative log-likelihood decreases toward a **finite infimum** (often zero for complete separation), not minus infinity. The coefficient vector has no finite unconstrained maximizer. This distinction matters: a finite objective and a small numerical gradient do not establish existence of a finite unconstrained MLE.

A finite box creates a constrained minimizer (the objective is continuous on a compact box), often on its boundary. That is a valid solution to the **constrained** optimization problem, not a recovered unconstrained MLE. Report it as boundary-limited and retain the box as part of the estimator definition.

### Quasi-separation wording

“Classifies all but a few points” is an intuition, not a sufficient mathematical definition. Complete and quasi-complete separation should be diagnosed using the formal sign/equality conditions for the design rows, preferably by a documented linear-programming implementation or a vetted separation routine. Large standard errors or extreme fitted probabilities are warning signs, not definitive separation tests.

### KKT sign convention

For minimization over l <= beta <= u, one consistent Lagrangian is
\[
\mathcal{J}=\ell(\beta)+\mu_L^T(l-\beta)+\mu_U^T(\beta-u),
\]
with \(\mu_L,\mu_U\ge0\). Stationarity is
\(g-\mu_L+\mu_U=0\), with complementary slackness. Thus at a lower bound the gradient must be nonnegative; at an upper bound it must be nonpositive. The V2 projected-gradient mapping is a useful box-KKT residual, but it should be validated against explicit active-set sign checks in tests.

### Rank deficiency is not automatically optimizer failure

A singular Hessian can prevent a unique coefficient solution while the fitted probability function remains identifiable on the observed support. V2 currently reports Hessian rank and condition number, and its rank-deficient test expects deterministic convergence. The correct interpretation is “coefficient non-identifiability / flat directions exposed,” not necessarily “fit failed.” Downstream predictions outside the observed support may nevertheless be unstable.

The dossier's \(10^{12}\) condition-number cutoff is a proposed heuristic, not a theorem. Use continuous condition diagnostics and calibrate any categorical threshold against perturbation sensitivity.

## 2. V2 algorithm classification

The inspected V2 source implements:
- ordinary Bernoulli logistic negative log-likelihood;
- analytic gradient and Hessian;
- a box-constrained damped Newton direction;
- projected Armijo backtracking;
- projected-gradient fallback;
- projected-gradient/KKT stopping;
- Hessian rank and condition diagnostics.

The diagonal term added to the Hessian is an **algorithmic damping device**, not an L2 penalty in the statistical objective. V2 is therefore not “specified with L2 regularization.” No Firth adjustment or Jeffreys-prior penalty is present.

Damping and trust regions can improve numerical steps; neither creates a finite unconstrained MLE under separation. Firth estimation changes the objective/score and is a distinct estimator. If compared, it must be an explicitly named third estimator arm with its own inferential interpretation—not silently substituted for V1/V2.

## 3. V1/V2 parity: what it can and cannot establish

The current test file contains same-data quadratic and linear parity tests. They compare fitted coefficients/probabilities and objective values on well-behaved synthetic data. These are useful tests of implementation agreement, but they do not prove robustness under separation, near-separation, ill-conditioning, or misspecification.

The ten-item checklist in the supplied dossier should be revised:
- Same observations, row order, feature basis, standardization, intercept, objective, and coefficient box are essential.
- Random seeds are irrelevant once the exact data arrays are fixed; both solvers are deterministic initialised at zero in the inspected implementations.
- Identical stopping tolerances are not required when algorithms use different stopping rules. Instead, compare independently computed objective values and a common KKT residual at a predeclared tolerance.
- Failure handling should preserve each solver's raw status and classify it independently; do not force identical exceptions.
- Serialization schemas need not be byte-identical if solver-specific diagnostics differ; define a common estimator-result contract with solver-specific metadata.
- Under separation, compare the same **bounded** optimization problem and active-set/KKT behavior. Do not compare against a nonexistent finite unconstrained MLE.

Required next parity fixtures: well-conditioned interior, complete separation, quasi-separation, near-separation, duplicated columns, near-collinearity, boundary-limited optimum, and extreme but finite predictor scaling. Record coefficient differences, prediction sup-norm differences on both observed support and full declared domain, objective difference, KKT residual, Hessian rank/condition, and status classification.

## 4. IPW diagnostics and positivity

For retained objects, \(w_i=1/\hat\pi(z_i)\). Kish ESS, entropy ESS, maximum normalized share, top-tail mass, probability quantiles, and validation-support diagnostics are useful complementary summaries. None alone certifies a reliable density estimate.

The current pre-posterior diagnostic implementation already computes:
- fitted probability min, 1st/5th percentiles, median and max;
- Kish ESS and ESS fraction;
- entropy ESS fraction;
- maximum normalized weight share;
- top 1% and 5% mass;
- maximum, 95th and 99th percentile inverse weights;
- catalogue fraction outside validation min/max and near validation edges.

The KDE rescales weights by their maximum before moments and normalization. This addresses common-scale overflow, not concentration, finite-sample variance, model misspecification, or lack of positivity.

Do not adopt “top 1% > 50% means unusable” as a universal law. It is a possible predeclared warning trigger for a particular design, but its consequences depend on catalogue size, redshift geometry, smoothing bandwidth, parent density, and downstream likelihood sensitivity. Calibrate decision thresholds by prospective simulations that measure density ISE, tail concentration, H0 error, interval width, coverage, and failure probability jointly.

Support diagnostics should not be limited to the validation sample's min/max. Also inspect local validation density and event/catalogue support over the full domain; min/max overlap can coexist with severe interior sparsity.

## 5. Correct IPW target and the mixture invariant

Let the parent galaxy density be f(z), retention probability be pi(z), and retained-catalogue density be proportional to f(z)pi(z). Under a valid validation design, positivity, and a sufficiently accurate fitted pi-hat, normalized IPW targets f(z) at the population level. If the merger-rate multiplier m(z) is applied afterward, the host target is proportional to f(z)m(z).

The statement “IPW requires f=0” is **not a universal theorem**. It is a model-specific invariant when a mixture parameter f denotes an additional missing-host component while the IPW density already reconstructs the full parent host population. In that particular factorization, retaining a nonzero missing component double counts missing population mass. Other mixture definitions may represent a different latent population partition and require a fresh derivation.

For a pure-IPW arm, use one reconstructed density p_IPW(z) in both:
\[
\int p_{IPW}(z) p(d_j\mid z,H_0)\,dz
\quad\text{and}\quad
\beta(H_0)=\int p_{IPW}(z)P(D=1\mid z,H_0)\,dz.
\]
If the numerator uses a mixture density, the selection normalizer must integrate that same mixture density. The numerator and denominator must describe the same generative host population.

Source inspection of the current estimated-inclusion pilot found no separate hidden-host mixture: it passes one normalized IPW/KDE host density to the selection-aware fixed-count toy likelihood. This is a finding about that inspected path only.

## 6. Uncertainty propagation: avoid the “Method C” shortcut

The dossier's claim that plug-in propensity uncertainty causes systematic undercoverage is too categorical. Ignoring first-stage estimation uncertainty often understates total repeated-sampling uncertainty, but the direction and magnitude of coverage error in a nonlinear posterior pipeline are not guaranteed without assumptions and empirical evaluation.

The current validation bootstrap refits the inclusion model and rebuilds the host density for each resampled validation dataset, conditional on a fixed observed catalogue and fixed GW distances. Its output is a distribution of posterior summaries across bootstrap datasets. This is a **frequentist sensitivity/uncertainty diagnostic**, not a Bayesian posterior over beta and not automatically a calibrated interval for H0.

A proposed “Method C likelihood marginalization” must not be implemented as simply adding a covariance matrix or averaging likelihoods. A coherent Bayesian joint model would specify a prior for inclusion parameters, the validation-data likelihood, the catalogue-generating mechanism, and the GW likelihood, then integrate over the shared parameters with the correct conditional dependencies. A modular/cut analysis is a different inferential target and must state which feedback is blocked. A bootstrap approximation is yet another procedure. These cannot be treated as interchangeable.

Before implementing a new method, specify:
1. the target probability statement (conditional posterior or repeated-sampling coverage);
2. the random quantities and conditioning set;
3. whether validation, catalogue, parent frame, and GW data are jointly generated or modularly separated;
4. the exact factorization and which data update which parameters;
5. a calibration experiment under correct specification and misspecification.

## 7. Confirmatory estimand and replication design

The paired contrast
\[
D_r=|\hat H_{0,cat}^{(r)}-H_{0,true}|-|\hat H_{0,IPW}^{(r)}-H_{0,true}|
\]
is interpretable: positive values favor the IPW arm in absolute posterior-mean error. It is not sufficient as a sole primary outcome. It should be accompanied by paired density ISE, bias, RMSE, posterior width, 68% coverage, calibration summaries, weight diagnostics, and all failure states.

“Interior stratum 55 <= H0 <= 85” is ambiguous if every replicate uses H0_true=70. Define whether this is:
- a grid-interior reporting region for posterior mass/coverage; or
- a stratum of different generating H0 truths sampled across [55,85].
These are different designs and must not be conflated. If truths vary, predefine the truth values and replicate allocation.

R=2,000 and Ndet=300 are not justified by assertion. Use a pilot to estimate the paired standard deviation of D, the Monte Carlo SE of coverage and failure rates, and the precision required for the smallest scientifically relevant effect. For a mean contrast, approximate planning uses
\[
R \approx (z_{1-\alpha/2}+z_{1-\beta})^2\sigma_D^2/\Delta^2,
\]
with design-specific simulation checks and adjustment for multiple primary cells. For coverage p-hat, MCSE is approximately \(\sqrt{p(1-p)/R}\), using the correct denominator and reporting failures separately. A fixed R must not be called adequately powered until these quantities are estimated.

The proposed four-state accounting is sensible, but define transitions and denominators: every planned replicate-arm cell is initially NOT_RUN; once attempted it becomes SUCCESS, NON_EVALUABLE, or FAILED under explicit rules. Never drop failed fits from the primary operational-success denominator. Report conditional-on-success performance separately from intention-to-simulate performance.

## 8. Release and research gates

| Gate | Current adjudication |
|---|---|
| Historical D-01 root cause | OPEN / UNKNOWN; no retained replay evidence recovered |
| V2 mathematical specification | Substantially specified; adversarial validation still required |
| V1/V2 same-data parity | Linear and quadratic tests present in source; exact-head execution must be verified |
| V2 in pilot | NOT INTEGRATED |
| Weight diagnostics | Rich observable diagnostics exist; acceptability thresholds not calibrated |
| IPW/mixture | No mixture in inspected pilot path; broader codebase not globally cleared |
| Uncertainty propagation | Bootstrap sensitivity exists; no validated joint Bayesian marginalization |
| Confirmatory design | NOT FROZEN / NO-GO |
| Manuscript | Methods/results skeleton possible; submission readiness not established |

## 9. Immediate implementation order

1. Recover historical D-01 raw inputs, logs, environment, and optimizer result if retained; otherwise formally close only the *replayability question* as unavailable, not the root-cause question.
2. Run the exact branch-head test suite, including linear/quadratic parity and V2 adversarial fixtures; capture commit SHA and CI provenance.
3. Add explicit active-set KKT sign tests and formal separation detection; distinguish rank deficiency from optimizer failure.
4. Stress fitted probabilities and weights over feasible coefficient boundaries and realistic fitted-data distributions; test KDE output perturbation, not just arithmetic finiteness.
5. Integrate V2 only behind an explicit solver option and preserve solver/version/failure metadata; do not replace V1 silently.
6. Specify a coherent uncertainty target before implementing any “Method C.”
7. Freeze primary and secondary estimands, truth strata, failure denominators, and MC precision before confirmatory simulation.

## Evidence statement

This document is a mathematical and source-level adjudication of the supplied dossier and inspected branch files. It reports no new test execution, simulation campaign, literature-systematic-review completion, historical D-01 replay, or confirmatory statistical result. It makes no change to main, no merge, and no release declaration.
