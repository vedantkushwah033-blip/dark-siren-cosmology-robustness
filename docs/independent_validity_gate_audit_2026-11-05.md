# Independent validity-gate audit — refreshed source and CI inspection

**Audit date:** 5 November 2026 (IST)  
**Repository:** `vedantkushwah033-blip/dark-siren-cosmology-robustness`  
**Branch inspected:** `research/logistic-inclusion-v2-newton-2026-11-05`  
**Inspected source revision:** `bdf3ccbf08108abd1e22db76bf923d041f9f3af7`  
**CI run:** [GitHub Actions run 42038139960](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/42038139960)

## Executive decision

**GO for source-level development and small, explicitly non-confirmatory diagnostics. NO-GO for a large confirmatory campaign and NO-GO for manuscript submission.**

This refresh supersedes status statements in earlier audit notes where the current source or exact-head CI has since changed. It does not erase historical uncertainty or convert smoke tests into scientific evidence.

The exact inspected branch head is `bdf3ccbf08108abd1e22db76bf923d041f9f3af7`. GitHub Actions run `42038139960` completed successfully. All three numerical-diagnostics matrix jobs (Python 3.10, 3.11 and 3.12) concluded success. The Python-3.12-only experiment/audit steps ran and passed on 3.12; those same conditional steps were skipped on 3.10 and 3.11 by design. This is evidence that the configured workflow passed on this exact SHA, not proof of scientific validity or real-data applicability.

## 1. Repository and workflow identity

The active research branch is 13 commits ahead of baseline `9195cf619aadd2e88e480bc6ffe4a7e42a606135`, zero commits behind. Its exact head is the source revision above. No merge or default-branch change was performed as part of this audit.

Open pull requests remain drafts or independent work streams. This audit did not merge, retarget, close or modify any pull request.

## 2. D-01 — historical optimizer defect

**Status: UNRESOLVED / HISTORICALLY UNIDENTIFIABLE.**

The current branch contains prospective tests for complete separation, near separation, boundary-limited fits and rank-deficient quadratic designs. Those fixtures test current software behavior. They do not recover the historical D-01 validation rows, traceback, optimizer result object, original convergence threshold or runtime environment.

No retrospective cause is assigned. In particular, it is not justified to label D-01 as separation, near-singular Hessian, line-search failure, floating-point failure or a coding defect without original records.

**Closure rule:** recover and inspect the original fixture, exact source revision, optimizer trace and environment; otherwise preserve the historically-unidentifiable classification permanently.

## 3. V1/V2 objective and implementation

### 3.1 What the current code actually optimizes

Files:
- `src/dark_siren_robustness/estimated_inclusion.py`
- `src/dark_siren_robustness/estimated_inclusion_v2.py`

V1 uses bounded L-BFGS-B for the ordinary logistic negative log-likelihood. V2 computes the same unpenalized logistic negative log-likelihood, analytic gradient and observed Hessian, then uses a damped Newton direction on free coordinates, projected-gradient fallback and monotone backtracking inside the same coefficient box.

**V2 does not add an L2 penalty to the objective.** The local variable named `ridge` is added to the Hessian system when constructing a numerical direction; it is not added to the reported objective. Describing this implementation as a ridge-logistic estimator is incorrect.

V2's box constraints are part of the estimator definition. A boundary-limited KKT solution is a constrained optimum, not an unconstrained finite MLE. Damping does not make a nonexistent unconstrained MLE exist.

V2 reports projected-gradient residual, Hessian numerical rank/condition, eigenvalue diagnostics and a coefficient-boundary flag. These are useful diagnostics but not a formal linear-programming separation detector. Hessian rank at the constrained fitted point is not a substitute for a formal separation classification.

### 3.2 Paired comparison evidence

The exact-head workflow ran the V1 and V2 pilot with:
- identical top-level seed and replicate count (2);
- identical parent and detected-event settings (1,200 and 8);
- validation sizes 100 and 300;
- identical redshift and H0 grid sizes (201 and 41);
- separate persisted V1 and V2 outputs, each passed through the estimated-inclusion artifact auditor.

The workflow then compared the persisted replicate/arm keys and explicitly checked equality of replicate seeds, catalogue seeds, audit seeds, catalogue counts, event counts, validation inputs, catalogue redshifts, redshift/H0 grids, observed distances, parent redshifts and true inclusion-probability scoring arrays.

This is genuine paired-input smoke evidence. It is **not** a large empirical solver benchmark. Two replicates cannot establish comparative stability, calibration, or superiority. The workflow's same-data tests also compare V1/V2 on ordinary well-conditioned linear and quadratic data, and on deterministic complete/near-separated fixtures; these are regression checks, not a comprehensive optimizer benchmark.

**Gate status: PARTIALLY CLOSED for smoke-level input parity; OPEN for empirical comparison and confirmatory inference.**

## 4. Separation, rank and KKT

V2's source computes
`g = X^T(pi-y)` and `H = X^T W X`, without a statistical penalty term. It checks a projected-gradient residual against a declared tolerance and identifies active box faces using gradient signs. The tests include:
- complete separation for linear and quadratic fits;
- one-discordant-observation quadratic near-separation;
- rank-deficient quadratic design;
- analytic gradient/Hessian finite-difference checks;
- active-bound KKT gradient-sign checks;
- same-input V1/V2 fitted-curve and objective comparisons.

No formal LP-based complete/quasi-separation detector is present in the inspected V2 module. The implementation therefore distinguishes some failure modes prospectively, but does not exhaustively classify all separation geometries. A near-separation fixture is not a general quasi-separation detector.

**Gate status: prospective constrained-solver regression coverage exists; formal separation classification remains OPEN.**

## 5. Positivity and inverse-weight stability

Files:
- `src/dark_siren_robustness/estimated_inclusion.py`
- `src/dark_siren_robustness/catalogue_density.py`
- `src/dark_siren_robustness/preposterior_diagnostics.py)

The estimator does not silently clip fitted probabilities or cap inverse weights. The KDE rescales weights by their maximum before weighted moments and normalization. This reduces avoidable common-scale overflow; it does not remove weight concentration or restore unsupported parent regions.

Pre-posterior diagnostics report fitted probability minimum and quantiles, inverse-weight maximum and quantiles, Kish ESS, entropy ESS, maximum normalized weight share, top 1%/5% weight mass, and validation-range/support-edge proxies. The support-gap metric is explicitly a min/max range proxy, not a formal local support or covariate-shift test.

The feasible coefficient boxes allow extremely small fitted probabilities. Those mathematical worst-case values must not be reported as observed pilot outcomes. No universal ESS threshold is established by the code.

**Gate status: diagnostic machinery exists; a design-calibrated positivity rule and confirmatory weight-stability evidence remain OPEN.**

## 6. Host-density construction and double-correction audit

The inspected estimated-inclusion path is:

1. Fit the inclusion model using external validation redshifts and retention indicators.
2. Predict inclusion probabilities at observed catalogue redshifts.
3. Form inverse-probability weights.
4. Build a reflected, weighted Gaussian KDE using the original observed sample's weighted spread and Kish ESS bandwidth rule.
5. Multiply by the grid-level merger-rate function.
6. Normalize the result to unit integral.
7. Pass that one density to `selection_aware_toy_likelihood`.

The likelihood uses that same normalized density in its event-host integrals and its detection-efficiency integral. It conditions on a fixed number of detected events and applies the toy hard distance threshold.

**No separate missing-host mixture term is present in this inspected estimated-inclusion pilot path.** Thus this path does not demonstrate simultaneous IPW and missing-host-mixture correction. It is not valid to generalize this finding to every likelihood module or historical experiment. The separate dataset-level mixture prototype supports an explicit missing-host component and documents its own contract; it is not the same production route.

The appropriate rule is not to impose a repository-wide `f=0` invariant blindly. Rather, each likelihood arm must declare whether its host density is already parent-reconstructed by IPW or instead represents a catalogue component requiring a missing-host mixture. Applying both corrections to the same missing population would be a model-definition error.

**Gate status: the inspected IPW pilot path is single-correction by construction; a repository-wide caller/factorization audit remains OPEN.**

## 7. Finite-validation uncertainty and bootstrap interpretation

Files:
- `src/dark_siren_robustness/estimated_inclusion.py`
- `src/dark_siren_robustness/posterior_validation_bootstrap.py`

Two distinct resampling implementations exist.

**Validation-density bootstrap:** resamples paired validation redshift/retention records with replacement, refits the inclusion model and rebuilds the host density while holding the observed catalogue fixed. The function supports V1/V2 solver selection. It is a frequentist conditional sensitivity diagnostic, not a posterior distribution over inclusion coefficients.

**Posterior-level validation bootstrap:** holds catalogue and GW distances fixed, resamples validation records, refits inclusion, rebuilds the density and recomputes the selected H0 posterior. Its paired-record bootstrap and multinomial-pseudo-population/SRSWOR branch are exploratory. The latter retains a pseudo-population posterior reference. The code explicitly labels its percentile summaries as resampling variability of posterior means, not H0 posterior credible intervals. Its scalar FPC-adjusted standard error is diagnostic only and has no claimed validity theorem for this nonlinear posterior functional.

Important integration limitation: `posterior_validation_bootstrap.py` calls `estimate_host_density_from_validation` without passing a solver selector. Consequently, this posterior-level bootstrap currently follows the default V1 path; it is not evidence of V2 bootstrap uncertainty propagation. The separate validation-density bootstrap does accept a solver version.

No joint Bayesian posterior over inclusion-model coefficients and H0 is implemented by these bootstrap routines. Averaging normalized posterior curves over arbitrary bootstrap or Wald draws must not be called Bayesian marginalization.

**Gate status: the bootstrap's computational meaning is documented and failure accounting is explicit; inferential calibration and V2 posterior-bootstrap integration remain OPEN.**

## 8. Prior domain and numerical grid

The estimated-inclusion pilot defaults to an H0 grid from 50 to 90 with 81 points (unless overridden), and requires the generating H0 to lie inside that interval. The exact-head smoke uses a much coarser 41-point H0 grid.

Other repository workflows contain grid-sensitivity experiments for other likelihood configurations. They do not automatically establish grid convergence for this estimated-inclusion/IPW host-density pipeline.

No evidence in this audit establishes a Track-A prior-domain sweep, a validated boundary-mass threshold, or convergence to a 3,201-point H0 grid. The proposal to lock every production grid to 3,201 points is not justified merely by choosing a large number; convergence should be assessed against a prespecified tolerance and computational cost.

**Gate status: OPEN for Track-A-specific prior-domain and grid-convergence analysis.**

## 9. Confirmatory protocol and manuscript readiness

The existing two-replicate V1/V2 runs are pipeline smoke tests. They do not estimate reliable Monte Carlo precision, coverage, bias or solver failure rates. The finite-validation bootstrap work is a separate estimand track and must not be pooled with the inclusion-misspecification campaign.

Before a large campaign, freeze:
- one primary estimand and its paired contrast;
- truth/scenario cells and validation-sampling design;
- H0 prior and grid-convergence tolerance;
- all failure states and denominators;
- positivity/weight diagnostic reporting;
- treatment of optimizer boundary/separation cases;
- whether bootstrap uncertainty is conditional sensitivity or a calibrated inferential interval;
- artifact provenance, source revision and checksums.

D-01 must remain unresolved unless historical evidence is recovered. Do not substitute a prospective fixture for a historical diagnosis.

## 10. Gate register

| Gate | Current decision | Evidence still required |
|---|---|---|
| Exact-head CI | **PASS** | Preserve run and SHA provenance |
| D-01 historical cause | **OPEN / UNIDENTIFIABLE** | Original logs/fixture/trace or permanent unknown classification |
| V1/V2 objective equivalence | **PASS at source level** | Continue adversarial and numerical parity checks |
| V1/V2 identical-input pilot | **PASS at smoke level only** | Adequately powered paired solver comparison |
| Formal separation detection | **OPEN** | LP-based detector and adversarial validation |
| Weight stability | **PARTIAL** | Design-calibrated positivity/instability criteria |
| IPW host-density path | **PASS for inspected pilot path** | Broader likelihood/caller audit |
| Double-correction | **PARTIAL** | Audit all relevant mixture callers and declare estimand per arm |
| Bootstrap interpretation | **PASS as conditional diagnostic description** | Calibration for any inferential interval; V2 integration if required |
| Prior/grid sensitivity | **OPEN** | Track-A-specific prior and grid refinement study |
| Confirmatory protocol | **OPEN** | Freeze estimands, cells, precision, failure policy and artifacts |
| Manuscript submission | **NO-GO** | Confirmatory evidence and independent final audit |

## 11. Explicit non-claims

This audit does not establish:
- a real LVK or galaxy-survey result;
- a measurement of observational H0;
- universal superiority of estimated IPW;
- nominal posterior coverage;
- validity of the scalar FPC adjustment for the nonlinear H0 functional;
- a new general statistical estimator;
- resolution of D-01;
- complete repository-wide absence of double correction;
- manuscript or journal readiness.

**Next engineering milestone:** add formal LP separation diagnostics as a separately tested diagnostic layer, without changing the V1/V2 estimator objective, and then add a deliberate V2 option to the posterior-level bootstrap only if that uncertainty track remains in scope.
