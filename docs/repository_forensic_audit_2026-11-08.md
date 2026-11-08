# Repository forensic audit — 8 November 2026

**Audited ref:** `research/independent-posterior-audit-claim-fix`  
**Audited head before this report:** `e012762e41899a2ae442b4a420b651bcc4492c00`  
**PR:** #26, open/draft/unmerged  
**Default branch:** `main`  
**PR #26 base:** `experiment/catalogue-density-estimators` (not `main`)  
**Scope:** source inventory, active inclusion/host-density/selection code, tests, documentation, workflows, history and open PR topology. This is a source audit, not a claim that all 481 files have been line-by-line re-reviewed in this session.

## 1. Repository inventory and branch topology

The recursive Git tree at the audited branch head contained 481 files:

| Area | Files |
|---|---:|
| Python source modules | 19 |
| Tests | 116 |
| Documentation | 167 |
| Experiment/audit scripts | 114 |
| GitHub workflows | 60 |
| Configurations | 1 |

The branch has 69 commits ahead of its PR base and 26 changed files in the current comparison. The repository default branch is `main`; PR #26 targets `experiment/catalogue-density-estimators`, so it is a stacked research PR, not a direct main-bound release. Keep it draft and unmerged until its intended base and stack dependencies are reconciled.

### Parallel open research PRs discovered

The open-PR inventory contains several separate, unmerged workstreams relevant to the same question:

| PR | Topic | Relationship / audit implication |
|---|---|---|
| #15 | Estimated-inclusion pilot and independent artifact audit | A fuller posterior-propagation pilot; its validation audit cohort can overlap the retained catalogue, unlike the disjoint design in #26. Its estimand and overlap handling must not be conflated with #26. |
| #16 | Finite-parent SRSWOR and bootstrap comparison | Adds finite-frame sampling and posterior-level resampling. Its own documentation marks nonlinear bootstrap validity as unresolved. |
| #20 | Independent mathematical audit | Identifies incompatible validation-design factorizations and an incomplete observer-facing likelihood. |
| #21 | Independent finite-parent optimizer audit | Adds a distinct analytic-gradient SLSQP reference to a bounded optimizer path. |
| #22 | Parent-law quadratic centering | Corrects a different ([-1,1]) parameterization's (c=1/3) assumption; reports (c_g\approx0.284500557016) for that parameterization. |
| #23 | Planned-attempt ledger | Adds `SUCCESS`, `FAILED`, `NON_EVALUABLE`, and `NOT_RUN` planning/checkpoint primitives. |
| #24 | Scenario lifecycle and durable checkpoints | Integrates failure accounting into the finite-parent runner; explicitly states it is not resume support. |
| #25 | Decimal active-set oracle for bounded logistic V2 | Adds a high-precision test-only optimizer oracle; several stress cases remain uncovered and exact-head CI is required. |

These PRs are not merged into #26. Their branches form a separate stack in places. Do not assume their code, tests or results exist on the #26 head. Do not cherry-pick or merge blindly: first reconcile the estimand, validation overlap, parameterization, solver revision, artifact schema and exact CI provenance.

## 2. File-by-file source review

| Module | Actual responsibility | Finding / remaining limitation |
|---|---|---|
| `catalogue_density.py` | Reflected, bounded weighted KDE; optional (1/\hat\pi) weights; merger-rate multiplier; finite-domain renormalization. | The reflected sample triples the rows and splits each weight by three. With SciPy's default Scott rule this also triples the augmented Kish effective size and can shrink bandwidth artificially. The current branch now supplies a callable Scott factor based on the original weighted sample when no explicit bandwidth is supplied. Explicit numeric bandwidth settings retain their previous semantics. Inverse-weight overflow now fails explicitly. New regression tests were added; exact-head CI is pending. |
| `generative_selection.py` | Synthetic parent intensity, monotone catalogue-retention curve, merger-weighted host draws, Gaussian noisy distances and hard GW detection threshold. | Separates truth tables from observable tables. This generator's default retention curve is not the quadratic-logistic DGP used by the finite-validation benchmark; do not treat them as the same experiment. The GW sample is fixed-count rejection sampling, not a rate model. |
| `estimated_inclusion.py` | Unpenalized BFGS logistic MLE, linear/quadratic redshift terms, validation-standardized predictors, observed-information covariance. | This is V1 only. Separation detection is heuristic (coefficient and Hessian-condition cutoffs); BFGS optimizer status is partly superseded by a score-residual acceptance rule. It is not Firth, Bayesian, profile-likelihood, sandwich or FPC-adjusted inference. The inverse-Hessian covariance is a local model-based approximation. |
| `validation_design.py` | SRSWOR validation draw from the full parent frame; validation units excluded from the retained analysis catalogue. | Good for the explicitly disjoint primary design. It is a different design from PR #15's overlapping audit cohort. It does not implement nested validation prefixes or a durable planned-attempt ledger. |
| `inclusion_diagnostics.py` | Held-out Brier/log-loss and inverse-weight ESS/concentration summaries. | ESS arithmetic now uses scaled inverse weights to avoid avoidable overflow. Calibration intercept/slope, reliability bins, probability quantiles, boundary frequencies, weight quantiles and support-gap diagnostics remain incomplete. |
| `estimated_inclusion_benchmark.py` | Probability-stage synthetic DGP and paired validation-size/model-degree experiment. | Uses (N=4000) by default, (z\in[0.01,0.40]), a non-uniform parent intensity, and (
\pi(z)=\operatorname{expit}(0.6-x+\gamma(x^2-c))), with (x\in[0,1]) and (c=E_f[x^2]) evaluated by quadrature. Curvature settings are now (-2,-1.2,0,1.2,2). Parent frames and retention uniforms are common across curvature scenarios. Validation objects are excluded from analysis. The new raw export records the complete parent partition and fitted predictions. It still stops at probability/weight diagnostics; it does not construct an H0 posterior. |
| `host_density_benchmark.py` | Three-arm host-density-to-H0 comparison: parent oracle, catalogue-only and exact-known-completeness. | Its `known_completeness` arm is an oracle-known probability correction, not a finite-validation fitted model. It does not use the new finite-validation fitter. |
| `selection_toy.py` | Fixed-count selected-event toy likelihood: Gaussian distance mixture and (-N_D\log\beta(H_0)). | Numerator and detection efficiency use the same normalized supplied host density and hard threshold. This is internally explicit for the declared toy model, not a detector pipeline or rate likelihood. |
| `selection.py` | Alternative finite-interval distance-threshold selection functions and missing-host toy mixture. | Separate module/model from `selection_toy.py`; avoid mixing their selection rules without an explicit derivation. |
| `cosmological_inference.py` | Weighted Gaussian distance-mixture posterior using low-z or flat-ΛCDM distance grids. | Toy posterior primitive; not a complete hierarchical event-population likelihood. |
| `inference.py` | Legacy simplified low-z weighted Gaussian-mixture posterior. | Explicitly omits GW detection selection. |
| `hardened_inference.py` | Validated reference implementation of the simplified low-z mixture, including optional missing-host and positive-distance truncation. | Separate from the selected likelihood; no general detector or survey model. |
| `dataset_likelihood.py` | Sums prior-free event log likelihoods and applies one uniform H0 prior; exposes separate selection helpers. | The dataset posterior function itself is explicitly not selection-aware and does not integrate shared nuisance parameters or event rates. |
| `completeness.py` | Generic distance-dependent inclusion and Bernoulli thinning utilities. | Illustrative mechanism only; not a survey calibration. |
| `cosmology.py` | Low-z distance and flat-ΛCDM numerical distance utility. | Flat-ΛCDM implementation loops over redshifts and uses fixed-step trapezoidal integration; convergence must be assessed for any chosen domain. |
| `uncertainty_likelihoods.py` | Three Gaussian distance-uncertainty likelihood parameterizations. | One-host low-z toy, separate from catalogue-inclusion uncertainty. |
| `uncertainty_sensitivity.py` | Paired one-host distance-uncertainty sensitivity experiment and contrasts. | Does not answer finite-validation catalogue inclusion. |
| `calibration.py` | Generic interval and repeated-sampling metrics. | Useful utility; downstream calibration needs dedicated simulation design and Monte Carlo uncertainty. |
| `__init__.py` | Package version and scope note. | No implementation logic. |

## 3. Finite-validation benchmark audit

### Current DGP

The parent redshift density is proportional to the declared two-overdensity intensity. The code's (x=(z-0.01)/0.39) parameterization is on ([0,1]), so its numerical centering constant is (c\approx0.3114962521) under that intensity. This is not the (c\approx0.2845005570) reported by PR #22, which uses a different (x\in[-1,1]) parameterization. These values are not contradictory; they must never be copied between parameterizations.

The curvature grid values are stress settings, not yet empirically justified labels for “moderate” or “severe.” The induced probability curves are asymmetric in effect because of the fixed intercept and slope. A future report should show each curve and its actual probability range rather than call the positive and negative values equally severe.

### Validation and pairing

- Validation is sampled uniformly without replacement from the entire finite parent frame before conditioning on retention.
- Validation objects are excluded from the analysis catalogue.
- The fitter receives only validation redshift and retention labels.
- The analysis catalogue contains only retained non-validation objects.
- Common parent redshifts, retention uniforms and validation IDs are now reused across curvature conditions within replicate.
- True inclusion probabilities are scoring/oracle information only; they are not passed to the fitter.

### Raw artifact and independent auditor

The new raw output now has one row per parent object per experiment/model arm, with explicit role `validation`, `analysis_catalogue`, or `nonretained_holdout`, object ID, redshift, retention label, simulation-truth probability and fitted probability. The summary also stores fitted coefficients, covariance, and standardization metadata.

The independent auditor:
- checks that all parent IDs (0,\ldots,N-1) occur exactly once;
- checks the validation/analysis/non-retained partition and the reported catalogue count;
- reconstructs held-out Brier score, log loss and probability RMSE;
- reconstructs estimated and oracle IPW diagnostics;
- reconstructs fitted probabilities from coefficients and standardization metadata;
- recomputes the logistic score, observed information and inverse-information covariance;
- rejects incomplete partitions, inconsistent predictions and tampered summaries.

This is a numerical/provenance audit of the probability-stage artifacts. It does not establish estimator superiority, probability calibration across all conditions, density reconstruction, or H0 posterior calibration.

## 4. Mathematical findings

For parent density (f(z)), true inclusion (pi_0(z)), and a fixed fitted function (	ilde\pi(z)>0), the retained density is
\[
f_{cat}(z)=\frac{f(z)\pi_0(z)}{E_f[\pi_0(Z)]}.
\]
The normalized IPW target is
\[
f_{IPW}(z;\tilde\pi)=
\frac{f(z)\pi_0(z)/\tilde\pi(z)}
{E_f[\pi_0(Z)/\tilde\pi(Z)]}.
\]
Thus misspecification creates a redshift-dependent distortion proportional to (pi_0/\tilde\pi), not merely a vague “bias.” The population linear-logistic pseudo-true coefficient solves
\[
E_f\{h(Z)[\pi_0(Z)-\pi_\beta(Z)]\}=0,
\]
for (h=(1,x)^T), under parent-frame SRS validation. The current math note records this derivation.

For the GW toy model, with parent host density (p_{host}(z)), Gaussian distance model and hard threshold (D=\{d\ge d_{th}\}),
\[
\beta(H_0)=\int p_{host}(z)
\Phi\!\left(\frac{c_{light}z/H_0-d_{th}}{\sigma}\right)dz,
\]
and the fixed-detected-count log likelihood subtracts (N_D\log\beta(H_0)). This is distinct from galaxy-catalogue inclusion. The host density in event numerators and the GW detection denominator must be the same. The current finite-validation stage has not yet propagated estimated inclusion into this likelihood.

## 5. Test and workflow coverage

Relevant tests reviewed include:
- `test_estimated_inclusion.py`: large-sample linear recovery, quadratic prediction, one-class sample, rank deficiency, invalid inputs and complete separation.
- `test_validation_design.py`: SRSWOR/disjointness, reproducibility, hidden-probability exclusion and input contracts.
- `test_inclusion_diagnostics.py`: ESS behavior, concentration and proper scores; extreme-probability stability test added on this branch.
- `test_estimated_inclusion_benchmark.py`: paired denominator and model-specification coverage.
- `test_audit_finite_validation_probability_benchmark.py`: raw reconstruction, tampered-summary rejection, full-parent partition, coefficient/prediction consistency and common-random-number checks.
- `test_catalogue_density.py`: normalization, known-probability weighting, merger weighting, invalid inputs, reflected-bandwidth reference and inverse-weight overflow.
- `test_generative_selection.py`: separated observable/truth generator and fixed-threshold data generation.
- `test_selection_toy.py`: independent quadrature, (-N\log\beta), invalid data and empirical detection-efficiency check.

The `tests.yml` workflow runs a broad cross-version suite and multiple smoke/full numerical diagnostics. The exact-head run for the latest changes has not yet completed successfully. Intermediate runs were cancelled as new commits arrived; cancellation is not a pass. No local execution is claimed.

## 6. Blocking scientific gaps

1. No finite-validation estimated-probability arm is connected to host KDE and H0 posterior on this branch.
2. No Firth, Bayesian logistic, profile-likelihood, sandwich or valid finite-population design covariance method is implemented here.
3. The current inverse-Hessian covariance is model-based; a scalar FPC multiplier is not justified for logistic MLE covariance.
4. No valid coefficient-uncertainty integration into the full selected likelihood is implemented. Bootstrap likelihood averaging is not authorized without a derived target.
5. No H0 coverage/SBC experiment exists for estimated-inclusion arms on this branch.
6. No complete planned-attempt ledger or durable checkpointing is integrated into this branch; PRs #23/#24 are separate unmerged work.
7. No empirical diagnostic threshold has been calibrated to predict downstream posterior failure.
8. Historical D-01 has not been recovered/reproduced on this branch.
9. The exact likelihood and observer-facing finite-frame observation model need a joint factorization before claiming a full Bayesian treatment.
10. The manuscript is not yet a results paper for finite-validation estimated inclusion; results must remain placeholders until experiments and audits complete.

## 7. File changes made on this branch during this audit

- `src/dark_siren_robustness/estimated_inclusion_benchmark.py`: full parent partition export, fitted parameter metadata, common random numbers across curvature scenarios, negative and positive curvature design.
- `scripts/run_finite_validation_probability_benchmark.py`: versioned JSON configuration support.
- `scripts/audit_finite_validation_probability_benchmark.py`: full parent-frame and partition audit; independent reconstruction of predictions, score, Hessian and covariance; stable ESS calculations.
- `tests/test_audit_finite_validation_probability_benchmark.py`: tamper, partition and pairing tests.
- `src/dark_siren_robustness/catalogue_density.py`: corrected default reflected-KDE bandwidth factor to use original weighted-sample effective size; inverse-weight overflow guard.
- `tests/test_catalogue_density.py`: reflected-bandwidth reference and overflow regression tests.
- `src/dark_siren_robustness/inclusion_diagnostics.py`: stable scaled inverse-weight ESS arithmetic.
- `tests/test_inclusion_diagnostics.py`: extreme-probability regression test.
- `docs/finite_validation_mathematical_framework_2026-11-08.md`: derivations for IPW misspecification, pseudo-true linear logistic fit, SRSWOR, FPC limits, and distinct GW selection.
- `configs/finite_validation_probability_pilot_v1.json`: reproducible probability-stage design.
- This forensic audit report.

## 8. Scientific decision

**Conditional GO for probability-stage software validation; NO-GO for downstream H0 claims.**

The current branch now has a clearer DGP, paired curvature design, full raw parent partition and an independent probability-stage auditor. Until exact-head CI passes and the configured pilot is executed and audited, even these changes are implementation work rather than validated numerical results. The central research question remains open: no evidence on this branch yet shows when finite-validation estimated inclusion improves, fails to improve, or harms H0 posterior calibration.
