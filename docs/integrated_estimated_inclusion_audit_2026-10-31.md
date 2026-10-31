# Integrated scientific audit and execution gates: estimated-inclusion study

**Date:** 31 October 2026  
**Branch:** `research/estimated-inclusion-audit`  
**Status:** design reconciliation; no new confirmatory results are claimed.  
**Scope:** the existing synthetic estimated-catalogue-inclusion pilot and its planned nonlinear/diagnostic extensions.

## Executive decision

Retain the estimated-inclusion project, but keep its claim narrow:

> Under a declared synthetic parent-frame validation design, how do finite validation samples and deliberate nonlinear misspecification of a fitted inclusion model affect inverse-probability-weighted host-density reconstruction and downstream fixed-truth (H_0) inference?

This is a calibration and identifiability study, not a new dark-siren likelihood, not a new completeness estimator, and not an observational (H_0) measurement. The current PR #15 implementation is a **non-confirmatory pilot**. Its existing nonlinear and pre-posterior documents are protocols, not results.

The two external audit notes are directionally useful, but their proposed four-arm matrix should not replace the repository's existing design. The repository already has the more interpretable practical comparator (catalogue-only), an oracle diagnostic, and estimated-inclusion arms at multiple validation sizes. Nonlinearity is a data-generating factor ((g)), not a reason to discard the validation-size contrast.

## 1. Reconciled estimand and arm structure

For replicate (r) and validation size (m), define

[
D_{r,m} =
left|widehat H_{0,mathrm{catalogue},r}-H_0^{mathrm{true}}ight|
-
left|widehat H_{0,mathrm{estimated}(m),r}-H_0^{mathrm{true}}ight|.
]

The primary contrast is the mean paired (D_{r,m}), with a replicate-level paired bootstrap interval. Positive values favor estimated-inclusion weighting over catalogue-only inference. Report the number attempted, successful in both arms, and failed in each arm; never hide differential failures by reporting complete-case contrasts alone.

Use these existing arms:

1. **Catalogue-only:** equal catalogue weights; the practical baseline.
2. **Oracle:** true inclusion probabilities; a simulation-only diagnostic, not a deployable method and not guaranteed to be an upper bound in every finite realization.
3. **Estimated-(m):** probabilities fitted from the independently sampled validation records, for each prespecified (m).

Vary the true quadratic coefficient (g) across the preregistered levels. The fitted model remains linear logistic. A permuted-validation-label arm may be added as a separately named negative control, but it must not replace catalogue-only or the oracle comparison.

Keep empirical 68% coverage, Wilson intervals, bias, RMSE, interval width, density error, fit failures, coefficient-boundary hits, and weight diagnostics as prespecified secondary outcomes. The primary comparison is estimated versus catalogue-only; the oracle contrast is diagnostic.

## 2. Separate the two selection mechanisms

Catalogue inclusion is a galaxy-survey process, represented in this study by

[
pi(z)=P(R=1mid z).
]

The event acceptance rule is a separate process. The current toy code implements

[
d_{m obs}mid z,H_0sim
mathcal N(cz/H_0,sigma^2),
qquad
mathrm{accepted}iff d_{m obs}ge d_{m cut}.
]

For the declared normalized host density (g(z)), its selection efficiency is

[
eta(H_0)=int g(z)
Phi!left(rac{cz/H_0-d_{m cut}}{sigma}ight),dz,
]

and, conditional on a fixed number (N) of accepted events, the likelihood is

[
mathcal L(H_0)=
prod_{k=1}^{N}
rac{int g(z),
mathcal N(d_kmid cz/H_0,sigma^2),dz}
{eta(H_0)}.
]

This is the normalization implemented by the current selection-aware toy likelihood. It is mathematically coherent for the stated right-tail truncation. However, **(d_{m obs}ge d_{m cut}) is not a realistic gravitational-wave detector selection rule**: a real GW detection threshold is normally expressed through signal-to-noise ratio and detector response, with detectability generally decreasing as luminosity distance increases. Therefore call the current mechanism a *synthetic upper-tail observed-distance acceptance rule*, not a physical GW detector model or LVK selection function. A physical detector model requires a detection probability over latent source and detector variables, normally calibrated with injections.

The generic expression (P(mathrm{det}=1mid d_{m obs})p(d_{m obs}mid S,H_0)) should not be inserted as an unexplained factor. Derive the selected likelihood from the actual data-generating acceptance rule. If selection depends on latent parameters (	heta), integrate (P(mathrm{det}mid	heta)) against the pre-selection population and normalize by the corresponding population efficiency. If acceptance is a deterministic region in observed data, include its indicator and the probability of that region under the measurement model.

Catalogue inclusion probabilities must never be substituted for GW-event detection efficiency.

## 3. Inclusion-model and positivity contract

The current data-generating law is

[
operatorname{logit}pi_{m true}(z)
=a+bx+g(x^2-1/3),
quad
x=rac{z-(z_{min}+z_{max})/2}{(z_{max}-z_{min})/2},
]

while the fitted model is linear in (x). The subtraction (1/3) centers (x^2) under a uniform reference distribution on ([-1,1]) only. The repository's parent-redshift distribution is non-uniform, so this term is not generally orthogonal to the intercept under the actual parent population. Describe (g) as a prespecified curvature perturbation; do not claim population-orthogonality.

Inverse-probability weighting requires positivity and support overlap. It cannot reconstruct a region with (pi(z)=0), nor can it invent hosts in a structural-zero region absent from the observed catalogue. Report the true probability range and low-probability mass for simulation explanation, but keep these oracle quantities out of operational pre-posterior predictors. Report fitted probability quantiles, maximum weights, ESS and ESS fraction for the estimator actually available to the analysis.

The current estimator treats fitted probabilities as plug-in values. It does not propagate the full sampling uncertainty of the fitted logistic coefficients into the (H_0) posterior. This is a material limitation, not a minor implementation detail. The paper must state it explicitly. A later, separately validated extension could use a validation-sample bootstrap or joint inclusion-model/posterior uncertainty propagation; do not claim that uncertainty has already been propagated.

## 4. What the existing pre-posterior diagnostics do—and do not do

The current feature builder computes observable-only fitted-weight and validation-support summaries before posterior construction. This is a useful feature-extraction layer, not yet a validated warning system.

A predictive study is a separate stage. It must:
- split by replicate seed, keeping every arm from a replicate in one partition;
- use development replicates to choose features, score, and thresholds;
- freeze the complete rule before evaluating an untouched test set;
- exclude true inclusion probabilities, true coefficients, posterior summaries, coverage, and (H_0) error from predictors;
- compare with simple baselines (validation size alone, ESS fraction alone, and a prespecified weight-tail rule);
- report prevalence, confusion matrix, sensitivity, specificity, false-alert rate, AUROC, precision-recall performance, and uncertainty; and
- include both difficult cases: extreme weights without large posterior error, and posterior error without extreme weights.

Until this held-out evaluation is implemented and audited, use wording such as **“pre-posterior diagnostic features are recorded”**, not “the diagnostics predict posterior miscalibration.”

## 5. Calibration and numerical validation

Keep these targets distinct:

- **Fixed-truth coverage:** repeated-sampling coverage at the injected (H_0).
- **SBC:** ranks when (H_0) is drawn from the declared prior and data are generated under the matching model.
- **Posterior predictive checks:** comparison of replicated selected data with observed data under posterior draws.
- **Frequentist conditional calibration:** repeated-sampling performance conditional on fixed generating parameters and declared strata.

A posterior CDF evaluated at one fixed true (H_0) is not automatically uniform across repeated datasets. Do not label such a fixed-truth PIT as SBC.

Report interior and boundary truths separately if the design later varies (H_0). The current pilot's default (H_0=70) is an interior fixed-truth experiment; it cannot establish boundary behavior. Grid-resolution checks, fit-failure accounting, and Monte Carlo uncertainty must be completed before interpreting small differences.

## 6. Data provenance: remaining independent-audit gap

The current PR #15 smoke workflow passed on commit `82ef8b0f751716a7f088ff98acd99ea359efb997` for Python 3.10, 3.11, and 3.12. Its estimated-inclusion pilot smoke experiment and independent artifact audit steps succeeded. This establishes the tested software and smoke-artifact checks at that exact commit; it does **not** establish confirmatory performance or a full independent regeneration of every simulated parent population.

The current auditor independently refits the validation logistic model and reconstructs weighted KDE, posterior calculations, and summaries from persisted primitive inputs. The next provenance gate should additionally regenerate the parent population, catalogue-retention indicators, and validation-frame sample from the recorded component seeds in a separately implemented audit path. Verify the immutable key

[
(	ext{replicate},g,m,	ext{arm},	ext{component seed})
]

and compare regenerated arrays or canonical hashes against the persisted raw inputs. The scoring-only truth payload must remain separate from the estimator-visible inputs. This closes a different gap from posterior reconstruction: a posterior can be reconstructed perfectly from a dataset whose generation or truth labels were recorded incorrectly.

## 7. Staged execution gates

**Gate 1 — protocol freeze.** Freeze the estimand, sign convention, (g) and (m) levels, seed derivation, primary contrast, failure policy, and multiplicity/reporting plan.

**Gate 2 — deterministic and provenance tests.** Test (g=0), positive and negative (g), invalid inputs, positivity diagnostics, component-seed regeneration, and the non-oracle information firewall.

**Gate 3 — pilot.** Run the prespecified 20-replicate smoke/pilot cells. Inspect retained-catalogue counts, validation outcome counts, fit failures, boundary hits, probability ranges, and weight tails. The pilot is for pipeline and precision planning, not effect selection.

**Gate 4 — confirmatory precision.** Choose replicate counts from the pilot SD of the paired primary contrast and a declared confidence-interval half-width target. Do not treat 200 replicates as automatically sufficient. Amend the protocol before confirmation if the runtime or precision target requires a change.

**Gate 5 — independent reconstruction.** Rebuild generation, fits, density estimates, posterior summaries, paired contrasts, and coverage counts from raw artifacts using independent calculations. Verify checksums and report all failures.

**Gate 6 — pre-posterior prediction.** Only after the diagnostic feature and outcome definitions are frozen should the held-out predictive study begin.

No result advances past a gate merely because the unit tests pass.

## 8. Claim boundary

Permitted after the relevant gates pass: “In a controlled synthetic low-redshift model with an independently sampled parent-frame validation subset, we quantified how specified inclusion-model curvature and validation-sample sizes affected inverse-weight diagnostics, host-density reconstruction, and fixed-truth posterior performance.”

Not permitted from this project as currently implemented: “new dark-siren method,” “realistic LVK validation,” “new (H_0) measurement,” “universal ESS threshold,” “uncertainty in estimated inclusion probabilities is fully propagated,” or “pre-posterior diagnostics predict real-data posterior failure.”

**Immediate priority:** implement seed-level generation reconstruction in the independent auditor, then complete the prespecified nonlinear pilot and overlap review. Do not add more simulation volume until that provenance gate is closed.


## Provenance-audit extension: detected-event RNG replay (2026-10-31)

The independent provenance auditor now also replays the detected GW distance observations generated by the low-redshift toy simulator. It reproduces the simulator's RNG consumption order: parent-redshift uniforms, the generator's internal catalogue-thinning uniforms, host draws, and Gaussian distance noise with rejection at the observed-distance threshold. The replay is implemented independently in the audit script rather than importing the production generator.

The persisted `observed_distances_mpc` in each successful estimated-inclusion inference payload is compared against the independently regenerated sequence. The regression fixture now obtains its parent and event observations from the production generator, then checks that the independent auditor accepts the intact artifact and rejects a deliberately modified event distance.

**Scope:** this extends the seed-level audit to the detected event-distance observations. It does not yet independently reconstruct GW host IDs or the full event-truth table, nor does it prove the physical adequacy of the simplified detector-selection model. It is a deterministic software/provenance check, not a scientific validation of real detector selection. Exact-head CI and a real persisted pilot-artifact audit remain required.


## Extended event-truth provenance (2026-10-31)

The estimated-inclusion study now stores GW host object IDs, host redshifts, and latent mean distances only in the separate scoring/provenance payload. These fields are not passed to the inclusion estimator, host-density estimator, or cosmological likelihood. The inference payload continues to contain only observed event distances.

The independent auditor now replays the generator's random stream and checks, in order, the detected host IDs, host redshifts, latent mean distances, and observed detected distances. A regression test deliberately changes a stored host identity and requires the audit to reject the artifact.

This closes a previously stated provenance gap for successful estimated-inclusion rows. It remains a deterministic reconstruction of the declared toy generator, not validation of detector realism. Failed-only replicates still require careful treatment because they may not contain successful-arm audit payloads.


## Failed-arm provenance retention (2026-10-31)

Raw parent/catalogue/validation/event provenance is now written for estimated-inclusion arms even when the catalogue-size gate or a later estimator/posterior step fails. The independent auditor checks all estimated-arm rows, not only successful fits, and reports the number of failed estimated rows reconstructed. This prevents a replicate from disappearing from provenance review merely because its inferential fit failed.

The audit still requires every replicate to contain its expected estimated-arm rows and their provenance payloads. It does not silently treat missing payloads as a passing result.


## Paired linear-versus-quadratic estimator integration (2026-10-31)

The pilot runner now creates two estimated-inclusion arms at every declared validation size: `estimated_linear_<m>` and `estimated_quadratic_<m>`. Both arms use the same replicate-level parent population, catalogue-retention realization, detected-event observations, and nested validation sample. Only the fitted inclusion-model degree differs. Catalogue-only and oracle arms remain as reference benchmarks.

Each estimated row records its fitted model family and degree. The independent posterior-artifact auditor now independently refits the linear or quadratic logistic model according to the arm, checks the fitted coefficients and predicted inclusion probabilities, and reconstructs the weighted KDE and selection-conditioned posterior. The summary also records paired linear-minus-quadratic absolute-H0-error contrasts at each validation size; positive values favor the quadratic arm. The auditor independently reconstructs those contrasts from replicate-level records.

The runner accepts a nonzero `true_quadratic` coefficient to create a controlled curvature stress condition. The default remains zero, so the quadratic arm is then a deliberately more flexible nested fit rather than a correctly specified nonlinear model. A meaningful nonlinear stress comparison must explicitly set and report nonzero curvature and should vary curvature and validation size under a prespecified design; one pilot setting cannot establish general performance.

**Verification status:** source and regression-test changes are committed, but the latest exact-head GitHub Actions result has not been retrieved and no new pilot artifact has been executed or audited in this session. Therefore no new numerical result or test pass is claimed. The next gate is exact-head CI, followed by a small paired pilot and its independent artifact audit.


## Protocol reconciliation after paired-model implementation (2026-10-31)

Updated `docs/estimated_inclusion_nonlinear_stress_protocol_2026-10-30.md` to match the implemented linear/quadratic estimator arms and the actual nonlinear data-generating law. The revised protocol:

- makes the paired linear-minus-quadratic absolute posterior-mean error the primary model-family contrast, while retaining each arm's catalogue-only comparison as secondary;
- explicitly states that positive contrast values favor quadratic only for absolute point-estimation error, not coverage or calibration;
- marks curvature and validation-size levels as candidates pending overlap review and precision planning;
- records that paired contrasts for coverage, interval width, inclusion-probability MAE, and density ISE are not yet implemented and must not be inferred from the current absolute-error contrast;
- adds an explicit gate for independent regeneration of parent/catalogue/GW/validation inputs and distinguishes it from posterior reconstruction;
- documents the plug-in inclusion-probability uncertainty limitation and the synthetic toy detector-selection boundary.

**Verification update:** CLI inspection confirms `--true-quadratic` is parsed and forwarded to the runner. The nonlinear CI smoke assertions match the current summary schema (`true_inclusion_model.family`, `quadratic`, and the paired arm labels). The latest documentation commit is `ca3a4401863fa8b572ec2e9d0467674967a26213`; the workflow-run connector returned no run for that SHA. No CI pass, local test execution, pilot artifact, or independent artifact-audit pass is claimed.

**Next gate:** obtain a completed exact-head CI run. Then execute the small nonlinear smoke/pilot artifact and independently audit it before increasing replicate count or launching the factorial screen.
