# Central scientific and execution audit — 2026-11-18

## Purpose and scope

This is a coordination memo for the synthetic dark-siren catalogue-inclusion project. It records a read-only audit of the live repository and the current literature landscape. It does not report new simulation results, does not establish novelty, and does not authorize an H0 experiment.

The project remains limited to synthetic low-redshift methodology. It is not an LVK/GWTC analysis, a real galaxy-survey validation, or an observational H0 measurement.

## Repository snapshot inspected

Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`

Default branch: `main`. The audit did not modify `main`, merge a PR, or alter any existing branch.

Relevant live PR states at audit time:

| PR | Head branch | Head SHA | State | Audit interpretation |
|---:|---|---|---|---|
| #44 | `research/five-arm-comparison-preregistration-20261113` | `2abb4c5275235fd97f87c0818f59b4dcc1e3223c` | open, draft | Core five-arm design and implementation; engineering validation is not equivalent to confirmatory evidence. |
| #54 | `research/five-arm-kde-kernel-reuse-v2-20261114` | `da3b3970775997593f47a22a6dbf114816eb5d74` | open, draft | Density/KDE engineering and bandwidth work; exact-head CI was not returned by the inspected workflow-run/status endpoints. |
| #55 | `research/curvature-known-probability-fix-20261114` | `95e8164bd1c17756e5416566ec3a28475b15f1b9` | open, draft | Main curvature-density implementation branch; exact-head CI not verified; seven-gamma workflow not dispatched. |
| #56 | `fix/curvature-audit-json-output-20261115` | `a4aa454225fa59f6ed1e51f16374021f9623ccd2` | open, draft | A related curvature audit-output fix; its PR-triggered test workflow completed successfully (run 43634413804). |
| #57 | `research/density-bandwidth-selection-design-20261115` | `085efa5c13528295180f496cc560cba8aac85720` | open, draft | Design-only bandwidth protocol; its PR-triggered test workflow failed (run 43634450957). |

Important CI distinctions:

- PR #55's current head SHA returned no associated PR workflow runs and no combined status entries. Therefore its exact-head CI is **NOT VERIFIED**.
- PR #56's related test workflow run 43634413804 completed successfully on Python 3.10, 3.11 and 3.12. It exercised the curvature density smoke, pilot and replicate test modules on that PR's tested revision. This is useful evidence for that revision, but it is not proof that PR #55's later head passes.
- PR #57's test workflow run 43634450957 failed on Python 3.10, 3.11 and 3.12 during collection of `tests/test_quadratic_logit_curvature_density_replicates.py`. The logged failure is a Python `SyntaxError` at line 148 in a `child_manifest.write_text(...)` string literal (over-escaped backslashes). The preceding smoke and single-replicate pilot test groups passed; the replicate module did not collect, so later steps were skipped. This is a test-source syntax failure, not a scientific result.
- The seven-gamma manual workflow in PR #55 is configured with 7 gamma cells × 20 replicates (140 planned replicates), but the workflow is manual-only and has **NOT RUN** according to the live PR description and the calibration note. No density-performance result exists for this curvature family.

The PR bodies are not perfectly synchronized with their current heads; this memo uses the fetched live head SHAs above and does not treat stale narrative status as current CI evidence.

## What is implemented versus what is scientifically established

### Implemented infrastructure (not equivalent to scientific completion)

- Five distinct arms are defined: catalogue-only; full-parent oracle; known-probability IPW; MAP plug-in estimated IPW; and coefficient-uncertainty-integrated estimated IPW.
- The explicit true inclusion-probability input is used for Arm C in the non-default quadratic-logit DGP. D/E are fitted from observable validation redshifts and retention labels, rather than receiving the simulation truth vector.
- The absolute-bandwidth density runner and independent raw-input auditor reconstruct weighted KDE densities and diagnostics, including a separate reconstruction of the degree-one logistic posterior/quadrature in the absolute-bandwidth smoke auditor.
- The curvature design uses a parent-rate-calibrated quadratic-logit inclusion mechanism and a cross-gamma audit requiring complete replicate cells and matched keys for gamma-versus-zero contrasts.
- The workflow configuration includes a density-only experiment, with 11 absolute bandwidth candidates, 3 validation sizes, 7 gamma values, and 20 paired replicate seeds per gamma.

### Results that are not established

- No seven-gamma curvature run or cross-gamma raw-artifact audit has been verified.
- No confirmatory bandwidth has been selected. The lower-boundary signal from an earlier incomplete 19/20 pilot is exploratory and cannot freeze a bandwidth.
- No H0 likelihood, H0 posterior, posterior calibration, PIT, SBC, coverage, or event-level inference result is produced by the curvature density workflow.
- No empirical superiority claim for D or E is justified.
- No real survey, LVK/GWTC, or observational validation is included.

## Deterministic curvature calibration audit

The separate calibration note on the PR #55 branch reports:

- Parent measure proportional to
  `1 + 2.2 exp(-0.5*((z-0.105)/0.025)^2) + 1.6 exp(-0.5*((z-0.275)/0.040)^2)`
  on `z ∈ [0.01, 0.40]`.
- Coordinate `x=(z-0.205)/0.195`.
- Parent-measure centering `c_g=E_g[x^2]=0.2845005579191241`.
- Linear reference slope `beta=-3.9`.
- Reference retention `r0=E_g[expit(-3.9*x)]=0.5101919962291073`.
- For each gamma, the intercept is calibrated so `E_g[pi_gamma]=r0`.

The calibrated gamma grid is `{-2,-1,-0.5,0,0.5,1,2}`. The calibration arithmetic is a deterministic numerical check, not a stochastic simulation result.

A critical interpretation point: this gamma grid is symmetric in the coefficient value but **not** symmetric in practical IPW stress. Under the fixed negative linear slope and parent measure, negative curvature creates smaller minimum inclusion probabilities and more variable ideal inverse weights than positive curvature. The paper must report the actual inclusion-probability and weight distributions, not call the cells equally difficult. Holding the parent-average retention rate fixed does not hold positivity, minimum inclusion probability, or weight variance fixed.

## Literature positioning: novelty must be narrowed

The broad statement "we study galaxy-catalogue incompleteness in dark-siren cosmology" is not novel and must not be used as the contribution claim.

Relevant work located in this audit:

1. **Borghi et al. (2026), A&A 706, A199, "Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology"**, DOI 10.1051/0004-6361/202557354; arXiv:2509.18243. This presents a self-consistent incompleteness and host-weighting framework in CHIMERA, including mock BBH catalogues and varying completeness. It is a close astrophysical precedent and must be discussed directly.
2. **Leyde, Baker & Enzi (2026), JCAP 01 (2026) 013, "Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology"**, DOI 10.1088/1475-7516/2026/01/013. This reconstructs the underlying galaxy field and magnitude distribution while respecting large-scale spatial correlations, validated on Millennium-based simulations.
3. **Tagliazucchi et al. (2026), "Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs"**, arXiv:2608.24992; published as a Physical Review D paper (DOI 10.1103/vmzf-1wxq). It develops a GP-based reconstruction of incomplete line-of-sight galaxy distributions with an explicit survey-selection likelihood and reports mock-catalogue H0 precision gains.
4. **Gray, Williams & Papadopoulos (2026), "Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies"**, arXiv:2609.04991. This jointly infers host-galaxy population properties and cosmological parameters, and is directly relevant to unknown host weighting.
5. **"The luminosity of the unseen: luminosity function in dark sirens"**, MNRAS 546 (2026), DOI 10.1093/mnras/stag218, studies the luminosity-function treatment of missing hosts.
6. **"Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant"**, MNRAS 549 (2026), DOI 10.1093/mnras/stag828, proposes a completeness-limit test in the gwcosmo context.

The plausible differentiating question is therefore much narrower:

> In a controlled synthetic benchmark with an independently sampled validation frame, how do finite validation size, parametric inclusion-model misspecification, and inverse-probability weight instability interact to affect reconstruction of a declared host-weighted redshift density—and which pre-inference diagnostics identify failure?

Even this narrower framing is **a candidate gap, not a verified novelty claim**. A systematic literature search must still cover survey-sampling propensity estimation, estimated inverse-probability weights, Bayesian propensity-score uncertainty, positivity violations, and calibration/coverage diagnostics, not only astrophysics papers.

## Methodological risks and decisions to resolve before confirmatory execution

1. **Target estimand must be written explicitly.** The density runner's target is proportional to latent galaxy intensity times the declared merger/host weighting function. State whether the scientific target is a finite-parent empirical distribution or a superpopulation density. The interpretation of oracle arm B and all error metrics depends on this choice.
2. **Validation design must be explicit.** The validation frame is generated independently and a subset is sampled. Document whether the inferential target is the superpopulation inclusion function or the fixed finite validation frame. Do not add a generic scalar finite-population correction to the logistic posterior covariance without deriving it for this sampling design.
3. **Arm E uncertainty is conditional.** Integrating the degree-one coefficient posterior captures uncertainty conditional on the fitted degree-one model. It does not integrate over omitted quadratic curvature or model selection. State this limitation in the estimand and conclusions.
4. **IPW positivity is a first-class factor.** Report minimum and quantiles of true/estimated inclusion probabilities, normalized-weight tails, Kish and entropy ESS, and support gaps. Do not silently clip weights in the primary analysis.
5. **Separate density recovery from H0 calibration.** A density estimator with low IAE is not automatically a calibrated cosmological estimator. No H0 experiment should begin until density-stage design, bandwidth, and audit gates are complete.
6. **Bandwidth selection needs disjoint data.** The incomplete earlier 19/20 run and boundary-minimum signal cannot select a bandwidth. Training/holdout replicate blocks must be disjoint, predeclared, failure-inclusive, and independently audited.
7. **Monte Carlo precision is not yet demonstrated.** Twenty replicates per gamma are an engineering/pilot scale, not automatically enough for stable method comparisons or interval estimates. Plan replicate counts from a precision target for the paired primary contrast, accounting for failures and multiplicity.
8. **Failure accounting must remain complete.** A cell with any failed replicate must not silently summarize only successes. The current cross-gamma completeness gate is directionally correct, but the handling of non-evaluable versus failed runs and the scientific reason for a failed run must remain explicit.

## Recommended next sequence

1. Repair and verify the malformed replicate-test source on the branch where the failure occurred; do not infer that PR #55 is repaired by a change on another branch.
2. Run exact-head CI on PR #55's current head (or a new, explicitly tracked follow-up head) across Python 3.10/3.11/3.12. Record the exact SHA and run IDs.
3. Inspect the full single-replicate and multi-replicate raw artifact schemas and test the child and outer auditors with adversarial corruption cases.
4. Freeze the density estimand, primary contrast, bandwidth-training/holdout plan, failure policy, and Monte Carlo precision target before the full run.
5. Execute the density-only gamma pilot only after exact-head CI passes; audit every child artifact and the cross-gamma paired outputs.
6. Use those results to decide whether the density-stage effect is large and stable enough to justify a separate H0-likelihood/calibration programme.
7. Conduct the full literature review and formulate a precise contribution before drafting a journal manuscript.

## Audit verdict on 2026-11-18

- Repository and relevant PR metadata inspected: **SUCCESS**.
- Curvature calibration arithmetic recorded in the existing branch note: **SUCCESS as deterministic arithmetic only**.
- Related PR #56 test workflow: **SUCCESS on its tested revision**.
- Related PR #57 test workflow: **FAILED at test collection due to a syntax error**.
- PR #55 exact-head CI: **NOT VERIFIED**.
- Seven-gamma curvature simulation: **NOT RUN / NOT VERIFIED**.
- Cross-gamma audit: **NOT RUN / NOT VERIFIED**.
- Density-performance conclusions: **NOT EVALUABLE**.
- H0 inference/calibration: **NOT IN SCOPE / NOT RUN**.
- Novelty: **NOT ESTABLISHED**.

This memo is an audit artifact, not a claim of scientific completion.


## Independent high-resolution check of the curvature weight-stress arithmetic

A separate numerical reconstruction was run on 20,001-, 200,001-, and 2,000,001-point redshift grids. The parent-measure centering and reference retention converge to approximately:

- `c_g = 0.2845005570074`
- `r_0 = 0.5101919962132`

The 20,001-point values already recorded above differ only at approximately the 1e-9 level. This supports the numerical stability of the deterministic calibration constants; it does not validate the stochastic generator or the experiment.

For an iid catalogue draw from `f_cat(z)=g(z) pi(z)/r0`, with ideal unnormalized IPW weight `w(z)=1/pi(z)`, the population large-sample Kish effective-sample-size fraction is

\[
\frac{ESS}{n}\;\longrightarrow\;
\frac{\{E_{cat}(w)\}^2}{E_{cat}(w^2)}
=\frac{1}{r_0 E_g(1/\pi)}.
\]

This follows because `E_cat(w)=1/r0` and `E_cat(w^2)=E_g(1/pi)/r0`. It is a population asymptotic diagnostic, not the realized ESS of a finite catalogue.

| gamma | min true inclusion | population CV of ideal IPW weights | asymptotic Kish ESS fraction |
|---:|---:|---:|---:|
| -2.0 | 0.003790 | 2.3889 | 0.1491 |
| -1.0 | 0.008661 | 1.8198 | 0.2319 |
| -0.5 | 0.013116 | 1.6014 | 0.2805 |
| 0.0 | 0.019840 | 1.4192 | 0.3318 |
| 0.5 | 0.029912 | 1.2682 | 0.3834 |
| 1.0 | 0.044822 | 1.1443 | 0.4330 |
| 2.0 | 0.097281 | 0.9649 | 0.5178 |

The negative-curvature extreme therefore has an asymptotic ideal-IPW ESS fraction of about 15%, compared with about 52% at positive-curvature extreme, despite equal parent-average retention. This is a major design feature, not a small nuisance. It means the gamma comparison jointly stresses misspecification direction and positivity/weight variance. The planned analysis should either (a) explicitly present it as a joint curvature-plus-weight-stress family, or (b) add a separately labelled exploratory severity-matched design; it must not reinterpret the current symmetric gamma grid as severity-matched after seeing outcomes.


### Adjacent statistical-methodology literature (essential, not optional)

The novelty review must also include the mature survey-sampling and propensity-score literature. The current dark-siren framing cannot claim a new general solution to estimated inclusion probabilities, misspecification, positivity, or uncertainty propagation.

- **Kim & Riddles (2012), "Some theory for propensity-score-adjustment estimators in survey sampling," Survey Methodology 38(2), 157–165.** Directly addresses estimated selection/response probabilities, nonresponse and undercoverage, asymptotic properties, and variance estimation in finite-population survey sampling.
- **Katsumata (2025), "How should we estimate inverse probability weights with possibly misspecified propensity score models?" Political Science Research and Methods 13, 315–336, DOI 10.1017/psrm.2024.23.** Explicitly studies IPW under propensity-model misspecification and alternatives based on weight-distance minimization. The project must distinguish its scientific target and diagnostics from this prior methodology.
- **"Uncertainty in the design stage of two-stage Bayesian propensity score analysis" (2020), Statistics in Medicine, DOI 10.1002/sim.8486.** Develops posterior propagation of design-stage uncertainty into later inference. This is directly relevant to the rationale for Arm E, so the manuscript must explain what is different about the validation/inclusion problem and the downstream dark-siren density target.
- **"Measuring Model Misspecification: Application to Propensity Score Methods with Complex Survey Data" (2018), Computational Statistics & Data Analysis 128, 48–57, DOI 10.1016/j.csda.2018.05.003.** Studies diagnostics and consequences of propensity-model misspecification, including weighting.
- **"A Tutorial for Propensity Score Weighting Methods Under Violations of the Positivity Assumption" (2025), Statistics in Medicine, DOI 10.1002/sim.70329.** Reviews overlap/positivity violations, unstable weights, estimand choice, diagnostics, and sensitivity analysis.
- **Leedy & Kim (2026), "Data integration using propensity score weighting under missing at random," arXiv:2611.03353.** Very recent work on non-probability sample inclusion, IPW, model misspecification and doubly robust adjustment. It should be screened in the full review because it is close to the validation-plus-target-population structure.

Implication: the defensible contribution, if supported by results, is likely an application-specific controlled benchmark and diagnostic study that connects these established statistical phenomena to a clearly defined host-weighted redshift reconstruction and (only in a later stage) dark-siren likelihood calibration. The novelty is not "Bayesian IPW" by itself, nor the fact that IPW weights can be unstable.


## Additional artifact-level inspection: existing finite-validation probability pilot

I downloaded and inspected the numerical-diagnostics artifact attached to workflow run 43634451023 (artifact `numerical-diagnostics-py3.12`). Its finite-validation probability pilot is a separate, earlier experiment, not the new seven-gamma curvature-density experiment.

Artifact provenance recorded in its own run manifest:

- Git ref: `refs/pull/57/merge`
- Git SHA: `b45a540b54f49a5fd50a4c9e7597f6fc5d458d24`
- Python 3.12.16, NumPy 2.5.4, SciPy 1.18.1, pandas 3.1.1
- 4,000 parent units; five replicates; gamma `{-2,-1.2,0,1.2,2}`; retention targets `{0.1,0.25,0.5}`; validation sizes `{100,300,1000}`; logistic model degrees 1 and 2.
- 450 model-fit rows and 1,800,000 raw unit records.
- Independent artifact audit: **PASS**, reconstructing all 450 experiments and all 1.8 million raw unit records at tolerance 1e-10.

This pilot uses an earlier DGP specification and is not interchangeable with the current planned seven-gamma design. Its curvature-centering constant is 0.3114962521, its linear coefficient is -1, and its parent/validation design differs from the new density pilot. It evaluates held-out inclusion-probability prediction and weight diagnostics, not five-arm host-density recovery or H0 inference.

### Exploratory bias-variance signal in the earlier pilot

The degree-two model is correctly specified for the quadratic-logit DGP, while degree one is deliberately misspecified when gamma is nonzero. Pairing the degree-one and degree-two held-out probability RMSE within each replicate gives:

| Validation size | Mean RMSE difference (degree 1 minus degree 2), all gamma | SD of paired differences | Paired cells |
|---:|---:|---:|---:|
| 100 | -0.01267 | 0.02395 | 75 |
| 300 | -0.00906 | 0.01125 | 75 |
| 1000 | +0.00229 | 0.00667 | 75 |

Restricting to nonzero gamma (thus excluding the correctly specified linear gamma=0 case), the mean paired differences are -0.01123, -0.00776, and +0.00357 at validation sizes 100, 300, and 1000 respectively (60 paired cells at each size). Negative means indicate lower held-out RMSE for the simpler degree-one fit in this small pilot; the sign reversal at n=1000 is consistent with a finite-sample bias-variance trade-off. This is an exploratory signal only: there are five replicates per design cell, and no inferential claim or estimator recommendation follows.

Across all 450 fits, 340 (75.6%) report ordinary MLE convergence and 110 (24.4%) are accepted by the implementation's score-tolerance rule despite the optimizer message "Desired error not necessarily achieved due to precision loss." All are labelled SUCCESS in the artifact. This deserves explicit convergence-sensitivity analysis in the confirmatory design: report these separately, test alternative optimizers/tolerances, and do not equate the SUCCESS label with identical numerical convergence quality.

The probability-stage results are useful for designing the larger study because they show why a correctly specified but more flexible quadratic model may not dominate at small validation sizes. They do not answer whether estimated IPW improves host-density reconstruction; prediction-score improvements and weighted-density performance are distinct targets.

### Workflow status clarification

The focused PR #57 curvature-density test workflow run 43634450957 failed on all three Python versions during collection of the replicate test module; its earlier smoke and single-replicate pilot test groups passed. A separate broad numerical-diagnostics workflow run 43634451023 completed successfully and produced artifacts on the PR #57 merge ref. Its 3.12 job log contains successful tests for the baseline estimated-inclusion, five-arm, KDE, selected-likelihood, and other existing modules. This broad run is useful regression evidence for the merge ref but does **not** retroactively turn the focused curvature workflow into a pass, nor does it establish exact-head PR #55 CI. The current PR #55 head SHA still has no associated PR workflow runs returned by the inspected endpoint.


## Formal target estimand and interpretation of the five arms

Inspection of the actual PR #55 runner and `five_arm_comparison.py` resolves the target-density question at the implementation level.

Let the superpopulation galaxy intensity on `[z_min,z_max]` be

\[
\lambda_g(z)=1+2.2\exp\!\left[-\frac12\left(\frac{z-0.105}{0.025}\right)^2\right]
+1.6\exp\!\left[-\frac12\left(\frac{z-0.275}{0.040}\right)^2\right],
\]

and let the declared merger/host-rate weight be

\[
h(z)=1+1.5\frac{z-z_{\min}}{z_{\max}-z_{\min}}.
\]

The superpopulation host-redshift target used by the density runner is

\[
f_{\mathrm{host}}(z)=
\frac{\lambda_g(z)h(z)}
{\int_{z_{\min}}^{z_{\max}}\lambda_g(u)h(u)\,du}.
\]

The parent sample is drawn from the normalized `lambda_g`; catalogue inclusion is Bernoulli with probability `pi_gamma(z)`; the validation frame is independently drawn from the same parent redshift law; and validation is simple random sampling without replacement from that frame. The validation observations expose redshift and retention label, not the true inclusion probability.

For a catalogue draw with density proportional to `lambda_g(z) pi_gamma(z)`, the target is recovered by weighting observed catalogue objects by `h(z)/pi_gamma(z)`, subject to the usual positivity and correct-probability conditions. The five arms therefore have the following interpretations:

- A: observed catalogue weighted by `h(z)`, without an inclusion correction.
- B: finite complete-parent sample weighted by `h(z)`; an oracle Monte Carlo benchmark, not the exact superpopulation target.
- C: observed catalogue weighted by `h(z)/pi_gamma(z)`, using simulation-truth probabilities; oracle IPW diagnostic.
- D: observed catalogue weighted by `h(z)/pi_hat_MAP(z)`.
- E: posterior mixture of node-specific weighted density estimates using the degree-one logistic coefficient posterior from validation data.

The code's primary Arm E density loss is the loss of the posterior-mean density, i.e. IAE of the quadrature-weighted mean density. The posterior expected node-wise IAE is separately recorded as a secondary quantity. These are different estimands and must not be interchanged.

This establishes a **superpopulation** density estimand. Arm B uses only a finite parent realization and so has its own finite-parent Monte Carlo error. Comparisons should preserve this distinction: the oracle arm is a finite-sample reference, not a zero-error ground truth. The analytic target is the appropriate reference for absolute IAE/ISE scoring.

## Artifact-level inspection of the separate 100-replicate parent-density bandwidth study

The `five-arm-kde-bandwidth-pilot` artifact from run 43634451023 contains a separate parent-density bandwidth study. Its raw audit reports:

- 100 parent-sample replicates, 4,000 galaxies per replicate (400,000 parent units total).
- 17 candidate Scott multipliers and 1,700 bandwidth attempts.
- Independent reconstruction of every attempt, with maximum absolute IAE discrepancy (7.5\times10^{-16}).
- 70 training and 30 holdout replicate IDs, separated by a declared split seed.
- Selected Scott multiplier: **0.5**, with holdout mean IAE 0.043997 and holdout MCSE 0.001677.

This is genuine evidence that the parent-only bandwidth selection pipeline and its holdout audit work for the stated analytic parent-host density. It is **not** a bandwidth selection for the incomplete-catalogue IPW estimators: the study uses parent samples, and does not reproduce the catalogue-thinning, estimated-probability, or weight-instability problem. It must not be used to declare 0.5 optimal for Arms C/D/E. The current curvature workflow instead tests a fixed grid of absolute redshift bandwidths and does not select a bandwidth.

The same artifact reports a batched-vs-scalar five-arm likelihood kernel equivalence check: maximum absolute log-likelihood discrepancy (1.14\times10^{-13}), with a hardware-specific observed speed ratio about 20.3×. This supports numerical equivalence and engineering performance for the isolated kernel, not scientific accuracy or calibration.

## Further methodological implication

The existing finite-validation probability pilot suggests a potentially important bias-variance interaction: a correctly specified quadratic logistic model did not uniformly outperform the simpler linear model in held-out probability prediction at smaller validation sizes. The new five-arm density experiment must therefore not treat "correctly specified model" as automatically optimal. It must measure both probability prediction and the downstream density estimand, and distinguish:
- probability-score performance on held-out validation records;
- stability of the induced inverse weights on the catalogue;
- density recovery against the analytic superpopulation host density;
- eventual H0 posterior bias and calibration.

These are related but non-equivalent endpoints.


## Proposed confirmatory protocol (draft; not yet frozen)

The existing 20-replicate-per-gamma workflow is appropriately labelled a development pilot. It should not be treated as the final confirmatory study or as enough information to estimate small method differences precisely.

### Primary scientific endpoints

1. **Primary density endpoint:** integrated absolute error of the Arm D MAP plug-in estimated-IPW density against the analytic superpopulation host density, evaluated at a common predeclared absolute bandwidth and on the fixed redshift grid.
2. **Primary method contrast:** paired `D - A` IAE difference within replicate, validation size, gamma, and bandwidth. Negative values favor estimated-IPW correction over the uncorrected catalogue-only baseline.
3. **Key diagnostic contrast:** paired `D - C` IAE difference. This quantifies the additional penalty from estimated/misspecified inclusion probabilities relative to the simulation-oracle known-probability IPW arm.
4. **Uncertainty-propagation contrast:** paired `E - D` IAE difference, with E defined as the posterior-mean density. The posterior expected node-wise loss remains secondary and must be reported separately.
5. **Oracle reference:** B is compared descriptively with the analytic target and other arms, but is not used as a real-world competitor.

All contrasts should be reported with paired replicate-level uncertainty intervals. A cell with failed/non-evaluable replicates must retain its full planned denominator and be labelled incomplete; never silently drop failed runs.

### Bandwidth handling

The current 11 absolute bandwidths should first be treated as a sensitivity surface, not as 11 opportunities to select whichever value makes a method look best.

A defensible two-stage option is:
- Generate at least 100 independent replicates per gamma (or enough to support a stable split), with a predeclared 70/30 training/holdout partition by replicate ID.
- Use training replicates only to select one common absolute bandwidth under a predeclared pooled criterion, shared by all five arms and all gamma values.
- Lock that bandwidth and report the untouched holdout performance as the confirmatory result.
- Also publish the full holdout IAE-versus-bandwidth curves and weight diagnostics. Do not reselect on holdout.
- If the pooled criterion or common-bandwidth policy is not acceptable scientifically, retain the entire sensitivity surface as the primary result and do not name a single winning bandwidth.

The earlier 100-replicate parent-only Scott-multiplier study is not a substitute for this catalogue/IPW-specific selection protocol.

### Replicate count and Monte Carlo precision

Twenty replicates per gamma can expose engineering failures and large effects, but do not justify narrow inferential intervals. A practical next stage is 100 replicates per gamma to provide a first independent train/holdout split; the final replicate count should then be increased if the holdout MCSE of the primary paired contrast exceeds a predeclared precision tolerance. The tolerance must be expressed in IAE units and justified relative to the smallest scientifically meaningful density difference. Estimate it from the pilot's paired contrasts, not from marginal standard deviations alone.

### Pre-inference diagnostic gate

Before any H0 propagation, require:
- all seven gamma cells complete and every raw child audit passing;
- independent cross-gamma paired-seed and validation-frame checks passing;
- no unreported failed/non-evaluable attempts;
- full weight-tail and positivity summaries by arm, gamma, and validation size;
- bandwidth selection/holdout audit passing (or no bandwidth selection claim);
- grid-refinement and KDE numerical stability within declared tolerances;
- a clear statement of whether the primary result is robust across the full bandwidth sensitivity surface.

Only after this gate should a separate H0 experiment be designed. That experiment must then have its own event-generation design, GW selection normalization audit, posterior-grid convergence, fixed-truth coverage, and SBC/PIT calibration protocol.

### Claims the paper must avoid

- Do not claim degree-one posterior integration repairs omitted quadratic curvature.
- Do not claim a correct inclusion model is always preferable at small validation sizes; the earlier probability-stage pilot suggests a variance penalty for the flexible model.
- Do not infer H0 calibration from density IAE alone.
- Do not interpret equal parent-average inclusion rates as equal positivity or IPW difficulty.
- Do not call the oracle arm a deployable estimator.
- Do not generalize from this Euclidean synthetic DGP to real surveys or LVK catalogues.

This protocol is a research recommendation, not a registered or frozen analysis plan. It should be versioned and approved before the confirmatory run.


## Critical estimand audit and corrective patch — 2026-11-18

A further source-level audit found a consequential double application of the host/merger weighting function in the absolute-bandwidth density runner.

### Finding

The declared target density is

\[
f_{\mathrm{host}}(z) \propto g(z)h(z),
\qquad
h(z)=1+1.5\frac{z-z_{\min}}{z_{\max}-z_{\min}}.
\]

Before the correction, `run_absolute_bandwidth_density_smoke.py` constructed
`parent_base_weights = h(z_parent)` and `catalogue_base_weights = h(z_catalogue)`
through `build_five_arm_weights`. It then called
`estimate_host_density(..., sample_weights=weights, merger_weight=h(grid), ...)`.
The production estimator multiplies the weighted KDE by `merger_weight` after evaluating the KDE. Therefore the sample-level h and grid-level h were both applied.

For the oracle parent arm B, the large-sample target was consequently proportional to
\[
g(z)h(z)^2,
\]
rather than the declared \(g(z)h(z)\). With correctly specified inclusion probabilities, the same double weighting affected arm C; the intended limiting target for D/E was also wrong whenever their fitted probabilities recovered the true inclusion function. The catalogue-only arm A had the additional catalogue-inclusion factor, but it too received the duplicated host weighting.

The independent smoke auditor reproduced the same construction: it expected h-weighted sample weights and multiplied its independent KDE by h again. Thus its PASS status established implementation/auditor agreement, but did **not** establish agreement with the declared estimand. This is an important distinction: an independent implementation can reproduce the same conceptual error unless the target is audited separately.

### Independent numerical size of the mismatch

On a 200,001-point grid over [0.01, 0.40], I normalized both candidate superpopulation targets:
- intended: \(g(z)h(z)\);
- accidentally implemented: \(g(z)h(z)^2\).

The resulting total-variation distance was approximately **0.09231** (equivalently, the L1 distance was **0.18462**). The mean redshift shifted from **0.22519** under the intended target to **0.24691** under the duplicated target, a difference of **0.02172**. This is large enough to invalidate density-recovery comparisons from the pre-fix implementation; it is not a negligible numerical correction.

### Corrective change made on the existing research branch

The correction preserves the explicit target convention that h(z) is applied once on the output grid:
- The runner now supplies unit parent and catalogue base weights to the five-arm weight builder.
- The density estimator continues to multiply the output-grid KDE by h(z), exactly once.
- The independent auditor now independently expects unit base weights and verifies the persisted base-weight arrays are all one.
- The smoke tests now assert unit base weights for both the legacy floor-logistic and quadratic-logit DGPs.

This choice also makes the weight-concentration diagnostics reflect catalogue-inclusion/IPW weights rather than the product of inclusion weights and h(z). The fixed absolute-bandwidth comparisons remain the primary density sensitivity surface; the ancillary Scott-width diagnostics are now computed from the pure inclusion weights.

Commits on `research/curvature-known-probability-fix-20261114`:
- `561177087ec894d5058eb9f087d39e1456fcb0d3` — runner uses unit base weights.
- `c3ce0b500551340ded56095a6df50fbd6a7b3253` — independent auditor checks unit base weights.
- `ecea98a77cabd5473dbb438f475d55f05ee59d2e` — regression tests for duplicate host-rate weighting.

The live PR #55 head was refreshed to `ecea98a77cabd5473dbb438f475d55f05ee59d2e` at the time of this addendum. No PR was merged and no change was made to `main`.

### Additional posterior-quadrature spot check

I independently reimplemented the degree-one logistic posterior and the corrected Laplace-proposal tensor Gauss-Hermite importance weights, using the declared parent intensity, calibrated gamma values, validation-frame construction, fixed coordinate system, and Gaussian prior. For one fixed validation frame of 1,500 records (seed 20261216), I compared order 5 with order 19 for all seven gamma values and validation sizes 100, 300, and 1,000.

Across these 21 cases, the largest order-5 versus order-19 differences were:
- posterior coefficient mean, maximum absolute component difference: about **0.00195** at n=100;
- posterior covariance, maximum absolute component difference: about **0.00391** at n=100;
- posterior-mean inclusion probability on a 301-point redshift grid: at most **8.3×10⁻⁵** at n=100, **7×10⁻⁶** at n=300, and **1×10⁻⁶** at n=1,000.

For selected extreme and reference cells (gamma = -2, 0, +2; n = 100 and 1,000), increasing quadrature order from 19 to 23 changed posterior predictive probabilities by at most approximately **1.1×10⁻¹⁰**; coefficient means, covariance, and log evidence were also stable to the displayed numerical precision. This is a targeted independent numerical spot check, not a test of every possible validation realization, not an execution of the repository's test suite, and not a validation of the degree-two model or any downstream H0 likelihood.

### Current execution status after the correction

- Double host-rate weighting identified from production and independent-auditor source: **CONFIRMED**.
- Corrective runner/auditor/test edits committed to PR #55 research branch: **DONE**.
- Current PR #55 head: `ecea98a77cabd5473dbb438f475d55f05ee59d2e`.
- Exact-head GitHub Actions CI for this head: **NOT VERIFIED**; the inspected commit-workflow/status endpoints returned no run/status entries.
- Local execution of the repository test suite: **NOT RUN**; a local checkout was not available in the connected runtime.
- Seven-gamma density workflow: **NOT DISPATCHED**; the connected GitHub tool set does not expose a workflow-dispatch action.
- Seven-gamma simulation and cross-gamma raw-artifact audit: **NOT RUN**.
- H0 inference/calibration: **NOT RUN**.

The pre-fix density artifacts, if any, must not be reused as evidence for the corrected estimand. After exact-head CI is verified, the engineering smoke must be regenerated and independently audited, followed by the full density-only pilot. No H0 work is authorized by this correction alone.


### Independent finite-sample sanity experiment after the weighting correction

I also ran a separate 50-replicate synthetic check for gamma = -2, 0, +2. Each replicate used N_parent=4,000, the same parent redshifts and retention uniforms across the three gamma cells, the calibrated parent-average inclusion rate, reflected Gaussian KDE, a fixed physical width of 0.016, and a 401-point redshift grid. The intended target was normalized g(z)h(z). The corrected construction used unit galaxy-level base weights (or 1/pi for Arm C) and applied h(z) once on the output grid. The pre-fix construction was reconstructed separately using h(z) in the galaxy-level weights and h(z) again on the output grid.

Mean IAE across the 50 replicates:

| gamma | mean catalogue size | Parent oracle B, corrected | Parent oracle B, duplicated h | Known-probability C, corrected | Known-probability C, duplicated h |
|---:|---:|---:|---:|---:|---:|
| -2 | 2043.46 | 0.0550 | 0.1813 | 0.1459 | 0.2544 |
| 0 | 2040.74 | 0.0550 | 0.1813 | 0.1006 | 0.2043 |
| +2 | 2041.46 | 0.0550 | 0.1813 | 0.0834 | 0.1957 |

This is a stand-alone Monte Carlo sanity check, not an execution of the patched repository runner or its tests. It shows that the corrected oracle constructions move toward the declared target under this finite design, while the duplicated-host-weight construction has substantially larger IAE. The known-probability IPW arm remains noisier than the parent oracle because it uses the thinned catalogue and can have unstable weights, especially at gamma=-2. These results are diagnostic only and must not be used as confirmatory evidence or as a claim about D/E.


## Exact-head CI and corrected density-pilot launch — 2026-11-18

### Corrected implementation CI

The final correction commit `ecea98a77cabd5473dbb438f475d55f05ee59d2e` received two successful GitHub Actions runs:

1. **Focused curvature regressions:** run `44362902183`, workflow `curvature-known-probability-fix`, event `push`, conclusion **success**. All three matrix jobs (Python 3.10, 3.11, 3.12) completed successfully. The logged test groups included:
   - five-arm weights, explicit truth probabilities, entropy ESS: 18 passed;
   - absolute-bandwidth density smoke and independent raw-input auditor: 8 passed;
   - multi-replicate density runner and auditor: 5 passed;
   - synthetic quadratic-logit DGP invariants: 8 passed;
   - cross-gamma completeness gate: 5 passed.
   The counts above are from the Python 3.10 log; the Python 3.11 and 3.12 jobs also completed with success and matching test-group success.

2. **Full repository numerical-diagnostics workflow:** run `44362902057`, workflow `tests`, event `push`, conclusion **success**. Python 3.10, 3.11 and 3.12 matrix jobs all completed successfully on the same code head. The Python 3.12 job completed the full configured sequence, including the relevant density, five-arm, posterior, calibration, and independent-audit test groups. Some explicitly marked large optional experiments were skipped by workflow conditions; this is not a claim that those skipped experiments ran.

This resolves the earlier exact-head CI uncertainty for the corrected runner/auditor/tests at `ecea98a`. The focused and full test runs are engineering validation, not evidence of estimator performance.

### Density-only pilot was launched

To execute the already configured seven-gamma pilot, I added a branch-specific `push` trigger to `.github/workflows/quadratic-logit-curvature-density-pilot.yml` while retaining `workflow_dispatch`. This change is confined to the research branch and does not alter `main`.

- Workflow commit: `35da06224585849839f3dc458f73a550f107a0a6`.
- New PR #55 head: `35da06224585849839f3dc458f73a550f107a0a6`.
- Pilot run: `44366440872`, URL: https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44366440872
- Design: seven gamma cells, 20 paired-seed replicates per cell (140 planned), N_parent=4,000, validation frame=1,500, validation sizes 100/300/1,000, 11 absolute bandwidths, 401-point redshift grid, quadrature order 5.
- At the latest status check, all seven gamma jobs had passed checkout, Python setup, and dependency installation and were executing the multi-replicate density-only runner. No gamma-cell artifact or cross-gamma result had yet been returned.

A second full test workflow was automatically queued for the workflow-trigger commit (`44366440973`); it was still running at the latest check. The prior full test workflow on the corrected code head had already passed.

### Interpretation and limits

The successful CI means the corrected density runner, independent auditor, and regression tests execute successfully on all supported Python versions in the configured test workflow. It does not mean the seven-gamma simulation has completed. The 20-replicate design remains a development pilot, not a confirmatory result or a bandwidth-selection dataset. No H0 likelihood, posterior, PIT, SBC, or coverage experiment is included in this workflow.
