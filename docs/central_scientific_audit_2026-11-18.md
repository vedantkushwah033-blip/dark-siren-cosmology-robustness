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
