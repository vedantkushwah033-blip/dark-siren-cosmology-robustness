# Scientific readiness audit — 2026-11-18

## Scope and evidence boundary

This is a focused audit of the active quadratic-logit / five-arm density-recovery line, not a claim that every repository branch, historical experiment, or local working-tree file has been audited. The inspected source is the branch `research/curvature-known-probability-fix-20261114`, whose head was verified identical to commit `95e8164bd1c17756e5416566ec3a28475b15f1b9` at the time of inspection.

Inspected:
- PR #55 description and its validation/status disclosures.
- PR #54 and #57 descriptions for bandwidth-pilot and holdout-protocol dependencies.
- The seven-cell manual workflow `.github/workflows/quadratic-logit-curvature-density-pilot.yml`.
- The multi-replicate runner, outer auditor, single-replicate runner/auditor, cross-gamma auditor, generator, five-arm weight construction, and focused orchestration tests.
- The deterministic calibration note committed as `docs/independent_quadratic_logit_calibration_check_2026-11-16.md`.

No repository Python tests or simulation workflows were executed in this audit. GitHub Actions dispatch is not exposed by the available connected GitHub actions. The seven-gamma workflow is manual-only and the inspected head has no associated PR workflow runs or combined status checks. Do not infer a CI pass from earlier commits.

## Verified implementation and design

The curvature workflow declares seven gamma cells `{-2,-1,-0.5,0,0.5,1,2}`, 20 replicates per cell, parent size 4,000, validation frame 1,500, nested validation sizes 100/300/1,000, degree-one logistic fitting for Arms D/E, quadrature order 5, a 401-point redshift grid, and 11 absolute KDE bandwidths from 0.006 to 0.080. The design uses common seed schedules across gamma values and an outer auditor that checks paired parent/validation redshifts, calibration constants, complete replicate cells, child audits, and paired metric keys.

The five arms in the density runner are:
A. catalogue-only;
B. full-parent oracle;
C. known-probability IPW (simulation truth);
D. degree-one logistic MAP plug-in IPW;
E. degree-one logistic posterior-integrated density.

The generator keeps catalogue inclusion separate from the toy GW event-detection selection model. The current curvature workflow is explicitly density-only. Its own configuration and PR description prohibit H0 posterior, PIT, coverage, bandwidth-selection, estimator-superiority, observational, or LVK/GWTC claims.

The code path uses the exact inclusion probabilities only for Arm C. Arms D/E are fitted from the observable validation redshifts and binary retention labels. The validation frame is generated independently of the analysis parent population and validation sample sizes are nested prefixes selected without replacement. Across gamma cells, common random-number streams intentionally pair redshifts and latent uniforms; retention outcomes can differ as the inclusion model changes.

## Independent deterministic calculation

I independently recomputed the parent-measure centering and rate calibration from the stated parent intensity on a 20,001-point redshift grid using trapezoidal integration and scalar root finding. Results reproduce the committed note:

- (c_g=E_g[x^2]=0.2845005579191241).
- Reference parent-average retention rate (r_0=0.5101919962291073).
- Calibrated intercepts: gamma -2: -0.24047847028208044; -1: -0.12472169347234588; -0.5: -0.06297071445790560; 0: 0; +0.5: 0.06312269530994806; +1: 0.12530550566919094; +2: 0.24244138342969593.

The practical IPW stress is asymmetric despite the symmetric gamma grid. The minimum true inclusion probability ranges from about 0.00379 at gamma -2 to 0.09728 at gamma +2. The population CV of ideal inverse-probability weights under the retained-catalogue law ranges from about 2.389 to 0.965. These are deterministic population calculations, not finite-sample ESS or pilot outcomes. Thus the experiment varies curvature direction and positivity/weight instability severity together; results must not call the gamma grid symmetric in difficulty.

## Current status matrix

| Component | Status | Evidence / limitation |
|---|---|---|
| Quadratic-logit parent-measure calibration | Independently recalculated | Deterministic arithmetic only |
| Five-arm density runner and raw-input audit | Implemented | Source inspection; no full execution in this audit |
| Multi-replicate lifecycle and outer audit | Implemented | Source/tests inspected; no full execution in this audit |
| Seven-gamma workflow | Configured, NOT RUN / NOT VERIFIED | Manual dispatch only; no run observed |
| Exact-head CI for commit 95e8164 | NOT VERIFIED | No associated PR workflow runs or status checks returned |
| Curvature density results | NOT AVAILABLE | No seven-cell artifacts |
| Bandwidth selection | NOT AUTHORIZED | Current curvature run sweeps candidate widths; no train/holdout selection |
| H0 likelihood/posterior/calibration | NOT RUN in this curvature line | Workflow is explicitly density-only |
| SBC / repeated-sampling H0 calibration for this line | NOT RUN | No H0 stage is authorized by current density artifacts |
| Real galaxy survey / LVK/GWTC validation | OUT OF SCOPE | Synthetic low-z methodology only |

PR #55 remains open/draft in the inspected search results. PR #56 is a serialization fix with CI requested but no run evidence in its description. PR #57 is a design-only bandwidth holdout protocol. PR #54 records earlier density-pilot problems and explicitly says the earlier incomplete pilot cannot authorize a bandwidth or H0 stage.

## Scientific interpretation and risks

1. **Arm E does not solve model misspecification.** It integrates uncertainty in coefficients conditional on the degree-one logistic model. At nonzero gamma, the true inclusion logit is quadratic while the fitted model remains linear. Posterior uncertainty around a misspecified pseudo-true projection is not uncertainty over the omitted curvature. Any interpretation must distinguish coefficient uncertainty from model uncertainty.

2. **Density recovery is not H0 calibration.** A low density IAE does not establish calibrated H0 posteriors. The downstream selection-conditioned event likelihood, event-level sampling, H0 grid/prior behavior, repeated-sampling coverage, and SBC require separate gates.

3. **Twenty replicates per gamma are an engineering pilot, not a confirmatory sample.** Paired t intervals from 20 replicates are exploratory and cannot support fine-grained ranking across 11 bandwidths, three validation sizes, multiple arms, and seven gamma cells. A confirmatory design needs predeclared primary estimands, Monte Carlo precision targets, multiplicity strategy, failure-inclusive bounds, and independent audit.

4. **Bandwidth is not frozen.** The current workflow evaluates an 11-point absolute-width grid on the same pilot replicates. The observed minimum cannot be treated as a selected bandwidth without disjoint training and holdout replicate blocks and a predeclared selection rule. A lower- or upper-grid-boundary minimum should trigger grid expansion under the protocol, not a post hoc choice.

5. **The curvature stress mixes two dimensions.** Negative gamma produces substantially more severe positivity/weight instability than positive gamma. Either interpret gamma contrasts as joint curvature-plus-overlap stress, or design a separately labelled severity-matched exploratory family. Do not change the current grid after looking at outcomes.

6. **The oracle arms are diagnostics.** Parent-oracle and known-probability IPW use information unavailable in practice. They benchmark attainable behavior under the simulation truth; they are not deployable methods.

7. **Finite validation sampling is a superpopulation approximation.** The validation frame is finite and the nested samples are drawn without replacement. The current model likelihood treats the sampled binary labels as Bernoulli observations conditional on redshift. This is coherent under the stated iid superpopulation construction, but the manuscript should state the sampling design and distinguish conditional-on-frame from superpopulation uncertainty. If a design-based finite-population estimand is intended, it needs separate treatment.

## Literature positioning update (not a completed systematic review)

Important adjacent work already published by 18 November 2026:

- Borghi et al., *Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology*, A&A 706, A199 (2026), develops a self-consistent incompleteness and host-weighting framework in CHIMERA and studies H0 implications under mock catalogues. DOI: https://doi.org/10.1051/0004-6361/202557354
- Tagliazucchi et al., *Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs*, Phys. Rev. D 114, 083523 (published 13 October 2026), fits a hierarchical selection-aware GP model to reconstruct the line-of-sight galaxy redshift distribution, including mock catalogues to z <= 0.4. DOI: https://doi.org/10.1103/vmzf-1wxq ; arXiv: https://arxiv.org/abs/2608.24992
- Leyde, Baker & Enzi, *Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology*, JCAP 01 (2026) 013, reconstructs the underlying galaxy field and magnitude distribution with spatial correlation information. DOI: https://doi.org/10.1088/1475-7516/2026/01/013

These works make broad claims about “correcting catalogue incompleteness” or “reconstructing missing galaxies” unsuitable as novelty claims. The narrower potential distinction here is a controlled statistical study of finite validation samples used to estimate a redshift-dependent catalogue-inclusion propensity, with explicit comparison of known-probability IPW, plug-in estimated IPW, and coefficient-posterior-integrated IPW under positivity stress and a deliberately misspecified inclusion model. This is a candidate gap, not a verified novelty claim; a systematic search must cover propensity-score, survey-sampling, missing-data, and gravitational-wave literatures.

## Publication-pathway assessment

- **Astrophysics pathway:** Current synthetic low-z Euclidean density-only evidence is not enough for a strong JCAP/MNRAS/ApJ/A&A claim. It would need the completed density audit, an explicitly validated H0 propagation stage, repeated-sampling calibration/SBC, a physically defensible cosmological distance relation and GW selection model, and careful comparison with existing completeness approaches. A&A/CHIMERA and the 2026 GP paper are especially close context.
- **Statistics pathway:** Current work applies established IPW and Bayesian logistic coefficient integration; no new estimator or general theory has yet been established. A strong statistics submission would need a general methodological contribution, such as a principled uncertainty procedure or diagnostic/decision rule with theoretical properties and broad simulation evidence, not only a dark-siren application.

## Recommended next gates (in order)

1. Verify exact-head CI on the current PR #55 head; if unavailable, run the documented test matrix in a reproducible environment and preserve logs.
2. Run the seven-gamma density workflow only after CI passes; preserve all raw artifacts and failure ledgers.
3. Independently audit every child, the outer manifests, the cross-gamma pairing, and all aggregate metrics. Do not report results until every planned gamma cell has 20/20 successful replicates and the independent audit passes.
4. Keep the 11-bandwidth sweep exploratory. Design a curvature-specific training/holdout protocol with disjoint replicate IDs and frozen equal-cell objectives before bandwidth selection.
5. Add a correctly specified quadratic-logit fitted-model diagnostic arm (or a separately prespecified model-degree comparison) to distinguish finite-sample uncertainty from structural model misspecification. Keep it separate from the five-arm pilot unless a new protocol and tests are approved.
6. Only after density-stage gates pass, implement H0 propagation with catalogue inclusion kept separate from GW detection selection, then run grid/prior convergence, selection-likelihood audits, fixed-truth coverage, and SBC.
7. Complete a systematic novelty review before choosing the final journal. Preserve one core research programme but do not assume one manuscript is simultaneously ready for a top statistics journal and a top astrophysics journal.

## Reproducibility note

The independent numerical calculation in this audit is a fresh deterministic calculation from the declared analytic DGP. It is not a simulation replicate, a repository test run, a GitHub Actions result, or an empirical density-recovery finding.
