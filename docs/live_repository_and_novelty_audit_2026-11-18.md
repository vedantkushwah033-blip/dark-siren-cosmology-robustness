# Live repository and novelty audit — 2026-11-18

## Purpose and scope

This note records a read-only inspection of the GitHub repository and selected current research literature on 18 November 2026. It is an audit checkpoint, not a claim that the seven-gamma experiment has run. The project remains a synthetic low-redshift Euclidean methodology study; no real survey, LVK/GWTC, or observational H0 result is claimed.

Repository: `vedantkushwah033-blip/dark-siren-cosmology-robustness`

## 1. Live repository state observed

- Default branch: `main`; repository is public.
- PR #55, “Implement audited quadratic-logit density pilot infrastructure”, is open and draft. Head: `research/curvature-known-probability-fix-20261114` at `95e8164bd1c17756e5416566ec3a28475b15f1b9`; base: `research/absolute-width-distribution-pilot-20261114` at `43797ff676a8e15080fbcec057960c0095dbaf48`. It changes 15 files.
- PR #54, five-arm KDE/bandwidth work, is open and draft. Its current head is `da3b3970775997593f47a22a6dbf114816eb5d74` on `research/five-arm-kde-kernel-reuse-v2-20261114`. Its PR description contains older status statements; later comments and branch documents should be reconciled before relying on the description alone.
- PR #56, the curvature-auditor JSON-output fix, is open and draft at `a4aa454225fa59f6ed1e51f16374021f9623ccd2`. Its description says CI was requested but does not claim a pass.
- PR #57, density bandwidth holdout protocol, is open and draft at `085efa5c13528295180f496cc560cba8aac85720`. It is design-only; it selects no bandwidth and opens no holdout.
- PR #44, the matched five-arm comparison, is open and draft at `2abb4c5275235fd97f87c0818f59b4dcc1e3223c`. The GitHub workflow-run lookup returned a completed successful test run (ID `43582953042`) associated with this commit. This validates the reported test workflow on that PR head; it does not validate the separate seven-gamma experiment.

The PRs are on multiple related but distinct branches. They must not be treated as one integrated, tested release. In particular, successful H0-likelihood tests on PR #44 do not establish that the curvature-density branch has passed exact-head integration tests.

## 2. Curvature experiment: what is implemented versus what is known

The PR #55 branch contains:
- a parent-measure-calibrated quadratic-logit inclusion DGP;
- a five-arm weight builder in which Arm C can receive the exact simulation inclusion probabilities while D/E are fitted from validation observations;
- single-replicate and multi-replicate density runners;
- child and outer raw-artifact auditors;
- a cross-gamma auditor checking common calibration, paired seeds and shared parent/validation redshift samples;
- paired same-arm curvature contrasts against gamma zero;
- a manual-only workflow at `.github/workflows/quadratic-logit-curvature-density-pilot.yml`.

The workflow specifies seven gamma cells, 20 replicates per cell, N=4,000, a validation frame of 1,500, nested validation sizes 100/300/1,000, and an 11-value absolute bandwidth grid. It runs density reconstruction and audits, not H0 inference.

**Execution status:** the workflow file is configured for `workflow_dispatch`; no run was found in the inspected repository records, and the PR description explicitly says it has not been dispatched. The exact PR #55 head had no PR-triggered workflow runs returned by the commit-run lookup. Therefore exact-head CI is **NOT VERIFIED** and the seven-gamma experiment is **NOT RUN / NOT VERIFIED**. No density-performance conclusion from this curvature experiment is currently evaluable.

The cross-gamma gate correctly refuses paired summaries unless all planned replicates succeed in all seven cells. This avoids silently changing the paired estimand after failures. Preserve every failed attempt and report it; do not relax the gate by dropping failed replicates.

## 3. Deterministic design audit

The stated DGP uses the parent intensity
[
g(z) propto 1+2.2e^{-\frac12((z-0.105)/0.025)^2}
+1.6e^{-\frac12((z-0.275)/0.040)^2}
]
on (z\in[0.01,0.40]), with (x=(z-0.205)/0.195), and
[
\pi_\gamma(z)=\operatorname{expit}\{\alpha_\gamma-3.9x+\gamma(x^2-c_g)\}.
]
The branch documentation reports (c_g=E_g[x^2]=0.2845005579) and a common parent-average inclusion rate (r_0=0.5101919962), with each intercept calibrated to that rate. A separate deterministic trapezoidal/root-solve calculation is documented in `docs/independent_quadratic_logit_calibration_check_2026-11-16.md`.

The same note correctly warns that symmetric gamma values are **not** symmetric in positivity or weight stress. The reported catalogue-law coefficient of variation of ideal (1/\pi) weights ranges from about 2.389 at gamma=-2 to 0.9649 at gamma=+2. These are population calculations, not observed finite-sample ESS or pilot outcomes.

Important interpretation:
- gamma=0 is correctly specified by the degree-one logistic analysis model;
- nonzero gamma values create structural functional-form misspecification;
- intercept calibration holds the marginal parent retention rate fixed, but does not hold the minimum inclusion probability, local slope, or inverse-weight tail behaviour fixed;
- Arm C is an oracle benchmark, not a practically available estimator;
- Arm E integrates coefficient uncertainty conditional on the degree-one model and does not integrate uncertainty about the omitted quadratic term.

Thus curvature direction and practical positivity stress vary together in this design. Report both the gamma value and the actual inclusion/weight profile. Do not describe the seven cells as equally difficult or as a pure one-dimensional “curvature severity” scale.

## 4. Experimental interpretation and gates

The current curvature workflow is a **density-only exploratory pilot**:
- no bandwidth is selected;
- no disjoint training/holdout bandwidth selection is performed;
- no H0 likelihood, posterior, PIT, coverage, or SBC result is produced;
- no claim of estimator superiority is supported by the workflow configuration alone.

Twenty replicates per gamma are useful for debugging and initial effect-size exploration, but may be too imprecise for fine ranking, tail diagnostics, rare failures, or small paired effects. Report paired Monte Carlo standard errors and intervals for every key contrast; avoid binary conclusions from point estimates. If a minimum lies at the edge of the width grid, expand the grid rather than freezing that endpoint.

The raw gamma-cell artifacts are configured for 14-day retention while the cross-gamma summary artifact is retained for 30 days. This is a reproducibility risk: the raw records required for future independent re-audit can expire before the summary. Archive the raw inputs, configs, logs, manifests and audit outputs in a durable versioned location before relying on the summaries.

The H0 stage must remain separate. The five-arm selected-likelihood implementation exists on related branches and PR #44 has a successful test run, but the curvature DGP has not been carried through a complete, exact-head, independently audited H0 experiment. Density improvement is not evidence of H0 calibration improvement. Before H0 claims, verify the selected event likelihood, the GW detection normalization, the common host-rate/base measure, and the complete shared-coefficient marginalization for Arm E on the exact integrated commit.

## 5. Literature positioning: material 2026 developments

The project must be positioned against a rapidly developing literature. The following are directly relevant and should be read in full before a novelty claim:

1. **Tagliazucchi et al. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” Physical Review D 114, 083523, published 13 October 2026.** It uses a selection-aware hierarchical Bayesian likelihood and Gaussian-process reconstruction of the line-of-sight galaxy redshift field, with mock catalogues to z≤0.4. This is a close astrophysical neighbour: the project cannot claim novelty for reconstructing incomplete galaxy information or for selection-aware Bayesian treatment in general. A narrower distinction may be finite validation-sample estimation of a parametric catalogue-inclusion mechanism, IPW instability and logistic functional-form misspecification. DOI: https://doi.org/10.1103/vmzf-1wxq ; article: https://journals.aps.org/prd/abstract/10.1103/vmzf-1wxq

2. **Leyde, Baker & Enzi (2026), “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology,” JCAP 01 (2026) 013, published 5 January 2026.** It reconstructs the underlying spatial galaxy field and magnitude distribution, validated on Millennium-based simulations. This is another direct catalogue-completion neighbour, though not the same finite-validation propensity-estimation question. DOI: https://doi.org/10.1088/1475-7516/2026/01/013 ; arXiv: https://arxiv.org/abs/2507.12171

3. **Borghi et al. (2026), “Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology,” A&A 706, A199, published 10 February 2026.** It studies catalogue incompleteness and host weighting in CHIMERA with O5-like simulations. This rules out broad novelty claims about studying catalogue incompleteness and host weighting in dark-siren H0 inference. DOI: https://doi.org/10.1051/0004-6361/202557354 ; arXiv: https://arxiv.org/abs/2509.18243

4. **Datrier & Hendry (2026), “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” MNRAS 549, stag828, published 1 May 2026.** It develops a robust apparent-magnitude completeness-limit test for gwcosmo. This is a related diagnostic contribution and should be distinguished from estimating a redshift-dependent inclusion model from a finite validation sample. DOI: https://doi.org/10.1093/mnras/stag828 ; arXiv: https://arxiv.org/abs/2502.14164

5. **Naveed, Turski & Ghosh (2026), “Dark standard siren cosmology with bright galaxy subsets,” MNRAS 548, stag605, published 28 March 2026.** It studies the use of bright galaxy subsets in dark-siren H0 inference. This is relevant context for catalogue selection and the information content of subsets. DOI: https://doi.org/10.1093/mnras/stag605

6. **Gray, Williams & Papadopoulos (2026), “Expanding the scope of dark siren cosmology: Inferring the population properties of gravitational wave-hosting galaxies,” arXiv:2609.04991, posted 4 September 2026.** This develops joint inference of host-galaxy population properties alongside cosmology and GW population parameters. It reinforces the need to keep catalogue inclusion distinct from host/merger weighting. arXiv: https://arxiv.org/abs/2609.04991

### Provisional novelty boundary

Do **not** claim a new general galaxy-catalogue completeness method, a new dark-siren likelihood, or the first Bayesian reconstruction of incomplete galaxy catalogues.

A defensible candidate contribution, subject to a full statistical-literature review, is:

> A controlled operating-characteristic study of estimated catalogue-inclusion probabilities under finite validation, parametric functional-form misspecification and inverse-weight instability, separating density-reconstruction error from downstream cosmological posterior calibration.

This is not yet an established novelty claim. A statistical-methodology paper would need a general insight beyond a single synthetic astrophysics example—e.g. a principled diagnostic or theoretically supported operating-characteristic result, plus careful uncertainty quantification. An astrophysics paper would need a clear connection to the dark-siren likelihood and selection model, with a fully validated H0 propagation stage. Neither pathway is publication-ready on the evidence inspected here.

## 6. Immediate next actions

1. Obtain exact-head CI for PR #55 and PR #56; do not infer success from older commits.
2. Preserve the raw artifacts for the curvature experiment beyond the short workflow retention window.
3. Run the seven-gamma workflow only after confirming the exact head, runtime budget, and immutable config/seed schedule; then independently audit every child and the cross-gamma artifact.
4. Report density-risk curves and paired contrasts across the full bandwidth grid, including Monte Carlo uncertainty and weight diagnostics. Do not select a bandwidth from this same exploratory pilot.
5. Reconcile the PR #44, #54 and #55 branch histories before treating the five-arm H0 implementation as integrated with the curvature DGP.
6. Design a separate fresh-seed training/holdout protocol if operational bandwidth selection is needed.
7. Only after density-stage gates pass, run an H0-stage experiment with independent likelihood/selection-normalization audits, posterior-grid and quadrature convergence, bias/RMSE, coverage and SBC.
8. Complete a statistical-methods literature search on estimated propensity scores, finite validation, model misspecification, positivity and Bayesian nuisance-parameter uncertainty before drafting novelty language.

## 7. Audit provenance

Files inspected on PR #55 head include:
- `.github/workflows/quadratic-logit-curvature-density-pilot.yml`
- `docs/quadratic_logit_curvature_stress_design_2026-11-14.md`
- `docs/independent_quadratic_logit_calibration_check_2026-11-16.md`
- `docs/curvature_known_probability_and_weight_diagnostic_audit_2026-11-14.md`
- `scripts/run_absolute_bandwidth_density_replicates.py`
- `scripts/audit_absolute_bandwidth_density_replicates.py`
- `scripts/audit_quadratic_logit_curvature_pilot.py`
- `src/dark_siren_robustness/five_arm_comparison.py`
- `src/dark_siren_robustness/generative_selection.py`
- `tests/test_curvature_cross_gamma_completeness.py`

This note was added on a separate research branch based on PR #55 head. No changes were made to `main`, no pull request was merged, and no simulation workflow was dispatched by this audit.


## 8. Static review of the five-arm H0 likelihood

I inspected `src/dark_siren_robustness/five_arm_kde_likelihood.py`, the five-arm smoke runner, the independent smoke auditor, and the likelihood/convergence tests on the PR #55 branch.

The implemented per-event selected likelihood has the form
[
\mathcal L_i(H_0\mid d_i,\mathrm{det})
=
\frac{\int p_k(z)\,\mathcal N(d_i;cz/H_0,\sigma_i^2)\,dz}
{\int p_k(z)\,\Phi((cz/H_0-d_{\rm th})/\sigma_i)\,dz},
]
and the dataset likelihood is the product over events (equivalently, the sum of event log-ratios). This is the appropriate conditional-on-detection form for a fixed number of independently detected events under the code's Gaussian distance-error model and threshold rule. A manuscript must state explicitly that the analysis conditions on the detected event count; a model for the event count/rate would require additional terms.

For Arm E, the implementation evaluates the full selected dataset log likelihood at each shared coefficient-posterior node and then computes a single `logsumexp(log L(theta_k) + log q_k)`. This is the correct order for integrating a shared nuisance parameter. It avoids the invalid alternatives of averaging host weights first or multiplying separately marginalized per-event likelihoods.

The likelihood accepts either one common absolute redshift bandwidth or a common Scott multiplier. The absolute-bandwidth option applies the same physical kernel SD to A-D and every E node; a common Scott multiplier does not generally imply a common realized physical bandwidth across differently weighted arms. Confirmatory comparisons should therefore use the absolute-bandwidth option or explicitly treat the Scott rule as a different estimand.

The code multiplies the normalized host-redshift density by the event-distance likelihood in the numerator and by the Gaussian-tail detection probability in the denominator. This is internally consistent with the documented threshold-on-measured-distance generator. The host/merger base factor is carried in the raw arm weights in the smoke and the grid-level multiplier is set to one, avoiding a second application of the same base factor.

**Limitations of this static review:** I did not execute the likelihood or independently recompute a complete H0 artifact in this session. The smoke auditor explicitly uses the persisted inclusion-posterior nodes; it reconstructs the KDE and selected likelihood independently but does not itself rebuild the logistic posterior quadrature from the validation records. Separate direct-integration/quadrature tests exist and passed in the PR #44 workflow run noted above, but that is not the same as an artifact-level posterior refit for every H0 replicate. The current H0 smoke uses small engineering settings and is not confirmatory evidence. The Euclidean distance relation, fixed Gaussian noise, simple threshold detector, fixed host-rate function and absence of real event-level likelihood data remain major toy-model assumptions.

**Decision:** the selected-likelihood algebra and shared-node marginalization look structurally coherent under the stated fixed-N toy assumptions, but the H0 stage remains NOT CLEARED for scientific claims. It requires exact integrated-head tests, independent artifact-level validation of posterior nodes, numerical grid/normalization checks, and a repeated-sampling calibration experiment under the same curvature DGP.


## 9. Statistical-methods prior art: important boundary conditions

The dark-siren papers above are only half of the novelty review. The statistical structure is also closely related to missing-data, nonresponse adjustment, propensity-score estimation, and two-phase sampling. A statistical-methodology contribution cannot treat estimated inclusion probabilities or Bayesian uncertainty propagation as new concepts.

Relevant foundations and adjacent work identified in an initial search:

- **Saegusa & Wellner (2013), “Weighted likelihood estimation under two-phase sampling,” Annals of Statistics 41, 269–295, DOI 10.1214/12-AOS1073.** Develops asymptotic theory for weighted likelihood under two-phase sampling, including estimated weights and calibration. It is directly relevant to the sampling/weighting framework and should inform the estimand and asymptotic discussion: https://doi.org/10.1214/12-AOS1073
- **Hasler (2024), “Inference from sampling with response probabilities estimated via calibration,” Survey Methodology 50(2), 235–259.** Studies estimated response probabilities under logistic response models and discusses asymptotic behaviour, efficiency, misspecification, existence/convergence problems and extreme weights. This is especially relevant to the distinction between maximum-likelihood propensity fitting and calibration-based weighting: https://publications.gc.ca/collections/collection_2024/statcan/12-001-x/12-001-x2024002-eng.pdf
- **An (2010), “Bayesian Propensity Score Estimators: Incorporating Uncertainties in Propensity Scores into Causal Inference,” Sociological Methodology 40(1).** Bayesian treatment of propensity-score uncertainty predates this project. It is not the same dark-siren estimand, but it rules out novelty claims for Bayesian propensity uncertainty integration in general: https://doi.org/10.1111/j.1467-9531.2010.01226.x
- **“Uncertainty in the design stage of two-stage Bayesian propensity score analysis” (Statistical Medicine, 2020).** Develops posterior inference that marginalizes over design-stage uncertainty in a two-stage propensity-score analysis: https://doi.org/10.1002/sim.8486
- **“Model misspecification and bias for inverse probability weighting estimators of average causal effects” (Biometrical Journal, 2022).** Analyses large-sample bias under propensity-score and outcome-model misspecification and compares normalized and unnormalized IPW: https://doi.org/10.1002/bimj.202100118
- **Katsumata (2025), “How should we estimate inverse probability weights with possibly misspecified propensity score models?”, Political Science Research and Methods 13(2), 315–336.** Proposes a weight-estimation strategy for possibly misspecified propensity models and evaluates finite-sample performance: https://doi.org/10.1017/psrm.2024.23
- **“A Tutorial for Propensity Score Weighting Methods Under Violations of the Positivity Assumption” (Statistics in Medicine, first published 27 November 2025).** Provides a current review of positivity violations and alternative target estimands/weighting methods: https://doi.org/10.1002/sim.70329
- **“Two-Phase Sampling Designs for Data Validation in Settings with Covariate Measurement Error and Continuous Outcome” (JRSS Series A, 2021).** Shows that validation-sample design itself can materially affect efficiency and that SRS is not necessarily optimal: https://doi.org/10.1111/rssa.12689

### Consequence for our statistical contribution

The project should not frame these as novel in themselves:
- estimating inclusion/response probabilities from validation data;
- using inverse estimated probabilities as weights;
- Bayesian integration over propensity/inclusion coefficients;
- studying positivity, extreme weights, or propensity-model misspecification;
- conducting a finite-sample simulation comparing IPW variants.

The possible contribution is narrower and must be demonstrated rather than asserted: a carefully designed, dark-siren-relevant operating-characteristic analysis that connects finite validation, misspecified inclusion functions, weight concentration, reconstructed host-redshift density, and a selection-conditioned cosmological likelihood. The added value must come from the way these pieces interact, the diagnostics/decision insight obtained, and a validated H0 calibration analysis—not simply from placing established IPW tools in an astrophysics setting.

A stronger statistical-journal case would likely require at least one of:
1. a new diagnostic with a defined target and independently evaluated ability to predict density or posterior-calibration failure;
2. a theoretical result for the specific weighted density or selected-likelihood functional under estimated/misspecified inclusion probabilities;
3. a principled estimator or uncertainty procedure that improves a clearly defined risk criterion and survives independent stress testing;
4. a generalizable design result for allocating a limited validation budget across redshift/support regions.

These are candidate directions, not recommendations to add all four. Choose one only after the baseline experiment establishes a concrete failure mode.


## 10. Additional statistical-computation audit: validation posterior quadrature

The degree-one logistic posterior implementation uses a proper independent Gaussian prior, a trust-region exact-Hessian mode fit with Newton refinement, and corrected tensor Gauss-Hermite importance weights. The Laplace Gaussian is a proposal, not the posterior itself; the code correctly reweights quadrature nodes against the declared logistic likelihood times prior. This is materially stronger than treating a Laplace approximation as the final posterior.

However, the curvature workflow currently calls the fitter with `--quadrature-order 5`. For a two-coefficient degree-one model, this is only (5^2=25) tensor nodes. Correct importance weights do not guarantee that a low-order quadrature resolves a skewed or heavy-tailed posterior, particularly when the finite validation sample has near-separation or the catalogue extrapolates beyond the observed validation redshift range.

The repository includes a quadrature-ladder utility that compares coefficient posterior means/covariances, posterior-mean inclusion probabilities, and log evidence over increasing orders. Its default predictive grid spans only the minimum to maximum observed validation redshift. For the scientific application, this default is insufficient by itself: catalogue redshifts outside that observed interval can receive the most influential extrapolated inverse weights. The ladder should be called with an explicit fixed prediction grid spanning the full analysis support [0.01, 0.40], and the downstream integrated H0 likelihood itself must also be checked for order stability. The current density workflow does not use the ladder as a precondition for the order-5 run.

**Required gate before confirmatory claims:** on each validation replicate, compare at least a prespecified ladder such as orders 5, 7, 9, 11 (and extend if needed) using the full catalogue support. Record posterior coefficient moments, posterior predictive inclusion curves, node-weight ESS, and the final density risk. For H0 work, compare the marginalized log-likelihood/posterior across quadrature orders, not only coefficient moments. Define tolerances before examining estimator rankings. If order 5 fails, the pilot output remains useful as a low-cost screen but cannot be used as the final uncertainty-integrated result.

A further diagnostic should report the validation redshift support versus the observed catalogue support for every replicate. A random validation sample from the same parent distribution can legitimately omit tails at small n; that makes the fitted logistic model extrapolate into regions that may dominate IPW. Report the fraction of catalogue objects outside the validation range, their total estimated inverse-weight mass, and posterior uncertainty in predicted inclusion probability there. Do not label this alone as a validation-design failure: it is a finite-sample support diagnostic.

## 11. Curvature geometry is not monotone in gamma

For the specified linear predictor,
[
\eta_\gamma(x)=\alpha_\gamma-3.9x+\gamma(x^2-c_g),
qquad
\frac{d\eta_\gamma}{dx}=-3.9+2\gamma x.
]
Since (x\in[-1,1]), at (gamma=+2) the derivative ranges from (-7.9) to (+0.1), so the high-redshift end can turn very slightly upward; at (gamma=-2), the derivative ranges from (-0.1) to (-7.9). Thus gamma sign changes not only curvature direction but also where the steepest decline and minimum inclusion probability occur. Intercept calibration fixes the parent-average retention rate, but it does not fix the inclusion curve's monotonicity, local slope, or overlap profile.

The primary seven-cell study should retain its frozen design and report the actual (pi(z)), derivative, minimum/maximum inclusion probability, parent-average rate, and population weight moments for each gamma. If a severity-matched comparison is later useful, treat it as a separate, explicitly exploratory design rather than retuning the existing gamma grid after observing results.

## 12. Further statistical prior art and implications

The first-pass statistical search also identified:
- Kim & Riddles (2012), “Some theory for propensity-score-adjustment estimators in survey sampling,” *Survey Methodology* 38(2), 157–165. It develops asymptotic properties and variance considerations for response-propensity adjustment in survey sampling. https://www150.statcan.gc.ca/n1/en/pub/12-001-x/2012002/article/11754-eng.pdf
- Chen & Kim (2014), “Two-phase sampling experiment for propensity score estimation in self-selected samples,” *Annals of Applied Statistics* 8(3), 1492–1515, DOI 10.1214/14-AOAS746. It uses a two-phase validation design to estimate response propensity parameters consistently. https://doi.org/10.1214/14-AOAS746
- Katsumata (2025), “How should we estimate inverse probability weights with possibly misspecified propensity score models?”, *Political Science Research and Methods* 13(2), 315–336, DOI 10.1017/psrm.2024.23. It directly studies the consequences of misspecified propensity models and methods targeting weight error rather than propensity-score error. https://doi.org/10.1017/psrm.2024.23
- “Uncertainty in the design stage of two-stage Bayesian propensity score analysis,” *Statistics in Medicine* 39 (2020), 2265–2290, DOI 10.1002/sim.8486. It explicitly marginalizes design-stage uncertainty in Bayesian propensity-score analysis. https://doi.org/10.1002/sim.8486

This makes it especially important that our primary estimand is not merely coefficient uncertainty or propensity prediction error. The relevant object is the downstream **host-redshift density functional and selection-conditioned cosmological likelihood**. A statistically meaningful comparison should therefore report:
1. calibration of inclusion probabilities under the known simulation truth;
2. error in the induced weighted host-density functional;
3. posterior uncertainty in that density functional;
4. downstream selected-likelihood and H0 calibration consequences.

A method may estimate coefficients poorly but still recover the relevant weighted density well; conversely, good average probability calibration may coexist with harmful tail-weight concentration. These should not be collapsed into one score.


## 13. Asymptotic target under a misspecified inclusion model

This is a central theoretical point that should guide the paper.

Let the parent redshift density be (g(z)), the true catalogue inclusion probability be (pi_0(z)), and the validation-fitted working model converge to a pseudo-true probability (pi_*(z)). Let (b(z)>0) be the declared galaxy-to-merger-rate factor. Ignoring finite-sample KDE smoothing for the moment, the observed catalogue has density proportional to (g(z)\pi_0(z)). Applying working-model inverse weights (1/\pi_*(z)) and the common rate factor (b(z)) yields the limiting normalized host density
[
p_*(z)=
\frac{g(z)\,\pi_0(z)\,b(z)/\pi_*(z)}
{\int g(u)\,\pi_0(u)\,b(u)/\pi_*(u)\,du}.
]
The target is
[
p_0(z)=
\frac{g(z)b(z)}{\int g(u)b(u)\,du}.
]
Therefore, consistency of the IPW density requires the ratio (\pi_0(z)/\pi_*(z)) to be constant over the target support (the constant cancels on normalization); exact inclusion-model correctness is sufficient, while a proportional misspecification can also cancel after normalization. In general, a degree-one logistic model fitted to a nonzero-curvature logistic truth will have a nonconstant ratio and hence a nonvanishing asymptotic density distortion.

This gives a clean separation:
- Increasing validation size reduces estimation uncertainty around the pseudo-true working model.
- It does not remove the pseudo-true functional-form bias.
- Consequently, the posterior over coefficients can become very concentrated while the IPW density remains biased.
- Arm E's coefficient integration can improve finite-validation uncertainty propagation, but cannot asymptotically repair omitted quadratic structure if the fitted model remains degree one.

The relevant validation target for the logistic fit is the Bernoulli log-score/KL projection under the validation sampling distribution. That projection is not generally the same as minimizing density IAE, H0 bias, or posterior miscalibration. This is why coefficient bias, probability calibration, density risk, and cosmological calibration must be reported as separate outcomes.

For the five-arm comparison, Arm C (known π) isolates finite-parent/KDE error under correct selection weights; Arm D exposes the combined finite-validation and model-projection error; Arm E adds posterior uncertainty propagation conditional on the same working model. Arm B remains a finite-parent oracle density reference rather than a noise-free analytic truth. Comparisons should include the analytic target density as well as Arm B, so that finite-parent oracle KDE error is not mistaken for irreducible truth.

This derivation is a theoretical consequence of the stated sampling model, not an empirical result from the unrun seven-gamma pilot. It should be independently checked and expanded into a formal proposition if the manuscript is to make an asymptotic claim.
