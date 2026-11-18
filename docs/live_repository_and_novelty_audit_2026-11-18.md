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
