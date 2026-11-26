# Research and conversation checkpoint — 2026-11-26

**Purpose.** Preserve the project-relevant decisions, mathematical clarifications, Consensus research output, unresolved issues, and next actions from the 2026-11-26 ChatGPT research continuation. This is a documentation-only checkpoint. It does not alter estimator code, launch simulations, validate a production H0 posterior, establish novelty, or authorize a confirmatory campaign.

**Repository:** `vedantkushwah033-blip/dark-siren-cosmology-robustness`  
**Saved on branch:** `research/prepared-bundle-independent-arm-auditor-20261122`  
**Related open draft PR:** [#181 — independently audit arm cells from preparation-only inputs](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/pull/181)  
**Default branch:** `main` (not modified by this checkpoint).

> Provenance distinction: the local worktree checkpoint described by the user on 2026-11-26 was at `C:\Users\vedant\AppData\Local\Temp\opencode\kde_robustness_phase12_impl_20261126`, branch `research/phase1-2-d3-arme-impl-20261126`, HEAD `7c7bf056fc91447cee49ff16f2a48c8f54ad88eb`. That branch/ref was not found on GitHub when checked. The present documentation commit is therefore saved on the verified remote PR #181 branch, not represented as a commit to that local-only worktree. Reconcile or cherry-pick deliberately when the local worktree is available.

---

## 1. Scientific scope and guardrails

The project is a **synthetic, low-redshift Euclidean dark-siren H0 methodology benchmark**. It is not an LVK/GWTC/public-data analysis, an observational H0 measurement, a validation against a real galaxy survey, or proof that one estimator is universally superior.

Keep three mechanisms distinct:

1. **Galaxy catalogue inclusion:** `pi_i = P(R_i=1 | x_i)`, the probability a galaxy appears in the catalogue.
2. **Astrophysical merger-host weighting:** `h_i`, the relative probability/rate that a galaxy hosts a merger.
3. **GW detection selection:** `P(D_GW=1 | theta,H0)`, the probability the detector pipeline detects an event.

Do not silently use one as a proxy for another. In particular, galaxy-inclusion IPW is not a substitute for a host-rate model or the GW detection-selection denominator. Avoid applying IPW and a missing-host mixture as two corrections for the same missing population unless the joint probability model explicitly derives both factors.

The current central research question is:

> Under what conditions does estimating galaxy-catalogue inclusion probabilities improve dark-siren inference, and under what conditions do misspecification, finite validation samples, and IPW instability eliminate or reverse the benefit?

The most defensible candidate contribution is **not** a new universal IPW estimator. It is a controlled study of uncertainty in per-galaxy inclusion probabilities estimated from finite validation data, and the consequences for a shared, selection-normalized multi-event H0 likelihood.

## 2. DGP and inclusion-model design to preserve

Latest proposed H0-stage design (keep separate from older density-only experiments):

- Parent population: approximately N=4,000.
- Detected GW events: 30.
- Independent validation frame: approximately M=1,500, drawn from the same latent redshift law as the parent population.
- Validation samples: nested SRSWOR prefixes n=100, 300, 1,000.
- True inclusion probability:
  `pi_0(x) = expit(alpha + beta*x + gamma*(x^2-c))`.
- Degree-1 logistic fit is misspecified when gamma is nonzero; degree-2 fit is correctly specified under this stated DGP.
- `c` must be the parent-law moment `E_P[x^2]`, not assumed to be 1/3 without deriving that from the actual parent distribution.
- Host-rate weighting is a separate model component.
- Use distinct RNG streams for parent generation, GW event generation, validation-frame generation, and validation-sample selection.

Do not mix this independent-frame DGP with older density-recovery experiments that split a fixed finite parent into catalogue/validation pieces. State the sampling design in every experiment report.

## 3. Four inclusion-coefficient estimands — do not conflate

The following targets are different and must be named explicitly:

- **Superpopulation pseudo-true coefficient** `beta_star`:
  `E_P[X {R - q_beta(X)}] = 0`, where `q_beta(x)=expit(x' beta)`.
- **Finite validation-frame empirical pseudo-true coefficient** `beta_F`:
  `sum_{i=1}^M X_i {R_i-q_beta(X_i)}=0`.
- **Finite-parent realized-label MLE** `beta_P`:
  `sum_{i=1}^N X_i {R_i-q_beta(X_i)}=0`.
- **Finite-parent covariate-only KL projection** `beta_P^KL`:
  `sum_{i=1}^N X_i {pi_0(X_i)-q_beta(X_i)}=0`.

The realized-label parent MLE and covariate-only KL projection are not interchangeable. The validation posterior must be evaluated against the target the H0 likelihood actually needs. Also state whether the frame and parent are fixed or random.

## 4. Finite validation sampling and FPC

For a fixed frame F_M and SRSWOR, the simple-mean design variance is

`Var_d(ybar_s | F_M) = (1-n/M) S_M^2/n.`

For a regular logistic M-estimator, the first-order conditional design approximation is

`Var_d(beta_hat | F_M) ~= A_M^{-1} [(1-n/M) S_{psi,M}/n] A_M^{-T},`

where `A_M=M^{-1} sum q_i(1-q_i) X_i X_i'` and `psi_i=X_i(R_i-q_i)`. It requires a well-behaved nonsingular information/sensitivity matrix, a root, controlled leverage/moments, adequate linearization, and the stated SRSWOR design. At n=M, conditional design variance is zero because the fixed frame has been fully enumerated.

**Important unconditional correction:** if F_M itself is an iid sample from P and n observations are selected by SRSWOR, the unconditional sampled records are iid from P. The between-frame component offsets the conditional FPC reduction to first order:

`E_F[Var_d(beta_hat|F)] + Var_F[E_d(beta_hat|F)] ~= A^{-1} B A^{-T} [(1-n/M)/n + 1/M] = A^{-1} B A^{-T}/n.`

This is a first-order result under the simple iid-frame/SRSWOR setup, not a universal formula for complex designs or nonlinear finite-sample behavior. For a misspecified model B need not equal A; the repeated-sampling covariance is sandwich `A^{-1} B A^{-T}/n`, not simply inverse information.

Nested prefixes share observations and have nonzero conditional covariance. Independent samples from the same fixed frame have zero conditional design covariance but can share between-frame covariance unconditionally. Independent frames have no shared frame component.

## 5. Bayesian posterior versus design-based variance — critical unresolved issue

The code paths represent different inferential objects:

- `inclusion_posterior.py`: proper Bayesian logistic likelihood with Gaussian prior and quadrature nodes.
- `shared_inclusion_likelihood.py`: integrates the H0 likelihood over shared inclusion-parameter posterior nodes.
- Older `estimated_inclusion.py`: MLE plus inverse observed-information covariance; not the same method.

Under misspecification, ordinary Bayesian posterior covariance generally follows inverse curvature `A^{-1}/n`, whereas repeated-sampling covariance for the pseudo-true parameter is sandwich `A^{-1} B A^{-T}/n`. They agree under the information identity (correct specification), but not generally under quadratic truth with a degree-1 fit.

Consequences:
- A model-based posterior does not automatically represent conditional-on-frame design uncertainty.
- Do not add D3 sandwich variance to posterior draws without a precise estimand and a derivation; it may double-count uncertainty already propagated through posterior draws.
- For a superpopulation target, compare ordinary posterior spread with repeated-sampling sandwich spread, or test a justified calibrated/generalized posterior (e.g. open-faced sandwich adjustment / Q-posterior) as a distinct method.
- For a finite-frame target, derive a coherent posterior conditional on the fixed frame; do not pretend a design covariance is an independent posterior component.
- For a finite-parent target, explicitly model transport from validation frame to parent or treat it as a different estimand.
- The current shared-node method is a modular/cut analysis if GW data do not update inclusion coefficients. Label it as such and justify it; a full joint posterior is a different model.

**Do not assume** a quasi-posterior or sandwich correction is automatically best. Compare methods against the declared estimand and evaluate both bias and calibration.

## 6. Shared nuisance integration and selection normalization

For a single shared inclusion nuisance vector gamma, the intended structure is

`p(D|H0,V) = integral [ product_{j=1}^J L_j(H0,gamma) ] p(gamma|V) dgamma.`

Do not replace it with

`product_{j=1}^J integral L_j(H0,gamma) p(gamma|V) dgamma,`

which effectively gives each event an independent gamma.

The event-host mixture/KDE and GW detection-selection efficiency must be computed from the same corrected galaxy measure and consistent normalization. If the numerator uses one measure and the selection denominator another, the selection-normalized likelihood is incoherent. Existing shared-node logic was previously described as summing event log-likelihoods at fixed nuisance nodes, subtracting J log beta(H0,gamma), and marginalizing once over nodes. That source-level assessment is historical and must be rechecked against the exact current branch before being treated as verified.

## 7. Normalized IPW derivative and diagnostics

Let `q_i=expit(x_i' beta)`, `g_i=grad_beta log q_i=(1-q_i)x_i`, `w_i=1/q_i`, and `p_i=w_i/sum_k w_k`. Then

`grad_beta p_i = p_i (sum_k p_k g_k - g_i).`

For a weighted KDE `f(z;beta)=sum_i p_i K_h(z-z_i)`,

`grad_beta f(z;beta)=sum_i p_i K_i (sum_k p_k g_k-g_i).`

This includes the normalization-denominator derivative. `J Sigma_beta J'` is a first-order delta approximation; do not add it if posterior draws already propagate the same beta uncertainty. Verify implementation with finite differences and Monte Carlo.

Required weight/support diagnostics:
- Kish ESS and entropy ESS.
- Maximum normalized weight.
- Top 1% and 5% total weight mass.
- Positivity and support-gap checks.
- Clipping/truncation sensitivity, clearly reported as a different estimator.
- Posterior tail/prior sensitivity, separation and optimizer failures.
- Quadrature convergence of H0 posterior summaries, not only node-level likelihoods.
- Full SUCCESS / NON_EVALUABLE / FAILED / NOT_RUN ledger, with every planned cell in the denominator.

Do not label ESS or other pre-posterior diagnostics as predictive of posterior failure unless their out-of-replicate predictive performance is tested.

## 8. User-reported local implementation checkpoint (not independently verified)

The user reported the following local state on 2026-11-26:

- Worktree: `C:\Users\vedant\AppData\Local\Temp\opencode\kde_robustness_phase12_impl_20261126`.
- Local branch: `research/phase1-2-d3-arme-impl-20261126`.
- Local HEAD: `7c7bf056fc91447cee49ff16f2a48c8f54ad88eb`.
- Modified: `scripts/run_inclusion_posterior_smoke.py` (manifest key path normalization using `.as_posix()`).
- Added `src/dark_siren_robustness/d3_finite_validation.py`: SRSWOR FPC covariance, nested-prefix covariance, zero conditional covariance for independent SRSWOR samples, census n=M zero design variance, normalized IPW Jacobian and delta propagation; ill-conditioned/non-PD information rejected by default, pseudoinverse only with explicit `allow_pseudoinverse=True`.
- Added tests: `test_d3_finite_validation.py`, `test_d3_downstream_propagation.py`, `test_d3_failure_guards.py`, `test_arm_e_certification.py`, `test_arm_e_eventwise_vs_joint_production.py`, `test_gaussian_shared_parameter_validation.py`.
- Added docs under `research/project_completion/final_audit/*`, `next_gate_20261126/*`, `continuation_20261126/*`, plus `research/project_completion/final_mission/INITIAL_REPOSITORY_AUDIT.md`, `IMPLEMENTATION_STATUS.md`, `FAILURE_LEDGER.csv`.
- User-reported focused tests:
  - 16 passed: Arm E certification, D3 finite validation, downstream propagation, Gaussian shared-parameter validation.
  - 61 passed: selection normalization, joint selection/missing hosts, selection quadrature audit, inclusion posterior smoke.
  - 3 passed: Arm E eventwise-vs-joint production and D3 failure guards.
  - 8 passed: D3 failure guards, finite validation, downstream propagation.
  - 35 passed: estimated inclusion, inclusion posterior, inclusion posterior smoke.
- Reported limitations: D3 production posterior integration incomplete; Arm E R5–R10 certification incomplete; adversarial selection-normalization certification incomplete; H0/z-grid and KDE convergence/campaign incomplete; no confirmatory calibration; no independent scientific audit.
- Reported ad hoc Arm E quadrature posterior-summary maximum delta vs q=21: q=5 0.0471; q=15 0.000126. This is not certification.
- Windows historical path mismatch was corrected; Python 3.15 SciPy `_flapack` DLL blocked by Application Control; Python 3.12 is the practical test environment.
- Production H0 posterior remains NOT VALIDATED; novelty NOT VERIFIED; no submission-ready ASA–SBSS / ASA–SLDS / IEOM package.

These are explicitly **user-reported local facts**, not independently reproduced by this documentation task. The branch was not found on GitHub; do not imply that this local checkpoint has been pushed.

## 9. Remote repository context observed while saving

Repository metadata confirmed that the repository is public, the default branch is `main`, and the connected account has push permission. The following open draft PR was inspected:

- PR #181, head branch `research/prepared-bundle-independent-arm-auditor-20261122`, head SHA `38547cbee384e7341aa7b45ba343936220397812`, base `research/preparation-only-shared-input-bundle-20261122`.
- PR #181 adds an independent raw-input per-arm auditor for preparation-only bundles; it rebuilds A/B/C weights and D/E logistic posterior nodes, reflected KDEs, selected multi-event likelihoods, detection efficiencies, posterior summaries, and weight diagnostics.
- PR #181's body reports focused exact-head workflows successful on Python 3.10/3.11/3.12 for its arm-audit and serial lifecycle tests. Its stated limits remain engineering/integrity only: no confirmatory calibration, no estimator-superiority claim, and aggregate lifecycle auditor not yet updated to independently reconstruct all raw-bundle arm-cell outputs.
- PR #180 was also listed as an open draft policy PR for non-evaluable accounting and freeze gates; it documents no frozen universal ESS cutoff and no confirmatory results.
- PR #177 provides preparation-only shared-input bundles with no inference outputs and distinct RNG streams. These PR descriptions are repository metadata, not an independent rerun of CI during this checkpoint.

The documentation file was saved on the PR #181 head branch only. No changes were made to `main`, no PR was merged, and no estimator code or existing file was overwritten.

## 10. Consensus research continuation — what was actually run

A new Consensus Pro research thread was created because no existing thread titled around “dark siren” was found:

- [Consensus research thread](https://consensus.app/search/continue-an-existing-research-programme-on-a-narro/RQNjVV9iSwCKvD6aXoZ5TA/)
- Thread ID: `RQNjVV9iSwCKvD6aXoZ5TA`.
- Interaction ID: `ayJnTH6CTy-RYB2V7VTwJA`.
- The returned research document was retrieved and archived in the prior ChatGPT exchange; a separate downloadable report was not created because the file-generation tool failed.

The Consensus memo covered misspecified Bayesian inference, quasi-posteriors, survey/validation designs, finite-population inference, and a few Bayesian sandwich papers. It returned the following references/metadata (all bibliographic details and transfer claims still require independent verification from publisher/arXiv full text):

1. Müller, U.K. (2013), “Risk of Bayesian Inference in Misspecified Models, and the Sandwich Covariance Matrix,” *Econometrica* 81, 1805–1849. DOI: 10.3982/ecta9097.
2. Vrugt, J.A. & Diks, C.G.H. (2025), “The Learning Rate Is Not a Constant: Sandwich-Adjusted Markov Chain Monte Carlo Simulation,” *Entropy* 27, 999. DOI as returned: 10.3390/e27100999.
3. Frazier, D.T., Kohn, R., Drovandi, C. & Gunawan, D. (2023), “Reliable Bayesian Inference in Misspecified Models,” arXiv:2302.06031.
4. Frazier, D.T., Drovandi, C. & Kohn, R. (2023), “Calibrated Generalized Bayesian Inference,” arXiv:2311.15485.
5. Huggins, J.H. & Miller, J.W. (2024), “Reproducible parameter inference using bagged posteriors,” *Electronic Journal of Statistics* 18, 1549–1585. DOI: 10.1214/24-EJS2237.
6. Agnoletto, D., Rigon, T. & Dunson, D.B. (2025), “Bayesian inference for generalized linear models via quasi-posteriors,” *Biometrika* 112, asaf022. DOI: 10.1093/biomet/asaf022.
7. Chiba, H. & Morikawa, K. (2026), “Semiparametric Efficient Inference under Non-Informative Complex Survey Designs,” arXiv:2609.17288.
8. Amorim, G., Tao, R., Lotspeich, S.C., Shaw, P.A., Lumley, T. & Shepherd, B.E. (2021), “Two-Phase Sampling Designs for Data Validation in Settings with Covariate Measurement Error and Continuous Outcome,” *JRSS A* 184, 1368–1389. DOI: 10.1111/rssa.12689.
9. Dagdoug, M. & Haziza, D. (2026), “Machine learning methods for finite population parameter estimation in survey sampling,” arXiv:2604.01160.
10. An, Z., Dagdoug, M., Haziza, D. & Tillé, Y. (2026), “Agnostic Model-Assisted Estimation with Machine Learning for Survey Data,” arXiv:2609.10904.
11. Hoff, P.D. & Wakefield, J.C. (2013), “Bayesian sandwich posteriors for pseudo-true parameters,” *Journal of Statistical Planning and Inference* 143, 1638–1642. DOI: 10.1016/j.jspi.2013.05.014.
12. Grünwald, P. & Van Ommen, T. (2017), “Inconsistency of Bayesian Inference for Misspecified Linear Models, and a Proposal for Repairing It,” *Bayesian Analysis* 12, 1069–1103. DOI: 10.1214/17-BA1085.
13. Giordano, R. & Broderick, T. (2023), “The Bayesian Infinitesimal Jackknife for Variance,” arXiv:2305.06466.

Consensus did **not** retrieve or substantiate the dark-siren papers specified in the research brief. It explicitly said the dark-siren prior-art list was not in its retrieved corpus. Therefore it did not establish a novelty gap or prove that dark-siren and survey-statistics literatures have not cross-cited.

### Critical corrections to the Consensus memo

- Its opening claim that inclusion-probability estimation “should help only when the inclusion model is correctly specified” is too strong. Misspecification can damage estimates and calibration, but whether IPW improves a particular H0 target depends on the estimand, support, transport, model error, and likelihood construction. A misspecified propensity fit may still reduce some bias in some settings; this is an empirical/theoretical question, not a blanket theorem.
- Its statement that degree-1 with gamma nonzero is “exactly” the posterior-misspecification case is valid as a working-model concern, but it does not alone prove that the H0 posterior is miscalibrated by a particular amount or that correction reverses benefit. That must be tested end-to-end.
- The correct generic sandwich form is `A^{-1} B A^{-T}/n` (or an equivalent convention based on definitions of A and B). The Consensus transcription for the Vrugt & Diks paper displayed a potentially erroneous `A^{-1} B^{-1} A^{-1}` form; verify against the original article before quoting.
- The FPC claim is only conditional on a fixed frame. For an iid random frame plus SRSWOR prefix, the between-frame term offsets the FPC reduction to first order in the simple setup. Do not repeat any claim that between-frame variability is negligible.
- Consensus gave only a compact paper matrix, not the requested exact estimands, sampling assumptions, detailed equations, or full-text dark-siren prior-art map. Its quotes and metadata are discovery leads, not publication-ready evidence.
- Some table/citation formatting in the returned memo was malformed or ambiguous. Check every title, DOI, volume/page, arXiv version, and quoted statement against publisher/arXiv sources before using it.
- The claim “no retrieved paper jointly models a shared inclusion nuisance inside a dark-siren likelihood” means only that this particular Consensus search did not retrieve one. It is **not** evidence that none exists.

## 11. Dark-siren prior art to verify directly

The following candidate prior art was identified in earlier project notes. The exact versions and metadata should be checked against the papers:

- Gair et al. (2023), “The Hitchhiker’s Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology,” *AJ* 166, 22, DOI 10.3847/1538-3881/acca78.
- Gray et al. (2023), “Joint Cosmological and Gravitational-Wave Population Inference Using Dark Sirens and Galaxy Catalogues,” *JCAP* 2023(12), 023, DOI 10.1088/1475-7516/2023/12/023.
- Gray, Williams & Papadopoulos (2026), “Expanding the Scope of Dark Siren Cosmology,” arXiv:2609.04991.
- Borghi et al. (2026), “Echoes from the Dark: Galaxy Catalog Incompleteness in Standard Siren Cosmology,” *A&A* 706, A199, DOI 10.1051/0004-6361/202557354; arXiv:2509.18243.
- Tagliazucchi et al. (2026), “Enhancing Dark Siren Cosmology via Gaussian Process Reconstruction of Incomplete Galaxy Catalogs,” arXiv:2608.24992.
- VanWyngarden et al. (2025/26), “How Low Can You Go: Constraining the Effects of Catalog Incompleteness on Dark Siren Cosmology,” arXiv:2511.04786.
- Dalang & Baker (2024), clustering-informed completion of invisible hosts, DOI 10.1088/1475-7516/2024/02/024.
- Dalang, Fiorini & Baker (2026), variance completion using LSS priors.
- Perna et al. (2025), astrophysical compact-binary host probability, *A&A* 698, A128, DOI 10.1051/0004-6361/202450840.
- Hanselman et al. (2025), “Gravitational-wave Dark Siren Cosmology Systematics from Galaxy Weighting,” arXiv:2405.14818.
- Leyde, Baker & Enzi (2026), “Cosmic Cartography II,” *JCAP* 01 (2026) 013, DOI 10.1088/1475-7516/2026/01/013.
- Datrier & Hendry (2026), MNRAS 549, stag828, DOI 10.1093/mnras/stag828.
- Li, Gray & Heng (2026), *ApJ* 1003(2), 176, DOI 10.3847/1538-4357/ae6100.

Classify each paper by whether it models a global completeness function, a missing-host/LOS mixture, clustering-based galaxy-field reconstruction, astrophysical host weighting, or per-galaxy inclusion propensities estimated from validation data. Similar-sounding “completeness” is not enough to establish overlap or novelty.

## 12. Historical pilot evidence — exploratory only

An earlier 200-replicate fixed-truth pilot under an older protocol (N=1,200, M=500, n=100/300/500, 30 events, gamma -1/0/+1, H0=70, prior [50,90], degree-1 fit) reported:

- Catalogue-only A median bias about -14 to -15; RMSE 13.8–15.5; coverage 0–1.5%.
- Parent-oracle B RMSE 3.793; coverage 69%.
- Known-probability IPW C RMSE 4.35–5.11; coverage 53–61%.
- Estimated-propensity MAP D RMSE 4.61–5.54; coverage 46.5–60.5%.
- Posterior-integrated E RMSE 4.50–5.37; coverage 49–63%.
- E was not uniformly better than C; at gamma=+1 the paired squared-error difference was +9.58, interval [5.93,13.23].
- A small 20-replicate degree-1/degree-2 check found degree-2 weight instability (node-level ESS minimum about 3.17), with no universal degree-2 advantage.

These are historical exploratory summaries, not confirmatory evidence. The protocol differs from the current independent-frame DGP. Do not combine these numbers with current branch results without matching definitions and provenance.

## 13. Required next research and engineering sequence

1. **Recover and reconcile the exact active worktree.** Confirm local branch/HEAD and whether the reported local branch has uncommitted files. Do not reset, clean, overwrite, rebase, or merge as part of recovery.
2. **Select the estimand.** Decide whether the H0 calculation needs `beta_star`, `beta_F`, `beta_P`, or `beta_P^KL`; define fixed/random frame and parent.
3. **Audit the inclusion posterior under misspecification.** Compare ordinary posterior spread against repeated-sampling sandwich spread in degree-1 misspecified and degree-2 correctly specified controls. If adding a calibrated posterior, define it as a separate method/arm.
4. **Resolve D3's role.** Is it a conditional design-variance diagnostic, a posterior calibration tool, or a new posterior construction? Do not double count uncertainty.
5. **Finite-difference test the normalized IPW Jacobian** and compare delta propagation with full posterior-node propagation.
6. **Re-audit measure consistency** across the KDE, event host mixture, and GW selection efficiency.
7. **Finish Arm E certification R5–R10:** posterior-level quadrature convergence orders 3–21; tail/ESS; prior sensitivity; near-separation; runtime/memory; independent end-to-end H0 posterior reconstruction.
8. **Adversarial selection-normalization tests:** missing/double correction, limiting cases, all arms, and shared nuisance.
9. **Numerical convergence:** H0 and redshift grids, KDE bandwidth/boundary, truncation/tails, log stability, MC SE.
10. **Calibration:** separate SBC/rank calibration under prior-predictive simulation from fixed-truth frequentist coverage; report truth-stratified results and Wilson/Monte Carlo intervals.
11. **Lifecycle/failure accounting:** freeze SUCCESS/NON_EVALUABLE/FAILED/NOT_RUN policy, preserve all planned denominators, independent audit and provenance.
12. **Prior-art and novelty audit:** read full text, check forward/backward citations, compare the exact statistical object rather than keywords.
13. **Independent scientific audit** before confirmatory simulations or manuscript submission.

Pass/fail tolerances must be justified and frozen before confirmatory runs. Do not invent universal ESS thresholds. Use exploratory pilots to characterize numerical error, then preregister tolerances based on the inferential target and decision relevance.

## 14. Current readiness verdict

- **Engineering:** substantial audit and reproducibility infrastructure exists in the remote PR stack; the user reports additional D3/Arm E implementation on a local-only branch. These are different evidence states.
- **Statistical uncertainty:** not settled. The most important open issue is whether the ordinary inclusion posterior has the right uncertainty for the chosen estimand under misspecification and validation-frame sampling.
- **Production H0 posterior:** NOT VALIDATED.
- **Confirmatory calibration:** NOT COMPLETED.
- **Novelty:** NOT VERIFIED.
- **Real-data / LVK / GWTC validation:** NOT PERFORMED.
- **Publication readiness:** NOT READY.

This checkpoint preserves the current reasoning and reported state. It is not a certification report and must not be used as evidence that code, tests, CI, calibration, or scientific conclusions were independently verified.
