# Conversation Research Checkpoint — 2026-11-26

> Purpose: preserve the research, implementation-status claims, mathematical cautions, literature leads, and next-step decisions discussed in the 2026-11-26 session. This is a research note, not a claim that local implementation has been merged, that the production posterior is validated, or that the literature review is complete.

## 1. Provenance and scope

- Repository: vedantkushwah033-blip/dark-siren-cosmology-robustness.
- Saved on a new branch, research/conversation-checkpoint-20261126, based on public main.
- The local worktree described below was not available through connected GitHub. Its files and test outputs are recorded as user-reported checkpoint information, not independently inspected or reproduced remote artifacts.
- User-reported local branch research/phase1-2-d3-arme-impl-20261126 at 7c7bf056fc91447cee49ff16f2a48c8f54ad88eb was not found as a GitHub ref. Reported local path: C:\Users\vedant\AppData\Local\Temp\opencode\kde_robustness_phase12_impl_20261126.
- No local files were uploaded, no local branch was pushed, and main was not changed. This branch contains this documentation checkpoint only.
- Consensus Deep Research could not be started because the connected service reported that included usage for the billing period had been exhausted. The literature leads below were gathered separately and remain preliminary, not a completed full-text systematic review.

### Scientific scope guardrails

The project is a synthetic low-redshift Euclidean dark-siren H0 methodology benchmark. It is not:
- an LVK/GWTC/public gravitational-wave data analysis;
- an observational measurement of H0;
- validation against a real galaxy survey;
- a claim of universal superiority of estimated-propensity IPW;
- a new general-purpose estimator unless a separate novelty audit establishes that claim.

Keep distinct:
1. pi_i = P(R_i=1 | x_i): probability galaxy i enters the observed catalogue;
2. h_i: astrophysical merger-host probability/rate weight;
3. P(D_GW=1 | theta,H0): GW detection probability and population-level selection normalization.

Oracle inclusion probabilities are simulation-only diagnostics, not an operational correction available for real surveys.

## 2. User-reported implementation checkpoint

### Reported worktree state

- Worktree: C:\Users\vedant\AppData\Local\Temp\opencode\kde_robustness_phase12_impl_20261126
- Branch: research/phase1-2-d3-arme-impl-20261126
- HEAD: 7c7bf056fc91447cee49ff16f2a48c8f54ad88eb
- User reported no reset, clean, delete, overwrite, main-branch change, commit, push, or merge.

### Reported modified / added files

- Modified: scripts/run_inclusion_posterior_smoke.py — manifest keys normalized using .as_posix() for Windows path portability.
- Added (untracked): src/dark_siren_robustness/d3_finite_validation.py — standalone SRSWOR FPC covariance, nested-prefix covariance, independent-SRSWOR zero cross-covariance, census n=M zero design variance, normalized IPW Jacobian, and delta-method propagation. Default guards reject ill-conditioned/non-positive-definite information matrices; pseudoinverse requires explicit allow_pseudoinverse=True.
- Added tests (untracked):
  - tests/test_d3_finite_validation.py
  - tests/test_d3_downstream_propagation.py
  - tests/test_d3_failure_guards.py
  - tests/test_arm_e_certification.py
  - tests/test_arm_e_eventwise_vs_joint_production.py
  - tests/test_gaussian_shared_parameter_validation.py
- Audit and continuation documents reported under:
  - research/project_completion/final_audit/
  - research/project_completion/next_gate_20261126/
  - research/project_completion/continuation_20261126/
  - research/project_completion/final_mission/INITIAL_REPOSITORY_AUDIT.md
  - research/project_completion/final_mission/IMPLEMENTATION_STATUS.md
  - research/project_completion/final_mission/FAILURE_LEDGER.csv

These local paths are a record of the user's report. They are not present in this remote branch unless separately uploaded or pushed.

### Reported tests (not independently rerun in this remote session)

Environment: Python 3.12.10, NumPy 2.5.5, SciPy 1.18.1, pytest 9.2.0, PYTHONPATH=src.

| Command group | Reported result |
|---|---:|
| Arm E certification + D3 finite validation + downstream propagation + Gaussian shared-parameter fixture | 16 passed, 12.67 s |
| Selection normalization + joint missing hosts + quadrature audit + inclusion posterior smoke | 61 passed, 13.94 s |
| Arm E eventwise-vs-joint production + D3 failure guards | 3 passed, 1.16 s |
| D3 failure guards + finite validation + downstream propagation | 8 passed, 3.45 s |
| Estimated inclusion + inclusion posterior + smoke | 35 passed, 8.01 s |

Reported ad-hoc Arm E quadrature probe: maximum posterior-summary delta relative to q=21 was 0.0471 at q=5 and 0.000126 at q=15 in the harder case. This is not a quadrature certification.

Reported independent in-session reconstructions (same agent, not external audit):
- Arm E R1 deterministic counterexample: joint integration 9 versus event-wise integration 25.
- Arm E marginalization reconstruction: logsumexp(node, axis=0, b=w[:,None]) matched arm_e_marginal_log_likelihood on a 25x21 case.
- Production joint vs event-wise Arm E marginal differed by max absolute difference 0.0147 in a small case.
- Separate Gaussian shared-parameter analytical fixture passed; it is not a reconstruction of the production model.

Reported environment issue: Python 3.15 SciPy _flapack DLL was blocked by Application Control; use Python 3.12. Historical Windows manifest path separator failure was reportedly fixed by .as_posix().

## 3. What is and is not established

### Reported as implemented / tested

- Standalone D3 SRSWOR, nested-prefix, independent-sample, census, FPC, normalized-IPW Jacobian, and delta-propagation behavior on targeted tests.
- Default D3 information-matrix guards reject ill-conditioned cases.
- Arm E small checks R1-R4; a production small-case joint-vs-eventwise difference was observed.
- Existing selection-normalization, joint-missing-host, quadrature-audit, and inclusion-posterior smoke suites passed in the reported environment.

### Explicitly not established

- D3 production posterior integration.
- A validated production H0 posterior incorporating finite-validation uncertainty, shared-parameter marginalization, and GW selection normalization consistently.
- Arm E R5-R10 certification.
- Adversarial selection-normalization certification.
- Quadrature sufficiency (q=5 vs q=15 remains unresolved).
- H0/redshift grid, KDE bandwidth/boundary, truncation/tail, log-stability, or Monte Carlo convergence certification.
- Confirmatory simulation campaign, SBC, or fixed-truth frequentist coverage campaign.
- Independent end-to-end production H0 posterior reconstruction.
- Independent scientific audit.
- Completed full-text prior-art/novelty matrix.
- Submission-ready manuscript or competition packages.

Current decision: production H0 posterior = NOT VALIDATED; novelty = NOT ESTABLISHED; confirmatory campaign = NOT AUTHORIZED; ASA–SBSS, ASA–SLDS, IEOM versions = NOT READY.

## 4. DGP and estimands to preserve

Latest H0 DGP as described by the user:
- Parent population size approximately N=4000.
- Approximately 30 detected events.
- Independent validation frame of size approximately M=1500, sampled from the same latent redshift distribution.
- Nested SRSWOR validation prefixes n in {100,300,1000}.
- Catalogue retention: pi_0(x) = expit(alpha + beta*x + gamma*(x^2-c)).
- x is standardized redshift, with approximate redshift support z in [0.01,0.40].
- Later stress design has beta around -3.9; degree-one logistic fit is misspecified for nonzero gamma; degree-two logistic fit is correctly specified under the stated DGP.
- c should be E_P[x^2] under the parent law; it must not be set to 1/3 without a distributional derivation.
- Host-rate weighting (example 1+1.5*x) is separate from catalogue inclusion.
- Parent, event, validation-frame, and validation-sample RNG streams are independent.
- Do not mix this independent-superpopulation validation-frame design with older finite-parent split designs in density-recovery code without explicitly distinguishing them.

Four coefficient targets must never be conflated:

1. Superpopulation pseudo-true root:
   E_P[X{R-q_beta(X)}]=0.
2. Finite validation-frame empirical pseudo-true root:
   sum_{i=1}^M X_i{R_i-q_beta(X_i)}=0.
3. Finite-parent realized-label MLE root:
   sum_{i=1}^N X_i{R_i-q_beta(X_i)}=0.
4. Finite-parent covariate-only KL projection:
   sum_{i=1}^N X_i{pi_0(X_i)-q_beta(X_i)}=0.

The realized-label MLE and covariate-only KL projection are not the same target.

## 5. D3 finite-population correction: mathematical decision and caveats

For a fixed validation frame F_M, SRSWOR sample size n, and a simple mean:
Var_d(Ybar_s | F_M) = (1 - n/M) S_M^2 / n.

For a regular logistic M-estimator, a first-order design-based approximation is:
Var_d(beta_hat | F_M) ~= A_M^{-1} [ (1 - n/M) S_{psi,M}/n ] A_M^{-T},
where A_M = M^{-1} sum q_i(1-q_i) X_i X_i^T, and S_{psi,M} is the finite-frame score covariance. This approximation requires a well-conditioned nonsingular information matrix, existence of the relevant root, controlled leverage/moments, adequate linearization, and SRSWOR independent of frame values. Census n=M has zero conditional design variance for a finite-frame target.

Critical target distinction: if F_M itself is iid from superpopulation P, and the validation sample is SRSWOR from that frame, then the unconditional sampled records are iid from P. By the law of total variance, to first order:
E_F[ Var_d(beta_hat | F_M) ] + Var_F( E_d[beta_hat | F_M] )
  ~= A^{-1} B A^{-T} [ (1-n/M)/n + 1/M ]
  = A^{-1} B A^{-T}/n.

Thus the between-frame component offsets the conditional FPC reduction for the superpopulation target. It is not generally negligible, particularly at high sampling fractions. Under misspecification, use the sandwich form with score covariance B; inverse Fisher information alone is not generally the correct repeated-sampling covariance.

This is a first-order result under the specified iid-frame/SRSWOR setup, not a universal identity for every validation design or estimand. Nested-prefix samples require their joint covariance structure; separate independent SRSWOR samples have zero design cross-covariance under the stated design, but shared frame variation can still induce unconditional dependence.

Practical consequence: retain D3 as a conditional finite-frame/design-uncertainty calculation if that is the intended target. Do not advertise the conditional FPC reduction as increased unconditional precision for estimating the superpopulation inclusion model.

## 6. D3–Arm E uncertainty reconciliation: blocking issue

Arm E's shared-parameter likelihood should have the structure:
L(H0) = integral { product_e L_e(H0, theta_I) } p(theta_I | V) d theta_I.

The same theta_I is shared across events when the inclusion mechanism is shared. In general:
integral { product_e L_e(H0,theta_I) } p(theta_I|V) dtheta_I
is not equal to
product_e integral L_e(H0,theta_I) p(theta_I|V) dtheta_I.
The latter behaves as if each event received an independently marginalized nuisance parameter and is a different model.

D3's delta-method covariance of fitted coefficients / normalized IPW weights may represent uncertainty already represented by p(theta_I|V) in Arm E. Do not add both uncertainty contributions without deriving that they correspond to distinct sources and a coherent joint generative model. Before production integration, specify:
- the estimand and target population/frame;
- what randomness D3 conditions on;
- what data enter p(theta_I|V);
- whether D3 covariance is design-based, superpopulation sandwich, or a posterior approximation;
- whether any D3 term is already integrated by Arm E;
- whether the cosmological likelihood depends on theta_I through both host weights and GW selection efficiency;
- how the resulting posterior avoids duplicate uncertainty.

A sandwich-calibrated posterior approximation or modular/cut construction is not automatically equivalent to the ordinary Bayesian posterior and requires separate justification and validation.

## 7. Selection normalization and the three mechanisms

If inclusion parameters affect the reconstructed population distribution used in GW detection efficiency, then the selection efficiency may depend on both H0 and theta_I:
xi_GW = xi_GW(H0, theta_I).

In that case the shared-parameter marginalization must consistently integrate the product of event likelihoods and selection-normalization terms under one shared theta_I, for the adopted fixed-detected-count likelihood. Do not assume xi_GW depends on theta_I without checking the actual generative model; derive it from the project code and DGP.

Keep separate:
- catalogue inclusion probability pi_i;
- astrophysical host-rate weight h_i;
- GW detection probability and its population-level normalization xi_GW.

Audit for missing correction, double correction, and limit behaviour. In particular, do not combine pure IPW correction with a missing-host mixture that corrects the same catalogue incompleteness a second time unless the likelihood derivation demonstrates why both terms are needed.

## 8. Preliminary literature leads from the session

These are leads and initial relevance assessments. Before citing them as definitive evidence in a manuscript, verify the final published version, read the full text, record exact equations/assumptions, and add contradictory work.

### Dark-siren catalogue incompleteness / reconstruction

1. Tagliazucchi et al. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs.”
   - arXiv: https://arxiv.org/abs/2608.24992
   - Relevant because it reconstructs incomplete galaxy redshift distributions with a GP and includes selection effects in hierarchical dark-siren inference.
   - The session's preliminary review noted a decomposition of the CBC redshift distribution into galaxy distribution and merger-rate weighting, plus a GW detection-selection normalization.
   - Reported precision gains in selected simulated completeness regimes are study-specific and must not be generalized to this project's low-z Euclidean DGP.
   - GP flexibility / calibration and overfitting are limitations to examine in the full text.
2. Borghi et al. (2026), “Echoes from the dark: Galaxy catalog incompleteness in standard siren cosmology,” A&A 706, A199.
   - DOI: https://doi.org/10.1051/0004-6361/202557354
   - Relevant to catalogue incompleteness, host weighting, and CHIMERA-based simulated dark-siren inference.
   - The session noted that incorrect host weighting need not cause a large H0 bias in every moderate-completeness configuration; consequences depend on completeness and modelling assumptions.
   - This is close prior art for broad claims about completeness correction. It does not by itself settle the narrower finite-validation / estimated-inclusion / IPW question.
3. Leyde et al. (2026), “Cosmic cartography. Part II. Completing galaxy catalogs for gravitational-wave cosmology.”
   - DOI: https://doi.org/10.1088/1475-7516/2026/01/013
   - Relevant to reconstructing galaxy fields and completing catalogues using spatial correlations and selection information. Assess full assumptions and applicability before drawing comparisons.
4. “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant” (2026).
   - MNRAS landing page: https://academic.oup.com/mnras/article/549/1/stag828/8666390
   - Relevant to completeness diagnostics in dark-siren analyses and the gwcosmo context; examine direction-dependent completeness and limiting-magnitude assumptions.

### Hierarchical selection / population inference

5. Mandel, Farr & Gair (2019), “Extracting distribution parameters from multiple uncertain observations with selection biases.”
   - MNRAS: https://academic.oup.com/mnras/article/486/1/1086/5421632
   - Foundational reference for hierarchical population inference with selection effects. Verify the precise conditioning and selection-factor convention against the project's likelihood.

### Misspecified likelihoods / weighting

6. White (1982), “Maximum Likelihood Estimation of Misspecified Models.”
   - Reference PDF used in the session: https://faculty.utrgv.edu/diego.escobari/teaching/Econ8370/Papers/White%281982%29Econometrica-Corey.pdf
   - Relevant to pseudo-true parameters and sandwich covariance under likelihood misspecification.
7. Zhou, Matsouaka & Thomas (2020), propensity-score weighting under limited overlap.
   - https://journals.sagepub.com/doi/10.1177/0962280220940334
   - The session's preliminary review noted comparisons in which overlap, matching, and entropy weighting can outperform conventional IPW in some limited-overlap / misspecification settings. Its estimands and application differ from dark-siren host-distribution reconstruction; treat as a warning and methodological lead, not direct evidence for a preferred dark-siren estimator.
8. Propensity-score weighting tutorial / uncertainty treatment (Wiley page surfaced in session).
   - https://onlinelibrary.wiley.com/doi/full/10.1002/sim.10078
   - Examine how estimated propensity uncertainty, overlap, and variance estimation are treated; confirm bibliographic details and exact claims from the full text.
9. Talts et al. (2018), “Validating Bayesian Inference Algorithms with Simulation-Based Calibration.”
   - arXiv: https://arxiv.org/abs/1804.06788
   - Relevant to SBC. SBC checks computational calibration under the assumed generative model; it does not establish physical realism or replace fixed-truth frequentist coverage studies.
10. Stan User's Guide, Simulation-Based Calibration.
    - https://mc-stan.org/docs/2_36/stan-users-guide/simulation-based-calibration.html
    - Practical reference for rank-based SBC design and interpretation.

## 9. Revised research question and safe novelty wording

Preferred research question:
> Under what conditions does estimating galaxy-catalogue inclusion probabilities improve dark-siren inference, and under what conditions do misspecification, finite validation samples, and IPW instability eliminate or reverse benefit?

More operational framing:
> In a controlled synthetic low-redshift dark-siren benchmark, how do finite validation samples, inclusion-model misspecification, and unstable estimated inverse-probability weights affect the calibration of shared-parameter, selection-normalized multi-event H0 inference?

Potential contribution, not yet established:
> A controlled synthetic benchmark of how finite validation samples, inclusion-model misspecification, and unstable estimated weights affect calibration in shared-parameter, selection-normalized multi-event dark-siren inference.

Do not claim:
- catalogue-incompleteness correction is novel;
- IPW is universally better;
- D3 is a new estimator;
- the production H0 posterior is validated;
- the current work is ready for publication or competition submission.

## 10. Recommended next gates

### Gate 1 — D3 estimand and uncertainty-target decision
Write a formal target/conditioning table for the four coefficient targets in Section 4. Derive conditional design covariance, between-frame variance, and unconditional sandwich variance for the exact validation design. State what D3 estimates and what it does not estimate.

### Gate 2 — D3–Arm E double-counting audit
Derive a single joint generative model for parent galaxies, catalogue inclusion, validation data, event hosts, GW observations, and GW detection. Identify whether D3 is a separate source of uncertainty or an approximation to uncertainty already represented in p(theta_I|V). No production integration until this is resolved.

### Gate 3 — Arm E R5–R10
Certify posterior-level quadrature convergence across q=3,5,9,15,21 (or a justified alternative); inspect posterior tails, ESS, difficult geometries / near-separation, prior sensitivity, runtime/memory, and compare against an independent end-to-end reconstruction. The q=5 vs q=15 discrepancy in the reported probe is not resolved by the q=15 vs q=21 result alone.

### Gate 4 — Selection-normalization audit
Independently derive the event likelihood and GW selection normalization from the exact DGP. Construct adversarial tests for missing correction, double correction, limiting cases, shared nuisance parameters, and event-wise vs joint marginalization.

### Gate 5 — Numerical certification
Predeclare acceptance criteria for H0 and redshift grids, KDE bandwidth and boundary treatment, support truncation/tails, log-domain stability, quadrature, Monte Carlo error, optimizer failures, and non-evaluable runs.

### Gate 6 — Calibration and diagnostic prediction
Run separate SBC and fixed-truth repeated-sampling coverage campaigns. Preserve a failure ledger with SUCCESS / NON_EVALUABLE / FAILED / NOT_RUN states. Train/evaluate any pre-posterior diagnostic prediction on held-out simulation replicates, never on the same replicates used to select thresholds.

### Gate 7 — Prior-art / novelty audit
Read full texts, snowball references, and construct a claim-by-claim matrix distinguishing established methods, adjacent work, true gaps, and unsupported novelty claims.

### Gate 8 — Reproducibility and independent review
Freeze configs; checksum inputs and outputs; capture environment; independently reconstruct posterior summaries from saved artifacts; commission an independent scientific audit before manuscript or competition submission.

## 11. Status snapshot

| Area | User-reported checklist completion | Scientific status |
|---|---:|---|
| Software implementation | 35% | Partial; D3 production integration absent |
| Numerical verification | 22% | Partial; convergence not certified |
| Scientific validation | 5% | Not validated |
| Literature / novelty | 3% | Not established |
| Manuscript / readiness | 10% | No submission-ready version |

These percentages are the user's reported checklist counts (20 items per area), not an independently designed maturity metric.

Bottom line: the reported work has useful software and mathematical infrastructure, but the production H0 posterior remains NOT VALIDATED. The confirmatory campaign should remain blocked until the D3 estimand and D3–Arm E uncertainty reconciliation are settled, followed by selection-normalization and numerical certification.
