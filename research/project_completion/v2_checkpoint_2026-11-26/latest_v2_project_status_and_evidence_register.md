# Latest v2 Project Status and Evidence Register — 2026-11-26

> This is a durable conversation-derived checkpoint. It preserves the latest user-supplied OpenCode/OpenHands project report and the subsequent audit decisions. It is not a substitute for the actual local worktree, raw artifacts, configs, or source code. All results below are labelled by provenance.

## 1. Provenance: do not merge incompatible checkpoints

### Latest reported local v2 state
- Worktree: `/workspace/v2`
- Branch: `research/v2-finite-frame-arm-e-2026-11-25`
- HEAD: `e06b88eb87fc4851105f34f9680de2f9deece7e6`
- Runtime reported: Python 3.13.15, NumPy 2.5.5, SciPy 1.18.1, pandas 3.1.1, pytest 9.2.0.
- Reported work was uncommitted; main untouched; no reset, clean, deletion, merge, or push in that local phase.

### Remote verification performed in this chat
- GitHub repository exists: `vedantkushwah033-blip/dark-siren-cosmology-robustness`.
- Default branch is `main`.
- Branch search did not return `research/v2-finite-frame-arm-e-2026-11-25`.
- Commit `e06b88eb87fc4851105f34f9680de2f9deece7e6` exists on the remote history. Its commit message is “Revert accidentally targeted default-branch pilot file; keep research work on draft branch”; its diff deletes `scripts/run_degree_matched_five_arm_quadrature_sensitivity.py`.
- Therefore, the reported local v2 worktree and its uncommitted files are NOT established as present on GitHub. The commit SHA alone does not prove that the local v2 work is saved remotely.
- This checkpoint is being stored on the separate documentation branch `research/conversation-checkpoint-20261126`. It does not alter main and does not upload local source files or simulation artifacts.
- An earlier checkpoint exists at `research/project_completion/conversation_checkpoint_2026-11-26.md`; this document supplements it and does not replace it.

## 2. Project scope and scientific guardrails

The project is a synthetic, low-redshift, Euclidean dark-siren H0 methodology benchmark only. It is not LVK/GWTC/public-data analysis, not an observational H0 measurement, not real galaxy-survey validation, and not a claim of universal estimator superiority.

Keep distinct:
1. Catalogue inclusion probability: pi_i = P(R_i=1 | x_i).
2. Astrophysical host/rate weighting: h_i.
3. GW detection probability and population selection normalization: P(D_GW=1 | theta,H0), alpha or xi_GW.

Oracle inclusion probabilities are simulation diagnostics only. Do not double-correct catalogue incompleteness by combining pure IPW with a missing-host mixture unless the likelihood derivation justifies both.

Preferred research question:
“Under what conditions does estimating galaxy-catalogue inclusion probabilities improve dark-siren inference, and under what conditions do misspecification, finite validation samples, and IPW instability eliminate or reverse benefit?”

More operational framing:
“In a controlled synthetic low-redshift dark-siren benchmark, how do finite validation samples, inclusion-model misspecification, and unstable estimated inverse-probability weights affect calibration of shared-parameter, selection-normalized multi-event H0 inference?”

No novelty claim is verified.

## 3. Latest v2 design reported by the user

Frozen config named `configs/finite_validation_probability_pilot_v1.json`:
- Curvature gamma: {-2, -1.2, 0, +1.2, +2}.
- Validation sizes n: {100, 300, 1000}.
- Parent size: 4000.
- Both degree-1 and degree-2 inclusion models are fitted; degree-1 is the preregistered primary arm in the reported pilot.
- Reported pilot: 20 replicates per curvature/validation-size condition, degree-1 primary arm.
- Five arms, 5 curvatures x 3 sizes x 20 replicates x 5 arms = 1,500 arm-cell outcomes.
- The config is reported as authoritative over an earlier brief that incorrectly listed {-2,-1,0,+1,+2}.

### Critical unresolved validation-frame/FPC conflict
Earlier project descriptions in this conversation specify an independent validation frame of approximately M=1500, sampled from the same latent redshift distribution. The latest D3 production report instead states M=4000 and FPC values 0.975, 0.925, 0.75 for n=100,300,1000.

These are not interchangeable:
- If the actual SRSWOR frame is M=1500, FPC = 1 - n/1500.
- If the actual frame is M=4000, FPC = 1 - n/4000.

Do not choose a denominator to match a report. Inspect frozen config, frame-generation code, sampling code, and artifacts. Determine whether 1500 is the validation-frame size, the parent size, a subset, or a different experiment. Until then, D3's reported FPC denominator and the claimed D3 production validation remain disputed.

## 4. Latest IPW stability and misspecification results (user-reported)

The user supplied a report titled `ipw_stability_and_misspecification.md`.

### Weight diagnostics (median over 20 replicates)
Selected reported rows:
| gamma | n | arm | n_cat | ESS | min ESS | max normalized weight | top-1% mass |
|---:|---:|---|---:|---:|---:|---:|---:|
| -2 | 100 | A | 1242 | 1242 | 1180 | 0.0008 | 0.0104 |
| -2 | 100 | C | 1242 | 880 | 824 | 0.0026 | 0.0322 |
| -2 | 100 | D | 1242 | 922 | 706 | 0.0022 | 0.0282 |
| -2 | 1000 | C | 960 | 682 | 625 | 0.0033 | 0.0319 |
| -2 | 1000 | D | 960 | 716 | 583 | 0.0028 | 0.0278 |
| -1.2 | 100 | C | 1201 | 971 | 920 | 0.0019 | 0.0241 |
| -1.2 | 1000 | D | 930 | 771 | 666 | 0.0022 | 0.0222 |
| 0 | 100 | C | 1162 | 1086 | 1023 | 0.0013 | 0.0154 |
| 0 | 1000 | D | 899 | 837 | 772 | 0.0017 | 0.0154 |
| +1.2 | 100 | C | 1146 | 1141 | 1055 | 0.0009 | 0.0109 |
| +1.2 | 1000 | D | 886 | 882 | 814 | 0.0012 | 0.0114 |
| +2 | 100 | C | 1142 | 1128 | 1066 | 0.0010 | 0.0120 |
| +2 | 1000 | D | 883 | 867 | 809 | 0.0017 | 0.0147 |

A = catalogue-only; C = known-probability IPW; D = plug-in estimated IPW in the user's arm naming. Preserve the exact arm map from the source config/code before reusing these labels.

### Curvature-level co-movement
| gamma | catalogue-only H0 bias (km/s) | ESS(C)/n_cat | max normalized weight |
|---:|---:|---:|---:|
| +2 | +0.64 | 0.99 | 0.0010 |
| +1.2 | -0.36 | 0.99 | 0.0009 |
| 0 | -2.25 | 0.93 | 0.0013 |
| -1.2 | -4.56 | 0.81 | 0.0019 |
| -2 | -6.29 | 0.71 | 0.0033 |

This is qualitative co-movement over five design points, not a fitted predictive relationship. No held-out conditions; no threshold calibrated.

### ESS does not diagnose the residual bias
At gamma=-2, n=1000, reported D-arm residual bias is -0.62 km/s with ESS=716/960 (~0.75). At gamma=0, bias is -0.22 km/s with nearly similar ESS. The bias differs by about a factor of three despite similar ESS. The proposed explanation is logistic degree-1 misspecification, which ESS does not encode. Do not infer that ESS is useless generally; the result only shows that it cannot discriminate these two cases.

### Misspecification measurements
- Parent second moment E_f[x^2] = 0.5923826935; closed-form and quadrature reportedly agree to 1e-15.
- Using wrong centering 1/3 changes max absolute pi_0 by 0.0968.
- Wrong centering changes IPW mean-redshift bias by 4.56e-3.
- Degree-1 fit to degree-2 truth has pseudo-true moment residual 1.46e-2 versus 2.5e-18 when correctly specified.
- Degree-1 pseudo-true coefficients: (2.05786, -3.94546).
- Host-density L1 distance from parent under misspecification: 0.17, reportedly insensitive to KDE bandwidth 0.4x–2x Scott.
- Interpretation: bandwidth tuning cannot repair inclusion-model misspecification in this setup.

### Failure accounting
- SUCCESS 1500
- NON_EVALUABLE 0
- FAILED 0
- NOT_RUN 0

Thus the frozen pilot did not exercise separation or positivity-failure regimes. Stronger curvature (e.g. |gamma|>=4), smaller n, or separating covariates are design changes requiring preregistration. Do not tune the existing confirmatory design post hoc.

## 5. Latest test execution log (user-reported)

Environment: `/workspace/v2`, `.venv/bin/python`, Python 3.13.15, NumPy 2.5.5, SciPy 1.18.1, pytest 9.2.0.

### Full suite (two pre-existing collection-error modules excluded)
Command:
```bash
.venv/bin/python -m pytest tests -q -p no:cacheprovider \
  --ignore=tests/test_audit_multiseed_grid_sensitivity.py \
  --ignore=tests/test_posterior_calibration_uncertainty.py
```
Reported result: 475 passed in 116.38s, exit 0. Earlier same command: 466 passed; +9 were `tests/v2/test_ledger.py`.

The two ignored modules fail at collection because they import missing files:
- `scripts/audit_multiseed_grid_sensitivity.py`
- `scripts/analyze_posterior_calibration_uncertainty.py`

User reports these collection errors were pre-existing on main. This has not been independently verified in this chat.

### v2 subset
`pytest tests/v2 -q` -> 72 passed in 52.01s:
- ledger 9
- parent 9
- finite-frame D3 9
- nested/independent 5
- IPW 7
- selection 9
- Arm E 11
- quadrature 6
- pipeline 7

### Targeted development runs
- ledger: 9 passed, 0.17s.
- quadrature + Arm E + pipeline after Newton-polish patch: 24 passed, 12.71s.

### Scripts
- `scripts/v2/run_certification_probes.py`: exit 0, about 47s, three JSON artifacts.
- `scripts/v2/run_calibration_pilot.py --reps 20 --gammas=-2,-1.2,0,1.2,2 --ns 100,300,1000 --n-parent 4000 --events 30 --h0-points 201 --order 5`: exit 0, 857.4s, 1500/1500 SUCCESS.

### Failures and fixes reported
1. Real package defect: Laplace mode search returned ABNORMAL and raised RuntimeError, killing first campaign after about 4.5 minutes. Newton polish stage added; large residual gradient becomes ArithmeticError; RuntimeError mapped/caught by pipeline; per-replicate exception handling. User reports code path and test now cover this.
2. CLI parser required `--gammas=-2,-1.2,...` syntax.
3. Probe script event-axis reduction shape bug: sum(axis=0) replaced sum(axis=1).
4. Probe Hessian sign error: information matrix should be X'WX + I/sigma^2, covariance inv(info).
5. Ad-hoc diagnostic extractor used wrong key; expected max_normalized_weight rather than max_weight.
Items 2–5 were reported as analysis-script issues, not package defects.

### Not run
- SBC not implemented.
- 100-replicate confirmatory campaign blocked.
- AIPW arm not implemented.
- Full-text literature screening not completed; report says web unavailable and arXiv API only.
- 230 MB CI artifact not downloaded/reanalyzed.

## 6. Latest scientific results and caveats (user-reported)

### Calibration pilot
- 20 replicates x 5 curvatures x 3 validation sizes x 5 arms = 1,500 cells.
- Catalogue-only H0 bias: -6.33 km/s at gamma=-2 to +0.64 at gamma=+2.
- Catalogue-only 68% coverage: 0.10–0.60.
- True-pi/oracle IPW bias: -0.10 to -0.28 km/s; 95% coverage 0.95–1.00.
- Estimated IPW residual bias at gamma=-2: approximately -0.62 to -1.61 km/s depending on validation size/arm; reported to shrink with n.
- Monte Carlo SE for bias approximately 0.55–0.73 km/s with only 20 replicates, so sub-1.5 km/s differences and D-vs-E ordering are not resolved.
- Critical finding: oracle arm's nominal 68% coverage is only approximately 0.60–0.65. Therefore estimated inclusion uncertainty is not the sole cause of undercoverage. No arm is calibrated at nominal 68%.
- Do not call the posterior calibrated.

### D3 integration report
User reports a new report `d3_production_integration_report.md` with:
- `finite_frame.census_coefficients` / `fit_logistic` for coefficient vectors.
- `finite_frame.design_covariance` for design covariance.
- `ipw.inclusion_probability`, `normalized_ipw_weights`, `normalized_weight_jacobian`, `propagate_weight_covariance`, `weighted_kde`.
- `KDE covariance = K Sigma_w K'`.
- `selection_likelihood.selected_log_likelihood` and posterior summary functions.
- Arm E uses quadrature over q(beta|V) in `arm_e.joint_and_product`.
- Delta-method vs Monte Carlo: weight covariance relative gap 3.5% (1500 draws); KDE covariance relative gap 2.8% (600 draws); Jacobian vs finite differences max discrepancy 1.8e-13.
- Three seeds, frame_size=2500, n=250, 12 events, quadrature order 10:
  - seed 20261125: A 65.07 (2.72), C 67.00 (2.04), D 66.74 (2.13), E joint 66.60 (2.20).
  - seed 20261126: A 65.69 (3.54), C 68.34 (2.50), D 68.21 (2.55), E joint 67.92 (2.73).
  - seed 20261127: A 70.46 (3.56), C 72.89 (2.55), D 72.62 (2.60), E joint 72.39 (2.71).
- E posterior sd exceeds D in all 3 seeds, consistent with marginalization but not evidence for universal ordering.
- D tracks C closely in these cases.
- D3 pilot reports FPC = 1 - n/M with M=4000; this conflicts with earlier M approximately 1500 validation-frame description and must be resolved from code/config/artifacts.
- No D3 covariance added to Arm E posterior, deliberately, to avoid double counting.
- D3 integration report verdict says Route A and Route B complete, but the FPC discrepancy and actual executable trace still need independent audit before accepting “production integration complete”.

### Selection and numerical certification reports
- Selection normalization limiting cases: constant pi and complete catalogue collapse to unweighted.
- Reported missing correction bias: -5.17 km/s.
- Reported double correction bias: +2.72 km/s and posterior sd narrows 5.95 to 4.44.
- Detection saturation: at d_th=250 Mpc, max spread in log P(detect) across coefficient nodes = 0.105 nats (selection term nearly inert); at 800 Mpc = 0.573 nats.
- Grid: 1601 points on H0 range 50–90 identical to 3201; wide same-spacing grid differs by 2.4e-4 (phase, not resolution).
- Quadrature hard geometry: q5 median error 0.0338 km/s; q15 0.0024 km/s relative to q21. Historical earlier q5=0.0471 and q15=0.000126 were not reproduced; new q15 error is 19x larger than historical report.
- Selection audit passes limiting cases but detection-limited validation remains a caveat.
- Quadrature certification is conditional, not universal.

## 7. Publication readiness (user-reported assessment)

One paper only, not three. The evidence currently supports a synthetic benchmark of catalogue-inclusion correction and its failure modes, not a calibrated cosmological inference method.

Not publication-ready:
- SBC absent.
- Novelty not verified; read arXiv:2609.04991 and arXiv:2507.12171 in full (the user's assessment says abstracts only were inspected).
- AIPW sixth arm absent.
- Confirmatory >=100 replicates/cell absent; 20-replicate pilot has wide MC uncertainty.
- Detection-limited selection regime needs stronger audit (d_th around 800 Mpc proposed as diagnostic).
- Ledger input hashing/resume not complete at campaign scale.
- Arm E R7 prior sensitivity open.
- Discrete-host vs reflected-KDE reconciliation open.
- Failure regimes (separation, positivity) not exercised by frozen pilot.
- No JCAP-level cosmological realism (Euclidean, flat H0 prior, simplified selection).
- No manuscript drafted.

Claims explicitly prohibited:
- “Calibrated posterior.”
- “IPW is superior.”
- “First to propagate inclusion uncertainty.”
- “Observational”, “LVK”, “GWTC”, or “H0 measurement.”
- “JCAP-ready.”

Venue fit reported:
- Frontiers in Astronomy and Space Sciences (Astrostatistics): potentially good fit after SBC, novelty audit, AIPW.
- Astronomy and Computing: infrastructure fit, but resume/input hashing and novelty remain blockers.
- JCAP: poor fit without cosmological distance relation, wCDM, magnitude-dependent completeness, realistic selection.

## 8. Mathematical decisions to preserve

### D3 versus Arm E
- D3: fixed-frame frequentist SRSWOR design uncertainty, conditional on the observed finite validation frame.
- Arm E: Bayesian posterior uncertainty q(beta|V,U,M), conditional on validation data, model, and prior.
- These are different statistical objects and conditioning sets.
- Do not combine as Sigma_D3 + Sigma_Bayesian or inject D3 covariance into Arm E without a formally derived non-overlapping decomposition and coherent joint generative model.

### Shared-parameter Arm E likelihood
The correct shared-beta structure is:
p(H0|D,V) proportional to p(H0) times
integral [ product_j { p(D_j|H0,beta) / alpha(H0,beta) } ] q(beta|V,U,M) d beta.

The same beta is shared across all events. Do not replace this by product_j integral L_j(H0,beta) q(beta|V) d beta.

If beta affects host distribution and therefore GW detection efficiency, beta-dependent alpha(H0,beta) remains inside the shared-beta integral/node likelihood.

### Four coefficient estimands
Keep distinct:
1. Superpopulation pseudo-true root E_P[X{R-q_beta(X)}]=0.
2. Finite validation-frame empirical pseudo-true root sum_{i=1}^M X_i{R_i-q_beta(X_i)}=0.
3. Finite-parent realized-label MLE root sum_{i=1}^N X_i{R_i-q_beta(X_i)}=0.
4. Finite-parent covariate-only KL projection root sum_{i=1}^N X_i{pi_0(X_i)-q_beta(X_i)}=0.

Finite-parent realized-label MLE and covariate-only KL projection are not identical.

## 9. Main unresolved scientific question: oracle undercoverage

Oracle inclusion probabilities still yield approximately 0.60–0.65 coverage for nominal 68% intervals. This is the top priority because it shows the deficit is not solely due to estimated inclusion probabilities.

Potential causes remain hypotheses until discriminated:
- likelihood/data-generating model mismatch;
- interval construction and posterior summarization;
- flat H0 prior or prior boundaries;
- event-generation and conditioning assumptions;
- host density / KDE approximation;
- GW detection selection and normalization;
- quadrature/grid numerical error;
- finite-event uncertainty and Monte Carlo precision.

Do not assume one cause. The next validation design should first reconcile the data-generating process and inferential model, then implement SBC (prior-drawn truths) and separate fixed-truth frequentist coverage experiments. SBC and fixed-truth coverage answer different questions.

## 10. Next-phase priorities (no campaign authorized)

1. Independent audit of exact validation-frame generation and FPC denominator.
2. Trace D3 uncertainty through the full executable beta -> pi -> weights -> density/KDE -> event likelihood -> selection normalization -> H0 posterior chain.
3. Reproduce oracle undercoverage from saved artifacts and construct a fault tree; distinguish demonstrated causes from hypotheses.
4. Reconcile inferential model and generative model before SBC.
5. Design SBC correctly: draw H0 and nuisance parameters from actual prior, simulate validation/catalogue/GW events and detection consistently, calculate posterior ranks with numerical CDF handling, report MC uncertainty.
6. Keep fixed-truth frequentist coverage as a separate experiment.
7. Complete Arm E R7 prior sensitivity and remaining certification only after mathematical review.
8. Audit detection-limited selection at the proposed higher threshold.
9. Only then revisit AIPW, preregister stress conditions for separation/positivity, and consider a confirmatory campaign.
10. Complete full-text novelty/prior-art review.

**Current decision: NO-GO for confirmatory campaign.** No new large campaign should run until FPC and oracle undercoverage are resolved.

## 11. Evidence labels

Every future report should label claims:
- USER-REPORTED: stated in user-supplied OpenCode/OpenHands output, not reproduced here.
- REMOTE-VERIFIED: checked against GitHub in this conversation.
- CODE-VERIFIED: inspected actual source in the active worktree.
- TEST-VERIFIED: exact test command run and result captured.
- ARTIFACT-VERIFIED: underlying output/config/manifest inspected.
- INDEPENDENTLY RECONSTRUCTED: recomputed independently from raw inputs/artifacts.
- SCIENTIFICALLY VALIDATED: assumptions, DGP, inference, uncertainty and repeated-sampling behavior all assessed.

In this chat, repository existence/default branch, the existence and diff of commit e06b88e, and the existing conversation checkpoint file were REMOTE-VERIFIED. The 475 tests, 1500-cell campaign, D3 Monte Carlo comparisons, coverage figures, and v2 source implementation remain USER-REPORTED until the local worktree/artifacts are made available and inspected.
