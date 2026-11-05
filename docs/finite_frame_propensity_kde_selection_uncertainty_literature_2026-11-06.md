# Literature evidence extension: finite-frame propensity uncertainty, weighted KDE, and selection normalization

**Date:** 2026-11-06  
**Status:** Research synthesis and design recommendations; not a completed numerical experiment.  
**Scope:** Synthetic dark-siren host reconstruction with a finite parent frame, Bernoulli catalogue retention conditional on parent redshift, and an SRSWOR validation sample drawn from that same frame. Catalogue/validation overlap is identified by persistent parent-unit IDs.

## Executive adjudication

The literature supports three distinct conclusions.

1. **Finite-population propensity uncertainty must include both sampling phases and their dependence.** Work on probability-reference/non-probability sample integration explicitly recommends paired, design-respecting replication and notes that overlap/membership information matters. This is close to the present design, but it is not identical: the present validation sample is SRSWOR from the same finite parent frame as the Bernoulli-thinned catalogue, with validation revealing redshift and retention status for sampled units. Do not import an independent-samples bootstrap unchanged.
2. **Extreme inverse-probability weights are an identification and estimand problem as well as a numerical problem.** Truncation, trimming, and overlap-focused weights are not interchangeable. In the present one-sample inclusion problem, standard causal overlap weights cannot simply be transplanted: the causal literature's treatment/control overlap estimands are not automatically estimators of the full parent redshift density.
3. **The selection denominator must match the event-count model.** A likelihood conditioned on a fixed detected-event count uses a per-event detection-normalization factor. A rate-aware point-process likelihood instead models the count through the expected number of detections. These are alternative factorizations, not corrections to stack on top of one another.

The directly relevant dark-siren literature has also moved recently: Tagliazucchi et al. (2026) use a hierarchical Bayesian likelihood and a Gaussian-process reconstruction of incomplete galaxy redshift structure, rather than the present logistic inclusion/IPW/KDE route. This should be added to the related-work comparison before making any novelty statement.

---

## A. Finite-population overlap and propensity-score propagation

### A1. Evidence table

| Source | Setting and supported result | Transfer to this project | Limitation / caution |
|---|---|---|---|
| Conti & Marella (first published 2025; *International Statistical Review*, 94(2), 264–290; DOI 10.1111/insr.12614) | Probability reference sample A plus non-probability sample B. Pseudo-inclusion probabilities for B are estimated from A. The paper recommends bootstrap variance estimation that accounts for both samples' sampling designs and dependence among units; it explicitly discusses observing B-membership indicators in A and record linkage for overlap. | Strong methodological analogue for propagating inclusion-model uncertainty and tracking overlap. It supports paired replication that repeats propensity fitting and downstream estimation, rather than conditioning on one fitted propensity model. | Its baseline setup assumes independent probability and non-probability samples; the present validation and catalogue are generated from the same finite frame and are not independent samples. The bootstrap must reproduce this project's actual joint design, not copy their independent-sample setup. |
| Mashreghi, Haziza & Léger (2016), “A survey of bootstrap methods in finite population sampling,” *Statistical Surveys* 10, 1–52. | Reviews bootstrap methods for finite-population sampling. | Supports design-aware replication as a family of methods rather than assuming ordinary iid row resampling. | Select a method only after specifying the sampling design and target variance; a generic bootstrap label is insufficient. |
| Kostouraki et al. (2024), *Statistics in Medicine* 43(13), 2672–2694; DOI 10.1002/sim.10078 | Propensity weighting is a two-step estimator. Variance methods must account for estimated propensity scores as well as the weighted estimator. The article discusses analytic and bootstrap approaches and the limitations of treating estimated scores as known. | Supports rerunning the entire inclusion-fit → IPW → density reconstruction → downstream inference chain within each replicate. | Causal treatment weighting is not the same estimand as estimating a catalogue inclusion mechanism. Its variance formulas should not be copied as if they were derived for this finite-frame design. |
| Ridgeway et al. (2015; published issue 2015; *Sociological Methods & Research*), “Propensity Score Analysis with Survey Weighted Data.” | Shows that survey design weights can matter in both propensity estimation and downstream analysis. | Relevant if future versions add unequal-probability parent-frame or validation sampling. | Current validation design is SRSWOR, so survey-weighted logistic fitting is not automatically required. |
| Tagliazucchi et al. (2026), *Physical Review D* 114, 083523; DOI 10.1103/vmzf-1wxq | Reconstructs incomplete galaxy line-of-sight redshift structure with a GP, using a hierarchical Bayesian likelihood that accounts for survey selection; evaluated on mock catalogs. | Direct dark-siren comparator: incomplete-catalogue uncertainty and reconstructed host-redshift structure propagate into H0 inference. | It is not the same finite-validation logistic inclusion experiment; compare assumptions, information sources, target density, and calibration rather than claiming equivalence. |

### A2. Correct joint design for this repository

Let the finite parent frame be fixed as \(Z=(z_1,\ldots,z_N)\). For each parent unit, catalogue retention is Bernoulli conditional on redshift, \(R_i\sim\mathrm{Bernoulli}(\pi_i(\gamma))\). Independently, validation IDs \(S\) are selected by SRSWOR with fixed size \(m\). Validation observes \((z_i,R_i)\) for every \(i\in S\), including units not retained in the catalogue. The catalogue observes \((z_i,R_i=1)\) for retained units. The overlap is \(C\cap S\), represented by IDs.

Conditional on the fixed parent redshifts, the joint design is

\[
p(R,S\mid Z,\gamma)=
\left\{\prod_{i=1}^{N}\pi_i(\gamma)^{R_i}
[1-\pi_i(\gamma)]^{1-R_i}\right\}
\frac{\mathbf 1\{|S|=m\}}{\binom Nm}.
\]

For the observed union \(U=C\cup S\), each unit's retention outcome contributes once:

\[
L_R(\gamma;R_U,Z)=
\prod_{i\in U}\pi_i(\gamma)^{R_i}
[1-\pi_i(\gamma)]^{1-R_i}.
\]

This is a **simulation-oracle conditional finite-frame likelihood** because it conditions on the complete parent redshift frame. It is not yet the observer-facing likelihood: the latter must integrate over parent redshifts not revealed by either the catalogue or validation sample, and must state exactly what is observed for each ID.

Important implementation implications:

- The validation sample's SRSWOR probability is a combinatorial design term, not a product of independent Bernoulli validation indicators.
- The same ID's redshift and retention flag must not be counted twice when forming a joint likelihood.
- A modular/cut analysis may intentionally prevent catalogue information from feeding back into the inclusion model, but that is a declared inferential target, not the full joint posterior.
- Because validation is drawn without replacement, the variance of fitted logistic coefficients is not generally obtained by multiplying an iid/SRSWR covariance by a scalar finite-population correction. Derive it for this nonlinear estimator and this joint observation map, or estimate it by design-respecting replication.

### A3. Recommended uncertainty propagation experiment

Use common parent frame, retention vector, catalogue IDs, GW events, and seeds across methods. For each validation replicate:

1. Draw validation IDs by SRSWOR, retaining the exact overlap with catalogue IDs.
2. Fit the inclusion model using the validation observations only for the declared modular baseline.
3. Predict inclusion probabilities for catalogue units; calculate inverse weights.
4. Rebuild the weighted host-redshift KDE, including its bandwidth rule and boundary correction.
5. Apply merger-rate weighting and normalize the host density.
6. Recompute all event integrals and the selection normalization using that same density realization.
7. Save fit status, all diagnostics, density error, posterior summaries, and failure reason.

Report separately:
- the distribution over repeated SRSWOR validation draws conditional on one fixed parent frame and catalogue;
- the simulation-oracle finite-frame target;
- any additional superpopulation uncertainty only if the parent frame itself is regenerated.

Do not call the spread of posterior means across validation replicates a Bayesian posterior interval. It is an outer repeated-sampling distribution.

---

## B. Extreme/irregular propensities and weighted-KDE uncertainty

### B1. Evidence table

| Source | Supported result | Project implication | Limitation |
|---|---|---|---|
| Liu et al. (2025), *Statistics in Medicine* 44, e70329; DOI 10.1002/sim.70329 | Reviews positivity violations and weighting estimands. Truncation can reduce variance while introducing bias for the original target; trimming changes the target population. Structural lack of positivity can make the original estimand unidentifiable. | Treat weight caps and redshift exclusions as explicit sensitivity estimands, not harmless numerical fixes. Separate random finite-sample extremes from structural regions with no support. | Causal treatment propensity is a two-arm probability; our inclusion probability is a one-arm catalogue-retention probability. |
| Kostouraki et al. (2024), DOI 10.1002/sim.10078 | Ignoring propensity estimation can misstate uncertainty; full bootstrap refitting is one approach, while analytic estimators can account for the first step. | Every replicate must refit the logistic inclusion model and recompute weights, not hold fitted probabilities fixed. | Their analytic variance results are not directly validated for this KDE-plus-selection-normalization pipeline. |
| Austin (2022), *Statistics in Medicine* 41(22), 4426–4443; DOI 10.1002/sim.9519 | Simulation comparison found bootstrap standard errors more accurate than asymptotic estimators for ATE IPW in smaller samples (in the studied settings); each bootstrap sample re-estimated the propensity model and weights. | Supports comparing design-aware replication against linearization, not assuming one is always superior. | The study's iid/causal simulation settings do not establish that an ordinary Efron bootstrap is valid for SRSWOR validation. |
| Conti & Marella (2025/2026 issue), DOI 10.1111/insr.12614 | Recommends finite-population bootstrap methods that respect the sampling designs of both samples and their dependence. | Strong reason to implement a paired, overlap-aware validation resampling design. | Their independent-sample assumptions differ from the shared-parent design here. |

### B2. What the remedies actually mean here

1. **Uncapped IPW baseline:** retain \(w_i=1/\hat\pi_i\) with no floor/cap. This preserves the intended inverse-probability construction but may have severe finite-sample variance and instability. Keep it as the reference arm.
2. **Weight truncation/capping:** cap \(w_i\) at a predeclared threshold. This can reduce variance but generally introduces bias relative to the full parent-frame target. Report threshold sensitivity and the induced density/H0 bias.
3. **Redshift trimming:** exclude catalogue support regions under a predeclared rule. This changes the target redshift population. It must be labeled as a restricted-support estimand, not a corrected estimate of the original full parent population.
4. **Smooth capping:** replace a hard cap with a smooth weight transformation. This can reduce discontinuity in sensitivity analyses but still changes the effective weighting target; it does not restore information in regions with no validation support.
5. **Regularized or bias-reduced logistic fit:** a possible numerical/model sensitivity arm if separation or near-separation occurs. It changes the fitted inclusion model and must be compared on held-out inclusion prediction, weight stability, density error, and downstream calibration. Do not label the current damped Newton solver as regularization: damping is an optimization device, not an added penalty.
6. **Overlap weighting:** do not directly import the causal formula \(e(X)[1-e(X)]\) as a replacement for \(1/\pi(z)\). In causal studies it targets a population with treatment-group equipoise. Here it would define a different redshift tilting of the parent frame and would not, without a separate derivation, reconstruct the full parent density.

No weighting remedy can identify a region where the inclusion probability is structurally zero from the observed catalogue alone. Validation must cover that region or an explicit model/extrapolation assumption must be introduced and tested.

### B3. Weighted KDE uncertainty is a full-pipeline quantity

For a weighted KDE,

\[
\hat f_h(z)=
\frac{\sum_{i\in C}\hat w_i K_h(z-z_i)}
{\sum_{i\in C}\hat w_i},
\qquad \hat w_i=1/\hat\pi(z_i),
\]

uncertainty enters through more than the weights. In this repository the bandwidth depends on weight-derived quantities (Kish ESS and weighted spread). Therefore each validation replicate must recompute:

- logistic coefficients and predicted inclusion probabilities;
- inverse weights and normalized weights;
- Kish ESS and weighted variance/spread;
- bandwidth;
- reflected-kernel contributions and boundary normalization;
- merger-rate weighting and final density normalization.

A delta-method calculation that treats bandwidth as fixed omits a pathway of uncertainty. A first-order derivative could be derived for the complete estimator, but the adaptive bandwidth rule, finite validation sample, boundary reflection, and occasional solver failures make repeated design-based simulation the more transparent first calibration tool.

**Bootstrap warning:** an ordinary iid row bootstrap of validation records does not reproduce SRSWOR. If used as an approximation, it must be justified against the actual finite-population design. The primary reference should be repeated SRSWOR draws from the fixed parent frame, preserving the catalogue-validation ID overlap. Keep failed fits in the operational attempt denominator; also report calibration conditional on successful fits.

### B4. Diagnostics to retain

- validation class counts and logistic separation/rank diagnostics;
- predicted inclusion-probability quantiles on catalogue support;
- inverse-weight quantiles, maximum normalized share, top 1%/5% weight mass;
- Kish ESS and entropy ESS (distinct concentration summaries);
- validation redshift support and boundary proximity;
- weighted spread and bandwidth per replicate;
- fit convergence, boundary flags, and failure category;
- density ISE against the declared finite-frame target;
- posterior bias, interval coverage, and interval width.

No single ESS threshold certifies validity. Diagnostics describe different failure modes and must not be collapsed into one pass/fail score.

---

## C. Selection-normalization semantics: fixed event count versus event rate

Let \(g_\theta(z)\) be the normalized source/host redshift density implied by an inclusion-model draw \(\theta\). Let \(p_{\rm det}(z,H_0)\) be the GW detection probability and define

\[
\beta(H_0,\theta)=
\int g_\theta(z)p_{\rm det}(z,H_0)\,dz.
\]

Keep this GW detection efficiency distinct from the galaxy-catalogue inclusion probability \(\pi(z;\gamma)\). They represent different selection stages.

### C1. Conditional-on-count likelihood

If the analysis conditions on the observed number \(n\) of detected GW events, each event is drawn from the detected population. Schematically,

\[
L_{\rm cond}(H_0,\theta)
=
\frac{\prod_{j=1}^{n} I_j(H_0,\theta)}
{\beta(H_0,\theta)^n},
\qquad
I_j(H_0,\theta)=
\int g_\theta(z)p(d_j\mid z,H_0)\,dz.
\]

There is no separate Poisson count likelihood in this conditional formulation. The denominator is required to normalize the event distribution conditional on detection.

### C2. Rate-aware point-process likelihood

If the number of detections is modeled, introduce an event-rate parameter \(R\) and observing time \(T\). A schematic inhomogeneous-Poisson likelihood is

\[
L_{\rm rate}(H_0,\theta,R)
\propto
\exp[-RT\beta(H_0,\theta)]
(RT)^n\prod_{j=1}^{n} I_j(H_0,\theta),
\]

up to factors independent of the parameters and the exact convention for whether \(I_j\) includes population normalization. Here the count informs the expected number of detections \(RT\beta\). If the rate is marginalized, its prior and integration must be specified.

The conditional-count and rate-aware forms are related by conditioning/factorization, but they are not to be multiplied together. Adding both a \(\beta^{-n}\) correction and a Poisson count term without deriving the joint observation model can double-count selection.

### C3. Propagating inclusion uncertainty consistently

For the modular uncertainty distribution \(p_{\rm mod}(\theta\mid V)\), conditional-count marginalization is

\[
\bar L_{\rm cond}(H_0)=
\int
\frac{\prod_j I_j(H_0,\theta)}
{\beta(H_0,\theta)^n}
p_{\rm mod}(\theta\mid V)\,d\theta.
\]

Average the **complete per-draw likelihood ratio**. Do not average the event numerators and selection denominators separately and then form their ratio. The same host-density draw \(g_\theta\) must be used in every event integral and in \(\beta\). This preserves their shared uncertainty and covariance.

For a rate-aware model, integrate the complete Poisson likelihood over the joint uncertainty in \((\theta,R)\), not a conditional likelihood with an extra count correction appended afterward.

### C4. Project-specific boundary

The current benchmark is fixed-count: it conditions on a chosen number of detected events and uses a selection-normalized event likelihood. It is not an event-rate analysis. A rate-aware extension would require a declared merger-rate model, observing-time convention, detection process, and rate prior or fixed-rate assumption. Do not describe the current code as modeling event counts.

---

## D. Recommended research sequence

**D1 — Finish the overlap-aware finite-frame benchmark first.**
- Confirm v2 artifact schema and independent audit on a small smoke artifact.
- Add repeated validation SRSWOR draws with paired seeds and overlap IDs.
- Compare plug-in modular estimator, repeated-validation distribution, and oracle finite-frame diagnostic.
- Keep parent frame, retention vector, catalogue, and GW data fixed across validation replicates.

**D2 — Add weight-stability sensitivity arms.**
- Uncapped IPW reference.
- Predeclared weight-cap levels.
- Predeclared redshift-support trimming levels.
- Optional regularized/bias-reduced inclusion fit only as a separate model arm.
- Report changed target and bias explicitly; do not rank methods by ESS alone.

**D3 — Audit adaptive KDE propagation.**
- Verify that every replicate recomputes ESS, weighted spread, bandwidth, reflection, and normalization.
- Compare adaptive-bandwidth propagation with a fixed-bandwidth diagnostic to isolate the bandwidth pathway.
- Keep density ISE and downstream posterior calibration as separate outcomes.

**D4 — Keep count semantics fixed for this study.**
- Retain the conditional-count likelihood as the primary declared target.
- Add a rate-aware likelihood only as a separately specified extension with explicit rate assumptions.
- Ensure one coherent selection factorization per analysis.

**D5 — Revisit Bayesian validation modeling only after the observation map is explicit.**
- State simulation-oracle versus observer-facing estimand.
- Define the observed-data likelihood for catalogue IDs, validation IDs, overlap, and redshifts not observed.
- Do not insert an ad hoc finite-population correction or treat validation units as independent Bernoulli draws.

---

## References

1. Conti, P. L., & Marella, D. (2025, first published online April 14). “An Uncertainty Based Approach for Dealing With Selection Bias in Non-Probability Samples.” *International Statistical Review*, 94(2), 264–290. DOI: 10.1111/insr.12614.
2. Mashreghi, Z., Haziza, D., & Léger, C. (2016). “A survey of bootstrap methods in finite population sampling.” *Statistical Surveys*, 10, 1–52.
3. Kostouraki, A., et al. (2024). “On variance estimation of the inverse probability-of-treatment weighting estimator: A tutorial for different types of propensity score weights.” *Statistics in Medicine*, 43(13), 2672–2694. DOI: 10.1002/sim.10078.
4. Austin, P. C. (2022). “Bootstrap vs asymptotic variance estimation when using propensity score weighting with continuous and binary outcomes.” *Statistics in Medicine*, 41(22), 4426–4443. DOI: 10.1002/sim.9519.
5. Liu, Y., et al. (2025). “A Tutorial for Propensity Score Weighting Methods Under Violations of the Positivity Assumption.” *Statistics in Medicine*, 44, e70329. DOI: 10.1002/sim.70329.
6. Ridgeway, G., et al. (2015). “Propensity Score Analysis with Survey Weighted Data.” *Sociological Methods & Research*. DOI: 10.1177/0049124113500475.
7. Tagliazucchi, M., Gair, J., Barbieri, R., & Moresco, M. (2026). “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs.” *Physical Review D*, 114, 083523. DOI: 10.1103/vmzf-1wxq.

**Evidence boundary:** The sources above support design-aware propagation, propensity-estimation uncertainty, positivity/weighting trade-offs, and the event-count versus rate distinction. They do not supply a ready-made variance formula for this exact finite-frame logistic-IPW-adaptive-KDE dark-siren pipeline. That exact combination remains an empirical/methodological question for the project's own controlled benchmark; no novelty claim is made here.
