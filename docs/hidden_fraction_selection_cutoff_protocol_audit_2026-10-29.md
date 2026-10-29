# Protocol audit and implementation contract: hidden-host fraction × event-selection cutoff

**Date:** 2026-10-29  
**Status:** design refinement plus a tested likelihood primitive; no experiment results are claimed  
**Source prompt:** user-supplied “Comprehensive Independent Scientific Audit & Experimental Protocol” (Mission 1 and Mission 2).  
**Base commit:** `78f3b5419569a21ed4ef0395d69e2c1a56157856`.

## Executive decision

The proposed question is worth investigating, but the submitted protocol should **not be implemented literally without a likelihood contract first**. The repository already contains separate hidden-host-fraction sensitivity studies, host-mixture/threshold studies, and a selection-aware likelihood. The defensible incremental question is whether **hidden-host mixture misspecification interacts with event-selection conditioning**, with matched simulations and diagnostics that are evaluated for both false alarms and detection power.

This document defines the minimum statistical contract before coding. It does not report numerical findings, and it does not claim that the proposed experiment has run.

## 1. Corrections and safeguards to the submitted protocol

### 1.1 Define what the hidden-host fraction means

The parameter (f) must be defined as one of the following, not left ambiguous:

- **Pre-selection host-population mixture fraction:** the fraction of the source population assigned to the missing-host component before event detection; or
- **Post-selection detected-event fraction:** the missing-host fraction among the events that survived the selection rule.

These are generally not equal if the two components have different detection efficiencies. The proposed likelihood below uses the **pre-selection mixture fraction**. The data generator must apply detection to the full mixture so the detected mixture proportion changes naturally.

If the simulator instead fixes the number of missing-host events after selection, that is a different experiment and needs a different likelihood derivation.

### 1.2 Keep catalogue inclusion and GW-event detection distinct

Catalogue inclusion (C_g) and GW-event detection (A) are separate stochastic processes. The model must not use a galaxy-catalogue inclusion probability as if it were a GW detection efficiency. If the experiment treats (f) as an abstract mixture parameter, document that it is a simplified host-population mixture, not a complete survey-generation model.

### 1.3 Make the event-selection region explicit

The submitted protocol calls ({500,1000,2000},mathrm{Mpc}) “distance cutoffs” but does not specify whether each value is:
- a maximum accepted observed distance (d_{max});
- a minimum detection threshold (d_{min}); or
- a proxy for a detector horizon that is not itself a hard observed-distance cut.

For the initial experiment, use a simple, explicitly synthetic observed-distance acceptance region
(R=[d_{min},d_{max}]), with (d_{min}=0) and (d_{max}) taking the prespecified values. Do not call this a realistic detector selection function. Generate the full mixture, then retain events satisfying the region. A later extension can use a more realistic detection model only as a separately versioned experiment.

For a Gaussian distance model (dmid z,Hsimmathcal N(mu(z,H),sigma^2)), the acceptance probability for a fixed redshift is

[
a(z,H;R)=Phi!left(rac{d_{max}-mu(z,H)}{sigma}ight)
-Phi!left(rac{d_{min}-mu(z,H)}{sigma}ight).
]

The implementation must derive the likelihood and its normalization from the exact simulator, not merely insert a generic (eta(H)) term.

### 1.4 Use a mixture-consistent selection normalization

Let (g_C(z)) and (g_M(z)) be normalized redshift densities for the catalogued and missing-host components, respectively. For a pre-selection missing-host fraction (f), define

[
q_C(dmid H)=int p(dmid z,H)g_C(z),dz,
qquad
q_M(dmid H)=int p(dmid z,H)g_M(z),dz,
]

[
alpha_C(H;R)=int a(z,H;R)g_C(z),dz,
qquad
alpha_M(H;R)=int a(z,H;R)g_M(z),dz.
]

For an accepted observation (din R), the fixed-count, selection-conditioned mixture likelihood is

[
p(dmid H,f,A=1,R)=
rac{(1-f)q_C(dmid H)+f q_M(dmid H)}
{(1-f)alpha_C(H;R)+falpha_M(H;R)}.
]

For conditionally independent detected events, sum the log of this per-event likelihood. The assumed fraction (f_{m model}) must appear in both the mixture numerator and the mixture selection denominator. The true fraction (f_{m true}) is used only by the generator and independent scoring/audit code. A non-oracle inference arm must not read (f_{m true}), the missing component's realized host labels, or any parent-population values unavailable under that arm.

This expression is a contract for the stipulated toy model; it is not a full LVK hierarchical likelihood and does not include event-rate information.

### 1.5 Do not assume the bias direction

The submitted H1 proposes that underestimating (f) biases (H_0) downward. That direction is **not guaranteed**: it depends on the redshift distribution of the missing component, the distance/noise model, and selection conditioning. Preserve the directional prediction as a falsifiable hypothesis, not as an expected fact built into code or reporting.

### 1.6 PIT/SBC needs a correct target distribution

A KS test of “PIT values are uniform” is not automatically valid for a fixed-truth experiment. The posterior CDF evaluated at a fixed true (H_0) over repeated data sets is not generally uniform in a frequentist fixed-parameter design.

Use distinct diagnostics for distinct questions:

- **Fixed-truth calibration:** empirical 68%/95% interval coverage, posterior bias/RMSE, interval width, and stratification by (H_0) and cutoff. Report Monte Carlo intervals.
- **Simulation-based calibration (SBC):** draw (H_0) from the declared prior for each replicated data set, simulate the data under that draw, compute posterior ranks using the matching generative/selection model, and assess rank uniformity with a prespecified diagnostic. SBC assesses computational calibration under that prior-predictive simulator; it does not prove real-world model adequacy.
- **Posterior predictive checks (PPCs):** simulate replicated *selected* data from posterior draws under the model and compare prespecified discrepancy statistics. Candidate statistics include selected-distance quantiles, tail fractions, and redshift/component summaries available to the simulator. State which statistics would be observable in a practical analysis versus known only in simulation.
- **Detection-rate checks:** only meaningful if the experiment models the number of attempted sources and the number detected. With a fixed number of accepted events by construction, the detected count itself cannot diagnose misspecification.

Do not use an uncorrected (p<0.01) rule as a universal detector. Estimate false-positive rate under correct specification and power under each misspecification, with Monte Carlo uncertainty and multiplicity handling stated in advance.

## 2. Recommended staged design

A full cross of every factor may be expensive and can encourage post-hoc selection. Run the stages below in order and preserve all outcomes.

### Stage A — deterministic likelihood contract tests

Before simulation, add tests for:

1. (f=0): likelihood and selection normalization reduce to the catalogue component.
2. (f=1): they reduce to the missing-host component.
3. (g_C=g_M): the likelihood is invariant to (f).
4. Correctly specified mixture: numerator and selection denominator use the same component weights.
5. (d_{min}=0), finite (d_{max}): analytic acceptance probabilities match numerical quadrature.
6. Probabilities are finite and in ([0,1]); (alpha(H)>0) throughout the inference grid.
7. The likelihood equals a separately coded direct calculation on a tiny known-answer case.
8. Non-oracle inference has no code path that reads true (f) or latent missing-host identities.
9. Invalid (f
otin[0,1]), invalid grids, and zero acceptance probabilities fail closed.
10. Numerical stability is tested for very small acceptance probabilities and extreme but valid log-likelihoods.

### Stage B — small pilot

Use (f_{m true}=0.3), (f_{m model}in{0.1,0.3,0.5}), (d_{max}in{500,1000,2000},mathrm{Mpc}), and 50 paired replicates per cell. Use a deterministic seed derivation from a stable hash of the full configuration so cell order cannot change results.

This is a 3 × 3 = 9-cell pilot, 450 data sets, before adding a second dimension for explicit cutoff misspecification. The pilot should estimate runtime and Monte Carlo precision; it is not enough for definitive coverage claims.

**Important distinction:** varying the true cutoff while using the same correct cutoff in inference tests sensitivity to the selection regime. It does **not** test cutoff misspecification. To test cutoff mismatch, add a separately labelled factor for (d_{max,m model}) (correct, tighter, looser) only after Stage B passes. Do not conflate the two questions.

### Stage C — confirmatory fraction-mismatch experiment

After inspecting pilot runtime but before looking at confirmatory outcomes, freeze:
- true fractions (f_{m true}in{0.1,0.3,0.5});
- model fractions (f_{m model}in{0.1,0.3,0.5});
- true selection cutoffs (d_{max,m true}in{500,1000,2000},mathrm{Mpc});
- the rule for whether (d_{max,m model}) is correct in this phase or intentionally varied in a distinct phase;
- truth values of (H_0), noise levels, event count, prior, redshift grid, host-component densities, and Monte Carlo precision targets;
- replicate count and independent seed ranges.

A fully crossed 3 × 3 × 3 design is 27 cells; at 500 replicates it requires 13,500 data sets. Adding a three-level model-cutoff factor increases this to 81 cells and 40,500 data sets. Choose the design based on a pilot and a stated precision target, not on which outcomes look most interesting. If compute is limited, preserve a smaller balanced design and explicitly narrow the question.

### Primary and secondary outcomes

Primary:
- posterior median error in (H_0);
- mean signed bias and RMSE across replicates;
- paired difference in absolute error between correctly specified and misspecified (f) models, with a paired bootstrap or a justified alternative and a deterministic analysis seed.

Secondary:
- 68% and 95% interval coverage with numerator, denominator, and Wilson intervals;
- interval width;
- fit failure/non-finite rate;
- selection efficiency (alpha(H)) and minimum numerical stability margin;
- PPC false-positive rate under the correct model and detection power under each misspecification.

Report each cutoff and truth stratum separately before any pooled summary. Do not pool cells with materially different selection regimes as if they were exchangeable.

## 3. Independent artifact audit contract

Every replicate should preserve at minimum:
- replicate and configuration IDs, seed, true (H_0), (f_{m true}), true cutoff, model (f), model cutoff, and event count;
- generated data sufficient to independently reconstruct the likelihood (or a clearly justified immutable seed/configuration replay);
- posterior grid and normalized posterior density, log-likelihood components, (alpha_C(H)), (alpha_M(H)), and mixture (eta(H));
- posterior mean/median, 68%/95% endpoints, coverage flags, and all PPC discrepancy statistics;
- failure records rather than silently dropping failed fits;
- environment/version metadata and SHA-256 hashes.

The auditor must independently reconstruct:
1. the selected-data acceptance rule and its acceptance probabilities;
2. the mixture numerator and denominator;
3. posterior normalization and quantiles;
4. bias, RMSE, coverage, interval width, and paired contrasts;
5. PPC/SBC diagnostic summaries, with their design-specific interpretation;
6. raw-to-summary provenance and checksums.

Add fault-injection tests for a wrong mixture weight, a missing selection-normalization factor, a wrong cutoff in the denominator, a corrupted posterior grid, a missing arm, duplicate replicate/configuration keys, and a modified summary with a recomputed manifest hash.

## 4. Stop/go gates

- **Gate 0 — protocol:** equations, parameter meanings, data-generation mechanism, and selection region reviewed.
- **Gate 1 — deterministic tests:** all limiting cases and independent known-answer tests pass.
- **Gate 2 — pilot:** all cells complete, raw artifacts exist, independent auditor passes, runtime and Monte Carlo precision are reported.
- **Gate 3 — confirmatory freeze:** factors, seeds, primary estimand, and decision rules are locked before outcomes are inspected.
- **Gate 4 — scientific interpretation:** report null and adverse outcomes; no claim that weighting always helps; separate model-based uncertainty from Monte Carlo uncertainty.
- **Gate 5 — release:** tests pass on the exact commit; artifacts and scripts are archived; limitations and AI-assistance disclosures follow the chosen venue's policy.

No numerical result should be written into a manuscript until the corresponding run and independent artifact audit pass.

## 5. Relationship to existing repository work

The repository already has:
- a hidden-host fraction/redshift assumption sensitivity study;
- host-mixture and detection-threshold sensitivity experiments;
- a selection-aware toy likelihood with a fixed-count detection-conditioning term;
- a catalogue-selection weighting pilot and independent raw-to-summary auditor;
- a derivation note explicitly distinguishing catalogue missingness from event detection.

The new experiment should reuse these validated components where appropriate, but **must not splice together separate likelihoods without a fresh derivation**. The contribution, if any, is the controlled interaction between the hidden-host mixture and selection conditioning, plus measured diagnostic sensitivity—not the generic discovery that host weighting or selection matters.

## 6. Literature and terminology notes

- Gray, Messenger & Veitch (2022), “A Pixelated Approach to Galaxy Catalogue Incompleteness,” *MNRAS* 512, 1127–1140, DOI [10.1093/mnras/stac366](https://doi.org/10.1093/mnras/stac366), establishes a pixelated approach to directional catalogue incompleteness.
- Gair et al. (2023), “The Hitchhiker’s Guide to the Galaxy Catalog Approach for Dark Siren Gravitational-wave Cosmology,” *The Astronomical Journal* 166(1), 22, DOI [10.3847/1538-3881/acca78](https://doi.org/10.3847/1538-3881/acca78), is a methodological guide to the galaxy-catalogue approach.
- Datrier & Hendry (2026), “Implementing a robust test of galaxy catalogue completeness for dark siren measurements of the Hubble constant,” *MNRAS* 549(1), stag828, DOI [10.1093/mnras/stag828](https://doi.org/10.1093/mnras/stag828), develops a magnitude-completeness test in `gwcosmo`. It is relevant context but not the same experiment as a synthetic hidden-host-fraction × event-selection-cutoff interaction.
- Tagliazucchi et al. (2026), “Enhancing dark siren cosmology via Gaussian process reconstruction of incomplete galaxy catalogs,” *Physical Review D*, DOI [10.1103/vmzf-1wxq](https://doi.org/10.1103/vmzf-1wxq), is a recent related reconstruction approach.

These sources establish close related work; they do not establish whether the proposed interaction experiment is novel. A more systematic search and direct method comparison are required before any priority claim.

## 7. Claim boundary

If implemented and independently audited, the strongest safe claim would be limited to the explicitly simulated generative process, event-selection rule, host-component distributions, prior, and tested parameter settings. Do not claim:
- a real-data (H_0) measurement;
- detector-level or LVK validation;
- universal correction of galaxy-catalogue incompleteness;
- calibrated observational robustness merely because SBC or CI passes;
- that PIT is uniform at a fixed truth by definition;
- that a nonsignificant result establishes equivalence;
- a novel method without direct comparison to the closest published methods.

**Current status:** protocol refinement and deterministic likelihood primitive implemented. Module: `src/dark_siren_robustness/selection_conditioned_hidden_host.py`; tests: `tests/test_selection_conditioned_hidden_host.py`. The dedicated Python 3.10/3.11/3.12 contract workflow passed at [run #40859955488](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40859955488). This validates the primitive's coded limiting cases and numerical contract only. The full fraction × cutoff experiment, independent artifact auditor, PPC/SBC study, and all scientific results remain unimplemented/unrun.
