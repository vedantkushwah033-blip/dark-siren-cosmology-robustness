# Confirmatory design lock: unresolved decisions and release gates

**Status:** decision register; no confirmatory design is locked by this document.  
**Scope:** the linear finite-population arithmetic-mean benchmark only.  
**Version date:** 2026-10-31.

This register prevents a planning utility, pilot result, or software test from being mistaken for a preregistered scientific design. It is deliberately conservative: unresolved items below are blockers to confirmatory execution, not invitations to choose settings after seeing outcomes.

## A. Decisions that must be made before confirmatory outcomes are generated

| Decision | Required prespecification | Current status |
|---|---|---|
| Primary scientific contrast | Name the reference method, comparator, estimand, sign convention, and the smallest practically relevant difference (or justify a purely precision-based objective). | OPEN |
| Primary design cell or multiplicity family | Identify the population regime(s), sampling fraction(s), and confidence level that define the primary claim. If several cells are co-primary, define family-wise error control (e.g. Holm) in advance. | OPEN |
| Frame-generating regimes | Fix distribution families, parameters, truncation/support rules, population size, and independent frame seed schedule. | PARTLY SPECIFIED; final grid not locked |
| Independent frame count | Use a justified between-frame variability estimate and a target interval half-width; report sensitivity to SD inflation and pilot degrees of freedom. | OPEN |
| Within-frame repetitions | Set from a prespecified Monte Carlo precision target for coverage and operational success; do not treat more within-frame draws as a substitute for independent frames. | OPEN |
| Inner bootstrap budget | Fix attempts per outer sample and rules for failed/degenerate bootstrap replicates. | OPEN |
| Method implementation | Freeze exact source commit, method definitions, numerical tolerances, and dependency environment. | OPEN |
| Stopping and exclusions | Fix stopping rule and deterministic handling of invalid draws, numerical failures, and incomplete artifacts. No outcome-driven exclusions or optional stopping. | OPEN |
| Inference and multiplicity | Fix frame-level aggregation, confidence-interval construction, and correction for the declared family of primary tests. | OPEN |
| Independent reproduction | Name a reviewer or independent execution procedure and define the artifact/audit checks required for release. | OPEN |

## B. Estimand and analysis hierarchy

For each fixed design cell, frame f, method m, and outer sampling attempt r, record whether an interval was constructed, whether it covered the exact realized finite-frame mean, and its width when defined.

Report separately:

1. **Construction probability:** successful interval constructions divided by all outer attempts.
2. **Conditional coverage:** covering intervals divided by successful constructions; undefined if there are no successful constructions.
3. **Successful-and-covering probability:** covering intervals divided by all outer attempts.
4. **Conditional width:** mean width among successful constructions, with its denominator stated.
5. **Failure distribution:** counts by prespecified canonical failure category.

For method contrasts, first calculate paired contrasts within each frame using common outer samples. The independent inferential units for generalization across generated frames are the frames. Summarize those frame-level contrasts with an explicitly justified interval; do not calculate a naive standard error by pooling all outer attempts as independent observations.

Do not pool heterogeneous population regimes into a single average unless the target mixture and its weights are defined before the experiment. Present regime-specific results as the default. If a pooled estimand is scientifically justified, define the mixture weights in the manifest before generation.

## C. Precision planning rules

The existing `frame_precision_planning.py` utility is a planning aid, not a design oracle. Its t-based frame-count calculation assumes independent, approximately normal frame-level contrasts and uses a chi-square upper confidence bound for the frame-level standard deviation. The sensitivity grid should be reported, not used to conceal dependence on the planning assumptions.

Before locking a frame count:

- Document the source and number of independent pilot frames used to estimate the frame-level SD.
- Do not represent the current three-frame pilot as a reliable estimate of a broad frame-generating distribution.
- Show frame-count sensitivity to at least the existing SD multipliers and variance-confidence levels.
- Evaluate whether the frame-level contrast distribution is plausibly compatible with the t-based planning approximation; if not, revise the planning method before locking.
- Independently plan within-frame repetitions from Monte Carlo precision for the coverage and operational-success quantities.
- State what happens if the target precision is not reached. The default is an inconclusive result, not a post hoc increase in sample size based on observed effects.

## D. Confirmatory execution gates

All gates are required, in order:

1. **Protocol gate:** all OPEN items in section A resolved, with protocol and manifest committed.
2. **Software gate:** full test suite passes on every supported Python version at the exact frozen commit.
3. **Pilot gate:** v3 smoke artifact is generated and the independent auditor passes, including checksum-refreshed tamper tests.
4. **Reproduction gate:** an independent execution reconstructs the manifest, frames, samples, method outputs, failures, and summaries from the frozen seeds.
5. **Execution gate:** confirmatory run completes without unplanned changes to design, methods, or stopping rules.
6. **Analysis gate:** all denominators, failures, Monte Carlo precision, frame-level uncertainty, and multiplicity adjustments are reported.
7. **Manuscript gate:** claims are limited to the prespecified estimands and synthetic finite-population regimes; code, manifests, raw artifacts, audit output, and environment are archived.

A failure at any gate blocks a confirmatory claim. Fixes that alter the analysis or experiment require a new versioned protocol and, where relevant, a fresh run.

## E. Claim boundary

A successful linear finite-population benchmark supports claims only about the tested interval procedures for the arithmetic mean under the specified SRSWOR designs. It does not establish calibration of an estimated galaxy-catalogue inclusion model, a weighted KDE, a nonlinear gravitational-wave likelihood, or a dark-siren posterior for (H_0). Those require a separate simulation design, estimands, precision plan, validation data, and independent audit.

## F. Publication and admissions integrity

The work should be presented as a reproducible student-led methodological simulation study, with the level of evidence stated plainly. Do not describe the pilot as confirmatory, the implementation as a new theorem, or synthetic data as observational gravitational-wave data. External faculty review and independent reproduction are meaningful next milestones; institutional admissions outcomes cannot be inferred from this project.
