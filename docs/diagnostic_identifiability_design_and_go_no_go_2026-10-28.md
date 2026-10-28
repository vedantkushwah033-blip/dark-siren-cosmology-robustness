# Diagnostic-identifiability study: operational design and implementation gates

**Status:** prespecified design proposal only; no diagnostic-identifiability experiment has been run.  
**Date:** 2026-10-28  
**Target branch:** `audit/research-readiness-2026-10-28`  
**Scientific scope:** controlled synthetic low-redshift model; no real GW inference claim.

## 1. Decision this study is meant to answer

The current repository contains separately audited sensitivity experiments for hidden-host assumptions and distance-uncertainty likelihoods. Those results do **not** by themselves establish that an analyst can tell which misspecification generated a discrepancy. This proposed study asks:

> Under a declared synthetic generator, how accurately can a prespecified diagnostic procedure distinguish host-population misspecification from distance-noise misspecification, and how often must it return “not distinguishable”?

The target is diagnostic identifiability under this simulator—not a universal misspecification detector, a new H0 estimator, or a validated dark-siren method. Novelty remains unestablished until the full methods/results of the closest literature are compared.

## 2. Critical feasibility gate: event-level versus dataset-level inference

The inspected `posterior_h0` API accepts a single scalar `gw_distance` and `gw_sigma` per call, with a vector of candidate host redshifts inside that event's mixture. That is an event-level posterior routine, not, by itself, a multi-event hierarchical likelihood or a dataset-level posterior-predictive diagnostic framework.

**Do not launch a diagnostic-classification experiment by treating individual events as independent replicates if the proposed diagnostic is defined on a whole multi-event dataset.** First implement and test an explicitly documented composition rule for independent events (or define a scientifically coherent single-event estimand and abandon dataset-level claims). The rule must specify how host mixtures, event-specific uncertainties, shared parameters, priors, and any selection conditioning combine. The existing scalar API must not be silently repurposed to imply a full population likelihood.

Minimum feasibility tests before simulation:
- A one-event dataset gives the same posterior as the current event-level routine under matching assumptions.
- For two independent events, the joint log likelihood equals the sum of their event log likelihoods, up to a parameter-independent constant.
- Reordering events leaves the posterior unchanged.
- Duplicating an event changes the posterior as expected from multiplying its likelihood twice.
- Zero hidden-host mass reduces exactly to the catalogue-only model.
- Generator and likelihood use the same positive-distance convention; if the generator redraws non-positive distances, the likelihood must use the matching truncated/conditional density.
- Prior normalization, grid stability, and extreme-likelihood behavior are tested before large runs.

If these checks cannot be passed, stop at a software/model audit. Do not manufacture a diagnostic study from incompatible abstractions.

## 3. Proposed four-condition factorial

Only after the feasibility gate passes, use a fixed, version-controlled design with four principal data-generating conditions:

| Condition | Host population | Distance-noise model | Purpose |
|---|---|---|---|
| M0: matched control | correctly specified | correctly specified | Negative-control/reference distribution |
| MH: host-only mismatch | misspecified | correctly specified | Isolate host-population discrepancy |
| MN: noise-only mismatch | correctly specified | misspecified | Isolate noise-model discrepancy |
| MJ: joint mismatch | misspecified | misspecified | Assess interaction and ambiguity |

Add two essential controls:
- **Zero-hidden-population control:** true hidden fraction is zero; adding a hidden component must not create a spurious benefit when its mass is zero.
- **Matched-generator negative control:** generate and fit under the same declared model and settings. This is a computational/calibration control, not evidence of robustness under misspecification.

The factor levels must be selected from already documented, scientifically interpretable values and frozen before the confirmatory run. At minimum, vary the hidden fraction (including exactly zero), hidden-host redshift assumption, and distance-noise scale. Include interior H0 truths as the primary coverage set and report prior-boundary truths separately. Do not choose levels after inspecting which values produce the most dramatic result.

Use paired common observations for inference-arm contrasts: generate each dataset once, then analyze that same dataset with each prespecified inference arm. If the data generator changes between arms, label the comparison as a generator-plus-inference contrast, not a likelihood-only contrast.

## 4. Define what “distinguish” means before coding

A plot of coverage, bias, or RMSE across conditions does not establish that the source of misspecification can be identified from observed data. Repeated-sampling coverage is a property evaluated across replicates; it is not automatically an individual-dataset classifier.

Before the confirmatory run, choose **one** of these paths and document why it answers the scientific question:

### Path A — diagnostic-distribution separation (preferred first pass)

Define a small vector of dataset-level discrepancy statistics with explicit reference distributions under M0. Candidate statistics may include a predeclared posterior-predictive discrepancy for standardized distance residuals and a separately defined discrepancy sensitive to host-mixture composition. The exact statistics must be derived from the actual implemented model, not selected after seeing results.

For each condition, estimate the distribution of each statistic and the joint diagnostic vector. Report:
- distributions and uncertainty intervals;
- pairwise separability using a prespecified distance or two-sample procedure;
- false-alarm rate under M0 and detection rate under each alternative at a threshold calibrated only on an independent M0 calibration set;
- sensitivity to event count and mismatch magnitude;
- whether MH and MN remain practically indistinguishable under the chosen diagnostics.

Use a held-out evaluation seed set. Calibration and evaluation datasets must not overlap. If diagnostic thresholds are tuned, use a separate training/calibration split and freeze them before evaluation.

### Path B — source-classification benchmark (only if scientifically justified)

If a classifier is used, prespecify the input statistics, model family, class balance, train/calibration/test split by independent simulation seed, tuning protocol, and primary metric. Report balanced accuracy, per-class sensitivity, confusion matrix, calibration, and uncertainty intervals against a simple baseline. Prevent leakage: replicates from one simulated dataset or seed family must not be split across training and test sets.

A classifier is not inherently a scientific contribution. Do not use one merely to make the project appear more advanced. If the result is that MH and MN cannot be reliably distinguished at realistic simulated sample sizes, report that limitation as a substantive finding only if the design and uncertainty support it.

## 5. Primary estimands and reporting

Use a small hierarchy of outcomes:

1. **Primary:** held-out pairwise diagnostic separability between MH and MN at the prespecified main sample size, with uncertainty interval and a prespecified practical-equivalence/indistinguishability criterion.
2. **Key controls:** false-alarm rate under M0; detection rate for MH, MN, and MJ relative to M0.
3. **Inference consequences:** H0 posterior bias, RMSE, 68% and 95% interval coverage, and interval width for each condition/arm.
4. **Robustness:** sensitivity to event count, prior support, grid resolution, hidden fraction, noise mismatch, and random seed.
5. **Failures:** numerical failures, invalid posteriors, rejected datasets, and any exclusion rule with counts and reasons.

Do not rank arms on coverage alone. A wider interval may cover more often while providing less precision. Report Monte Carlo uncertainty for all primary estimates. For a binomial proportion near 0.5, the approximate standard error is √[0.25/R]; therefore select the replicate count from a target precision calculation rather than habit. State the exact number of independent simulated datasets separately from the number of posterior fits.

The practical-equivalence margin and all decision thresholds must be set before looking at confirmatory results. If no defensible threshold can be justified, report estimates and uncertainty without a categorical “distinguishable/not distinguishable” claim.

## 6. Independent audit and reproducibility requirements

The runner and auditor must be separate implementations at the summary/reconstruction layer. The audit should reconstruct from raw synthetic records:
- condition labels and truth/assumed parameter values;
- the data-generation protocol and pairing identifiers;
- posterior summaries and diagnostic statistics;
- cell-level coverage, bias, RMSE, widths, and Monte Carlo intervals;
- paired contrasts and their resampling units;
- threshold calibration and held-out evaluation metrics;
- failure counts, exclusions, and source revision metadata.

Save one record per independent simulated dataset and inference arm. Store deterministic seed identifiers, source revision, Python and dependency versions, exact equations, prior, grid, event count, and all factor levels. Include tamper tests for raw observations, condition labels, pairing, summaries, thresholds, and held-out split assignments. Publish checksums and clean-checkout reproduction commands.

Use an exploratory pilot only for runtime and numerical feasibility. Do not reuse pilot data in the confirmatory evaluation if the pilot influenced model or threshold choices.

## 7. Literature and external-validity gates

Before describing this as novel:
- read methods and results, not just abstracts, for existing dark-siren completeness and catalogue-reconstruction studies;
- compare the exact proposed estimand with Bayesian calibration, posterior predictive checking, simulation-based calibration, and misspecified-model diagnostics;
- record what is mathematically different, what is already known, and what the toy model cannot answer;
- update the comparison for relevant 2026 publications.

The exact pinned CHIMERA v2.2 benchmark remains a separate open gate because its documented input fixtures have not been verified as available. Do not silently substitute a different version or claim real-data validation. If external validation remains infeasible, the manuscript must explicitly be a synthetic-model study and discuss the consequences.

## 8. Go/no-go criteria

**GO to a confirmatory simulation only if all are true:**
- [ ] The event/dataset likelihood composition is defined and its limiting-case tests pass.
- [ ] The four generator conditions differ only in the factors explicitly named; paired arms use identical observations.
- [ ] Diagnostic statistics and reference distributions are defined mathematically.
- [ ] The primary estimand, practical-equivalence margin, replicate count, seeds, and exclusions are frozen in version control.
- [ ] Training/calibration/evaluation data are independent by seed, with no leakage.
- [ ] An independent audit can reconstruct all primary results from raw records.
- [ ] The literature comparison identifies a precise unresolved gap rather than assuming novelty.

**NO-GO / narrow the claim if any fail.** A successful CI run demonstrates that code executed and tests passed; it does not prove the statistical question is identifiable, astrophysically realistic, novel, or publication-ready.

## 9. Current status

This file is a design and decision gate, not a completed experiment. No new diagnostic result, source-classification performance, external CHIMERA validation, real-data inference, or publication claim is asserted. The most important next technical task is to resolve the event-level versus dataset-level inference boundary before writing diagnostic code.
