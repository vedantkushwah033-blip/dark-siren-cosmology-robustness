# Preregistration protocol for a confirmatory multi-frame finite-population bootstrap study

**Status:** prospective protocol; not an executed experiment, not a result, and not a theorem.
**Applies to:** the linear finite-population arithmetic-mean benchmark only.
**Does not validate:** nonlinear dark-siren posterior means, credible intervals, catalogue-selection corrections, or cosmological measurements.

## 1. Purpose and confirmatory question

The confirmatory study will assess interval procedures for the arithmetic mean of a realized finite population under simple random sampling without replacement (SRSWOR), while separating variation conditional on a finite frame from variation across independently generated frames.

The question is not whether one procedure is universally best. It is whether prespecified procedures exhibit practically important differences in (i) conditional coverage, (ii) probability of successfully returning a covering interval, and (iii) interval width across a prespecified collection of finite-population regimes.

## 2. Experimental units and hierarchy

There are two stochastic levels:

1. **Frame level:** independently generated finite populations within each prespecified distribution regime.
2. **Sampling level:** repeated SRSWOR samples and interval constructions within each realized frame and sampling fraction.

Frames—not individual sample repetitions—are the independent units for generalization across generated finite populations. Repetitions within a frame estimate conditional operating characteristics and must not be treated as independent frame replicates.

Each frame receives a stable identifier and deterministic seed. The frame-generation seed sequence, distribution parameters, population size, sample sizes, and all analysis settings must be committed before confirmatory outcomes are inspected.

## 3. Estimands

For frame f, method m, and design cell c, let R be the number of outer sampling attempts, S_fmc the number of successfully constructed intervals, and C_fmc the number of successful intervals covering the exact finite-frame mean.

Report all of the following:

- **Interval-construction success:** S_fmc / R.
- **Conditional coverage:** C_fmc / S_fmc, undefined when S_fmc = 0.
- **Unconditional successful-and-covering probability:** C_fmc / R, counting an unsuccessful interval construction as not delivering a covering interval.
- **Mean interval width among successful constructions:** the arithmetic mean of upper minus lower endpoint.
- **Failure frequency and canonical failure categories:** separately from coverage.

These are distinct operating characteristics. The unconditional quantity is an operational delivery metric; it must not be described as a universal replacement for conditional coverage.

For method comparisons, compute paired frame-level contrasts using the same frame, sampling design and outer replicate seed schedule wherever the methods permit. Define the sign convention in every table so that positive values consistently favour the named reference method (or explicitly state otherwise).

## 4. Scenarios and method set

The final scenario grid must be frozen in a versioned design manifest before confirmatory runs. It should include distributional regimes motivated by the pilot (symmetric light-tailed, skewed positive-valued, and heavy/right-tailed finite populations), but parameter values and frame-generation rules must be explicit.

The method set is limited to procedures implemented and independently tested before the design lock. The primary comparison should be the studentized pseudo-population bootstrap against the prespecified SRSWOR Wald interval. The ordinary percentile bootstrap may be retained as a secondary descriptive comparator. No method may be removed after observing its confirmatory performance; any numerical failures remain in the denominator for operational success metrics.

A separate pilot or planning run may be used to estimate computational cost and between-frame variability. Planning data must not be silently pooled into confirmatory estimates.

## 5. Precision-driven sample-size selection

Do not select the number of frames or repetitions merely because a round number is conventional. Use a documented planning stage:

1. Estimate between-frame variability from pilot frames that are excluded from confirmatory inference.
2. Choose a target Monte Carlo precision for the primary frame-averaged coverage contrast and for the unconditional success-and-covering probability.
3. Simulate the proposed hierarchical design using planning-only estimates, including plausible low-success scenarios.
4. Increase the number of independent frames until the planned confidence-interval half-width for the primary frame-level contrast meets the prespecified tolerance under conservative planning assumptions.
5. Select within-frame repetitions to make conditional operating-characteristic estimates sufficiently precise without pretending they compensate for too few independent frames.
6. Freeze frame count, repetitions, inner bootstrap replicates, and all seeds in the manifest before confirmatory outcomes are examined.

Report the achieved Monte Carlo standard errors and confidence-interval widths. If the planned precision is not achieved, label the study inconclusive for that estimand rather than changing the target after seeing results.

The repository also provides `src/dark_siren_robustness/frame_precision_planning.py`, a small planning utility that takes a pilot estimate of the standard deviation of independent frame-level contrasts and a target confidence-interval half-width. It inflates the planning SD using a one-sided chi-square upper confidence bound based on the supplied pilot degrees of freedom, then finds the smallest frame count whose two-sided Student-t half-width meets the target. Its output is explicitly labelled `PLANNING_ONLY_NOT_CONFIRMATORY`. The utility separates `pilot_frame_count` (which determines the chi-square degrees of freedom for the upper SD bound) from `minimum_frames` (the lower bound for the planned confirmatory frame count). If `pilot_frame_count` is omitted it defaults to `minimum_frames` for backward compatibility. The frame-count search uses a monotone binary search rather than scanning every candidate count. The utility does not estimate the frame-level SD from the existing three-frame pilot, does not select within-frame repetitions, and relies on independent approximately normal frame-level contrasts. Those assumptions must be assessed and the planning inputs justified before a confirmatory design is locked.

## 6. Uncertainty and inference

The primary generalization target is the mean operating characteristic across the prespecified frame-generating process, not merely across repeated samples from one convenient frame.

Use frame-level summaries as the units for between-frame uncertainty. Report the number of frames and the distribution of frame-level results. A frame-cluster bootstrap or another prespecified cluster-aware procedure may be used for uncertainty in the frame-averaged estimand, with the exact resampling scheme and interval method locked in advance.

Within-frame binomial Monte Carlo uncertainty and between-frame uncertainty must be distinguished. A pooled binomial interval treating every outer repetition as an independent population draw is not an acceptable substitute for frame-aware uncertainty.

For multiple scenarios and secondary outcomes, define a hierarchy:
- one primary method contrast and one primary estimand;
- secondary estimands and scenario interactions;
- exploratory subgroup or diagnostic analyses.

Use a prespecified multiplicity procedure for confirmatory secondary claims (for example, Holm adjustment across a declared family). Label unadjusted exploratory intervals as descriptive.

## 7. Failure accounting and artifact audit

The v2 artifact schema must retain every outer attempt, including unsuccessful interval constructions. For studentized pseudo-population resampling it must also retain the number of inner attempts, successes, failures, and canonical failure-code counts, or an explicit NA state when the procedure fails before inner draws.

The independent auditor must reconstruct, from the locked seed schedule and source configuration:
- frame-generation provenance and exact finite-frame means;
- sample indices and sample estimates;
- interval endpoints and success/failure classification;
- inner resampling attempt accounting and failure categories;
- conditional and unconditional estimands;
- summary counts, widths, and Monte Carlo standard errors;
- all expected scenario × frame × design × repetition × method keys.

Checksum verification alone is insufficient. At least one adversarial test per major artifact field must demonstrate that a checksum-refreshed but scientifically inconsistent artifact is rejected.

## 8. Stopping, exclusions, and deviations

No optional stopping based on favourable coverage or method ranking is permitted. Run the frozen number of frames and repetitions unless a documented infrastructure or safety failure makes continuation impossible.

Do not discard a frame because its distributional realization is inconvenient, unusually skewed, or produces numerical failures. Exclusions are permitted only for prespecified, outcome-independent provenance failures (for example, corrupted input or irreproducible seed) and must be listed with reasons before outcome summaries are interpreted.

Any deviation from the locked design must be recorded in a deviation log, with the original and revised settings, date, reason, affected outputs, and whether the affected analysis remains confirmatory. A modified analysis after outcome inspection is exploratory unless independently rerun under a new frozen design.

## 9. Reporting and release gates

A confirmatory release requires all of the following:

1. Exact-head CI passes on the supported Python versions.
2. The frozen manifest and environment metadata are archived.
3. The complete v2 run finishes with all attempts represented.
4. The independent artifact auditor returns PASS.
5. An independent reconstruction reproduces primary summary values within declared numerical tolerances.
6. The manuscript reports all prespecified estimands, failure rates, frame-level variability, uncertainty, deviations, and limitations.
7. The release states clearly that this linear benchmark does not establish validity for nonlinear dark-siren inference.

If any gate fails, the artifact may be released as a diagnostic or pilot with the failed gate prominently stated, but it must not be described as a confirmatory result.

## 10. Interpretation boundary

This protocol can support a careful finite-population simulation comparison. It cannot establish a new bootstrap theorem, guarantee performance for arbitrary finite populations, or validate the nonlinear hierarchical likelihood and posterior calibration used in dark-siren cosmology. A separate nonlinear study must define its own estimands, simulation truth, calibration tests, matched controls, and independent audit.

## 11. References

Wang, R., Peng, L., & Kim, J. (2022). Bootstrap inference for the finite population mean under complex sampling designs. *Journal of the Royal Statistical Society: Series B*, 84(4), 1150–1174. https://doi.org/10.1111/rssb.12506

This reference motivates the pseudo-population bootstrap family; the present software benchmark is not a proof or independent verification of the theorem's assumptions or conclusions.

### Planning-input sensitivity analysis

The executable planner also exposes a planning-only sensitivity grid over:
- multiplicative stresses to the externally supplied pilot frame-level SD; and
- alternative confidence levels for the upper SD bound.

The grid reports the resulting minimum frame count for each combination and retains the `PLANNING_ONLY_NOT_CONFIRMATORY` label. Regression tests check that, holding other inputs fixed, the planned frame count is non-decreasing as the SD stress multiplier or variance-confidence level increases. This is a software consistency check, not evidence that the pilot SD is representative or that the normal-theory bound is valid for a particular scientific estimand.

The sensitivity grid must be specified before examining confirmatory outcomes. It is intended to expose how design recommendations depend on uncertain planning inputs, not to select whichever row yields a preferred or affordable sample size.
