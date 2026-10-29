# Preregistered validation protocol: corrected host-density and selection pilot

**Status:** protocol and staged validation record. Stage A and the locked Stage B synthetic campaigns have completed with independent artifact audits; the sensitivity and broader-model gates below remain outstanding.  
**Scope:** synthetic methodological study only. This is not an LVK analysis, detector validation, or measurement of the real Hubble constant.

## 1. Research question

Under the repository's declared toy generative model, how much do catalogue incompleteness and host-density estimation change H0 posterior error, and can those effects be separated from KDE approximation error and the selection-conditioned likelihood?

The experiment must distinguish software correctness, numerical convergence, and statistical calibration. Passing one category does not imply the others.

## 2. Model and estimands

Use the declared latent galaxy intensity, merger-weight function, catalogue inclusion law, Gaussian distance-noise model, and fixed distance threshold. Hold the true H0 and other parameters fixed for the primary experiment; conclusions therefore apply only to this controlled setting.

The six paired arms are:

- `analytic_truth`: exact declared normalized target density on the redshift grid; diagnostic control only.
- `oracle`: KDE from the complete simulated parent population; not deployable in a real analysis.
- `catalogue_only`: observed thinned catalogue without inverse-inclusion correction.
- `correct_weighting`: observed catalogue weighted by the true declared inclusion probability.
- `shifted_transition`: deliberately incorrect inclusion-transition midpoint.
- `shallow_slope`: deliberately incorrect inclusion-transition width.

Use identical replicate seeds across arms. Truth-only columns may be used for scoring and the oracle arm, never for non-oracle fitting. Record the precise software commit, configuration, seeds, runtime, dependency versions, raw replicate-level metrics, and SHA-256 manifest.

### Primary estimand

For each replicate, calculate absolute posterior-mean error for `catalogue_only` and `correct_weighting`. Define the paired contrast

`D_r = |mean_H0(catalogue_only,r) - H0_true| - |mean_H0(correct_weighting,r) - H0_true|`.

Positive mean D indicates lower absolute error under correct weighting in this toy setting. Report the mean paired contrast, a 95% paired-bootstrap interval, its Monte Carlo standard error, the number of replicates, and the full replicate-level distribution. Do not treat individual events or grid points as independent replicates.

### Prespecified diagnostic contrasts

1. `oracle` versus `analytic_truth`: estimates the combined finite-parent/KDE approximation effect.
2. `catalogue_only` versus `analytic_truth`: describes the total error under catalogue incompleteness plus estimation error.
3. `correct_weighting` versus `analytic_truth`: checks whether known-probability correction recovers the diagnostic reference in finite samples.
4. `shifted_transition` and `shallow_slope` versus `correct_weighting`: assess the two declared forms of selection-law misspecification.

These are contrasts within the declared simulation, not universal statements about galaxy catalogues.

## 3. Two-stage simulation plan

### Stage A: exploratory pilot

- Use at least 200 independent replicate seeds, subject to practical runtime.
- Use the corrected Scott bandwidth by default.
- Preserve every arm for every replicate; fail the campaign if any replicate-arm record is missing or duplicated.
- Estimate the standard deviation of the primary paired contrast.
- Do not promote pilot point estimates to confirmatory evidence.

### Stage B: precision-based confirmation

Before generating confirmatory results, lock the code commit, configuration, analysis script, and seeds in a manifest. Choose the confirmatory replicate count using the Stage A paired standard deviation and a target 95% confidence-interval half-width of 0.5 km s^-1 Mpc^-1 for the primary contrast:

`N = ceil((1.96 * s_D / 0.5)^2)`.

Use at least 500 and at most 5,000 confirmatory replicates. If the calculated N exceeds 5,000, use 5,000 and explicitly report that the target precision was not achieved. Confirmatory seeds must not overlap pilot seeds. Do not stop early because a p-value or effect estimate looks favorable.

The target half-width is a design goal, not a guarantee. The normal approximation is only used for planning; report the empirical paired-bootstrap interval as the primary uncertainty interval.

## 4. Numerical and estimator sensitivity

Run these checks before making a scientific performance claim:

1. **H0-grid convergence:** compare 161, 321, and 641 H0 points on the same selected pilot seeds. Compare posterior means, 68% interval endpoints, and interval widths.
2. **Redshift-grid convergence:** compare 801 and 1,601 redshift points on the same seeds, with all other settings fixed.
3. **Bandwidth sensitivity:** compare Scott default against explicitly declared factors 0.15, 0.25, and 0.35 under the documented convention h = factor × weighted SD of the original sample. Do not describe these factors as SciPy covariance multipliers.
4. **Boundary/normalization checks:** independently integrate every estimated density and every saved posterior; reject nonfinite or nonpositive densities, invalid normalization, and posterior summaries inconsistent with the saved grid.
5. **Selection checks:** independently reconstruct the detection efficiency by direct quadrature and verify that every observed event meets the threshold.
6. **Artifact integrity:** verify checksums, replicate/seed pairing, all six arms, row counts, and independently reconstructed summaries. Rehashing a modified manifest must not bypass structural checks.

Set numerical tolerances before viewing results. If grid changes alter posterior means by more than 0.1 km s^-1 Mpc^-1 or any 68% endpoint by more than 0.2 km s^-1 Mpc^-1 on the sensitivity subset, increase resolution or label the result numerically unresolved. These are practical software acceptance thresholds, not physical accuracy standards.

## 5. Calibration and uncertainty reporting

For every arm, report bias, RMSE, empirical 68% and 95% interval coverage, mean interval width, and Wilson intervals for coverage. Report the number of independent replicates beside all coverage estimates. Include paired contrasts with confidence intervals, not only separate arm averages.

A nominal 68% interval's coverage estimate is noisy at small N. Do not call a difference meaningful merely because point estimates differ. For a fixed single H0 truth, coverage is only a repeated-sampling check at that truth; it does not establish calibration across the prior or parameter space. Later work must vary H0 and the relevant nuisance parameters using a predeclared design.

## 6. Stop/go gates

**GO to confirmatory interpretation only if all are true:**

- CI passes on the exact commit being evaluated.
- The smoke experiment completes and the independent artifact auditor passes.
- All six arms are present for every replicate and the truth-only data boundary is respected.
- Posterior grids reconstruct saved summaries within declared tolerances.
- Independent likelihood and detection-efficiency calculations agree within numerical tolerances.
- Grid and bandwidth sensitivity are reported, including any unresolved numerical dependence.
- The confirmatory run uses locked, non-overlapping seeds and a fixed analysis plan.

**NO-GO for a performance claim if any are true:**

- Exact-head CI is pending or failing.
- The corrected campaign has not been regenerated after estimator changes.
- Analytic-truth bias remains unexplained.
- Artifact validation fails or replicate pairing is broken.
- Numerical sensitivity changes the ranking of arms or materially changes the primary contrast.
- The result depends on post-hoc selection of seeds, bandwidths, or reported metrics.

## 7. Interpretation boundaries

Even a successful result supports only a claim about this synthetic model and its specified assumptions. It does not validate realistic gravitational-wave detector selection, galaxy-catalogue completeness in the real Universe, event-rate marginalization, peculiar-velocity treatment, uncertain redshifts, or an empirical H0 measurement. Those require separate models and independent validation.

## 8. Reporting template

Publish, at minimum:

- Commit SHA and environment versions.
- Locked configuration and seed ranges.
- Number of planned and completed replicates.
- Primary paired contrast with 95% interval and Monte Carlo uncertainty.
- Per-arm bias, RMSE, coverage with Wilson intervals, and interval widths.
- Grid and bandwidth sensitivity tables.
- Independent audit status and checksums.
- Deviations from this protocol, including failed runs.
- A concise statement of what the toy experiment cannot establish.

**Historical status note (superseded by Section 9):** at the time this protocol was first drafted, it was a prespecification proposal and did not imply that a pilot or confirmation had run. Section 9 records the subsequent locked design and its completed Stage B run. The result remains scoped to the fixed-H0 synthetic model.

## 9. Stage A calibration and locked Stage B design (2026-10-29)

Stage A completed with 200 paired replicates. See [the audited Stage A report](selection_misspecification_stage_a_report_2026-10-29.md) and [workflow run #1](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40793766654).

For the primary paired estimand define `D_r = |posterior_mean(catalogue_only,r)-H0_true| - |posterior_mean(correct_weighting,r)-H0_true|`. Stage A gives mean D = 6.535691 km s^-1 Mpc^-1, sample SD = 4.321762, and SE = 0.305595 (n=200). The precision formula gives `N = ceil((1.96 * 4.321762 / 0.5)^2) = 288`.

The minimum therefore governs: **Stage B is locked at 500 paired replicates**, with planning half-width 1.96 × 4.321762 / sqrt(500) ≈ 0.379 km s^-1 Mpc^-1, conditional on the Stage A SD being useful for planning. Confirmatory seeds are fixed at **20262030–20262529**, disjoint from Stage A seeds **20261030–20261229**.

Locked Stage B configuration: 4,000 parent galaxies per replicate; 30 detected events; fixed true H0 = 70 km s^-1 Mpc^-1; 81-point H0 grid over [50, 90]; 401-point redshift grid over [0.01, 0.40]; same selection law and six arms; corrected Scott bandwidth default. The primary result is mean D with a percentile-bootstrap 95% interval using 10,000 resamples; also report paired SE and normal interval. Deterministic bootstrap seed is `seed_start + 104729`.

This is prospective planning, not evidence of favorable confirmatory results. Do not change the seed range, replicate count, primary contrast, or interval method after observing Stage B output. Code changes affecting the estimand/audit invalidate this run and require a new locked commit and disjoint seed stream.


## 10. Current release gate (2026-10-29)

The locked Stage B primary result is complete and independently audited as documented in `selection_misspecification_stage_b_report_2026-10-29.md`. Do not rerun or relabel it as a fresh confirmation. The next scientific gate is a new, explicitly versioned robustness campaign over multiple fixed H0 values, event counts, completeness strengths, and host-population shapes, with grid/bandwidth sensitivity and independent artifact reconstruction. Any changes to the likelihood, estimator, primary endpoint, or simulation design require a new preregistration and disjoint seeds. The latest repository-head test workflow must also reach a terminal passing state before release; as of the current check, its run is pending. CHIMERA v2.2 remains NO-GO until the pinned, version-compatible input bundle is independently verified.
