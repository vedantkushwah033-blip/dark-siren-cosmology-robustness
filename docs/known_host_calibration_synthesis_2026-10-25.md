# Known-host calibration synthesis and next decision gate (2026-10-25)

## Purpose

This note consolidates the completed known-host calibration surface and the later two-host posterior-level calibration experiment. It is a decision note, not a new analysis result. It does not replace either experiment's raw artifact or independent audit.

## Evidence reviewed

### A. Single-known-host calibration surface

- Prespecified design: 9 fixed H0 truths × 4 fractional distance-noise levels × 300 replicates = 10,800 raw records.
- Prior support: uniform H0 on [50, 90] km s^-1 Mpc^-1; posterior grid: 1,601 points.
- Simulator revision pinned to `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`.
- Raw-record audit independently reconstructed posterior summaries and checked the fixed design and hashes; paired contrast audit independently reconstructed 135 comparisons and bootstrap intervals.
- Report: [known_host_calibration_run_20261025.md](known_host_calibration_run_20261025.md).

The seven interior truth cells (55–85) show increasing interval widths and larger median error as distance noise increases. The report's equally weighted interior-truth summaries show 68% coverage from 0.695 at 5% noise to 0.770 at 20% noise; 95% coverage from 0.954 to 0.983. This is compatible with a tendency toward conservative intervals in this matched synthetic setup, but cell-level Wilson intervals should govern any stronger statement.

The exact prior-boundary truths (50 and 90) have zero coverage for the reported central equal-tail intervals. This is a structural boundary/prior-support issue, not evidence that the same procedure has zero coverage for interior truths. Keep boundary truths separate from interior calibration summaries.

### B. Two-host posterior-level calibration

- 2,400 synthetic observations across four injected truths (55, 65, 75, 85), three noise levels (3%, 10%, 20%), and 200 replicates per truth/noise cell.
- Eight inference arms cross prior support ([50,90] or [20,220]), prior measure (uniform or log-uniform H0), and catalogue weights (correct or reversed).
- Independent audit reconstructed posterior medians and 68%/95% interval endpoints from raw observations and verified input hashes.
- Report: [posterior_calibration_results_20261025.md](posterior_calibration_results_20261025.md).

Within that specific two-host mixture, broad support plus reversed weights produced pooled median bias +48.36 km s^-1 Mpc^-1 (stratified-bootstrap 95% interval +46.98 to +49.72) and 68% coverage 30.3% (95% interval 28.4%–32.1%). This is strong evidence of failure **under this specified synthetic mixture and misspecification**, not a general estimate of the effect in real gravitational-wave data. The broad prior is a diagnostic arm, not a recommended astrophysical prior.

## What the two experiments do and do not establish

They provide useful, independently checked software and repeated-sampling evidence for two different synthetic models:
1. the single-host, matched low-redshift likelihood behaves in a measurable way across fixed truths and noise levels; and
2. in the chosen two-host mixture, prior support, prior measure, and catalogue-weight misspecification can materially change posterior summaries and coverage.

Do **not** compare the numerical coverage values across A and B as if they were a head-to-head method comparison: the data-generating models, truth grids, noise levels, prior arms, and estimands differ. An audit PASS establishes that the stored outputs are reproducible under the audited definitions; it does not establish astrophysical realism or validate a real-event H0 analysis.

## Priority next step: calibration by design, not another unplanned sweep

Before expanding the scenario matrix again, produce one compact calibration scorecard from already-audited artifacts, with the following rules:

1. **Keep estimands separate.** Report fixed-truth frequentist coverage, point-error metrics, and simulation-based calibration (SBC) rank diagnostics in separate sections. Do not use the term “SBC” for a fixed-truth coverage experiment.
2. **Stratify before pooling.** Show truth × noise coverage and Wilson intervals first. Then give any pooled summary with the weighting rule stated explicitly. Do not let high-coverage cells conceal failures in another cell.
3. **Boundary-aware interpretation.** Mark truth-at-prior-boundary cells separately and explain why central equal-tail intervals are not a useful nominal-coverage check for an exact boundary truth under this prior.
4. **Prior-predictive calibration.** Use the existing SBC implementation only after verifying that the generative prior, likelihood, and inference prior are the same model. Uniform ranks should be assessed with uncertainty and dependence checks; rank uniformity does not validate astrophysical realism.
5. **Model-checking gate.** Document which assumptions are deliberately idealized (linear low-redshift distance relation, known redshift/host, Gaussian distance error, no detector selection or realistic event likelihood). A later realistic extension must add these ingredients explicitly rather than implying the toy model already covers them.
6. **No new grand claim from repeated simulation alone.** The next result should answer a prespecified question, include raw artifact hashes, independent reconstruction, and Monte Carlo uncertainty. If it only repeats that a deliberately wrong weighting scheme performs poorly in the same toy setup, treat it as a robustness demonstration rather than a novel cosmology result.

## Go/no-go recommendation

**Continue, but consolidate before adding more factorial dimensions.** The current work has passed meaningful reproducibility checks, but the scientific contribution remains a simulation-methods result. The next milestone should be a compact, audited synthesis and a validated prior-predictive SBC diagnostic, followed by a targeted extension with a clearly motivated selection mechanism. Do not describe this project as a real H0 measurement or as validating catalogue-incompleteness correction for LVK data.
