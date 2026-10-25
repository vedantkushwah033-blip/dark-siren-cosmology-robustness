# Known-host H₀ calibration surface: prespecified protocol

**Status:** protocol and implementation documented; calibration outputs are not yet claimed here.  
**Purpose:** test the numerical and statistical behaviour of the repository's simplified one-host likelihood before interpreting more complicated catalogue-mixture experiments.

## Scope and estimand

This is a controlled synthetic calibration experiment, not an LVK analysis, not a measurement of the real Hubble constant, and not a validation of a full cosmological distance-redshift relation. It asks how posterior summaries behave when a single simulated host is known to have generated an event and the same host redshift is supplied to the likelihood.

The simulator generates luminosity distance with the low-redshift relation (d_L = c z/H_0). The inference likelihood uses that same relation and a Gaussian distance likelihood. This matched-model design is useful for debugging and measuring conditional calibration, but cannot quantify model misspecification from realistic cosmology, peculiar velocities, selection effects, detector calibration, or uncertain host identity.

## Prespecified design

| Dimension | Values |
| --- | --- |
| True H₀ | 50, 55, 60, 65, 70, 75, 80, 85, 90 km s⁻¹ Mpc⁻¹ |
| Fractional distance uncertainty | 0.05, 0.10, 0.15, 0.20 |
| Replicates per truth × uncertainty cell | 300 by default; workflow accepts 100–2,000 |
| Number of cells | 9 × 4 = 36 |
| Raw observations at default | 36 × 300 = 10,800 |
| Prior support | Uniform on [50, 90] km s⁻¹ Mpc⁻¹ |
| Posterior grid | 1,601 equally spaced points |
| Simulator revision | `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f` |

The two boundary truths, 50 and 90, are deliberate stress cases: a posterior is truncated by the prior at these values, so their behaviour should not be pooled indiscriminately with interior truths.

## Pairing and reproducibility

Within each true-H₀ value, the catalogue seed and host-selection seed are held fixed across the four distance-uncertainty levels. The event-noise seed is also held fixed, so the same standardized Gaussian draw is reused when the uncertainty changes. This common-random-number design reduces Monte Carlo noise in comparisons between uncertainty levels. It does **not** make the four observations independent; they are paired by replicate. Analysis should compare matched replicate-level differences where possible.

Seed blocks are disjoint across the nine true-H₀ strata: the seed for truth index j and replicate r is `base_seed + j × replicates + r`. This avoids accidentally reusing identical redshift/luminosity catalogues and host draws across different truths. The pairing unit remains the replicate within one fixed truth across the four uncertainty levels; no cross-truth pairing is presumed. Results remain conditional on this particular simulator and host-generation rule.

## Quantities to report

For every cell, the implementation records:

- empirical coverage of central 68% and 95% posterior intervals, with Wilson binomial intervals and Monte Carlo standard errors;
- median bias and mean posterior-median error, reported as distinct quantities;
- RMSE of the posterior median and of the posterior mean;
- mean interval widths and posterior standard deviation;
- median-pull diagnostics and the fraction of posterior medians within 0.5 H₀ units of either prior boundary.

Coverage is a repeated-sampling property of the procedure under this simulation protocol. It is not the probability that a particular observed interval contains the true value. A coverage estimate near its nominal target should be considered together with its uncertainty interval and the per-cell sample size.

## Audit and acceptance checks

The workflow is designed to preserve the raw per-run CSV, cell summary, metadata, an independent audit report, and a Markdown report. The independent auditor recomputes cell statistics from raw records instead of trusting the experiment's stored summary implementation. The report renderer checks the prespecified 36-cell grid, replicate counts, finite values, interval nesting/width relationships, coverage/Wilson consistency, metadata counts, and input hashes.

A result should be interpreted only after all of the following are true:

1. the workflow run has completed successfully;
2. the artifact is available and contains the raw records, summaries, metadata, and audit output;
3. the independent audit explicitly reports PASS and confirms posterior summaries were reconstructed from raw observations;
4. the report's checksums match the audited inputs;
5. the report is read alongside boundary-specific results and Monte Carlo uncertainty.

A successful unit-test suite validates code paths and fixtures; it does not substitute for a completed experiment or prove that the scientific model is realistic.

## How to run

In GitHub, open the repository's **Actions** tab, choose **known-host-calibration-surface**, click **Run workflow**, and retain the defaults for the prespecified run:

- `replicates = 300`
- `seed = 20261025`

The workflow checks out the pinned simulator commit, executes the surface, independently audits it, renders the report, and uploads the result directory as the `known-host-calibration-surface` artifact. If changing the replicate count or seed for exploratory work, document the change and do not silently compare it as though it were the default protocol.

## What this experiment cannot establish

Even a successful audit does not establish realistic dark-siren H₀ coverage. The model has a known host, Gaussian distance noise, a low-redshift linear distance law, a finite uniform prior, and a matched generation/inference model. The results should be described as a calibration and numerical diagnostic for that model only. Subsequent work should test realistic cosmological distance-redshift relations, uncertain host identity, selection effects, catalogue incompleteness, and model misspecification separately, with clear validation targets for each.
