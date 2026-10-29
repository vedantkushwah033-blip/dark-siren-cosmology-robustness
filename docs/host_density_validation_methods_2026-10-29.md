# Host-density validation protocol and reproducibility status

Date: 2026-10-29. Scope: synthetic low-redshift (H_0) toy only.

## What changed

The benchmark now accepts a physical KDE bandwidth factor, defaulting to weighted Scott's rule, and deliberately misspecified catalogue-inclusion midpoint/width values for the inference correction. The generator still uses its fixed true selection law (floor 0.1, midpoint 0.22, width 0.05). The scoring-only target density is computed from the known toy galaxy intensity multiplied by the declared merger weight.

Each replicate records L1 host-density error, maximum-to-mean inverse-weight ratio, and effective-sample-size fraction alongside posterior mean/median, equal-tail 68% interval, and coverage. The command-line runner writes a manifest with SHA-256 hashes for its CSV and summary. The independent audit re-simulates selected replicate seeds, uses a separately implemented reflected-kernel KDE and Simpson quadrature on finer grids, and checks posterior mean/median/interval endpoints against fixed tolerances.

## Paired sensitivity design

The expedited screening campaign uses 40 replicates per setting, initial seed 20261029, 6000 parent galaxies, 40 detected events, Gaussian distance noise (sigma=140) Mpc, threshold 250 Mpc, and fixed truth (H_0=70). Settings are Scott plus physical bandwidth factors 0.05, 0.1, 0.2, 0.5, 0.8, and 1.2; coarse and fine numerical grids; and selection-law assumptions that shift the midpoint or widen/narrow the transition. The script computes arm bias/SE, RMSE, Wilson intervals for empirical 68% coverage, mean interval width, host-density L1 error, weight diagnostics, and paired absolute-error differences versus the oracle arm.

The workflow runs the campaign and uploads raw tables, aggregate summary, and SHA-256 manifest as a 90-day artifact. A separate smoke run is re-simulated and reconstructed by the independent auditor. The supplied patch contains an earlier locally generated campaign snapshot, but this workflow is being used to reproduce it on the current code; do not treat the snapshot as CI verified until the new job completes.

## Interpretation guardrails

- A correctly specified inclusion law is a strong toy-model control; real surveys rarely provide the exact true inclusion law.
- A finite weight ratio and reasonable effective sample size in this toy do not guarantee stable inverse-probability weighting under other selection laws.
- Oversmoothing can flatten the inferred host density and change (H_0) error without correcting catalogue incompleteness.
- Prior-averaged Bayesian calibration and fixed-truth estimator bias answer different questions. A posterior mean can be biased at fixed truth while remaining calibrated when truth is drawn from the prior.
- The study does not model a realistic GW detector, event-rate uncertainty, sky/inclination marginalization, peculiar velocities, redshift uncertainty, or full galaxy-survey hierarchical selection.

## Scientific gate

GO for exploratory synthetic-methods reporting with explicit limitations. NO-GO for claims that this result validates real LVK inference, universally corrects incompleteness, or measures real (H_0). No human expert review or CHIMERA v2.2 benchmark-data recheck is claimed by this protocol.


## Expedited execution and decision rule

The first 100-per-setting CI run remained in progress longer than useful for rapid iteration. The current workflow therefore uses 40 paired replicates per configuration and 300 replicates per oracle truth regime on coarser 401-redshift/81-H0 grids. This is explicitly a screening campaign. Do not use its empirical coverage as a high-precision estimate or make a confirmatory claim from the exploratory sensitivity surface. For any configuration whose conclusion appears decision-relevant, extend it using a new, disjoint, preregistered seed range and a precision target based on paired variability. The independent fine-grid smoke audit is separate.
