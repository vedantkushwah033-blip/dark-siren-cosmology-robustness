# Prior-support × detection-threshold factorial audit

## Research question

How do prior boundaries and misspecification of the detection threshold jointly affect posterior point estimation, credible-interval width, and repeated-sampling coverage in the existing synthetic dark-siren model?

## Design

- Data-generating detection interval: 0 < observed distance < 120 Mpc.
- Inference thresholds: 120, 150, 180, 210 Mpc.
- Prior intervals: [50,90], [40,100], [30,110].
- Truths: H0 = 52, 70, 88.
- Host redshifts: 0.02 and 0.04, mixture weights 0.7 and 0.3.
- Gaussian distance noise: sigma = 20 Mpc.
- Each detected observation is reused across all 3 priors × 4 thresholds × 2 methods.
- Methods: uncorrected likelihood and selection-corrected likelihood.
- Default: 300 repetitions per truth, seed 20261028, 4,001 grid points.

The default design produces 21,600 per-run rows (3 truths × 300 replicates × 3 priors × 4 thresholds × 2 methods) and 72 summary rows. The paired contrast table compares every configuration against the reference configuration (prior [50,90], assumed threshold 120 Mpc), separately by truth and method, across six metrics. The compact CI job uses 100 repetitions and a 1,001-point grid for structural validation.

## Metrics and uncertainty

Summaries include bias, MAE, RMSE, 68% and 95% coverage, Wilson 95% intervals for coverage, and mean interval widths. Paired contrasts report mean differences and percentile-bootstrap intervals. The reference configuration is a fixed comparison baseline; it is not asserted to be universally optimal.

## Reproduction

```bash
python scripts/run_prior_threshold_factorial.py --replicates 300 --seed 20261028 --grid-size 4001 --output-dir results/prior-threshold-factorial
python -m pytest tests/test_prior_threshold_factorial.py -q
```

GitHub Actions runs focused tests and a compact experiment on Python 3.10, 3.11, and 3.12, uploading the raw rows, summary table, and paired contrasts.

## Interpretation guardrails

This is a simplified Gaussian distance-threshold model, not a realistic gravitational-wave detector selection pipeline. It omits SNR-dependent detection efficiency, detector duty cycles, waveform dependence, multiple-event population inference, redshift uncertainty, and calibrated galaxy-catalogue incompleteness. Truths 52 and 88 lie close to the baseline prior boundaries, so their coverage is strongly affected by prior truncation. The factorial design is intended to expose interactions, not to tune the prior or assumed threshold to whichever gives the most favorable results. Results must be interpreted with Monte Carlo uncertainty and must not be described as real LVK measurements.


## Full-run CI workflow

The workflow `.github/workflows/full-prior-threshold-factorial.yml` can be started from the repository's Actions tab. For the protocol run, select 300 repetitions, grid size 4001, and seed 20261028. It runs focused tests and validates expected output row counts and numeric finiteness before uploading `per_run.csv`, `summary.csv`, and `paired_contrasts.csv` for 90 days. This workflow is manually dispatched and does not run merely because this protocol was updated. Treat the experiment as complete only after checking the actual run conclusion and inspecting its artifacts.
