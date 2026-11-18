# Failure-inclusive uncertainty readout for the 100-replicate H0 pilot

## Purpose
The 600-cell follow-up reports fixed-truth H0 coverage among SUCCESS cells and separately records FAILED, NON_EVALUABLE, and NOT_RUN attempts. Reporting success-only coverage alone can hide operating failures. This post-hoc analysis therefore preserves both views.

## Outputs
Run `scripts/analyze_h0_100rep_calibration_uncertainty.py <artifact_dir> --output-dir <output_dir>` after the 600-cell artifact has passed its independent orchestration audit.

It writes:
- `cell_uncertainty_summary.csv`: planned/success/non-success counts, mean posterior-median bias and t interval, posterior-median RMSE and bootstrap interval, interval width, weight diagnostics, success-only Wilson coverage intervals, and all-attempt coverage bounds.
- `paired_gamma_contrasts.csv`: paired gamma-versus-zero differences in signed bias, absolute median error, and 68% interval width, with t and bootstrap intervals.
- `paired_validation_size_contrasts.csv`: paired largest-versus-smallest validation-size contrasts for the same metrics.
- `posthoc_report.json`: status counts and interpretation limits.

## Coverage estimands
For a nominal interval level, let (k) of (s) successful cells cover the fixed truth, with (m) planned cells and (u=m-s) non-successes.

- Success-only coverage: (k/s), with a Wilson interval. This describes calibration conditional on the run producing an evaluable posterior.
- All-attempt lower bound: (k/m), treating every non-success as uncovered.
- All-attempt upper bound: ((k+u)/m), treating every non-success as covered.

The lower/upper bounds are intentionally conservative operating-characteristic bounds, not estimates of latent coverage for failed cells. They make the effect of missing outcomes visible without silently dropping them.

## Pairing
Contrasts are computed only on replicates with SUCCESS in both compared cells. The script verifies that matched replicate IDs share the same analysis seed, validation seed, and H0 truth. The paired count is reported; if failures remove pairs, the contrast remains conditional on jointly successful runs and must be read alongside all-attempt bounds.

## Limits
This is fixed-truth repeated-sampling coverage in a synthetic low-redshift Euclidean model, conditional on exactly four detected events. It is not joint SBC, event-rate calibration, real galaxy-survey validation, LVK/GWTC analysis, or an observational H0 measurement. One hundred replicates per cell reduce but do not eliminate Monte Carlo uncertainty. No inference runner or production likelihood code is changed.
