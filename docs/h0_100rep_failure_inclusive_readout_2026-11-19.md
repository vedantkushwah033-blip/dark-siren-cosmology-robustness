# Failure-inclusive uncertainty analysis for the 100-replicate H0 pilot

## Purpose

Success-only posterior coverage can hide operational failures. This post-hoc analysis reports both calibration conditional on an evaluable posterior and conservative all-attempt bounds that preserve every preregistered cell.

## Inputs and integrity gates

The analysis requires:
- `planned_cells.csv` — the preregistered cell inventory and seeds;
- `replicate_ledger.csv` — every planned attempt, including SUCCESS, NON_EVALUABLE, FAILED, RUNNING, and NOT_RUN;
- `manifest.json` — SHA-256 hashes for the aggregate files.

Before computing metrics, the script verifies every manifest hash and exact equality between the planned inventory and the result ledger. It rejects missing/duplicate cells, changed seeds, unknown statuses, non-finite metrics in SUCCESS rows, and malformed posterior quantiles.

## Outputs

Run:

```bash
python scripts/analyze_h0_100rep_calibration_uncertainty.py <artifact_dir> --output-dir <output_dir>
```

Outputs:
- `cell_uncertainty_summary.csv`: per-cell planned/success/non-success counts, bias t intervals, RMSE bootstrap intervals, interval width, weight diagnostics, success-only Wilson coverage intervals, and all-attempt coverage bounds.
- `paired_gamma_contrasts.csv`: paired gamma-versus-zero differences in signed median bias, absolute median error, and 68% interval width.
- `paired_validation_size_contrasts.csv`: paired largest-versus-smallest validation-size differences for those metrics.
- `posthoc_report.json`: manifest status, ledger status counts, completion status, and interpretation limits.

## Coverage estimands

For a nominal interval level, let (k) of (s) SUCCESS cells cover the fixed truth, with (m) planned cells and (u=m-s) non-successes.

- **Success-only coverage:** (k/s), with a Wilson interval. This describes coverage conditional on producing an evaluable posterior.
- **All-attempt lower bound:** (k/m), treating every non-success as uncovered.
- **All-attempt upper bound:** ((k+u)/m), treating every non-success as covered.

The bounds are operating-characteristic bounds, not estimates of latent coverage for failed cells. They expose the impact of non-evaluable, failed, running, or unstarted attempts without silently dropping them.

## Pairing

Contrasts use only replicates successful in both compared cells. Before calculating a paired contrast, the script checks equality of analysis seed, validation seed, and H0 truth. Paired sample counts are reported. If failures remove pairs, contrasts remain conditional on joint success and must be interpreted alongside the all-attempt bounds.

## Scope

This is fixed-truth repeated-sampling coverage in a synthetic low-redshift Euclidean model, conditional on exactly four detected events. It is not joint SBC, event-rate calibration, real galaxy-survey validation, LVK/GWTC analysis, or an observational H0 measurement. One hundred replicates per cell improve Monte Carlo precision but do not make this a confirmatory calibration study.
