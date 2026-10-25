# Paired host-weight sweep with prior-support factor — results (2026-10-25)

**Status:** completed and independently audited.  
**Workflow:** [run 40168837084](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40168837084).  
**Artifact:** [paired-weight-sweep-results ZIP](https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12589821737/zip).  
**Audit:** PASS; 9,000 generated observations, 18 prior-support/design summaries, 36,000 posterior fits, and 4,000 paired bootstrap resamples per metric/cell. The independent auditor reconstructed all summary estimates and bootstrap intervals and verified the raw/summary SHA-256 hashes.  
**Protocol:** [paired-weight sweep protocol](paired_weight_sweep_protocol.md).

## Design

The injected value is $H_0=70$ km s$^{-1}$ Mpc$^{-1}$. The data generator selects between hosts at $z=0.02$ and $z_2\in\{0.025,0.04,0.06\}$ with weights $(0.7,0.3)$. For each observation, the same measured distance and noise scale are analyzed using correct weights $(0.7,0.3)$ and reversed weights $(0.3,0.7)$. Each paired comparison is run under both a narrow $H_0$ prior $[50,90]$ and a deliberately broad prior $[20,220]$. The grid spacing is 0.025 in both arms.

The reported RMSE and coverage differences are **reversed minus correct weights**. RMSE confidence intervals and coverage/absolute-error intervals are deterministic percentile-bootstrap 95% intervals based on paired replicate resampling.

## Results

Positive RMSE differences indicate worse median point estimation under reversed weights. Coverage differences are percentage-point changes in 68% interval inclusion.

| Prior support | $z_2$ | Distance error | RMSE difference | 95% bootstrap CI | Coverage difference | 95% bootstrap CI | Mean absolute-error difference |
|---|---:|---:|---:|---:|---:|---:|---:|
| Narrow [50,90] | 0.025 | 3% | +3.41 | [+2.90, +3.90] | -16.0 pp | [-19.1, -12.8] | +3.89 |
| Broad [20,220] | 0.025 | 3% | +6.30 | [+5.75, +6.86] | -21.9 pp | [-25.3, -18.6] | +5.90 |
| Narrow [50,90] | 0.025 | 10% | +0.98 | [+0.78, +1.17] | -7.7 pp | [-10.2, -5.4] | +1.10 |
| Broad [20,220] | 0.025 | 10% | +5.03 | [+4.75, +5.31] | -15.9 pp | [-18.9, -13.0] | +4.10 |
| Narrow [50,90] | 0.025 | 20% | +0.41 | [+0.31, +0.50] | -2.6 pp | [-3.8, -1.4] | +0.51 |
| Broad [20,220] | 0.025 | 20% | +4.24 | [+4.06, +4.42] | -9.5 pp | [-11.6, -7.3] | +3.81 |
| Narrow [50,90] | 0.040 | 3% | +0.00 | [-0.00, +0.00] | +0.0 pp | [+0.0, +0.0] | -0.00 |
| Broad [20,220] | 0.040 | 3% | +38.34 | [+36.24, +40.37] | -39.2 pp | [-42.8, -35.8] | +34.40 |
| Narrow [50,90] | 0.040 | 10% | -0.08 | [-0.13, -0.04] | +2.2 pp | [+1.3, +3.2] | -0.05 |
| Broad [20,220] | 0.040 | 10% | +31.29 | [+29.99, +32.59] | -31.3 pp | [-34.8, -27.9] | +23.49 |
| Narrow [50,90] | 0.040 | 20% | +0.67 | [+0.56, +0.77] | +0.0 pp | [+0.0, +0.0] | +0.74 |
| Broad [20,220] | 0.040 | 20% | +16.50 | [+15.81, +17.18] | -25.3 pp | [-28.6, -22.0] | +13.22 |
| Narrow [50,90] | 0.060 | 3% | +0.00 | [-0.00, +0.00] | +0.0 pp | [+0.0, +0.0] | +0.00 |
| Broad [20,220] | 0.060 | 3% | +38.57 | [+35.71, +41.59] | -44.4 pp | [-47.8, -41.0] | +46.17 |
| Narrow [50,90] | 0.060 | 10% | -0.00 | [-0.00, -0.00] | +0.0 pp | [+0.0, +0.0] | -0.00 |
| Broad [20,220] | 0.060 | 10% | +56.37 | [+53.71, +59.03] | -43.9 pp | [-47.3, -40.4] | +49.85 |
| Narrow [50,90] | 0.060 | 20% | -0.06 | [-0.12, +0.00] | +1.2 pp | [+0.4, +2.0] | -0.04 |
| Broad [20,220] | 0.060 | 20% | +35.55 | [+34.09, +36.92] | -37.6 pp | [-40.8, -34.2] | +29.69 |

RMSE and absolute-error differences are in km s$^{-1}$ Mpc$^{-1}$. “pp” means percentage points. Displayed zero contrasts are numerical zeros to the reporting precision.

## Interpretation

### 1. Prior support changes whether weights appear to matter

Under the narrow prior, the weighting contrasts are exactly or nearly zero in several $z_2=0.04$ and $0.06$ cells. The implied alternative-host $H_0$ modes fall outside $[50,90]$ for these redshift pairs. Once that likelihood component is cut off by the prior, reversing its weight cannot meaningfully change the normalized posterior. This is a support-truncation effect, not evidence that host weights are generally unimportant.

The broad-prior arm retains the alternative-host modes and shows much larger differences. For example, at $z_2=0.06$ and 10% distance error, the reversed-weight RMSE is 90.12 versus 33.75 under correct weights; the paired RMSE difference is +56.37 (95% CI [+53.71,+59.03]), and 68% coverage falls by 43.9 percentage points (95% CI [-47.3,-40.4]).

### 2. The broad prior is a diagnostic, not a recommended prior

Even the correct-weight broad-prior arm has high RMSE in some cells. A broad uniform prior allows distinct host-redshift explanations to map to widely separated $H_0$ modes. The posterior median and central interval then depend strongly on mixture mass and prior support. This is a demonstration of multimodality and prior sensitivity in the toy model—not a claim that a broad prior or the resulting numerical errors are appropriate for real dark-siren analyses.

### 3. The narrow-prior results are not universally favorable to correct weights

A few narrow-prior cells have small negative RMSE differences under reversed weights, with intervals excluding zero. This is possible in a finite, nonlinear posterior mixture: a misspecified weighting can occasionally move the posterior median closer to the injected value. It does not make the reversed model generally superior. The correct-weight model is the data-generating reference, but individual loss functions and prior truncation can produce non-monotone paired outcomes.

## Audit and reproducibility

- Raw rows: `per_run.csv` (9,000 records; each observation reused across all four inference variants).
- Summary rows: `summary.csv` (18 prior-support × redshift/noise cells).
- Provenance: `run_metadata.json`, including seed 20261028 and SHA-256 hashes.
- Independent audit: `audit_report.json`, status PASS. It validates design completeness, per-run median-error/interval-width/coverage identities, paired absolute-error differences, summary statistics, deterministic bootstrap intervals, and CSV hashes.
- Script: `scripts/run_paired_weight_sweep.py`; independent audit: `scripts/audit_paired_weight_sweep.py`.
- Workflow: `.github/workflows/paired-weight-sweep.yml`.

## Limitations and next step

This remains a one-truth, two-host, low-redshift synthetic experiment with Gaussian distance noise. It does not include galaxy-catalogue selection, realistic redshift uncertainty, detector selection, or real GW data. The next analysis should separate host-mixture ambiguity from prior-induced posterior mass by auditing the per-host marginal likelihood under the chosen $H_0$ prior, then repeat calibration across multiple injected $H_0$ values and prespecified prior supports before any astrophysical interpretation.
