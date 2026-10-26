# Posterior-level H₀ calibration: corrected and independently audited results

**Status: corrected-model workflow PASS; generation and posterior-reconstruction audits PASS.** This report supersedes the pre-correction numerical tables. Results are from a synthetic two-host model only, not a real gravitational-wave measurement.

> **Tail-stable artifact audit (26 October 2026):** I inspected the latest artifact from [workflow run 40337117092](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40337117092), generated from revision `7943eecf3147d9e119cc5584705d8b141da830c4`. The artifact ([ID 12658030412](https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12658030412/zip)) passed the raw-generation, posterior-reconstruction, and paired-contrast audits. I verified the SHA-256 hashes recorded in `run_metadata.json`, independently aggregated the 96 summary rows into the eight pooled arms, and confirmed that the displayed pooled estimates below match the tail-stable artifact to rounding. The tail-stable rerun's raw and summary CSV hashes are `1714ccdc34697f831cf7777e1076552bd4506590740248c6505377af127153ab` and `4cdb79913cf1b2a14cf1851d54cb322ca8a9c5defff025a36fb4bbf23ff8e30c`, respectively. Identical hashes to the immediately preceding artifact show that this stability correction did not change these particular experiment outputs at the recorded precision; it does not make the correction unnecessary for extreme-tail cases. All results remain conditional on the stated synthetic two-host model and are not astrophysical validation.

## Provenance and design

- Workflow: [Posterior-level H₀ calibration, run 40337117092](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40337117092)
- Source revision used by that run: `7943eecf3147d9e119cc5584705d8b141da830c4` (includes the stable positive- and negative-tail interval-probability implementation and corrected tail regression tests).
- Artifact: [posterior-calibration-results.zip, artifact 12658030412](https://api.github.com/repos/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/artifacts/12658030412/zip)
- Seed: `20261030`; 200 replicates per each of 12 truth/noise cells; 2,400 raw observations and 96 summary rows.
- Injected H₀ truths: 55, 65, 75, 85 km s⁻¹ Mpc⁻¹.
- Reference noise settings: 3%, 10%, 20% of the expected distance at z = 0.04. The same absolute sigma is used for both hosts, so the fractional noise for the z = 0.02 host is twice the labelled reference setting.
- Inference arms: narrow support [50, 90] and broad diagnostic support [20, 220], uniform-H₀ and log-uniform-H₀ priors, and correct or reversed host weights.
- Posterior grid: as configured in the experiment runner; all arms use the same generated observations.

### Scientific correction

The generator redraws nonpositive Gaussian distance observations until positive. The prior likelihood treated each host component as an ordinary Gaussian, which did not exactly match that generator. The corrected likelihood uses the host-specific positive-truncated density

\[
p(d\mid H_0,j,d>0)=
\frac{\mathcal N(d;\mu_j(H_0),\sigma)}
{\Phi(\mu_j(H_0)/\sigma)},\qquad
\mu_j(H_0)=\frac{cz_j}{H_0}.
\]

Because the generator selects the host first and then redraws distance for that host until positive, the host-mixture weights are not changed by this truncation; the normalization is applied separately to each host component. The corrected runner and independent auditor both implement this likelihood, and a regression test independently reconstructs a posterior median by quadrature.

### Audit results

The archived workflow completed successfully. Its artifact records:

- Generation-provenance audit: **PASS** — all 2,400 raw observations reproduced from the declared seed and process; input hash verified.
- Posterior reconstruction audit: **PASS** — all posterior medians and 68%/95% interval endpoints reconstructed from raw observations; summary metrics reconstructed; file hashes verified.
- Paired-contrast audit: **PASS** — 48 truth/noise-cell contrast rows and 24 aggregate contrast rows independently reconstructed; 4,000 bootstrap resamples and input hashes verified.

The artifact metadata reports these checks explicitly. This is evidence of reproducibility and internal consistency for the synthetic experiment, not astrophysical validation.

## Pooled results

Each inference arm is evaluated on the same 2,400 observations, equally allocated across the 12 truth/noise cells. Bias is the mean posterior-median error (median minus injected H₀); MAE and RMSE are in km s⁻¹ Mpc⁻¹. Coverage is the fraction of intervals containing the injected truth.

| Prior support | Prior measure | Host weights | Bias | MAE | RMSE | Mean 68% width | 68% coverage | Mean 95% width | 95% coverage |
|---|---|---|---:|---:|---:|---:|---:|---:|---:|
| [50, 90] | Uniform H₀ | Correct | +0.63 | 6.24 | 8.33 | 16.68 | 68.9% | 26.87 | 98.0% |
| [50, 90] | Uniform H₀ | Reversed | +2.28 | 6.52 | 9.03 | 16.34 | 70.9% | 26.74 | 97.8% |
| [50, 90] | Log-uniform H₀ | Correct | −0.69 | 6.13 | 8.31 | 16.63 | 71.7% | 26.91 | 97.8% |
| [50, 90] | Log-uniform H₀ | Reversed | +1.01 | 6.31 | 8.52 | 16.68 | 69.9% | 26.92 | 97.8% |
| [20, 220] | Uniform H₀ | Correct | +14.28 | 30.09 | 37.88 | 74.84 | 59.7% | 109.92 | 93.1% |
| [20, 220] | Uniform H₀ | Reversed | +48.62 | 50.22 | 60.52 | 68.23 | 30.2% | 110.13 | 86.7% |
| [20, 220] | Log-uniform H₀ | Correct | +1.50 | 22.40 | 28.94 | 67.65 | 66.2% | 104.42 | 95.0% |
| [20, 220] | Log-uniform H₀ | Reversed | +40.07 | 42.57 | 51.85 | 73.50 | 52.4% | 109.08 | 93.7% |

## Monte Carlo uncertainty

The workflow reports 95% percentile-bootstrap intervals from 4,000 resamples, stratified within each of the 12 truth/noise cells. These intervals quantify Monte Carlo uncertainty conditional on this fixed simulation design; they are not posterior credible intervals for H₀ and do not capture uncertainty about whether the toy model represents nature.

| Support | Prior measure | Weights | Bias (95% CI) | RMSE (95% CI) | 68% coverage (95% CI) | 95% coverage (95% CI) |
|---|---|---|---:|---:|---:|---:|
| [50,90] | Uniform H₀ | Correct | +0.63 [0.43, 0.82] | 8.33 [8.12, 8.54] | 68.9% [67.4%, 70.5%] | 98.0% [97.4%, 98.5%] |
| [50,90] | Uniform H₀ | Reversed | +2.28 [2.08, 2.47] | 9.03 [8.83, 9.23] | 70.9% [69.4%, 72.5%] | 97.8% [97.2%, 98.3%] |
| [50,90] | Log-uniform H₀ | Correct | −0.69 [−0.89, −0.48] | 8.31 [8.08, 8.54] | 71.7% [70.0%, 73.4%] | 97.8% [97.2%, 98.3%] |
| [50,90] | Log-uniform H₀ | Reversed | +1.01 [0.82, 1.21] | 8.52 [8.31, 8.71] | 69.9% [68.3%, 71.5%] | 97.8% [97.2%, 98.4%] |
| [20,220] | Uniform H₀ | Correct | +14.28 [13.01, 15.53] | 37.88 [36.92, 38.80] | 59.7% [57.7%, 61.6%] | 93.1% [92.1%, 94.1%] |
| [20,220] | Uniform H₀ | Reversed | +48.62 [47.24, 49.98] | 60.52 [59.46, 61.54] | 30.2% [28.3%, 32.0%] | 86.7% [85.3%, 88.0%] |
| [20,220] | Log-uniform H₀ | Correct | +1.50 [0.42, 2.65] | 28.94 [28.19, 29.70] | 66.2% [64.3%, 67.9%] | 95.0% [94.2%, 95.9%] |
| [20,220] | Log-uniform H₀ | Reversed | +40.07 [38.81, 41.37] | 51.85 [50.90, 52.77] | 52.4% [50.3%, 54.4%] | 93.7% [92.7%, 94.6%] |

The artifact also includes `cell_uncertainty.csv` and `cell_uncertainty_report.md`, with separate intervals for all 96 truth/noise × inference-arm cells, and paired host-weight contrasts calculated on shared observations.

## Interpretation

1. **Prior support strongly affects the apparent error.** Under broad uniform support, correct weights yield RMSE 37.88 and 68% coverage 59.7%; reversed weights yield RMSE 60.52 and 68% coverage 30.2%. The narrow support looks much better on this truth grid partly because it is more restrictive; this does not establish that [50, 90] is a generally appropriate astrophysical prior.
2. **Prior measure changes the result.** Under broad support with correct weights, log-uniform H₀ reduces pooled bias from +14.28 to +1.50 and RMSE from 37.88 to 28.94. That is a result of this specific mixture and design, not evidence that log-uniform priors are universally superior.
3. **Host-weight misspecification can be consequential.** Under broad support, reversed weights reduce 68% coverage to 30.2% for the uniform prior and 52.4% for the log-uniform prior. The narrow support masks some alternate-host explanations and can make the weighting error appear smaller.
4. **Pooled metrics are not enough.** Truth/noise-specific reports show heterogeneous behavior across the design. Results should be presented with the cell-level metrics and their Monte Carlo intervals, not only a single pooled score.

## Scope, limitations, and next gate

This is a synthetic low-redshift two-host Gaussian-distance experiment. It does not implement a realistic gravitational-wave event likelihood, a calibrated detector-network selection function, redshift uncertainty, a survey-calibrated galaxy-catalogue completeness model, or hierarchical multi-event inference. The broad [20, 220] prior is a diagnostic arm, not a recommended astrophysical prior. Neither a passing audit nor reasonable PIT/coverage behavior in a toy simulation establishes validity for real LVK data.

The full test suite passed on the current branch revision across Python 3.10, 3.11, and 3.12 in [workflow run 40337117199](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40337117199). Numerical diagnostics also passed across all three Python versions in [run 40337105606](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40337105606). The next scientific gate is to extend the joint catalogue-incompleteness × event-selection model only after its matched generator/likelihood controls and independent numerical checks pass. Do not label a full joint model implemented or validated until those tests and experiments actually exist.
