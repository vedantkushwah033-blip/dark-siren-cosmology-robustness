# Noise-scale likelihood consistency audit (23 October 2026)

**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/39974976415  
**Artifact:** `noise-scale-consistency-audit` (14,000 posterior fits)  
**Test-suite status:** passed.

## Question

The toy simulator generates distance noise with standard deviation proportional to the true distance, `sigma_true = 0.15 * d_true`. Does treating the reported uncertainty as fixed across the candidate H0 grid materially change the single-host inference compared with a likelihood where the standard deviation depends on predicted distance?

This is a controlled single-host diagnostic. It does not reproduce the unknown-host galaxy-catalogue likelihood and is not an LVK analysis.

## Design

Seven injected values of H0 (55, 60, 65, 70, 75, 80, 85 km/s/Mpc), 1,000 replicates each, same observation for both modes, two likelihood specifications. Total: 14,000 fits.

- `fixed_reported_sigma`: reported measurement uncertainty is held constant over H0.
- `sigma_depends_on_h0`: sigma(H0) = 0.15 * c*z/H0, with the Gaussian `-log(sigma)` normalization term included.

The posterior grid spans 30–120 km/s/Mpc. The synthetic distance is generated from the low-redshift relation d = c*z/H0, with Gaussian noise whose standard deviation is 15% of the true distance.

## Results

Coverage is empirical across 1,000 replicates per cell. Approximate MC standard errors near 68% coverage are 1.5 percentage points and near 95% coverage about 0.7 percentage points.

| H0 truth | Likelihood | Mean signed median error | RMSE | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---:|---|---:|---:|---:|---:|---:|---:|
| 55 | Fixed reported sigma | +4.45 | 10.97 | 65.7% | 92.9% | 20.28 | 41.61 |
| 55 | H0-dependent sigma | +2.74 | 9.55 | 67.9% | 94.6% | 16.60 | 32.56 |
| 60 | Fixed reported sigma | +4.74 | 11.93 | 63.6% | 92.5% | 21.58 | 43.39 |
| 60 | H0-dependent sigma | +3.05 | 10.79 | 65.9% | 92.7% | 18.12 | 35.58 |
| 65 | Fixed reported sigma | +4.25 | 11.64 | 68.4% | 94.0% | 22.37 | 44.45 |
| 65 | H0-dependent sigma | +2.68 | 10.89 | 69.4% | 93.9% | 19.41 | 38.06 |
| 70 | Fixed reported sigma | +4.53 | 11.95 | 67.2% | 94.9% | 23.42 | 45.50 |
| 70 | H0-dependent sigma | +3.07 | 11.44 | 68.0% | 94.9% | 20.84 | 40.76 |
| 75 | Fixed reported sigma | +3.95 | 11.92 | 66.2% | 95.0% | 23.75 | 45.31 |
| 75 | H0-dependent sigma | +2.74 | 11.81 | 66.9% | 95.1% | 21.93 | 42.66 |
| 80 | Fixed reported sigma | +3.72 | 11.93 | 66.9% | 95.1% | 23.98 | 44.75 |
| 80 | H0-dependent sigma | +2.77 | 12.10 | 67.1% | 94.6% | 22.90 | 44.21 |
| 85 | Fixed reported sigma | +3.00 | 11.46 | 68.9% | 95.5% | 23.84 | 43.69 |
| 85 | H0-dependent sigma | +2.30 | 11.81 | 68.4% | 96.1% | 23.59 | 45.09 |

RMSE and widths are in km/s/Mpc.

## Interpretation

1. **The likelihood specification matters, but it does not explain all observed problems.** The H0-dependent-noise model reduced mean signed median error in all seven truth cells. The reduction ranged from about 0.7 to 1.7 km/s/Mpc.
2. **Coverage is broadly near nominal for this single-host diagnostic.** The H0-dependent-noise model's 68% coverage ranges from 65.9% to 69.4%; 95% coverage ranges from 92.7% to 96.1%. These results are compatible with ordinary Monte Carlo variation and finite-grid/model effects; they do not prove exact calibration.
3. **Bias remains.** Mean signed median errors remain positive in all cells. A likely contributor in this toy setup is the nonlinear inverse relationship between distance and H0, combined with noisy distance measurements and the posterior summary; this is a plausible interpretation, not an independently isolated causal result.
4. **The full unknown-host catalogue model remains the key unresolved component.** This single-host experiment does not marginalize over possible hosts and therefore cannot explain the earlier 56% coverage at H0=70 in the wide-support unknown-host experiment.
5. **Do not infer real-world performance.** These experiments use a simplified low-redshift relation, a single fixed redshift in this audit, Gaussian noise, and a flat bounded H0 prior. They are not end-to-end analyses of LVK events or observed galaxy catalogues.

## Next experiment

Build a minimal matched unknown-host test with a deliberately transparent generative model: draw a finite catalogue, choose a host using a specified rule, generate distance data, and perform inference by summing over exactly that catalogue with the same weights and noise model. Verify the model in a small hand-checkable case before scaling the number of replicates. Then use posterior-rank/SBC diagnostics under matched prior-predictive simulations.

## Reproducibility

The successful GitHub Actions run validates 14 summary rows and 14,000 per-run records. The uploaded artifact contains `summary.csv` and `per_run.csv`.