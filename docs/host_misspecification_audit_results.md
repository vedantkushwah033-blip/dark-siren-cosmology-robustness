# Host-weight and omitted-host misspecification audit

**Status:** completed; GitHub Actions passed 4,000 synthetic fits (1,000 replicates × 4 scenarios).  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40066102259  
**Artifact:** `host-misspecification-audit-results` (CSV outputs, subject to Actions artifact retention).  
**Scope:** controlled low-redshift synthetic stress test; not real gravitational-wave or galaxy-catalogue validation.

## Design

The injected value is fixed at $H_0=70$ km s$^{-1}$ Mpc$^{-1}$. Each synthetic observation is generated from a latent host at one of the specified redshifts, with Gaussian distance noise whose fixed standard deviation is shared with the inference likelihood. The likelihood uses the low-redshift approximation $d_L \approx cz/H_0$ and a uniform prior on $H_0 \in [50,90]$.

Four scenarios are compared:

1. **Correct weights:** observed hosts at $z=0.02,0.04$; data-generation weights 70/30 and inference weights 70/30.
2. **Reversed weights:** same data-generation weights 70/30, but inference weights 30/70.
3. **Hidden host ignored:** the data generator includes a third host at $z=0.03$ with 30% weight, but inference omits it and uses equal weights for the two listed hosts.
4. **Hidden host modelled:** the generator includes the same third host and 30% weight; inference includes a missing-host component with 30% weight and uses equal relative weights for the two listed hosts.

Each scenario has 1,000 replicates. The run reports posterior-median error, RMSE, empirical 68% and 95% interval coverage, and average interval widths.

## Results

| Scenario | Mean median error | RMSE of median | 68% coverage | 95% coverage | Mean 68% width | Mean 95% width |
|---|---:|---:|---:|---:|---:|---:|
| Correct weights | +1.258 | 6.651 | 73.0% | 98.7% | 15.366 | 27.849 |
| Reversed weights | +1.253 | 6.656 | 74.3% | 99.0% | 15.456 | 28.138 |
| Hidden host ignored | +1.285 | 8.565 | 66.4% | 99.3% | 17.003 | 29.255 |
| Hidden host modelled | +1.529 | 4.721 | 88.5% | 98.8% | 22.187 | 33.066 |

All distance and $H_0$ quantities are in the model's corresponding units; interval widths are in km s$^{-1}$ Mpc$^{-1}$.

## Interpretation

### 1. Reversing the two observed-host weights barely changes this experiment

The correct-weight and reversed-weight summaries are very similar. This does **not** demonstrate that host weights are generally unimportant. Under the current noise level, the two-host likelihood components overlap substantially, so this experiment has limited power to distinguish the two weighting choices. A follow-up should vary measurement precision and host redshift separation, and use paired datasets so each weighting choice is evaluated on identical observations.

### 2. Omitting the third host worsens median RMSE in this setup

When the data generator assigns 30% of observations to a host at $z=0.03$ but inference excludes that host, median RMSE increases from 6.65 to 8.57 compared with the correct two-host model. The empirical 68% coverage falls to 66.4%. This is consistent with the omitted-host support affecting inference in this particular synthetic setup.

### 3. Modelling the hidden-host component reduces RMSE but widens intervals

Including the third-host component lowers median RMSE to 4.72, while the average 68% interval widens to 22.19 and empirical coverage rises to 88.5%. That is not an unqualified win: the model becomes more conservative and produces substantially wider intervals. The 88.5% coverage is above the nominal 68% target, suggesting the intervals are over-covering in this scenario.

### 4. The 95% intervals over-cover in all four scenarios

Observed 95% coverage ranges from 98.7% to 99.3%, substantially above the nominal 95% level. The current grid-based posterior intervals and the simple distance-mixture model therefore need additional calibration work before this model is treated as a reliable uncertainty-quantification engine.

## Limitations and next steps

- These are synthetic results at one fixed $H_0$, two observed redshifts plus one deliberately hidden redshift, and one noise scale. They should not be generalized to real catalogues or LVK analyses.
- The experiment is a stress test, not a full prior-predictive simulation-based calibration analysis.
- The 68% and 95% coverage estimates each use 1,000 replicates; approximate binomial Monte Carlo standard errors are about 1.5 percentage points near 68% coverage and 0.7 percentage points near 99% coverage.
- Weight-misspecification effects should next be tested across a grid of redshift separations and noise levels using paired observations.
- Interval coverage should be audited across several injected $H_0$ values, especially near prior boundaries, and with finer numerical grids.
- Before using the results in a paper, independently inspect the likelihood, interval construction, and all generated CSV outputs.

## Reproducibility

- Script: `scripts/run_host_misspecification_audit.py`
- Tests: `tests/test_host_misspecification_audit.py`
- Workflow: `.github/workflows/host-misspecification-audit.yml`
