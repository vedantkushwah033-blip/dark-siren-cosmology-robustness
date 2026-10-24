# Two-host mixture simulation-based calibration audit

**Status:** completed; 20,000 posterior fits validated by GitHub Actions.  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40065075227  
**Scope:** synthetic low-redshift two-host model only. This is a computational calibration check, not a measurement of the real Hubble constant and not validation against real gravitational-wave events or galaxy catalogues.

## Question

Does the weighted two-host Gaussian-mixture likelihood recover calibrated posterior ranks when the data are generated from the same model used for inference?

## Design

- Draw injected $H_0$ uniformly from the inference prior $[50,90]$ km s$^{-1}$ Mpc$^{-1}$.
- Two candidate hosts have fixed redshifts $z_1=0.02$ and $z_2=0.04$.
- Generate one latent host using the scenario's catalogue weights, then generate one observed luminosity-distance proxy from a Gaussian with fixed standard deviation.
- Infer $H_0$ using the matching weighted mixture likelihood and the same prior.
- Run 5,000 prior-predictive replicates for each of two weight settings: $(0.7,0.3)$ and $(0.5,0.5)$.
- Repeat posterior CDF-rank calculation on grids of 1,001 and 4,001 points. Total: 2 scenarios × 5,000 replicates × 2 grid resolutions = 20,000 posterior fits.
- Assess ranks using the Kolmogorov–Smirnov (KS) test against Uniform(0,1), a ten-decile chi-square diagnostic, and grid-resolution differences.

The observation uncertainty is fixed and shared by the generator and likelihood. This avoids the earlier generative/likelihood mismatch found in an initial single-host calibration attempt.

## Results

| Host weights | H0 grid points | Replicates | Mean rank | Rank variance | KS statistic | KS p-value | Decile-test p-value |
|---|---:|---:|---:|---:|---:|---:|---:|
| 70/30 | 1,001 | 5,000 | 0.50395 | 0.08377 | 0.00976 | 0.7243 | 0.9843 |
| 50/50 | 1,001 | 5,000 | 0.50161 | 0.08424 | 0.00797 | 0.9059 | 0.9461 |
| 70/30 | 4,001 | 5,000 | 0.50396 | 0.08408 | 0.01073 | 0.6083 | 0.9449 |
| 50/50 | 4,001 | 5,000 | 0.50163 | 0.08457 | 0.00790 | 0.9114 | 0.8826 |

Neither diagnostic rejects uniform rank calibration at the 5% level for these two scenarios and grid resolutions. The rank means are close to 0.5. This is evidence consistent with calibration under this specific matched toy model, not proof that every implementation detail or the astrophysical model is correct.

### Grid-resolution sensitivity

| Host weights | Median absolute rank difference (1,001 vs 4,001) | 95th percentile | Maximum |
|---|---:|---:|---:|
| 70/30 | 0.000353 | 0.002778 | 0.020204 |
| 50/50 | 0.000305 | 0.002613 | 0.032255 |

Typical differences are small, while a few ranks differ more noticeably. A finer grid is not identical to a continuous exact posterior; tail ranks and sharp/multimodal shapes warrant further numerical checks before using this as a high-precision inference engine.

## An initially concerning result and why the larger run matters

An earlier 1,000-replicate run using seed 20261026 returned KS p-values of 0.023 and 0.022 for the equal-weight scenario on the two grids, while the 70/30 scenario did not reject. The 5,000-replicate rerun with the same seed yielded KS p-values 0.906 and 0.911 for equal weights and decile-test p-values 0.946 and 0.883. The earlier signal did not persist with the larger sample and is consistent with a finite-sample fluctuation; it is retained here rather than silently omitted.

## Interpretation and limitations

1. The data generator and inference likelihood use the same redshifts, host weights, fixed Gaussian uncertainty, and low-redshift relation $d_L \approx cz/H_0$.
2. This audit checks whether posterior CDF ranks behave as expected in this particular prior-predictive setup. It does not establish real-data validity, accurate catalogue selection, or realistic astrophysical host probabilities.
3. The two hosts are fixed and the model omits a realistic detector selection function, cosmological distance-redshift integration, measurement-systematic uncertainties, and real galaxy-catalogue structure.
4. The decile chi-square test is a diagnostic, not an independent proof; the KS and decile tests are based on the same ranks.
5. The next model-validation step is to test deliberately mismatched host weights and omitted-host support, and to add deterministic unit tests showing that swapping or perturbing weights changes the likelihood in the expected direction.

## Reproducibility

- Script: `scripts/run_two_host_sbc_audit.py`
- Tests: `tests/test_two_host_sbc_audit.py`
- Workflow: `.github/workflows/two-host-sbc-audit.yml`
- Artifact: `two-host-sbc-audit-results` (contains per-run ranks and CSV summaries; GitHub Actions retention is time-limited).
