# Posterior grid-refinement audit

**Status:** experiment completed successfully: 21,600 posterior fits across 36 design cells and three grid resolutions.  
**Workflow:** https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40066736877  
**Scope:** numerical sensitivity audit for the synthetic low-redshift two-host model.

## Design

For injected $H_0=52,70,88$ km s$^{-1}$ Mpc$^{-1}$, second-host redshift $z_2=0.025,0.04,0.06$, noise fractions 3% and 10%, and correct/reversed host weights, 200 synthetic observations were generated per design cell. Each same observation was inferred on $H_0$ grids of 1,001, 2,001 and 4,001 points, spanning 25 to 125. This gives 36 design cells × 200 observations × 2 weightings × 3 grids = 21,600 posterior fits.

We compared the posterior median and the endpoints of the 68% and 95% equal-tail credible intervals. Differences are absolute changes in km s$^{-1}$ Mpc$^{-1}$ between grid resolutions.

## Results

| Quantity compared: 1,001 vs 4,001 grid points | Median of cell-level median absolute differences | Largest cell-level 95th percentile | Largest observed difference |
|---|---:|---:|---:|
| Posterior median | 0.0123 | 0.1324 | 0.3462 |
| Lower 68% endpoint | 0.0102 | 0.1463 | 0.2820 |
| Upper 68% endpoint | 0.0150 | 0.0813 | 0.7972 |
| Lower 95% endpoint | 0.0065 | 0.0483 | 0.5871 |
| Upper 95% endpoint | 0.0305 | 0.0587 | 1.4726 |

For posterior medians, the median of the cell-level median differences between 2,001 and 4,001 points was 0.00411; the largest cell-level 95th percentile was 0.04554, and the largest observed difference was 0.12652.

## Interpretation

- Typical posterior-median and interval-endpoint changes are small when moving from 1,001 to 4,001 points.
- Some individual cases show materially larger endpoint changes, especially the upper 95% endpoint (maximum observed difference about 1.47 km s$^{-1}$ Mpc$^{-1}$). The coarser grid should not be assumed accurate for every posterior, particularly for tail quantiles and complex mixture shapes.
- Refinement from 2,001 to 4,001 points reduced typical median variation, but the maximum differences show that a single grid-size rule is not sufficient as a universal guarantee.
- A practical next step is adaptive-grid or interpolation validation against a high-resolution reference, with explicit tolerances for posterior quantiles and normalization. Coverage estimates should then be recomputed with the validated numerical method.

## Limitations

This only checks numerical grid sensitivity for selected synthetic configurations. It does not establish astrophysical realism, detector-selection validity, or calibration against real gravitational-wave events. The experiment uses a low-redshift distance approximation, two candidate hosts, a bounded uniform prior, and Gaussian noise.

## Reproducibility

- Script: `scripts/run_grid_refinement_audit.py`
- Tests: `tests/test_grid_refinement_audit.py`
- Workflow: `.github/workflows/grid-refinement-audit.yml`
- Artifact: `grid-refinement-audit-results` (per-run and cell-level CSVs; retention is time-limited).
