# Joint posterior-grid resolution audit (24 October 2026)

## Status

The GitHub Actions workflow completed successfully: [run 40117473534](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117473534). The run executed the focused test selection successfully (12 passed) and generated the numerical sensitivity artifact.

This is a **synthetic numerical sensitivity check**, not an observational cosmology result.

## Design

- Seed: `20261025`
- 27 host/redshift/weight configurations
- 30 replicates per configuration (810 simulated observations)
- Four inference arms: `lcdm_two_host`, `lowz_two_host`, `lcdm_collapsed`, and `lowz_collapsed`
- Uniform prior support: [40, 100] km s⁻¹ Mpc⁻¹
- Compared H0 grids: 2,001, 4,001, and 8,001 points
- Reference: the 8,001-point grid

The same simulated observations were reused across grid resolutions, isolating discretization sensitivity from simulation noise.

## Verified output counts

- Posterior summary rows: 9,720 (= 810 observations × 4 arms × 3 grid resolutions)
- Grid-difference rows: 9,720
- All reported posterior medians were finite.
- The workflow's row-count assertions and artifact upload succeeded.

Artifact: [joint-grid-resolution-audit (ID 12570341647)](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40117473534/artifacts/12570341647)

## Maximum absolute differences from the 8,001-point reference

Across all observations, arms, and configurations:

| Posterior quantity | Maximum absolute difference |
|---|---:|
| Median H0 | 0.00894 km s⁻¹ Mpc⁻¹ |
| Lower endpoint of 68% interval | 0.00392 km s⁻¹ Mpc⁻¹ |
| Upper endpoint of 68% interval | 0.00744 km s⁻¹ Mpc⁻¹ |
| Lower endpoint of 95% interval | 0.00130 km s⁻¹ Mpc⁻¹ |
| Upper endpoint of 95% interval | 0.00879 km s⁻¹ Mpc⁻¹ |

These maxima are small relative to the 60 km s⁻¹ Mpc⁻¹ prior width and support the practical numerical adequacy of the 4,001-point grid **for this tested synthetic setup and these posterior summaries**. This does not establish convergence for every possible dataset, prior, likelihood, or model; it is not a substitute for checking more difficult or narrower posteriors.

## Interpretation and limitations

The audit tests discretization error in posterior summaries. It does not validate the scientific assumptions of the host-galaxy model, the Gaussian distance likelihood, the redshift-distance relation, catalogue completeness, or detection-selection treatment. Those assumptions require separate validation. The run uses 30 replicates per configuration, appropriate for this numerical sensitivity diagnostic but not intended to estimate final scientific coverage precisely.

## Reproducibility

The workflow runs:

```bash
python -m pytest -q tests/test_joint_grid_resolution_audit.py tests/test_joint_distance_model_host_mixture.py tests/test_cosmological_inference.py
python scripts/run_joint_grid_resolution_audit.py --replicates 30 --seed 20261025 --output-dir results/grid-resolution-audit
```

The artifact contains the manifest and CSV outputs for further inspection.
