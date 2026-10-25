# Posterior-level H₀ calibration: protocol

## Objective
Measure posterior-level bias, RMSE, interval width, and truth coverage in a synthetic two-host low-redshift mixture. This extends component-mass diagnostics to cosmological posterior summaries. It is a simulation study, not a real gravitational-wave result.

## Prespecified design
- Injected H₀ values: 55, 65, 75, 85 km s⁻¹ Mpc⁻¹.
- Candidate host redshifts: z₁ = 0.02 and z₂ = 0.04.
- Generative host weights: (0.7, 0.3); inference weights: correct (0.7, 0.3) and reversed (0.3, 0.7).
- Distance-noise fractions: 3%, 10%, 20%, with sigma defined as noise_fraction × c × z₂ / H₀,true. This is a fixed reference scale for each truth/noise cell; generated distances use the selected host's luminosity-distance proxy c z_host/H₀,true.
- Prior support: narrow [50,90] and broad [20,220] km s⁻¹ Mpc⁻¹.
- Prior measure: uniform in H₀ or log-uniform in H₀ (density proportional to 1/H₀).
- Posterior grid: 10–250 km s⁻¹ Mpc⁻¹, 2,401 equally spaced points.
- Default: 200 replicates per truth/noise cell; seed 20261030; 2,400 raw observations and 96 summary cells.

## Paired design
Each synthetic observation and its sigma is evaluated under all 2 × 2 × 2 inference combinations: two supports, two prior measures, and two catalogue-weight choices. This enables paired comparison without changing the generated data between analysis arms.

## Posterior and estimands
For each observation, the unnormalized posterior is
$$
p(H_0\mid d)\propto\left[\sum_{i=1}^2 w_i\,p(d\mid H_0,z_i)\right]\pi(H_0),
$$
truncated to the prespecified support and normalized numerically on the grid. Report median error (posterior median minus injected truth), median absolute error, RMSE, mean central 68% and 95% interval widths, and empirical truth coverage at 68% and 95%. Coverage is estimated separately for every truth/noise/support/measure/weight cell.

## Validation
The independent auditor reconstructs each posterior median and interval from raw distances and the fixed protocol; checks truth, host identity, sigma, unique complete replicate blocks and design counts; reconstructs all aggregate calibration metrics; and verifies SHA-256 hashes. The run is marked PASS only after these checks.

## Interpretation cautions
- 200 replicates provide useful pilot resolution, not high-precision coverage estimates. For nominal 68% coverage, the binomial standard error is roughly 3.3 percentage points at n=200; for nominal 95%, roughly 1.5 points.
- Truth coverage under a narrow support is only interpretable for truths inside that support; the chosen truths are inside [50,90].
- This two-host Gaussian-distance toy omits catalogue selection/incompleteness, redshift uncertainty, detector selection, realistic event likelihoods, and population-level multi-event inference.
- A prior/weight contrast in this simulation is not by itself evidence that a method is calibrated for real observations. Report support sensitivity and prior-measure sensitivity separately.
