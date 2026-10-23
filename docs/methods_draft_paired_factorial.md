# Methods draft: paired synthetic factorial experiment

## Study objective

This experiment evaluates how specified galaxy-catalogue selection mechanisms and host-galaxy weighting choices affect single-event synthetic inference of the Hubble constant, $H_0$. It is a controlled simulation study, not a reanalysis of gravitational-wave observations and not a measurement of the real Hubble constant.

## Simulation implementation

The experiment imports the simulator from the pinned source revision `b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f` of `vedantkushwah033-blip/bayesian-dark-siren-cosmology`. Each replicate creates a synthetic parent catalogue of 500 galaxies with the simulator's catalogue generator, fixing the data-generating value to $H_0=70$ km s$^{-1}$ Mpc$^{-1}$. The posterior is evaluated on a 1,601-point grid spanning $H_0=50$ to $90$ km s$^{-1}$ Mpc$^{-1}$.

For each simulated parent universe, the experiment generates a host under a uniform or luminosity-proxy weighting rule and generates one mock gravitational-wave distance observation for that host. Within a host-generation setting, that host and distance observation are reused across the relevant catalogue-selection scenarios. Catalogue masks use scenario-specific seeds. For weighting comparisons, the same event and mask are used while the inference weighting is changed.

## Factorial scenarios

Nine prespecified scenarios combine:
1. Catalogue state: complete, random 70% target completeness, or redshift-dependent 70% target completeness.
2. Host generation: uniform or luminosity-proxy weighting.
3. Inference weighting: uniform or luminosity-proxy weighting.

The uniform-host-generation block includes complete, random-incomplete, and redshift-dependent-incomplete catalogues, all inferred with uniform weights. The luminosity-host-generation block includes the same three catalogue states, each inferred with either matched luminosity-proxy weights or mismatched uniform weights.

The redshift-dependent mask is the simulator's specified toy selection rule; it is not a calibrated observational survey selection function. Likewise, luminosity is a synthetic proxy rather than an observed galaxy property or an empirically calibrated galaxy luminosity function.

## Outcomes

For each replicate, the analysis records the posterior median and mean, posterior standard deviation, central 68% and 95% credible intervals, interval widths, signed posterior error relative to the known generating value, realized catalogue completeness, and whether the true host remains in the observed catalogue.

Scenario summaries include bias summaries, RMSE, empirical 68% and 95% interval coverage, Monte Carlo standard errors for coverage, mean interval widths, pull summaries, and host-observation fraction. Paired contrasts are computed by universe identifier for prespecified scenario pairs. The contrast set covers seven comparison families and seven metrics per family: posterior-median error, posterior-mean error, 68% and 95% interval widths, 68% and 95% coverage indicators, and host visibility. The analysis reports mean and median paired differences and percentile-bootstrap 95% intervals from 5,000 bootstrap resamples.

## Pairing and uncertainty

Pairing improves comparability by holding the parent catalogue and, within a host-generation block, the host and mock distance observation fixed across relevant catalogue scenarios. Mask seeds are scenario-specific. For matched-versus-mismatched inference-weight contrasts, the same mask and event are used. The bootstrap intervals quantify Monte Carlo uncertainty in paired contrasts conditional on the simulator and scenario definitions; they do not quantify uncertainty about model realism.

Coverage is evaluated against the fixed synthetic truth $H_0=70$. A nominal Bayesian credible interval does not necessarily achieve the same frequentist coverage at each fixed parameter value. Coverage estimates must be interpreted together with bias, RMSE, and interval width.

## Limitations and scope of inference

The study inherits the pinned simulator's simplified low-redshift distance relation, Gaussian mock distance likelihood, synthetic galaxy properties, and simplified single-event inference. It does not model the full gravitational-wave detection selection process or a survey-calibrated galaxy catalogue selection function. Results should therefore be stated as properties of this specified simulation design only.

Before making observational claims, future work would need to validate a more realistic cosmological distance-redshift model, specify a justified host-population model, incorporate relevant detection and catalogue selection effects, and test against suitable public gravitational-wave and galaxy-catalogue data.

## Reproducibility

The main experiment is defined in `scripts/run_paired_factorial_sensitivity.py`; its GitHub Actions workflow pins the original simulator commit, fixes the base seed, validates the scenario and contrast counts, and stores per-run records and summaries as an artifact. Diagnostic figures are generated by `scripts/plot_paired_factorial_results.py`.

**Status note:** This is a methods draft for review. Final sample size, results, figures, and conclusions must be updated from the successful 300-replicate run before manuscript submission.
