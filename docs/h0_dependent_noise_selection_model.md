# H0-dependent distance noise and selection-conditioned dataset likelihood

## Purpose and scope

This module implements a transparent synthetic model step: heteroscedastic distance uncertainty whose scale depends on candidate H0, plus a fixed-sample selection correction using the same noise law. It is a methodological prototype, not an LVK likelihood, detector-response model, real-data analysis, or estimate of the real Hubble constant.

## Generative model

For a discrete host at redshift z and candidate H0,

`mu(z,H0) = c*z/H0`, `sigma(z,H0) = sqrt(sigma_floor_mpc^2 + (fractional_sigma*mu)^2)`.

Here `fractional_sigma > 0` and `sigma_floor_mpc >= 0` (Mpc). For event i with measured distance d_i, host redshifts z_ij, and normalized discrete host masses w_ij,

`L_i(H0) = sum_j w_ij Normal(d_i | mu(z_ij,H0), sigma(z_ij,H0)^2)`.

An optional positive-truncation setting conditions each host-wise distance density on d>0, dividing by Phi(mu/sigma). The same setting must be applied to the selection efficiency.

For a population represented by discrete redshifts z_k and normalized masses q_k, with detection defined by d >= d_min > 0,

`beta(H0) = sum_k q_k * SF((d_min-mu(z_k,H0))/sigma(z_k,H0))`.

For positive-truncated noise, divide the survival function by P(d>0)=SF(-mu/sigma). For a redshift grid approximating a continuous density, supply quadrature masses (density times quadrature weights), not raw density values.

For a fixed-size set of N detected events,

`logL_selected(H0) = sum_i logL_i(H0) - N*log(beta(H0))`.

This conditions on the number of detections. It does not infer event rates or event counts. The population mixture supplied for beta must actually represent the selection-normalization population; an event-specific candidate-host list is not automatically that population.

## Software interface

The module exposes `event_log_likelihood_h0_dependent_noise`, `population_detection_efficiency_h0_dependent_noise`, `selection_conditioned_dataset_log_likelihood`, and `dataset_posterior_h0_selection_conditioned`. Inputs are validated fail-closed. The posterior helper inserts exact prior endpoints into the normalization quadrature even when they fall between the supplied H0 grid nodes.

## Validation and limitations

Tests include direct Normal and weighted-mixture calculations, the positive-truncation normalizer, analytic detection probabilities, the fixed-N correction identity, endpoint-exact posterior normalization, invalid-input checks, zero-weight equivalence, and seeded Monte Carlo calibration of detection efficiency.

This prototype still omits realistic detector response, sky and inclination marginalization, waveform/reconstruction uncertainty, calibration uncertainty, population-rate inference, realistic cosmological distance-redshift relations, peculiar velocities, redshift uncertainty, and galaxy-survey hierarchical selection. Passing these tests verifies software against the stated toy model; it does not establish physical realism or support a real-data H0 claim. The scientific gate remains NO-GO for those claims until relevant ingredients and benchmarks are independently validated.
