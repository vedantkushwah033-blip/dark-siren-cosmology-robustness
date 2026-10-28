# Synthetic population generator with separated selection mechanisms

**Status:** implemented simulation primitive; tests added; CI validation pending.  
**Scope:** synthetic data generation only. This is not real survey data, a realistic detector model, or an H0 inference result.

## Why this stage exists

Catalogue incompleteness and gravitational-wave detection are different selection mechanisms. This generator represents them separately and stores hidden simulation truth apart from observable tables.

## Generative process

1. Draw a fixed number of parent galaxies from a positive redshift intensity: a constant baseline plus two Gaussian overdensity features.
2. Retain each parent galaxy independently with a smooth redshift-dependent catalogue inclusion probability.
3. Draw GW host galaxies from the complete parent population using a declared increasing toy merger-rate weight. Hosts are not sampled from the thinned observed catalogue.
4. Generate a distance-like measurement using the low-redshift relation d_obs ~ Normal(c*z/H0_true, sigma_mpc^2).
5. Apply a separate GW detection rule: d_obs >= detection_threshold_mpc. Rejection sampling produces the requested number of detected events.

## Output contract

- catalogue_observed.csv: observable catalogue only (object_id, redshift).
- gw_events_observed.csv: observable events only (event_id, distance_mpc, detected).
- parent_truth.csv: truth-only parent population, inclusion probabilities, and retention flags.
- gw_event_truth.csv: truth-only host IDs, host redshifts, and latent mean distances.
- config.json: seed, units, distributions, selection definitions, parameters, and limitations.

Truth tables are for simulation scoring and oracle diagnostics only. They must not be passed to non-oracle inference. The observable event table deliberately excludes true host redshift and host object ID.

## Correctness gates and limitations

Tests cover intensity positivity and clustering, monotone catalogue incompleteness, reproducibility, observable/truth schema separation, detection-threshold compliance, host-ID consistency, and invalid parameters.

This stage does not yet implement an H0 posterior, a catalogue-selection likelihood, uncertainty propagation for the latent galaxy intensity, or a homogeneous/GP comparison. The distance relation is a low-redshift Euclidean approximation with fixed true H0 and fixed Gaussian noise. These limitations must remain explicit.
