# Posterior normalization audit: scope and artifact schema

## Current implementation

Artifact schema version 2 stores, for every replicate and estimator arm:

- the explicit H0 grid as a JSON array;
- the corresponding normalized posterior density vector as a JSON array;
- the density vector's SHA-256 checksum;
- the posterior mean, median, 16th/84th percentiles, interval width, and recorded integral.

The summary metadata declares schema version 2 and the JSON encoding. The independent auditor verifies the density checksum, finite/nonnegative values, grid shape and strict monotonicity, agreement with the declared grid, and the independently recomputed trapezoidal integral. It then recomputes the posterior mean and quantiles from the stored density and compares them with the recorded summaries. Aggregate metrics are reconstructed from replicate rows.

## Interpretation and limitations

The audit now independently reconstructs the posterior quantities represented by the saved grid and density vector. It does **not** independently reconstruct the likelihood from the event-level inputs or prove that the producer's likelihood model is scientifically correct. The density is an artifact produced by the benchmark; the auditor checks its internal numerical consistency and consistency with its reported summaries.

The checksum detects accidental or uncoordinated modification of the density field. It is not a digital signature and does not establish provenance against an actor able to rewrite both data and checksum.

## Adversarial tests

Tests must alter the density while leaving stored summaries unchanged, update the density checksum to avoid a trivial checksum-only failure, and verify that independently recomputed posterior summaries disagree. Separate tests cover missing schema metadata, interval tampering, aggregate tampering, and missing paired arms.

## Scientific scope

The benchmark remains a synthetic exploratory toy study. It is not a real H0 measurement, does not use LVK event data, and does not validate a real detector or event-selection pipeline.
