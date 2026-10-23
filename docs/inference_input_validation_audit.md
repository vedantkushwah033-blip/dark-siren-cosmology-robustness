# Inference input-validation audit of the pinned simulator

## Scope and reproducibility

This audit inspects, but does not modify, the original inference implementation pinned by the robustness experiment:
`vedantkushwah033-blip/bayesian-dark-siren-cosmology`, commit
`b9277e7bf9ad5e053c40c1c671b70aa7f3d0708f`, file `src/bayesian_inference.py`.

The findings below are code-level risks inferred from the function's validation logic. They are not claims that the default experiment necessarily triggers these cases. The pinned source should remain unchanged for reproducibility; any hardened implementation should be versioned separately and compared against the pinned implementation.

## Findings

### 1. H0 grid validation is incomplete

`posterior_h0` converts `H0_grid` to a float array but does not explicitly require that it be one-dimensional, finite, contain at least two points, or be strictly increasing. Trapezoidal normalization assumes a sensible ordered grid. A malformed grid can produce invalid or misleading integration results.

**Required guard:** one-dimensional, at least two finite points, strictly increasing; also validate that prior bounds are finite and `prior_min < prior_max`.

### 2. Distance uncertainty and GW distance are not validated explicitly

The function passes `gw_distance` and `gw_sigma` to `scipy.stats.norm.pdf` without explicit checks for finite values and strictly positive uncertainty. Invalid uncertainty can generate warnings, non-finite likelihood values, or a late normalization failure rather than a clear input error.

**Required guard:** finite scalar observed distance and finite positive `gw_sigma`.

### 3. Redshift and weight validation can miss non-finite values

The function checks catalogue weights for negative entries and non-positive sums, but those comparisons alone do not reliably reject NaN or infinity. Redshifts are converted to arrays but are not explicitly checked for finite, non-negative values. Similar checks are needed for missing-host support and weights.

**Required guard:** finite, non-negative redshifts and weights; compatible array shapes; strictly positive total weight for every active component.

### 4. Posterior-summary inputs are assumed valid

`posterior_summary` computes moments and a cumulative distribution without first checking that the grid is one-dimensional and strictly increasing, the posterior has matching shape, values are finite and non-negative, and the total area is positive. Its CDF normalization divides by the final cumulative value.

**Required guard:** validate grid and posterior shapes and values, verify positive finite integrated area, then normalize before computing summaries.

## Recommended hardening plan

1. Add tests for malformed grids, NaN/infinite redshifts and weights, non-positive uncertainty, invalid prior bounds, and mismatched shapes.
2. Implement explicit validation in a separately versioned hardened inference module rather than silently changing the pinned baseline.
3. Compare the hardened implementation with the pinned baseline on valid default inputs; verify posterior and summary agreement within a declared tolerance.
4. Run the posterior-grid stability audit on fixed simulated observations.
5. Only then rerun the full synthetic factorial experiment and report failures/exclusions explicitly.

## Scientific interpretation

Input validation improves failure handling; it does not prove posterior calibration or model validity. The current experiment remains synthetic, not an observational LVK result. The main factorial output must be generated and independently checked before drawing conclusions about catalogue incompleteness or host-galaxy weighting.
