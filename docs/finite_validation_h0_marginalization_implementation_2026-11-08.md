# Finite-validation selected-likelihood marginalization: implementation note

**Status:** computational plumbing implemented; no confirmatory H0 benchmark, posterior calibration, or coverage result is claimed.

## Implemented operation

`finite_validation_h0_marginalization.py` connects the fixed-basis weighted inclusion-posterior approximation to the existing synthetic host-density and selection-aware toy likelihood primitives.

For each coefficient proposal \(\theta_b\):

1. Evaluate \(\pi(z_i;\theta_b)\) for the analysis-catalogue redshifts in the fixed domain basis.
2. Construct one merger-weighted IPW host density using the frozen Scott multiplier.
3. Evaluate the unconditional event-mixture log likelihood over the complete H0 grid.
4. Evaluate \(\beta_{GW}(H_0;\theta_b)\) from that same density and same grid.
5. Form the selected log likelihood \(\ell_b(H_0)=\ell_{uncond,b}(H_0)-N_{det}\log\beta_{GW,b}(H_0)\).
6. Marginalize the selected likelihood with the proposal's **normalized importance weights** in log space.

The same coefficient proposal is used across every H0 point. The selection correction is applied within each proposal before marginalization. The implementation does not average beta values, does not average unconditional likelihoods before selection correction, and does not silently drop failed draws.

## APIs

- `selection_correct_each_draw`: validates and applies the fixed-count GW selection term to arrays of shape (coefficient draw, H0 grid).
- `marginalize_selected_draws`: calls the weighted log-marginal-likelihood implementation using the proposal log weights.
- `propagate_inclusion_posterior_to_h0`: builds each draw's IPW density and selected likelihood, then returns the draw-level arrays, marginal likelihood vector, importance ESS and maximum normalized weight.

## Validation

Focused tests cover:
- draw-specific selection correction versus an incorrect average-beta shortcut;
- equal-weight reduction to ordinary log-mean-exp;
- rejection of invalid efficiencies and malformed/missing weights;
- end-to-end smoke propagation with two coefficient proposals and distinct beta(H0) curves.

The importance sampler's weights are self-normalized. Zero-valued normalized weights caused by floating-point underflow are permitted only when their corresponding log weights remain finite; the log weights, not rounded normalized weights, are used in marginalization.

## Limitations and next gates

- This module returns a marginalized selected log-likelihood vector, not a normalized H0 posterior.
- It does not yet provide posterior summaries, prior application, boundary-mass accounting, raw artifact manifests, or an independent full-artifact auditor.
- The proposal's importance ESS and maximum normalized weight must gate downstream use. Importance Monte Carlo error must be checked by larger draw counts and independent proposal seeds.
- Prior sensitivity (candidate scales 0.75, 1.0, 1.5) remains mandatory.
- Redshift-grid and H0-grid convergence must be established before fixed-truth coverage or SBC.
- The fixed-count toy selection rule is not a real GW detector model.


## Additional resource and support guards

The propagation function now rejects redshift grids that do not span the complete fixed prior domain. It also estimates the memory needed for the returned density, unconditional-likelihood, beta, and selected-likelihood draw arrays before allocation. The default one-gigabyte guard fails early with an actionable message rather than risking an uncontrolled allocation; callers can reduce proposal draws or grid resolution and must still establish Monte Carlo/grid convergence separately.

A regression test covers both the full-domain requirement and the pre-allocation memory guard.
