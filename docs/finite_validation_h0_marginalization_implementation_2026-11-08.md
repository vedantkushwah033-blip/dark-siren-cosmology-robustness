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


## Relationship to the vectorized H0 implementation

The primary implementation is now `finite_validation_h0.py`, which reconstructs density draws, evaluates event mixtures and per-draw selection efficiencies in H0 batches, performs weighted marginalization, and normalizes a declared-prior H0 posterior. This module remains a deliberately slower loop-based reference path: it calls the established `selection_aware_toy_likelihood` once per coefficient proposal. A cross-implementation regression test compares its unconditional likelihoods, beta curves, and selected likelihoods against the vectorized implementation using the same host-density draws.

The reference path is useful for small independent numerical checks, not for large proposal counts. Its one-gigabyte guard estimates returned draw-array storage; it is not a strict peak-memory bound. The vectorized implementation has separate output/temporary allocation guards.


## H0-specific likelihood-contribution effective sample size

The coefficient-posterior importance ESS alone is not sufficient for likelihood marginalization. At a given H0, a small subset of coefficient proposals may dominate the integrand even when the proposal weights themselves have a healthy ESS. The propagation output now reports, at every H0 grid point,

\[
a_b(H_0)=\frac{\exp\{\log w_b+\ell_b(H_0)\}}
{\sum_j\exp\{\log w_j+\ell_j(H_0)\}},
\qquad
ESS_L(H_0)=\frac{1}{\sum_b a_b(H_0)^2}.
\]

It also reports the maximum normalized likelihood-contribution weight and the same contribution ESS within each disjoint Monte Carlo batch. This is distinct from (i) the validation-posterior proposal ESS and (ii) the batch ESS of the proposal weights alone. The default diagnostic status is LOW_CONTRIBUTION_ESS if any H0 point has contribution ESS below 5; callers may set a stricter threshold before a confirmatory run. A high coefficient-posterior ESS must not override a low H0-specific contribution ESS.

The memory guard now accounts for the additional full draw-by-H0 log-contribution and normalized-contribution arrays, as well as the returned batch diagnostics. This remains a conservative allocation estimate, not a strict process peak-memory guarantee.


---

## Updated implementation checkpoint — 2026-11-08

The earlier status statement is retained as a snapshot from the first implementation stage. Since then:

- The vectorized primary path in `finite_validation_h0.py` has been connected to the weighted fixed-basis Bayesian inclusion importance sampler.
- The returned H0 posterior now includes declared-prior normalization, 68%/95% summaries, and boundary-mass diagnostics.
- The reference path in this module is tested against the vectorized per-draw unconditional likelihood, GW detection efficiency, and selected likelihood.
- The integration suite includes a real importance-sampler posterior flowing through host-density reconstruction and selected-likelihood marginalization, plus a parent-frame-to-normalized-H0-posterior smoke test.
- Marginalization diagnostics now include coefficient-proposal ESS by batch, H0-specific likelihood-contribution ESS, maximum contribution weight, batch-level contribution ESS, and log-likelihood MCSE. Low contribution ESS is separately flagged; a high proposal ESS does not override it.

The focused sampler workflow and full numerical-diagnostics matrix passed on Python 3.10, 3.11, and 3.12 at commit `dae374c17a462b65aec94c06ca05f730ac094e37` (runs `42482040468` and `42482040471`).

These are software and numerical cross-checks, not repeated-sampling calibration. The production/reference paths share SciPy's KDE primitive, and the reference path's returned-array memory guard is not a strict process peak-memory guarantee. The repeated-sampling experiment and independent H0 artifact auditor remain future work; see `docs/finite_validation_h0_repeated_sampling_protocol_2026-11-08.md`.
