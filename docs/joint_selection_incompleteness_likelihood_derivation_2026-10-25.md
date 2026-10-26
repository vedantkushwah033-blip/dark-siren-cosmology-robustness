# Joint-model derivation checkpoint: catalogue missingness and event selection

**Status: derivation checkpoint; no combined-model code or numerical results are claimed here.** This document records the likelihood issue that must be resolved before implementing the proposed factorial.

## 1. Separate the two selection processes

Let H denote H0, G the latent host galaxy, D the observed GW data, A=1 the event-detection indicator, and C_g=1 indicate that galaxy g is present in the external catalogue.

- Event detection is governed by P(A=1 | G, H, source/event properties).
- Catalogue inclusion is governed by P(C_g=1 | galaxy properties, survey/selection information).

These are not the same conditioning event. Catalogue missingness should not be represented by the GW detection efficiency, and GW detection correction should not be estimated from the observed catalogue alone unless the required population assumptions justify that calculation.

## 2. Event-level likelihood with a complete host population

For a discrete, known candidate population with host weights pi_g(H) (normalized over the full relevant population), an illustrative detected-event likelihood has the form

L_det(D | H, A=1) = I(D in R) * [sum_g pi_g(H) p(D | g,H)] / alpha(H),

where R is the accepted-data region, I(D in R) is its indicator, and alpha(H) = P(D in R | H) is the population-averaged detection probability under the exact event-generation and detection rule. For a realised detected datum inside R, the indicator equals one. In the repository's toy distance-threshold model, the observed distance itself determines detection, so the numerator and alpha must be derived consistently from that truncated-data process. This schematic expression is not a substitute for deriving the precise likelihood for the simulator.

If selection is defined on noisy observed distance d and the retained sample requires d_min < d < d_max, then for fixed host g with Gaussian distance model,

alpha_g(H) = Phi((d_max - mu_g(H))/sigma) - Phi((d_min - mu_g(H))/sigma).

For a complete known host mixture, alpha(H) = sum_g pi_g alpha_g(H), matching the structure of `selection.py`. The weights used in alpha(H) must correspond to the population that generates detected events, not automatically to whatever subset happens to be listed in a catalogue.

## 3. What changes when the catalogue is incomplete?

Let O be the observed galaxy subset and M the unobserved/missing subset. In a finite-population simulation where the full catalogue is available to the simulator, a diagnostic oracle can compute

L_full(D | H) = sum_{g in O} pi_g p(D | g,H) + sum_{g in M} pi_g p(D | g,H),

with the relevant detection conditioning then applied consistently. The observed-only likelihood drops the second term and renormalizes the remaining observed weights; that is a different model, not merely a lower-precision version of the same likelihood.

A practical inference procedure cannot generally know the identities, redshifts, or true weights of missing galaxies. To make a practical correction, the model must instead integrate over a specified latent missing population, for example

L_catalogue(D | H, O, survey data) = sum_{g in O} pi_g p(D | g,H) + Pi_M(H) * integral p(D | x,H) p(x | missing, survey data,H) dx,

where Pi_M is the modeled prior mass assigned to the missing population and x contains the latent host properties needed by the likelihood. This is only a schematic decomposition: the conditional distribution for missing hosts and the normalization of observed/missing population mass must be derived from a joint population and catalogue-inclusion model.

Crucially, if event detection and catalogue inclusion share latent dependencies, the joint model must represent those dependencies. One must not multiply two marginal correction factors by default.

## 4. Required compatibility tests

Before adding code, write down and test these invariants:

1. **Complete-catalogue limit:** when every galaxy is observed, the joint likelihood reduces to the existing complete-population selection likelihood.
2. **No-event-selection limit:** when every event is retained, the detection normalization tends to one and the inference reduces to the catalogue-only model.
3. **No-missingness limit:** catalogue completeness of 100% removes the latent missing-population term.
4. **Weight conservation:** observed and missing population masses sum to one under the declared generative model; no double counting and no silent renormalization over only observed hosts.
5. **Selection denominator:** alpha(H) is computed over the same event/source/host population and detection rule used by the numerator.
6. **Matched generator/inference:** prior-predictive SBC is tested only after all prior, host-population, catalogue-inclusion, noise, and detection steps match.
7. **Independent numerical calculation:** compare selected-likelihood normalization and posterior quantiles against a separately implemented quadrature or enumeration on small finite examples. The current focused tests now include an independent direct Gaussian-mixture/quadrature reconstruction and a high-threshold no-selection limit check; these are regression gates, not evidence that the full joint catalogue-selection model is validated.
8. **Prior-predictive selected-likelihood calibration:** simulate H₀ from the declared prior, draw hosts and noisy distances from the same finite mixture, reject events outside the declared detection interval, and compute posterior PIT values using the selected likelihood. A deterministic 300-draw regression test checks for gross PIT non-uniformity. The generator draws H₀ once from the declared prior per replicate, then resamples host and distance conditional on that fixed truth until the distance falls inside the detection window. It must not redraw H₀ after rejected events: doing so weights the accepted truth distribution by detection efficiency and changes the prior-predictive design. This is a smoke test for a matched toy likelihood, not a definitive SBC study or astrophysical validation. The earlier implementation reported mean PIT 0.493 and a maximum absolute deviation of 0.030 from the expected 0.10 frequency across ten equal-width bins, but those historical diagnostics apply to the earlier generator and must not be cited as validation of the corrected generator. The corrected generator requires a fresh workflow run and fresh results before any numerical calibration claim is made. A separate 1,000-replicate runner, `scripts/run_selected_likelihood_sbc.py`, has now been added to the numerical-diagnostics workflow. It archives replicate-level truth/host/distance/PIT records and reports PIT mean Monte Carlo standard error, decile Wilson intervals, and a descriptive KS statistic. The first completed 1,000-replicate run and independent audit are documented in `docs/selected_likelihood_sbc_results_2026-10-25.md`. Its PIT behavior is consistent with calibration for that matched synthetic design only; it is not astrophysical validation. Subsequent numerical-diagnostics commits rerun the audit.

Each test should include a simple finite two- or three-host example where the answer can be computed by direct enumeration. Only then scale to repeated simulation.

## 5. Interpretation of existing missing-host work

The repository's existing `docs/missing_host_model_audit.md` identifies the current oracle missing-host comparison as a hybrid diagnostic: the hidden component can use oracle simulation information while the observed component retains uniform inference weights. This is useful for probing the hidden term but is not a fully oracle-weighted posterior.

Therefore, do not treat that existing oracle arm as a validated practical correction. A clean follow-up should separately label:
- observed-only practical baseline;
- missing-population model using explicitly assumed inference weights;
- current hybrid oracle diagnostic;
- full-population oracle benchmark using consistent weights for observed and missing galaxies.

The full-population oracle is a simulation-only information benchmark, not an estimator available to real catalogue analyses.

## Go/no-go decision

**Go** to a small finite-population derivation and unit-test prototype. **No-go** to the large combined factorial until the denominator, latent missing-population model, and complete/no-selection/no-missingness limits are independently verified. A mathematically coherent small model is the next milestone; a large simulation matrix is not.
