# Independent likelihood and selection-normalization audit

**Date:** 2026-10-29  
**Scope:** equation-level review of the discrete-host, H0-dependent-noise synthetic model in `selection_conditioned_h0_dependent_noise.py`. This is a documentation audit, not an independent rerun of the artifact or a physical validation of gravitational-wave inference.

## 1. Declared generative model

For a candidate host at redshift z_j and trial Hubble parameter H, the toy model defines

\[
\mu_j(H)=c z_j/H, \qquad
\sigma_j(H)=\sqrt{\sigma_{\rm floor}^2+[f\,\mu_j(H)]^2}.
\]

With normalized discrete candidate-host masses w_j, the distance-mixture likelihood for one event is

\[
L_i(d_i\mid H)=\sum_j w_{ij}\,
\mathcal{N}\!\left(d_i;\mu_j(H),\sigma_j(H)^2\right).
\]

The implementation evaluates the sum in log space. When the optional positive-truncated noise model is enabled, each Gaussian component is divided by P(d>0 | z_j,H), so the component density integrates to one on positive distances.

## 2. Selection conditioning

For the declared population redshifts and normalized population masses v_j, the detection efficiency is

\[
\beta(H)=\sum_j v_j\,
P(d\ge d_{\min}\mid z_j,H).
\]

For untruncated Gaussian noise, each survival probability is the Gaussian upper tail at d_min. For the positive-truncated option, the tail is divided by P(d>0 | z_j,H). Given a fixed number N of detected events, the implemented conditional log likelihood is

\[
\log L_{\rm sel}(H)=\sum_{i=1}^{N}\log L_i(d_i\mid H)-N\log\beta(H).
\]

This is appropriate for the stipulated fixed-count, independently detected toy sample when the supplied population mixture represents the selection-generating population. It does not model event-rate information, count likelihoods, detector response, realistic GW distance posteriors, or a joint latent galaxy-catalogue model.

## 3. Prior and numerical integration

The prior is uniform on the declared finite interval [H_min,H_max]. The implementation inserts exact prior endpoints into a normalization grid and linearly interpolates the sampled log likelihood to those endpoints before trapezoidal integration. The returned `normalization_density` is normalized on that endpoint-inclusive grid.

The separate `posterior_density` array is evaluated only at the requested H grid nodes and is zero outside prior support. If the prior endpoints do not coincide with H-grid nodes, integrating that sampled array directly on the original grid need not reproduce the endpoint-inclusive integral exactly; consumers requiring normalized summaries should use the endpoint-inclusive normalization arrays or explicitly renormalize on their integration grid. This distinction should remain explicit in downstream interfaces and reports.

## 4. Checks already present in the repository

The test suite includes checks for:
- direct equality of the selection-conditioned likelihood to the sum of event log likelihoods minus N log beta(H);
- positive-truncation survival normalization;
- Monte Carlo agreement between analytic detection efficiency and simulated detection frequency;
- exact prior endpoints in the normalization grid when they fall between H-grid nodes;
- deterministic simulation and the rejection of invalid or below-threshold inputs;
- a prior-predictive rank/coverage diagnostic and a negative control omitting selection normalization.

These tests check internal consistency under the toy assumptions. They are not evidence of realistic detector-level calibration.

## 5. Interpretation and release boundary

1. Treat prior-predictive SBC as a check under the stipulated generative model; fixed-truth experiments answer a different conditional-calibration question.
2. Report Monte Carlo uncertainty and avoid interpreting a non-significant contrast as equivalence.
3. Preserve the distinction between a known-selection oracle and a realistic analysis where selection and host-population uncertainty must be inferred.
4. Do not claim a CHIMERA v2.2 benchmark until the pinned notebook's exact inputs are obtained, checksummed, and executed in a version-compatible environment.
5. Do not claim a real-data H0 result or detector-level validation from these synthetic studies.

**Status:** equation-level documentation only. It does not add new simulation results, change the scientific estimand, or certify publication readiness.
