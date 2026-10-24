# Protocol: separate catalogue incompleteness from event-detection selection

**Status:** protocol plus an initial analytic implementation and unit tests; no selection-normalized simulation result is claimed.  
**Purpose:** define a minimal, testable selection model before implementing it.  
**Scope:** controlled synthetic experiments, not a realistic gravitational-wave detection pipeline.

## Why this extension is needed

Catalogue incompleteness and gravitational-wave event-detection selection are different mechanisms:

- **Catalogue inclusion** determines whether a galaxy is present in the redshift catalogue available to the analysis.
- **Event detection** determines whether a gravitational-wave event enters the analysed event sample.

Changing which galaxies are listed does not by itself implement the selection normalization required when inference conditions on a detected event. Conversely, a detection correction does not restore redshifts for galaxies missing from a catalogue.

The current baseline likelihood is a weighted host mixture using a simplified low-redshift relation. It does not yet provide a validated, realistic event-detection selection function. This document defines a toy extension that can be tested analytically before any attempt to use real observations.

## Toy generative model

For a fixed finite set of candidate hosts with redshifts $z_j$, prior host weights $w_j$ satisfying $w_j\geq0$ and $\sum_j w_j=1$, and fixed measurement uncertainty $\sigma$, define

$$
\mu_j(H_0)=\frac{c z_j}{H_0},\qquad
d_{\rm obs}\mid(H_0,j)\sim\mathcal{N}(\mu_j(H_0),\sigma^2).
$$

For this first selection experiment only, define detection by a transparent finite noisy-distance interval:

$
D=1\quad\text{if}\quad 0<d_{\rm obs}<d_{\rm th}.
$

This is a deliberately simplified pedagogical selection rule. It is **not** a model of a real detector's matched-filter signal-to-noise ratio or search pipeline.

For host $j$, the detection probability under this toy rule is

$$
P(D=1\mid H_0,j)=\Phi\!\left(\frac{d_{\rm th}-\mu_j(H_0)}{\sigma}\right)
-\Phi\!\left(\frac{-\mu_j(H_0)}{\sigma}\right),
$$

where $\Phi$ is the standard normal CDF. The population detection efficiency is

$$
\beta(H_0)=\sum_j w_j P(D=1\mid H_0,j).
$$

For an observed event conditioned on detection, the toy likelihood is proportional to

$$
p(d_{\rm obs}\mid H_0,D=1,\{z_j,w_j\})
=\frac{\sum_j w_j\,\mathcal{N}(d_{\rm obs}\mid\mu_j(H_0),\sigma^2)}
{\beta(H_0)}
$$

for observations in the detected region. With a prior $\pi(H_0)$, the posterior is proportional to this conditional likelihood times the prior. The implementation must be careful not to normalize away the $H_0$ dependence in $\beta(H_0)$ before forming the posterior.

For multiple independently detected events under a shared population model, the corresponding selection normalization generally enters once per event in the conditional likelihood, or as $\beta(H_0)^{-N}$ in the product, subject to the exact event/population conditioning assumptions. Do not implement a multi-event expression until its conditioning assumptions are written down.

## Required experiment matrix

1. **No-selection baseline:** generate and infer without a detection cut. This must reproduce the existing baseline within Monte Carlo uncertainty.
2. **Selected data, uncorrected inference:** generate only observations satisfying $D=1$, then infer without dividing by $\beta(H_0)$. This is a deliberate misspecification control.
3. **Selected data, correctly normalized inference:** use the same detection rule in the generator and $\beta(H_0)$ in inference.
4. **Selection mismatch:** generate with threshold $d_{\rm th}$ but infer with a different threshold. Label this as deliberate misspecification.
5. **Catalogue incompleteness crossed with detection:** compare complete, random-incomplete and redshift-dependent-incomplete catalogues separately from the detection rule. If a host is missing from the observed catalogue, do not silently use its true redshift in the practical inference branch.
6. **Oracle benchmark:** a branch allowed to use the full simulated host population may be included only as an explicitly labelled simulation-only upper-information control.

Use paired random seeds across scenarios wherever possible so the same underlying mock universe can be compared under different analysis choices. Record both the latent host population and the observed catalogue mask for reproducibility, but ensure the practical inference code receives only information it would have access to in that scenario.

## Required verification tests before interpreting results

### Analytic identities

- Verify $0<\beta(H_0)\leq1$ for finite, valid inputs.
- Verify that the detection probability approaches 1 as $d_{\rm th}-\mu_j$ becomes large and negative/positive in the mathematically expected direction.
- For a single host, verify $\beta(H_0)=\Phi((d_{\rm th}-c z/H_0)/\sigma)$.
- For a mixture, verify $\beta(H_0)$ equals the weighted sum of host-specific detection probabilities.
- Verify that equal host redshifts make the mixture likelihood invariant to redistributing normalized host weights.
- Verify the no-selection limit as $d_{\rm th}\to\infty$: $\beta(H_0)\to1$ and the selection-corrected posterior converges to the unselected posterior.

### Simulation checks

- Prior-predictive SBC under the correctly specified selected-data model.
- Fixed-truth repeated-sampling coverage at several $H_0$ values, including prior interiors and boundaries.
- Bias, RMSE, interval width and coverage, with Monte Carlo standard errors or intervals.
- A deliberately omitted-$\beta(H_0)$ negative control that demonstrates whether the selection mechanism actually matters under the chosen setup.
- Grid-refinement comparison for $\beta(H_0)$ and posterior quantiles.
- Checks that no latent host or hidden catalogue information leaks into the practical inference path.

### Statistical interpretation

Do not require frequentist fixed-truth coverage to equal nominal Bayesian credible mass at every parameter value. Report prior-predictive calibration and fixed-truth repeated-sampling performance as distinct diagnostics. A non-significant SBC uniformity test is not proof of calibration.

## Implementation sequence

1. **Completed:** add a separate selection module; do not silently alter the legacy simulator. The initial implementation is `src/dark_siren_robustness/selection.py`.
2. **Completed:** unit-test the detection probability and $\beta(H_0)$ against the analytic formula. See `tests/test_selection_normalization.py` and [the successful focused workflow](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40068498540).
3. Compare the selection-corrected likelihood to an independently coded reference calculation.
4. Add a one-host selected-data benchmark before a multi-host model.
5. Run prior-predictive SBC with selection.
6. Run the crossed catalogue-completeness experiment only after the selected-data baseline passes.
7. Document the toy threshold and explicitly prohibit claims that it approximates real LVK selection without a separately justified detector model.

## Out of scope for this protocol

This toy threshold model does not include real GW strain data, detector duty cycle, antenna patterns, waveform/template-bank effects, network signal-to-noise ratio, observing-run sensitivity, redshift uncertainty, lensing, peculiar velocities, or the actual GW detection pipeline. It is a controlled test of selection-normalization logic, not a physical detector model.

The analytic unit tests have passed, but no selected-data simulation result should be reported until the reference comparison and prior-predictive checks pass and the resulting artifacts have been inspected.
