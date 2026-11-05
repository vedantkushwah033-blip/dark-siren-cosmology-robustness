# Statistical model and implemented likelihood boundaries

**Status:** source-audited description of implemented synthetic models. This document is not a claim that the repository contains a complete realistic gravitational-wave cosmology likelihood.

## 1. Shared low-redshift distance model

The core toy relation is

[
d_L(z,H_0)\simeq \frac{c z}{H_0},
]

with a Gaussian observed-distance model

[
d_{obs}\mid z,H_0\sim\mathcal N\!\left(\frac{cz}{H_0},\sigma^2\right).
]

This is a deliberately simplified low-redshift approximation. It is not a full luminosity-distance calculation in a specified cosmology and does not include peculiar velocities, lensing, redshift uncertainty, or a detector strain likelihood.

## 2. Finite catalogue/missing-host mixture

The dataset-level event likelihood in `src/dark_siren_robustness/dataset_likelihood.py` uses a catalogue component and a missing-host component. Let (f\in[0,1)) be the **pre-detection population mass** assigned to the missing component. Let (	ilde w_{c,j}) and (	ilde w_{m,k}) be non-negative within-component weights normalized separately to sum to one. Then

[
L_e(H_0)=
(1-f)\sum_j\tilde w_{c,j}\,
\mathcal N\!\left(d_e\mid cz_{c,j}/H_0,\sigma_e^2\right)
+
f\sum_k\tilde w_{m,k}\,
\mathcal N\!\left(d_e\mid cz_{m,k}/H_0,\sigma_e^2\right).
]

The two component masses sum to one before detection. The missing fraction is not an after-the-fact fraction of a normalized likelihood. When (f>0), missing-host redshifts are required. The dedicated selected-posterior implementation in `selection.py` follows the same convention by combining the full host support and full-population weights before evaluating the event likelihood and selection efficiency.

For conditionally independent events sharing (H_0), the unconditional fixed-sample log likelihood is

[
\ell_{uncond}(H_0)=\sum_{e=1}^{N}\log L_e(H_0).
]

This assumes conditional independence given the shared parameters and fixed supplied host mixture. It is not an event-rate model.

## 3. Detection selection

For the finite-interval toy model in `selection.py`, detection is

[
d_{min}<d_{obs}<d_{max}.
]

For host (j), with (mu_j(H_0)=cz_j/H_0),

[
P(D=1\mid H_0,j)=
\Phi\!\left(\frac{d_{max}-\mu_j(H_0)}{\sigma}\right)
-
\Phi\!\left(\frac{d_{min}-\mu_j(H_0)}{\sigma}\right).
]

The population efficiency must use the same complete host mixture and weights as the event numerator:

[
\beta(H_0)=
(1-f)\sum_j\tilde w_{c,j}P(D=1\mid H_0,j)
+f\sum_k\tilde w_{m,k}P(D=1\mid H_0,k).
]

For a fixed number (N) of conditionally independent detected events, when the event numerator is unconditional,

[
\ell_{det}(H_0)=
\sum_{e=1}^{N}\log L_e(H_0)-N\log\beta(H_0).
]

The one-event selected likelihood is the corresponding (L_e(H_0)/\beta(H_0)). Do not apply this correction to event likelihoods that already condition on detection. Do not mix an event-specific candidate-host distribution with a population-level selection efficiency unless the generative model explicitly justifies that identification.

A separate one-sided-threshold toy helper exists in `dataset_likelihood.py`; a continuous-density, one-sided toy likelihood exists in `selection_toy.py`. These are related but not interchangeable without matching their detection domains and host-population assumptions.

## 4. Prior and numerical posterior

For prior density (pi(H_0)), the posterior is proportional to

[
p(H_0\mid\{d_e\},D=1)
\propto \pi(H_0)\exp\{\ell_{det}(H_0)\}.
]

The current baseline utilities generally use a uniform prior on an explicitly declared finite interval. Grid-based posterior summaries use numerical quadrature. The prior support, grid endpoints, grid spacing, and numerical convergence must be reported; a normalized grid posterior is not evidence that the grid is sufficiently resolved.

## 5. Catalogue inclusion and host-density estimation

Catalogue inclusion (R\) is distinct from GW event detection (D). The estimated-inclusion study models

[
P(R=1\mid x)=\operatorname{logit}^{-1}(\alpha+\beta x+\gamma(x^2-c_{parent})).
]

The parent-law centering constant is (c_{parent}=E_{p(parent)}[x^2]), not automatically (1/3). The code computes it under the declared non-uniform synthetic parent distribution. The analyst fits linear or quadratic logistic inclusion models to an independently sampled parent-frame validation subset and uses inverse fitted probabilities to weight the retained catalogue KDE.

This is an experiment about estimating catalogue retention from an external validation design. It is not a model for GW detector selection. Oracle probabilities are available only in simulation scoring/diagnostic arms.

## 6. Implemented versus unimplemented components

Implemented toy components include:
- low-redshift Gaussian distance likelihoods;
- finite weighted host mixtures;
- an explicit missing-host population-mass component;
- toy detection efficiencies and fixed-count selection normalization;
- catalogue inclusion/thinning and weighted host-density estimation;
- linear and quadratic inclusion-probability fitting from validation data;
- pre-posterior diagnostics and independent artifact reconstruction.

Not established as implemented or validated:
- a realistic cosmological luminosity-distance model;
- a GW strain-level or LVK detector-selection likelihood;
- an event-rate/population-count likelihood;
- joint marginalization over uncertain host-population, catalogue-inclusion, and detector-selection parameters;
- realistic galaxy-survey masks, clustering, luminosity functions, and redshift errors;
- validation on LVK/GWTC observations or real survey catalogues.

## 7. Interpretation requirements

Keep separate:
1. prior-predictive Bayesian calibration under the assumed model;
2. fixed-truth repeated-sampling coverage;
3. posterior interval width/sharpness;
4. Monte Carlo uncertainty in estimated performance;
5. numerical grid/quadrature error;
6. model-misspecification sensitivity.

None of these alone establishes the others. In particular, simulation-based calibration under a correctly specified model does not validate robustness to misspecification, and synthetic experiments do not constitute an observational (H_0) measurement.

See [the 5 November repository baseline and release audit](research_release_baseline_2026-11-05.md) for the source-level inventory and outstanding release gates.
