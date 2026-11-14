# Parent-calibrated quadratic-logit inclusion stress design (2026-11-14)

## Research purpose

The completed floor-family density pilot changes the shape of the inclusion curve through an additive floor. A separate quadratic-logit family is needed to isolate curvature misspecification while holding parent-averaged catalogue completeness fixed.

The experiment asks whether finite-sample estimation of a degree-one logistic inclusion model deteriorates as the true inclusion mechanism acquires quadratic curvature, even when the mean inclusion rate is held constant.

This is a synthetic statistical stress test, not a real survey model.

## Parent measure and coordinate

The parent redshift distribution is proportional to the existing latent intensity:

[
g(z)=1+2.2expleft[-rac{(z-0.105)^2}{2(0.025)^2}ight]
+1.6expleft[-rac{(z-0.275)^2}{2(0.040)^2}ight],
quad zin[0.01,0.40].
]

Use:

[
x=(z-0.205)/0.195.
]

Numerical integration on the same 20,001-point grid used by the implementation gives:

[
c=E_g[x^2]=0.2845005579.
]

The centering is with respect to the parent sampling measure (g(z)), not the host/merger-weighted target density (g(z)b(z)). The host/merger weight must not enter catalogue inclusion calibration.

## Inclusion mechanism

[
pi_gamma(z)=
operatorname{expit}left[
alpha_gamma-3.9x+gamma(x^2-c)
ight].
]

The reference model is gamma=0 with alpha_0=0. Its parent-average inclusion rate is:

[
r_0=E_g[operatorname{expit}(-3.9x)]
=0.5101919962.
]

**Correction to an earlier preliminary note:** 0.5487448435 was an incorrect preliminary rate and is superseded. It must not be used in configs or analysis.

For every gamma, alpha_gamma is deterministically calibrated so:

[
E_g[pi_gamma(z)]=r_0.
]

The calibrated intercepts on the implementation's 20,001-point quadrature grid are:

| gamma | alpha_gamma | achieved parent-average rate |
|---:|---:|---:|
| -2.0 | -0.2404784703 | 0.5101919962 |
| -1.0 | -0.1247216935 | 0.5101919962 |
| -0.5 | -0.0629707145 | 0.5101919962 |
| 0.0 | 0.0000000000 | 0.5101919962 |
| 0.5 | 0.0631226953 | 0.5101919962 |
| 1.0 | 0.1253055057 | 0.5101919962 |
| 2.0 | 0.2424413834 | 0.5101919962 |

Gamma zero is exactly degree-one logistic truth. Nonzero gamma values are structurally misspecified by a degree-one logistic analysis model. Equal parent-average inclusion prevents average completeness from being a confounder of the curvature comparison.

## Implemented DGP component

Added to `src/dark_siren_robustness/generative_selection.py`:

- `quadratic_logit_parent_measure`: computes the parent-measure centering constant.
- `calibrate_quadratic_logit_intercept`: calibrates alpha for a specified curvature.
- `quadratic_logit_inclusion_probability`: evaluates the calibrated probability.
- `generate_quadratic_logit_inclusion_samples`: creates the parent catalogue and independent validation frame.

The independent validation sample is generated with a separate seed. Its observable table contains only validation ID, redshift, and retention outcome. Inclusion probabilities and the full validation frame remain truth-only. The generated configuration records the model, calibration constants, seeds, and sample sizes.

Regression tests check:
- parent-measure centering;
- equal parent-average rates across all seven gamma values;
- the gamma-zero degree-one logistic identity;
- nonzero curvature changing the inclusion curve;
- separation of observable and truth-only fields;
- deterministic reproduction under fixed seeds.

## Implementation readiness update (2026-11-14)

The following infrastructure is now implemented on the curvature-fix draft PR:

- The single-replicate density runner supports the parent-rate-calibrated quadratic-logit DGP and supplies the exact truth probabilities only to Arm C.
- The independent raw-input auditor independently reconstructs the parent-measure centering constant, gamma-zero target rate, achieved parent-average inclusion rate, validation-frame probabilities, and degree-one logistic posterior/quadrature.
- The multi-replicate runner and outer auditor support one gamma cell per run, preserve all replicate attempts, and verify each child independently.
- The cross-gamma auditor checks the seven gamma artifacts, common configuration, shared seed schedule, paired parent/validation redshift samples, calibration constants, and produces density summaries plus within-gamma estimator contrasts.
- It additionally reports paired same-arm IAE contrasts of each nonzero-curvature cell against the correctly specified gamma=0 reference, requiring identical replicate/validation-size/bandwidth/arm keys rather than comparing unpaired aggregate means.
- A manual-only workflow is defined at `.github/workflows/quadratic-logit-curvature-density-pilot.yml`.

The focused regression workflow and full repository numerical-diagnostics matrix both passed on Python 3.10, 3.11, and 3.12 at code head `afb2b647fae28c46f06c18f554019f007b4f8b1b`. Subsequent additions of the gamma-versus-zero contrast helper and its tests are pending exact-head CI. The seven-gamma workflow has **not** been dispatched; no curvature-family experiment results exist yet. The cross-gamma gate now also requires every planned replicate to succeed in every gamma cell; a child-level audit may preserve and report failures, but the paired curvature summary will refuse an incomplete cell rather than silently analyze a reduced successful subset.

## Proposed full density experiment (not yet run)

The next stage should use:

- gamma values: -2, -1, -0.5, 0, 0.5, 1, 2;
- validation sizes: 100, 300, 1,000;
- independent validation frame: 1,500;
- parent population: 4,000;
- five arms A-E, preserving the existing definitions;
- primary Arm E estimand: loss of posterior-mean density;
- secondary Arm E estimand: posterior expected node-specific loss;
- common physical bandwidth grid, initially inherited from the expanded-grid review;
- paired replicate seeds across gamma values and validation sizes;
- explicit lifecycle ledger and preservation of all failures;
- independent child reconstruction and outer aggregate audit.

The analysis should report paired contrasts against catalogue-only, known-probability IPW, and parent oracle, along with ESS and weight-tail diagnostics. Because mean inclusion is fixed across gamma, differences can be interpreted as curvature stress rather than a shift in overall catalogue retention.

Do not use gamma values or validation sizes to tune the model after inspecting holdout results. Freeze the density estimand and bandwidth-training/holdout split before confirmatory runs.

## Scope boundaries

This experiment is density-only. It does not evaluate the GW detection-selection function, H0 likelihood, H0 posterior calibration, real galaxy surveys, or public LVK/GWTC data. The separate GW host-weighting model and GW detection-selection normalization remain distinct from catalogue inclusion.

No curvature-family performance claim is made until the seven-gamma workflow has completed, every per-gamma child and outer audit has passed, and the cross-gamma paired audit has verified the shared design.
