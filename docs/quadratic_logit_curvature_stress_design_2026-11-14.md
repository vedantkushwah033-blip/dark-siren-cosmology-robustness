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

## Density smoke implementation (added 2026-11-15)

Two separate scripts now exercise the five-arm density calculation under the quadratic-logit DGP:

- `scripts/run_quadratic_logit_curvature_density_smoke.py` creates one gamma-specific synthetic replicate, stores observable and truth-only records separately, passes true quadratic inclusion probabilities to Arm C as an oracle diagnostic, and fits the degree-one logistic analysis model for Arms D/E.
- `scripts/audit_quadratic_logit_curvature_density_smoke.py` verifies the manifest, regenerates the population and independent validation frame from recorded seeds, checks the observable/truth artifacts, then independently reconstructs density metrics, weight diagnostics, and Scott widths from raw inputs.
- `tests/test_quadratic_logit_curvature_density_smoke.py` checks the gamma-zero and curved cases, nested validation prefixes, observable-field separation, oracle-weight construction, and output structure.

This is a one-replicate engineering smoke path, not the seven-gamma experiment. Its CI result is pending at the time of this design update; no curvature-family performance conclusion is made from it.

## Paired multi-replicate orchestration (implemented 2026-11-15; CI pending)

Added:

- `scripts/run_quadratic_logit_curvature_density_replicates.py`: enumerates the replicate × gamma design, runs the audited single-cell density runner, preserves every cell attempt in a lifecycle ledger, and aggregates density, weight-diagnostic, and Scott-width records.
- `scripts/audit_quadratic_logit_curvature_density_replicates.py`: verifies the outer manifest and complete planned cell grid, checks paired seed schedules and failure preservation, invokes the independent child raw-input auditor for every successful cell, and reconstructs all aggregate tables from child outputs.
- `tests/test_quadratic_logit_curvature_density_replicates.py`: checks a compact paired 2-replicate × 3-gamma run and rejects aggregate corruption even when the outer manifest is refreshed.
- The focused curvature workflow now includes both single-cell and multi-replicate tests plus a compact Python 3.12 end-to-end run.

The runner uses common random numbers across gamma values within each replicate: all gamma cells share that replicate's population and independent-validation seeds. Replicate seed pairs remain disjoint. This improves precision of within-replicate gamma contrasts, but does not make the cells statistically independent. All gamma-specific DGPs still recalibrate the intercept to the same parent-average inclusion rate.

The new orchestration has not yet completed CI at this exact head. The full seven-gamma × 20-replicate experiment has not been run, and no curvature-family result is reported.

## Proposed full density experiment (not yet run)

The production-scale design should use:

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

No curvature-family performance claim is made until the multi-replicate runner, independent auditor, and artifacts are complete.
