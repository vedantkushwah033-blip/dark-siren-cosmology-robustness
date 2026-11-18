# Fixed-realization H0 grid-refinement audit

## Question
Does the discretization of the uniform H0 grid materially change the selected-likelihood posterior for one fixed synthetic realization?

## Paired design
The runner is executed three times with identical analysis and validation seeds, generator settings, observed-data process, inclusion model, prior, and quadrature tolerances. Only the H0 grid size changes:

| Grid points | H0 interval | Spacing |
|---:|---:|---:|
| 141 | [20, 300] | 2.0 |
| 281 | [20, 300] | 1.0 |
| 561 | [20, 300] | 0.5 |

The synthetic DGP is the quadratic-logit catalogue inclusion model with gamma=+1, beta=-3.9, H0_true=70, 1,200 parent galaxies, four detected events, and a separate validation frame of 500 with an SRS validation sample of 150. The inference model remains degree-one logistic and is therefore deliberately misspecified for gamma=+1.

Each run independently audits its own artifact, including selected-likelihood reconstruction, posterior normalization, truth firewall, and manifest integrity. The cross-grid auditor then requires byte-identical observable catalogue, GW-event, and validation-observation files before comparing normalized posterior densities on a common fine grid.

## Readouts
For each grid, the auditor reports:
- posterior median and 16/84/2.5/97.5 percentiles;
- maximum absolute density difference from the finest grid;
- total-variation distance from the finest grid;
- maximum absolute difference among the five posterior quantiles.

The finest grid (561 points) is the numerical reference, not mathematical ground truth. The comparison uses interpolated normalized densities on a 11,201-point common grid.

## Interpretation limits
This is a numerical discretization check for one fixed synthetic realization. It does not establish repeated-sampling calibration, SBC, coverage, estimator superiority, realistic GW inference, or an observational H0 measurement. Coefficient-quadrature convergence is separately required in each run; a run that fails that check is not included in the comparison.
