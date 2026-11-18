# Matched quadratic-logit H0 smoke readout

## Execution and artifact provenance
- Workflow: 44404473699 (completed SUCCESS)
- Commit: e9be5bc6cea9abad3b90dab0fd3dd41971ab1354
- Artifact ID: 14472608790
- Artifact archive SHA-256: 2f6a0a2dce9b792f3b7c061d0430af76cfeb6c2bf42b5c44c1cf9392affb15dd
- Downloaded archive SHA-256 independently recomputed and matched the GitHub artifact digest.
- Independent audit: PASS; all 16 manifest entries verified; selected-likelihood reconstruction maximum absolute density differences were exactly 0.0 for quadrature orders 3 and 5; posterior normalization checks passed at both orders; prior-scale sensitivity reconstruction differences were exactly 0.0 for multipliers 0.5, 1, and 2.
- Quadrature order comparison: order 3 to 5 passed configured thresholds; status CONVERGED at order 5.

## Synthetic design
- Quadratic-logit catalogue inclusion DGP: gamma=+1, beta=-3.9, alpha=0.12530550566919094.
- Degree-one logistic inclusion model fitted in inference (deliberately misspecified at gamma=+1).
- Parent N=1,200; observed catalogue 576.
- Independent validation frame 500; SRS validation sample 150.
- Four detected GW events; 1,000 generation attempts; 932 detected candidates; 68 rejected candidates.
- H0 truth 70 km s^-1 Mpc^-1; uniform prior [20,300]; 281-point H0 grid (spacing 1 km s^-1 Mpc^-1).
- Low-redshift Euclidean distance toy with fixed Gaussian noise sigma=140 Mpc and fixed detection threshold 250 Mpc.

## Posterior result
At quadrature order 5:
- H0 median: 130.03365 km s^-1 Mpc^-1.
- Central 68% interval: [108.71496, 153.39148].
- Central 95% interval: [86.36860, 181.06454].
- The generating value 70 is outside this single-realization 95% interval.
- Quadrature ESS: 6.93; maximum quadrature node mass: 0.281.
- Minimum predicted catalogue inclusion probability across quadrature nodes: 0.00376.
- Maximum normalized catalogue weight: 0.08110; minimum catalogue Kish ESS across nodes: 43.20.

Prior slope-scale sensitivity (same observations, grid and likelihood):
- multiplier 0.5: median 127.23544.
- multiplier 1.0: median 130.03365.
- multiplier 2.0: median 130.89954.
This is a modest shift relative to the large posterior offset; it does not identify the cause.

## Truth-only post-hoc diagnostic
The observed event distances were 780.48, 866.76, 322.96, and 610.12 Mpc. The truth-only latent mean distances were 1096.47, 779.26, 315.33, and 562.50 Mpc. Standardized residuals (observed minus latent mean, divided by sigma=140 Mpc) were approximately -2.26, +0.62, +0.05, and +0.34. These truth-only values are not inference inputs; they indicate that the first event is a comparatively low distance-noise realization. With only four events, this can materially move the posterior.

## Interpretation
This is a matched-DGP software and arithmetic smoke, not a repeated-sampling result. The posterior median is far above H0 truth and the truth is outside the central 95% interval for this realization. That is a warning requiring repeated-sampling calibration, not evidence that the method is biased or that the DGP is broken. The event-level noise realization, small event count, catalogue completeness, model misspecification and selection normalization all require assessment across replicates.

The H0 grid has 281 points over [20,300] (spacing 1). The separate PR #79 grid-refinement experiment compares 141, 281 and 561 points on identical observable inputs. No SBC, coverage, real-data, LVK, or observational H0 claim is supported by this smoke.
