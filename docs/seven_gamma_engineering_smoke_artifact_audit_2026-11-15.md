# Seven-gamma curvature engineering smoke: independent artifact inspection

Date: 2026-11-15  
Artifact: `quadratic-logit-curvature-seven-gamma-smoke`  
Workflow run: `43671806461`  
Artifact ID: `14144815968`  
Artifact ZIP SHA-256: `290d77d9ff07bb949b10895d56fc510c84a8d6957d302fd3c18ab2f9f393a7e0`

## What was inspected

The downloaded ZIP contains seven gamma directories:
`-2,-1,-0.5,0,0.5,1,2`. Each directory contains one engineering replicate with:

- Parent population size: 4,000.
- Independent validation frame: 1,500.
- Nested validation sizes: 100, 300, 1,000.
- Nine absolute bandwidth candidates: 0.008 through 0.080 as specified in the workflow.
- Quadrature order 5 and 401-point redshift grid.
- Common population seed 20261115 and validation seed 20261116 across gamma values (common random numbers for paired stress comparisons).
- Degree-one logistic fitted inclusion model; Arm C uses true quadratic probabilities only as an oracle diagnostic.

The archive has 98 ZIP members and seven child manifests. I independently checked every file listed by every child manifest against its byte count and SHA-256: **84/84 files matched**, with zero mismatches. Each child `independent_audit.json` reports PASS. Each child audit reconstructed 837 density metric rows, 90 weight-diagnostic rows and 87 Scott-width rows. Across seven cells this is 5,859 metric rows, 630 diagnostic rows and 609 Scott-width rows.

The archive contains no outer seven-gamma summary file. The child manifests and independent audits establish integrity and reconstruction for each child artifact; they do not establish a paired multi-replicate experiment.

## Workflow failure and correction

The originating workflow run `43671806461` completed with **failure**, but the seven gamma simulations and all seven child audits completed successfully and the artifact was uploaded. The final summary step failed because it asserted that each child should contain 42 weight-diagnostic rows. The independently audited child output actually contains 90 rows: four non-posterior arms plus 25 posterior quadrature nodes and one aggregate posterior diagnostic for each of the three validation sizes.

The workflow assertion was corrected from 42 to 90 in commit `180011a96790284ef01b93d8752924290c8f1653`. A new workflow run was triggered by that change; its final result must be checked separately. The old failed run must not be relabelled as successful.

## Descriptive smoke values (not inferential evidence)

For traceability, the table below reproduces the single-replicate IAE values at absolute bandwidth 0.024, validation size 1,000, aggregate node `-1`. These are descriptive outputs from one common-seed engineering smoke only.

| gamma | A catalogue-only | B parent oracle | C known-probability IPW | D MAP estimated IPW | E posterior-mean density |
|---:|---:|---:|---:|---:|---:|
| -2.0 | 0.6534 | 0.1816 | 0.2057 | 0.2574 | 0.2563 |
| -1.0 | 0.6520 | 0.1816 | 0.1885 | 0.1520 | 0.1513 |
| -0.5 | 0.6483 | 0.1816 | 0.1758 | 0.1565 | 0.1588 |
| 0.0 | 0.6514 | 0.1816 | 0.1552 | 0.1946 | 0.1979 |
| 0.5 | 0.6380 | 0.1816 | 0.1699 | 0.2664 | 0.2706 |
| 1.0 | 0.6232 | 0.1816 | 0.1898 | 0.3483 | 0.3534 |
| 2.0 | 0.5921 | 0.1816 | 0.1800 | 0.5231 | 0.5276 |

These values are **not** a result about average performance, calibration, estimator superiority, or bandwidth choice. There is one replicate per gamma, and no sampling uncertainty can be estimated from this artifact. The table must not be used to choose a bandwidth or to make the research claim.

## Decision

**Engineering artifact integrity: PASS for the seven child directories.**  
**Originating workflow: FAILED at its final summary assertion.**  
**Scientific benchmark: NOT EVALUABLE from this one-replicate smoke.**  
**Training/holdout and H0 stages: remain blocked.**

The next valid gate is the corrected workflow run followed by the full multi-replicate seven-gamma pilot and its outer auditor. The training/holdout protocol remains unapproved until historical source-register completeness and seed-overlap checks are independently established.
