# Fresh-seed expanded parent-density bandwidth pilot

**Status:** preregistered engineering / density-only selection pilot; not an H0 experiment.  
**Branch:** `research/five-arm-comparison-preregistration-20261113`  
**Master seed:** `20261116`  
**Selection split seed:** `20261117`

## Motivation

Earlier density-only pilots repeatedly placed their minimum at the upper candidate boundary of the original 0.05–0.40 grid. An expanded 100-replicate analysis found an interior candidate at 0.50, but its replicate IDs 0–99 and master seed overlapped the earlier 200-replicate experiment. It is therefore exploratory, not independent confirmation.

This pilot uses a new master seed and a fresh 200-replicate source schedule. The runner uses `SeedSequence([master_seed, replicate_id, 1])`; with master seed 20261116, its parent samples are distinct from all earlier 20261114 pilots. The selection split uses a separate seed.

## Frozen design

- 200 independent replicates; 4,000 parent objects per replicate.
- 801-point redshift grid on [0.01, 0.40].
- Same declared synthetic parent intensity and common merger weight as the earlier parent-density pilot.
- Candidate weighted-Scott multipliers: `0.25, 0.30, 0.35, 0.40, 0.45, 0.50, 0.55, 0.60, 0.70, 0.80, 1.00`.
- 140 training replicates and 60 holdout replicates, determined by split seed 20261117.
- Selection metric: training mean integrated absolute error against the analytic merger-weighted parent-host redshift density.
- No H0 likelihood, posterior error, coverage, PIT, or event-level metric enters selection.
- The source auditor independently reconstructs raw IAE values and hashes; the selection auditor independently reconstructs the split and scores.

## Interpretation rules

1. This pilot is independent of prior bandwidth pilots at the parent-sample seed level.
2. The holdout is untouched within this new pilot and is descriptive; it cannot be used to reselect.
3. If the training minimum is at either candidate endpoint, label the result boundary-limited and do not freeze a bandwidth.
4. An interior minimum is only a candidate for this synthetic DGP, N=4,000, grid, and candidate set. It still requires holdout reporting, grid-resolution checks, and a separately preregistered H0 evaluation on disjoint replicate seeds.
5. The 200 source replicates are density-selection data. They must not also be used as confirmatory H0 evaluation replicates.
6. Keep raw inputs, attempt ledger, manifest, raw audit, selection result, and selection audit immutable as one artifact.

## Claim boundary

This experiment concerns synthetic parent-density recovery only. It does not compare inclusion-correction arms, establish H0 calibration, demonstrate estimator superiority, validate observational galaxy catalogues, or use LVK/GWTC data.
