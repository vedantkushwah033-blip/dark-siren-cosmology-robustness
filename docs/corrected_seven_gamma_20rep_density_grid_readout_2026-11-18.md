# Corrected seven-gamma 20-replicate density-grid readout

**Date:** 2026-11-18  
**Status:** All seven gamma cells and the workflow-level cross-gamma raw-reconstruction audit completed successfully. The final paired-summary artifact was downloaded and independently checksum-verified.  
**Scope:** Synthetic low-redshift Euclidean catalogue-inclusion density recovery only.

## 1. Provenance and computational audit

Workflow: [44369630799](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44369630799)  
Workflow branch: `research/corrected-curvature-grid-20rep-20261118`  
Workflow head: `17ad47e3c4211abb7f0d3317c535d7fd4cbe6186`  
Base implementation: corrected PR #59 head `390ba806c42376b1a2d4b8d083b92b8d206c8a88`.

Design per gamma:
- 20 paired replicates.
- Parent population N=4,000; independent validation frame=1,500.
- Validation sizes n=100, 300, 1,000.
- Curvature gamma in {-2, -1, -0.5, 0, 0.5, 1, 2}; beta=-3.9.
- Parent-average inclusion target 0.5101919962291073, preserved across gamma by intercept calibration.
- Degree-one logistic inclusion fit for D and E; E integrates the degree-one validation-posterior quadrature.
- Quadrature order 5; redshift grid 401.
- Seven fixed absolute bandwidth SD candidates: 0.012, 0.016, 0.020, 0.024, 0.028, 0.032, 0.040.
- Scott multiplier 0.7 recorded for realized-Scott diagnostics.

All seven cell jobs completed successfully. Each cell's independent audit reported 20/20 successful replicates, zero failed replicates, 20/20 child raw-input audits passed, 244 manifest files verified, 13,020 aggregate metric rows reconstructed, and 1,740 realized-Scott-width rows reconstructed. Across the full grid: 140/140 successful replicates, zero failed replicates, 140/140 child audits, 1,708 manifest entries, 91,140 aggregate metric rows, and 12,180 realized-Scott-width rows.

I downloaded all seven raw ZIP artifacts and independently recomputed every manifest byte count and SHA-256 digest. All 244 entries matched in each of the seven archives (1,708/1,708 total). The downloaded ZIP digests matched the GitHub artifact digests:

| gamma | artifact ID | ZIP SHA-256 |
|---:|---:|---|
| -2.0 | 14458384420 | `dfa914e497be3d69efbaa989ff1fb66a98bea2a9b10c9b22a49e4fe5db8190ea` |
| -1.0 | 14459806525 | `3e0664234917414377d83578e0146ff02526a52d6d61fb2a23737fd6d1f37e22` |
| -0.5 | 14458623924 | `42a24e176a5570bf79b65c7963b86e8036db1e2c1330fbc9f20d62ea9d0885ba` |
| 0.0 | 14460640021 | `b4c22d261398d280cbb33b4e6d8b0790ce41fbc1a83f57bcdc030f5c833d367e` |
| 0.5 | 14457819695 | `22391674c846a8891bdb8ebe826e90fbe059cfb4da8ef80ba8a13b935e910dde` |
| 1.0 | 14458964091 | `158e1ad3627f02225baf2f75df7e7098be6be1c104d0c45e33a2947b7da559c7` |
| 2.0 | 14459496034 | `d3b29dcd15b951450abc70131235b3d69e4804701b3ce0949f551cbbc571d9df` |

Raw-weight spot checks at gamma=-2 and gamma=+2 independently confirmed A/B sample weights exactly equal one and C weights equal inverse true inclusion probabilities to floating-point precision (maximum absolute differences 8.9e-15 and 7.1e-14, respectively). The common host-rate factor is a downstream grid factor, not duplicated in A/B/C sample weights.

### Independent paired-design check

I compared all seven archives directly against the gamma=-2 reference for every replicate:
- Parent object IDs and redshifts.
- Validation-frame IDs, selected-for-validation flags, and redshifts.
- Observed validation IDs and redshifts.

All 7 x 20 = 140 matched replicate pairs agreed. I also independently checked the child calibration metadata: beta=-3.9, x2_center=0.2845005579191241, and parent-average inclusion target 0.5101919962291073 were common across gamma; achieved rates differed only at floating-point roundoff.

The workflow-level cross-gamma auditor subsequently completed successfully in run `44369630799`; its final raw-reconstruction report and paired-summary artifact are documented in Section 7. The manual paired-design checks above are supplementary and separate from that workflow audit.

## 2. Density IAE readout

For each replicate and validation size, I averaged IAE across the seven fixed bandwidth candidates, then summarized those replicate-level averages across the 20 paired replicates. This is a descriptive grid-average, **not** bandwidth selection. A/B/C are independent of validation size within each gamma because validation size affects the estimated-inclusion arms, not the fixed parent/catalogue realization.

| gamma | n_val | A catalogue | B oracle | C known pi | D MAP | E posterior-mean density |
|---:|---:|---:|---:|---:|---:|---:|
| -2 | 100 | 0.7174 | 0.1009 | 0.1700 | 0.2843 | 0.2644 |
| -2 | 300 | 0.7174 | 0.1009 | 0.1700 | 0.2557 | 0.2476 |
| -2 | 1000 | 0.7174 | 0.1009 | 0.1700 | 0.2557 | 0.2531 |
| -1 | 100 | 0.7240 | 0.1009 | 0.1538 | 0.2318 | 0.2161 |
| -1 | 300 | 0.7240 | 0.1009 | 0.1538 | 0.1827 | 0.1772 |
| -1 | 1000 | 0.7240 | 0.1009 | 0.1538 | 0.1769 | 0.1750 |
| -0.5 | 100 | 0.7234 | 0.1009 | 0.1486 | 0.2208 | 0.2150 |
| -0.5 | 300 | 0.7234 | 0.1009 | 0.1486 | 0.1641 | 0.1635 |
| -0.5 | 1000 | 0.7234 | 0.1009 | 0.1486 | 0.1532 | 0.1525 |
| 0 | 100 | 0.7224 | 0.1009 | 0.1359 | 0.2171 | 0.2191 |
| 0 | 300 | 0.7224 | 0.1009 | 0.1359 | 0.1667 | 0.1706 |
| 0 | 1000 | 0.7224 | 0.1009 | 0.1359 | 0.1477 | 0.1487 |
| 0.5 | 100 | 0.7191 | 0.1009 | 0.1281 | 0.2223 | 0.2354 |
| 0.5 | 300 | 0.7191 | 0.1009 | 0.1281 | 0.1903 | 0.1993 |
| 0.5 | 1000 | 0.7191 | 0.1009 | 0.1281 | 0.1715 | 0.1746 |
| 1 | 100 | 0.7147 | 0.1009 | 0.1234 | 0.2555 | 0.2753 |
| 1 | 300 | 0.7147 | 0.1009 | 0.1234 | 0.2369 | 0.2496 |
| 1 | 1000 | 0.7147 | 0.1009 | 0.1234 | 0.2264 | 0.2305 |
| 2 | 100 | 0.6984 | 0.1009 | 0.1201 | 0.3747 | 0.4058 |
| 2 | 300 | 0.6984 | 0.1009 | 0.1201 | 0.3783 | 0.3925 |
| 2 | 1000 | 0.6984 | 0.1009 | 0.1201 | 0.3754 | 0.3801 |

## 3. Paired estimated-IPW contrasts against known-probability IPW

Each contrast is computed within replicate after averaging across the seven bandwidth candidates. Intervals are ordinary paired t intervals over 20 replicates; they are exploratory, not confirmatory guarantees.

| gamma | n_val | D − C mean (95% paired t CI) | E − C mean (95% paired t CI) | D worse reps | E worse reps |
|---:|---:|---:|---:|---:|---:|
| -2 | 100 | 0.1143 [0.0786, 0.1500] | 0.0944 [0.0568, 0.1321] | 18/20 | 17/20 |
| -2 | 300 | 0.0857 [0.0507, 0.1207] | 0.0776 [0.0424, 0.1128] | 17/20 | 17/20 |
| -2 | 1000 | 0.0857 [0.0518, 0.1196] | 0.0831 [0.0492, 0.1171] | 17/20 | 17/20 |
| -1 | 100 | 0.0780 [0.0274, 0.1287] | 0.0623 [0.0110, 0.1136] | 16/20 | 16/20 |
| -1 | 300 | 0.0289 [-0.0053, 0.0631] | 0.0234 [-0.0087, 0.0556] | 13/20 | 13/20 |
| -1 | 1000 | 0.0231 [-0.0100, 0.0561] | 0.0212 [-0.0111, 0.0534] | 13/20 | 13/20 |
| -0.5 | 100 | 0.0722 [0.0173, 0.1271] | 0.0664 [0.0122, 0.1206] | 12/20 | 14/20 |
| -0.5 | 300 | 0.0155 [-0.0047, 0.0356] | 0.0148 [-0.0036, 0.0332] | 13/20 | 13/20 |
| -0.5 | 1000 | 0.0046 [-0.0123, 0.0215] | 0.0039 [-0.0121, 0.0199] | 9/20 | 8/20 |
| 0 | 100 | 0.0811 [0.0276, 0.1347] | 0.0832 [0.0288, 0.1376] | 15/20 | 17/20 |
| 0 | 300 | 0.0308 [0.0092, 0.0523] | 0.0347 [0.0110, 0.0584] | 14/20 | 15/20 |
| 0 | 1000 | 0.0117 [0.0035, 0.0200] | 0.0127 [0.0042, 0.0213] | 14/20 | 13/20 |
| 0.5 | 100 | 0.0943 [0.0394, 0.1492] | 0.1073 [0.0501, 0.1645] | 19/20 | 20/20 |
| 0.5 | 300 | 0.0622 [0.0313, 0.0932] | 0.0712 [0.0379, 0.1046] | 18/20 | 19/20 |
| 0.5 | 1000 | 0.0434 [0.0268, 0.0601] | 0.0465 [0.0293, 0.0638] | 19/20 | 19/20 |
| 1 | 100 | 0.1322 [0.0618, 0.2025] | 0.1519 [0.0793, 0.2245] | 20/20 | 20/20 |
| 1 | 300 | 0.1135 [0.0770, 0.1500] | 0.1262 [0.0880, 0.1645] | 20/20 | 20/20 |
| 1 | 1000 | 0.1030 [0.0814, 0.1245] | 0.1071 [0.0852, 0.1291] | 20/20 | 20/20 |
| 2 | 100 | 0.2546 [0.1888, 0.3204] | 0.2856 [0.2164, 0.3549] | 20/20 | 20/20 |
| 2 | 300 | 0.2582 [0.2160, 0.3004] | 0.2724 [0.2287, 0.3160] | 20/20 | 20/20 |
| 2 | 1000 | 0.2553 [0.2284, 0.2821] | 0.2599 [0.2329, 0.2870] | 20/20 | 20/20 |

## 4. Weight concentration at validation n=1000

Mean Kish ESS over the 20 replicates:

| gamma | C known-pi ESS | D MAP ESS | E posterior-node ESS |
|---:|---:|---:|---:|
| -2 | 381 | 1078 | 666 |
| -1 | 478 | 854 | 465 |
| -0.5 | 565 | 745 | 386 |
| 0 | 663 | 674 | 339 |
| 0.5 | 779 | 603 | 298 |
| 1 | 881 | 554 | 275 |
| 2 | 1051 | 512 | 267 |

The ESS pattern is not a monotone proxy for density accuracy. At negative curvature, fitted D weights can be less concentrated (higher ESS) than oracle C while still producing larger IAE; this is consistent with misspecification bias rather than a pure variance/weight-tail explanation. At positive curvature, D/E are both more concentrated than C and have much larger IAE. E's posterior-node ESS is particularly low across the grid.

## 5. Scientific interpretation

1. The known-probability IPW arm C is the strongest operational density reference among the non-oracle arms in this pilot. It is not a real-world method because its true inclusion probabilities are known only to the simulator.
2. Near mild negative curvature, the estimated arms approach C as validation size grows: at gamma=-0.5 and n=1000, D-C and E-C intervals overlap zero. At gamma=-1 and n=300/1000, the intervals also overlap zero.
3. At gamma=0, the degree-one model is correctly specified, but finite validation uncertainty still leaves a small positive D/E IAE gap at n=1000.
4. At gamma=+0.5, +1, and +2, D/E remain materially worse than C as validation size increases. The gap is particularly large at gamma=+2. At gamma=-2, the true oracle weights are highly concentrated; the fitted linear model has higher ESS but still higher IAE, again showing that ESS alone does not diagnose model bias.
5. Posterior integration over degree-one coefficients does not automatically improve density recovery. E is often similar to or worse than D; it is not a correction for omitted quadratic curvature.

This pattern supports the research question: finite validation can dominate under mild curvature, while model misspecification can create persistent error that does not disappear merely by increasing validation size. The positive/negative curvature conditions are not symmetric in difficulty because the calibrated inclusion-probability surfaces and weight tails differ.

## 6. Limitations and next gates

- Density recovery only: no H0 likelihood, H0 posterior, SBC, PIT, frequentist coverage, or posterior calibration is evaluated in this grid.
- No bandwidth selection: the seven absolute bandwidth candidates were fixed in advance and averaged only for descriptive readout.
- Twenty replicates per cell provide exploratory precision only.
- No observational galaxy catalogue or public LVK/GWTC data are used.
- Oracle C probabilities are simulation-only diagnostics.
- The child-level config still inherits stale wording that each child is “one engineering replicate only.” PR #62 proposes a metadata-only correction and regression test; do not treat that wording issue as a numerical error, but fix it before publication-ready provenance.
- The final workflow cross-gamma audit has passed, and its paired-summary artifact checksum was verified. This is still a synthetic density-only pilot, not an H0 inference result.

Next: finish the independent bandwidth-sensitivity/holdout gate before deciding whether density reconstruction is stable enough to propagate into the H0 likelihood. The H0 propagation design contract is `docs/finite_validation_h0_marginalization_design_contract_2026-11-08.md`. Density IAE must not be presented as evidence of H0 posterior calibration.


## 7. Final workflow-level cross-gamma audit (completed)

The workflow run `44369630799` finished with conclusion `success`. The `Cross-gamma paired audit and summaries` job completed all steps successfully, including full raw-input reconstruction and upload. The final cross-gamma artifact is ID `14462504250`, SHA-256 `c4e86b59dfa37f93b0fd0682b00c440615544fb0f4e10a4989120a7f257b88fc`; the downloaded archive digest matched the GitHub artifact digest. Its `audit.json` reports:

- `status=PASS` across gamma values -2, -1, -0.5, 0, 0.5, 1, 2.
- 20/20 successful replicates and 20/20 child raw audits for every gamma; zero failed replicates.
- 244 manifest files verified per gamma; 1,708 total.
- 91,140 reconstructed metric rows; 4,557 density-summary rows.
- 2,205 paired IAE contrast rows and 756 paired gamma-versus-zero contrast rows.
- Paired seed schedule and parent/validation redshifts verified.

The archive contains the audit JSON and three CSV summaries. The paired contrast CSV is at the individual absolute-bandwidth level (not merely the seven-bandwidth average used in the main text). These results complete the computational density-grid gate only. They do not establish H0 posterior performance, coefficient-posterior calibration, SBC, PIT, or fixed-truth coverage.

