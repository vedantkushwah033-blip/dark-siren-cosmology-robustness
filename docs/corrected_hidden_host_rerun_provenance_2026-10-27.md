# Corrected hidden-host rerun provenance — 2026-10-27

## Status

The four hidden-host experiments identified in the missing-host mixture likelihood impact audit have now been rerun on corrected branch revision `ecb260886716c54fdd2d5e47390701949f5a18d5`. Each workflow completed successfully and ran its independent reconstruction auditor before uploading the artifact. The repository's general numerical-diagnostics workflow also passed on Python 3.10, 3.11, and 3.12 at this revision.

This closes the **rerun-and-implementation-audit gate** for these four named experiments. I downloaded the historical and corrected artifacts and compared every CSV row and field. All categorical/string fields match exactly; summary and contrast tables agree to floating-point precision. The maximum absolute difference in any raw per-run numeric field was 8.74e-12 (hidden-host assumption), 6.11e-13 (sensitivity surface), 2.37e-12 (prior-support factorial), and 1.85e-13 (paired host-misspecification). The differences are negligible floating-point changes and the published tables agree at displayed precision. This does not establish that the toy model is astrophysically realistic, that its scientific assumptions are correct, or that the research question is novel.

## Corrected artifacts

All runs below use source revision `ecb260886716c54fdd2d5e47390701949f5a18d5`, branch `research/distance-uncertainty-likelihoods`.

| Experiment | GitHub Actions run | Artifact ID | Output design verified in workflow | Downloaded ZIP SHA-256 |
|---|---:|---:|---|---|
| Hidden-host fraction/redshift assumptions | [40417072147](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072147) | 12696226361 | 105,600 posterior rows; 528 summaries; 2,880 paired-contrast rows | `ee684f9ee3b59717f62c6149ea15e05cb0ea3bf2ff1703294e0b7006dc4b9ff8` |
| Host-misspecification sensitivity surface | [40417072115](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072115) | 12696905126 | 86,400 posterior rows; 432 summaries; 1,296 paired-contrast rows | `c2c2c485d6b5357ffe58c1540a013b0da3cbdbada48c653094d08f1ce1974974` |
| Prior-support × hidden-host factorial | [40417072086](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072086) | 12695694723 | 86,400 posterior rows; 432 summaries; 1,296 arm contrasts; 1,728 prior contrasts | `0e82705ebc2cb6808bb320a0b3e7e9b94a4d8fa2bcb9dd51cf567366fdcc4fab` |
| Paired host-misspecification experiment | [40417072077](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/40417072077) | 12695609236 | 4,000 posterior rows; 4 summaries; 12 paired-contrast rows | `ab680e0de5256397f76513f320a99b81f40fcc5ca11c9efa4dabf5b4aa779743` |

The row counts above are design/shape checks, not independent evidence by themselves. The workflow jobs also invoke the corresponding audit scripts, which regenerate seeded observations and reconstruct posterior summaries and contrasts. The run records show successful job conclusions.

## Historical-artifact comparison

The historical ZIP hashes were `94b3cb3fa5167a680eeeeb52c18180ed244b04bf7178947be6911aecbc3e4356` (hidden-host assumption), `452cc37bab97fb95605d38a9e1eff1643e8cb6d66a417ddb107af7beff56bd67` (sensitivity surface), `61f78104305fc7d3fb4992e0841b7aa8addc9b27864bd090a6bf719243279c16` (prior-support factorial), and `abc6f451311b921a0370e94356848a65734ecfb554f8443b43b740bf2b9f0d47` (paired host-misspecification). All were compared row-by-row with the corrected artifacts listed above. The maximum raw numeric differences are caused by floating-point evaluation and do not change any summary or contrast beyond the reported precision. The earlier source revisions for these runs contain the linear-space implementation (same inference source blob across the historical revisions), predating the transient extra-factor bug in the later log-space rewrite. The reruns therefore confirm the historical numerical tables rather than revealing material invalidation.

## Interpretation boundary

- The likelihood correction is a model-correctness change: the intended catalogue/missing-population mixture is
  `L(H0) = (1 - f_miss) L_cat(H0) + f_miss L_miss(H0)`.
- The explicit regression test compares the direct API against an explicitly concatenated full mixture at missing fractions 0.05, 0.35, and 0.80.
- Reconstructed agreement establishes that outputs match the declared implementation and seed. It does not establish that the generator or likelihood represents real gravitational-wave observations, detector selection, or galaxy-catalogue structure.
- The prior-support factorial is especially important to interpret cautiously: changing prior support can change posterior bias, RMSE, width, and empirical coverage simultaneously. A broader interval is not automatically a better or calibrated inference.
- External CHIMERA/reference-pipeline validation remains blocked until the required benchmark inputs and their provenance are verified. Do not describe the external benchmark as completed.

## Next release gate

1. Update the four experiment reports with corrected-revision provenance and the row-by-row artifact comparison; preserve historical artifact hashes for traceability.
2. Verify the downstream caller inventory, including wrappers that may call the likelihood indirectly.
3. Run the targeted novelty review against a documented completeness-aware baseline before making any originality or publication-readiness claim.
4. Keep the project framed as a synthetic robustness/audit study unless and until realistic benchmark validation supports a stronger claim.
