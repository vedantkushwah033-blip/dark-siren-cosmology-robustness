# Five-replicate prior-predictive five-arm smoke — 2026-11-20

## Execution and independent audit

Workflow `44723390611` completed successfully; artifact `14622399730` (five-arm-prior-predictive-smoke), SHA-256 `7580fad91233d15a09abd288ca64fba7cbdddbe7eb0d82c12bbd9e7c099f2797`.

- 225/225 lifecycle cells SUCCESS across 15 child runs.
- Independent aggregate audit PASS.
- All 293 aggregate manifest entries independently rechecked against byte counts and SHA-256 digests; no mismatch.
- All 15 child audit reports PASS; each verifies 14 child manifest entries and reconstructs 15 posterior summary rows.
- 45 aggregate metric rows and 360 paired contrast rows reconstructed.
- 15 cross-curvature paired-validation blocks verified.
- Five unique prior-predictive H0 truths were shared across curvature conditions within replicate.

This is a five-replicate engineering smoke, not a calibration study.

## Selected results at validation n=500

| Curvature | Method | Bias | RMSE | 68% coverage |
|---:|---|---:|---:|---:|
| -1 | Parent oracle | +3.41 | 6.04 | 80% |
| -1 | Known-probability IPW | +0.07 | 5.95 | 60% |
| -1 | MAP degree-one IPW | -1.78 | 5.39 | 40% |
| -1 | Integrated degree-one IPW | -1.71 | 5.27 | 60% |
| 0 | Parent oracle | +3.41 | 6.04 | 80% |
| 0 | Known-probability IPW | -0.03 | 5.35 | 40% |
| 0 | MAP degree-one IPW | +0.17 | 5.02 | 40% |
| 0 | Integrated degree-one IPW | +0.22 | 4.92 | 40% |
| +1 | Parent oracle | +3.41 | 6.04 | 80% |
| +1 | Known-probability IPW | +0.37 | 5.18 | 60% |
| +1 | MAP degree-one IPW | +2.78 | 6.05 | 20% |
| +1 | Integrated degree-one IPW | +2.73 | 5.88 | 20% |

Catalogue-only posterior medians are strongly biased in this design, as expected when catalogue incompleteness is not corrected. The estimated degree-one arms show curvature-dependent behavior. The five-replicate estimates are too noisy to rank the methods or establish calibration.

## Limits

- Five prior-predictive truths yield very coarse rank and coverage resolution.
- These results are not directly comparable to the 50-replicate run as a confirmatory update; this is a small engineering smoke with its own frozen seed schedule.
- The degree-one inclusion model remains misspecified at nonzero curvature.
- No observational H0 or real-survey inference is involved.

The appropriate next comparison is the prior-matched centered D1/D2 experiment now running on PR #133, with weight-instability diagnostics reported alongside posterior error and coverage.


## Same-seed workflow rerun and reproducibility check

A later PR-head rerun, workflow `44726607857`, produced artifact `14624961566` (ZIP SHA-256 `293cb6d2ec5215f87b863fd12f848543e5553dc36c5a29577c413ad3adc7afc8`). It used the **same** master seed, analysis seeds, validation seeds, and five prior-predictive truth draws as workflow `44723390611`; it is therefore a rerun, not an additional independent five-replicate stream.

I independently verified all 293 outer-manifest entries and all 15 nested child manifests (14 files each); zero hash/size mismatches. The independent aggregate audit reports 225/225 SUCCESS lifecycle cells, 15/15 child audits, 45 metric rows and 360 paired contrasts reconstructed. The maximum child likelihood reconstruction discrepancy was (2.27\\times10^{-13}).

Comparing the two same-seed artifacts:
- Maximum numeric difference in aggregate metrics: (2.49\\times10^{-14}).
- Maximum numeric difference in replicate posterior summaries: (1.42\\times10^{-13}).
- Maximum numeric difference in paired contrasts: (3.53\\times10^{-12}).
- Weight-diagnostic records were identical.

The results are numerically reproducible to floating-point precision under the rerun. They do not increase the effective replicate count beyond five and must not be pooled as ten independent replicates.
