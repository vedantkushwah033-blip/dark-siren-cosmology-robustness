# Independent wide-prior five-arm readout — 2026-11-20

## Scope

This is a descriptive readout of the 100-replicate, 30-detected-event synthetic five-arm experiment using the expanded (H_0) prior support ([30,110]) km s(^{-1}) Mpc(^{-1}). It is not an observational (H_0) analysis, a full SBC study, or evidence of universal estimator superiority.

Source: GitHub Actions run [44483923518](https://github.com/vedantkushwah033-blip/dark-siren-cosmology-robustness/actions/runs/44483923518), master seed 2026111940, artifact 14527343221, SHA-256 `aeae71bb5088cabdac75bce060b9afcd3964d7d73945aa7d3b4f6390d490d137`.

The run completed 4,500/4,500 lifecycle cells successfully. The independent orchestration audit reported PASS, 300/300 child runs, and 5,708 aggregate manifest entries verified. No planned cells were failed or non-evaluable.

## Results at validation size 500

Bias is the mean posterior-median error (posterior median minus the fixed truth); RMSE is across the 100 replicate posterior medians. Coverage is empirical 68% interval coverage.

| Curvature γ | Arm | Bias | RMSE | 68% coverage |
|---:|---|---:|---:|---:|
| -1 | Catalogue-only A | -17.047 | 17.503 | 0.00 |
| -1 | Full-parent oracle B | +0.818 | 3.408 | 0.71 |
| -1 | Known-π IPW C | +0.219 | 4.383 | 0.54 |
| -1 | Estimated MAP IPW D | -2.713 | 4.659 | 0.56 |
| -1 | Integrated estimated IPW E | -2.487 | 4.539 | 0.58 |
| 0 | Catalogue-only A | -16.585 | 17.203 | 0.00 |
| 0 | Full-parent oracle B | +0.818 | 3.408 | 0.71 |
| 0 | Known-π IPW C | +0.621 | 4.002 | 0.63 |
| 0 | Estimated MAP IPW D | +0.581 | 4.050 | 0.58 |
| 0 | Integrated estimated IPW E | +0.732 | 4.105 | 0.64 |
| +1 | Catalogue-only A | -15.665 | 16.534 | 0.00 |
| +1 | Full-parent oracle B | +0.818 | 3.408 | 0.71 |
| +1 | Known-π IPW C | +0.647 | 3.819 | 0.64 |
| +1 | Estimated MAP IPW D | +3.344 | 5.372 | 0.45 |
| +1 | Integrated estimated IPW E | +3.343 | 5.378 | 0.49 |

## Validation-size pattern

For A, B, and C, the reported bias/RMSE/coverage are identical across validation sizes 100/300/500 in this design; these arms do not estimate the inclusion model from the validation sample. D and E vary with validation size.

- At γ=-1, E RMSE decreases from 4.878 (n=100) to 4.487 (n=300) to 4.539 (n=500); its bias remains near -2.4 to -2.5.
- At γ=0, E RMSE is 4.626, 4.194, and 4.105; bias is +0.480, +0.783, and +0.732.
- At γ=+1, E RMSE is 5.453, 5.447, and 5.378; bias remains about +2.7 to +3.3.

The wider prior exposes a severe catalogue-only failure in this seed stream: all three curvature conditions have zero 68% coverage and large negative median bias. In contrast, the parent-oracle benchmark remains close to the fixed truth. At positive curvature, both estimated-IPW arms remain materially worse than known-π IPW C in RMSE and 68% coverage. Integrating over the fitted degree-one coefficient uncertainty does not remove the misspecification pattern.

## Interpretation limits

1. This is one 100-replicate seed stream; report Monte Carlo uncertainty before making fine distinctions.
2. These figures are the runner's success-only aggregate summary. The independent orchestration audit verifies lifecycle, child attestations, and artifact manifests; it is not an independent reconstruction of every posterior summary.
3. This stream (master seed 2026111940) is independent of the paired prior-support stream (master seed 2026111929). Its first five replicate seeds overlap the earlier five-replicate engineering smoke; do not count that smoke as additional independent replication.
4. The paired ([50,90]) versus ([30,110]) analysis must use the seed-2026111929 narrow/wide pair and verify identical raw inputs, validation samples, weights, bandwidths, and common-support likelihoods before attributing posterior changes to prior support.
5. Do not pool this wide-prior run with the earlier narrow-prior runs as if they shared a prior. No observational or universal estimator-ranking claim follows.

## Next audit gate

Run `scripts/analyze_paired_h0_prior_sensitivity.py` on the matching narrow/wide seed-2026111929 artifacts. Then independently reconstruct posterior summaries from the immutable wide-prior artifact and add paired Monte Carlo intervals. Keep all pull requests draft and unmerged until those audits pass.
