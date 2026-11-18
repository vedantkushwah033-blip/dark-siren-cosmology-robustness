# Full-parent fixed-count H0 calibration baseline

**Status:** exploratory repeated-sampling diagnostic; 2,000 fresh event realizations evaluated across 10 archived finite-parent populations.

## Question

The 90-cell estimated-inclusion pilot showed positive posterior-median bias even at \(\gamma=0\). Is that entire bias caused by estimating catalogue inclusion probabilities, or does it also appear when the exact finite parent population is known?

## Design

- Reuse the 10 distinct parent-galaxy realizations archived in the corrected 90-cell pilot (one parent per replicate; each parent has 1,200 galaxies).
- For each parent, draw 200 new independent event sets, each containing exactly four detected events: 2,000 event realizations total.
- Fix generating \(H_0=70\;\mathrm{km\,s^{-1}\,Mpc^{-1}}\).
- Use the exact parent host weights \(w_g\propto 1+1.5(z_g-z_{\min})/(z_{\max}-z_{\min})\), with no catalogue thinning or inclusion estimation in the inference.
- Use a uniform H0 prior on [20, 300] and a 281-point grid (spacing 1).
- Seed 990000, with a separate deterministic stream for each archived parent replicate.

The event sampler is exactly conditional on the one-sided detection rule. First draw a host from

\[
P(g\mid A=1,H_0)\propto
w_g\,\Phi\!\left(\frac{cz_g/H_0-d_{\rm th}}{\sigma}\right),
\]

then draw its measured distance from the corresponding Gaussian truncated below at \(d_{\rm th}\). This is distributionally equivalent to the original rejection sampler but avoids generating rejected candidates.

The posterior is independently reconstructed from the finite-parent selected-event likelihood, including the same detection-efficiency denominator. No catalogue-inclusion correction is present in this experiment.

## Results

| Quantity | Estimate |
|---|---:|
| Parent populations | 10 |
| Event realizations per parent | 200 |
| Total event realizations | 2,000 |
| Mean posterior-median bias | +6.42 km s\(^{-1}\) Mpc\(^{-1}\) |
| Median of replicate median-bias values | +2.11 km s\(^{-1}\) Mpc\(^{-1}\) |
| RMSE of posterior-median estimator | 19.58 km s\(^{-1}\) Mpc\(^{-1}\) |
| Mean absolute posterior-median error | 11.48 km s\(^{-1}\) Mpc\(^{-1}\) |
| Fraction of replicate posterior medians above truth | 57.9% |
| Monte Carlo SE of mean bias | 0.414 |
| Paired t 95% interval for mean bias | [+5.61, +7.23] |
| Mean 68% interval width | 24.20 km s\(^{-1}\) Mpc\(^{-1}\) |
| 68% interval coverage | 67.2% (Wilson 95%: 65.1%–69.2%) |
| 90% interval coverage | 88.4% (Wilson 95%: 86.9%–89.7%) |
| 95% interval coverage | 94.3% (Wilson 95%: 93.2%–95.2%) |

Mean posterior-median bias was positive within each of the 10 parent populations (approximately +4.9 to +9.9 km s\(^{-1}\) Mpc\(^{-1}\)). There are two uncertainty levels. The reported event-draw MCSE (0.414) and t interval [+5.61,+7.23] quantify event-realization Monte Carlo uncertainty **conditional on these 10 fixed parent populations**. I also recomputed uncertainty using the ten parent-level means as independent clusters: the between-parent SD is 1.453, the parent-level SE is 0.459, and the parent-cluster t 95% interval for mean bias is [+5.38,+7.46]. A parent-cluster percentile bootstrap (200,000 resamples; seed 12345) gives [+5.69,+7.37]. For coverage, parent-cluster t 95% intervals are [64.2%,70.2%] for 68% coverage, [86.4%,90.4%] for 90% coverage, and [92.7%,95.9%] for 95% coverage. The parent-level intervals are more appropriate when generalizing beyond the ten fixed parent populations, but remain uncertain with only ten clusters. The aggregate 68% and 95% coverage estimates are close to nominal at this conditional Monte Carlo precision; the 90% interval coverage is modestly below nominal.

## Interpretation

This establishes an important conditional baseline for the project: positive posterior-median bias is present across the ten archived finite parent populations even with the complete finite parent population and exact host weights. Therefore the positive bias in the 90-cell pilot cannot be assigned wholly to finite-validation inclusion estimation or logistic curvature misspecification.

The result does **not** establish that the likelihood is wrong. With only four detected events, posterior medians can have finite-sample bias even when the selected-event likelihood is matched to the generator; fixed-truth credible-interval coverage is also not identical to joint SBC. The observed 90% undercoverage is modest and should be rechecked with more parent realizations and a separate calibration design.

The follow-up 600-cell experiment in PR #89 will estimate fixed-truth coverage under estimated inclusion. A subsequent oracle-arm reconstruction on that artifact should compare estimated-IPW, known-\(\pi\) IPW, and full-parent posteriors under the same event realizations. This is necessary to distinguish the event-sample baseline from catalogue finite-sample and estimated-inclusion effects.

## Reproduction

The independent runner is `scripts/run_full_parent_h0_fixed_count_baseline.py`. It verifies the source artifact's aggregate manifest, generates conditional detected events, computes full-parent posteriors, and writes replicate-level and parent-level CSVs plus a checksum manifest. The 2,000-draw numerical result above was obtained with the same conditional sampler and likelihood calculation. The committed CLI was smoke-tested at 10 event sets per parent; its first 100 posterior medians matched the corresponding rows from the 2,000-draw analysis to a maximum absolute difference below \(5\times10^{-14}\). The full 200-draw-per-parent CLI rerun is now complete in GitHub Actions workflow `44417807075` on commit `38a5c94f272f711d4a499ba3cc71a0e8e0a17912`. The workflow downloaded the corrected 90-cell source artifact `14476042847`, reran all 2,000 event realizations, independently reconstructed the aggregate summary and output manifest, and completed successfully. The uploaded result is artifact `14479644393` (`full-parent-fixed-count-h0-baseline-2000`), with archive SHA-256 `e0c702607ff632c71c6268f34fe247519c6c2c8eb8b92c37f35a6a8fd672f7e6`; the downloaded archive independently matched this digest. Its summary exactly reproduces the reported mean median bias (+6.4203), event-draw MCSE (0.41374), 68/90/95% coverage (0.672/0.884/0.943), and mean 68% interval width (24.1983). Recomputing from the ten parent-level means gives cluster SE 0.45938 and t 95% interval [+5.3811,+7.4595]; parent-cluster coverage intervals are [0.6424,0.7016], [0.8637,0.9043], and [0.9268,0.9592]. These cluster intervals are based on only ten parent populations. The workflow-side audit verified 2,000 rows, 10 parent populations, 200 event draws per parent, quantile ordering, summary arithmetic, and every output-manifest hash.
