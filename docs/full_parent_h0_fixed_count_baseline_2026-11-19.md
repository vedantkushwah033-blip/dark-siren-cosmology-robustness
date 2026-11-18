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
| Monte Carlo SE of mean bias | 0.414 |
| Paired t 95% interval for mean bias | [+5.61, +7.23] |
| Mean 68% interval width | 24.20 km s\(^{-1}\) Mpc\(^{-1}\) |
| 68% interval coverage | 67.2% (Wilson 95%: 65.1%–69.2%) |
| 90% interval coverage | 88.4% (Wilson 95%: 86.9%–89.7%) |
| 95% interval coverage | 94.3% (Wilson 95%: 93.2%–95.2%) |

Mean posterior-median bias was positive within each of the 10 parent populations (approximately +4.9 to +9.9 km s\(^{-1}\) Mpc\(^{-1}\)). The aggregate 68% and 95% coverage estimates are close to nominal at this Monte Carlo precision; the 90% interval coverage is modestly below nominal.

## Interpretation

This establishes an important baseline for the project: positive posterior-median bias is present even with the complete finite parent population and exact host weights. Therefore the positive bias in the 90-cell pilot cannot be assigned wholly to finite-validation inclusion estimation or logistic curvature misspecification.

The result does **not** establish that the likelihood is wrong. With only four detected events, posterior medians can have finite-sample bias even when the selected-event likelihood is matched to the generator; fixed-truth credible-interval coverage is also not identical to joint SBC. The observed 90% undercoverage is modest and should be rechecked with more parent realizations and a separate calibration design.

The follow-up 600-cell experiment in PR #89 will estimate fixed-truth coverage under estimated inclusion. A subsequent oracle-arm reconstruction on that artifact should compare estimated-IPW, known-\(\pi\) IPW, and full-parent posteriors under the same event realizations. This is necessary to distinguish the event-sample baseline from catalogue finite-sample and estimated-inclusion effects.

## Reproduction

The independent runner is `scripts/run_full_parent_h0_fixed_count_baseline.py`. It verifies the source artifact's aggregate manifest, generates conditional detected events, computes full-parent posteriors, and writes replicate-level and parent-level CSVs plus a checksum manifest. The 2,000-draw numerical result above was obtained with the same conditional sampler and likelihood calculation; a small one-event-set-per-parent smoke run of the committed CLI completed, while a full CLI rerun remains a longer computational task.
