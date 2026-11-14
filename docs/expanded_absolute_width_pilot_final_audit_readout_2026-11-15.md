# Expanded absolute-bandwidth density pilot: independent audit readout

**Date:** 2026-11-15  
**Branch:** `research/absolute-width-distribution-pilot-20261114`  
**Experiment run:** GitHub Actions workflow `43617052462`, exact run head `d94fdd62e514c81c20444196e166d2ebbb341a40`  
**Scope:** density-only synthetic benchmark; no bandwidth freeze and no H0 inference.

## Audit completion

The expanded 11-width run is now complete for all three catalogue-floor conditions. Each independent outer/child audit reports:

| Catalogue floor | Replicates successful | Child audits passed | Manifest files verified | Aggregate metric rows reconstructed | Scott-width rows reconstructed | Audit |
|---|---:|---:|---:|---:|---:|---|
| 0.00 | 20/20 | 20 | 244 | 20,460 | 1,740 | PASS |
| 0.10 | 20/20 | 20 | 244 | 20,460 | 1,740 | PASS |
| 0.25 | 20/20 | 20 | 244 | 20,460 | 1,740 | PASS |

No failed replicates were recorded in these expanded-grid artifacts. The audit scope includes seed/ledger/manifest/aggregation checks and independent raw-input reconstruction for every successful child. This is an internal synthetic-artifact audit, not an external validation.

The grid is **h = 0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080** in the physical redshift coordinate. The audit reconstructed realized Scott widths from 0.01217–0.03302 (floor 0.00), 0.01325–0.02741 (floor 0.10), and 0.01356–0.02444 (floor 0.25).

## Fixed-bandwidth comparison (n=1000, h=0.016)

Mean integrated absolute error (IAE), over 20 paired replicates:

| Arm | Floor 0.00 | Floor 0.10 | Floor 0.25 |
|---|---:|---:|---:|
| A catalogue-only | 0.53049 | 0.37809 | 0.21349 |
| B parent oracle | 0.18378 | 0.18378 | 0.18378 |
| C known-probability IPW | 0.20391 | 0.18653 | 0.18065 |
| D MAP estimated IPW | 0.20720 | 0.29618 | 0.27331 |
| E posterior-mean density | 0.21032 | 0.30006 | 0.27605 |

Paired mean contrasts (first arm minus second; MCSE of paired mean):

| Floor | Contrast | Mean | MCSE |
|---|---|---:|---:|
| 0.00 | D − A | −0.32330 | 0.01145 |
| 0.00 | D − C | +0.00329 | 0.01189 |
| 0.00 | C − B | +0.02013 | 0.00654 |
| 0.00 | E − D | +0.00312 | 0.00035 |
| 0.10 | D − A | −0.08192 | 0.01470 |
| 0.10 | D − C | +0.10964 | 0.01359 |
| 0.10 | C − B | +0.00275 | 0.00401 |
| 0.10 | E − D | +0.00388 | 0.00012 |
| 0.25 | D − A | +0.05983 | 0.00907 |
| 0.25 | D − C | +0.09266 | 0.01012 |
| 0.25 | C − B | −0.00313 | 0.00256 |
| 0.25 | E − D | +0.00274 | 0.00008 |

This fixed h=0.016 comparison reproduces the corrected five-width readout. It is useful as a common-bandwidth comparison but is not itself a bandwidth-selection result.

## What the expanded grid says—and does not say

The 11-width surface shows materially different smoothing behavior across arms. At n=1000, the lowest *observed mean across the same 20 pilot replicates* occurred at:

| Floor | Arm A | Arm B | Arm C | Arm D | Arm E |
|---|---:|---:|---:|---:|---:|
| 0.00 | 0.080 (0.4517) | 0.016 (0.1838) | 0.016 (0.2039) | 0.016 (0.2072) | 0.016 (0.2103) |
| 0.10 | 0.080 (0.3643) | 0.016 (0.1838) | 0.016 (0.1865) | 0.016 (0.2962) | 0.012 (0.3000) |
| 0.25 | 0.012 (0.2127) | 0.016 (0.1838) | 0.016 (0.1806) | 0.012 (0.2732) | 0.012 (0.2759) |

Each cell is bandwidth (mean IAE). These minima are **descriptive only**: selecting the minimum after examining the same pilot replicates would be optimistic and would leak evaluation information into tuning. In particular, the catalogue-only arm's broad-bandwidth minima at 0.080 are not evidence that 0.080 is a valid common bandwidth for all arms. A common-bandwidth training/holdout protocol must be frozen and run on independent seeds before choosing a production bandwidth.

The broad-bandwidth catalogue-only improvement also demonstrates why a single bandwidth comparison can mislead: arm-specific density roughness and bias/variance tradeoffs differ. The fixed h=0.016 results remain useful for paired comparisons, but should not be promoted to a final estimator ranking without a prespecified bandwidth protocol.

## Interpretation

1. Under the correctly specified floor-0.00 condition, estimated IPW is substantially better than catalogue-only at h=0.016 and close to known-probability IPW in this pilot.
2. Under floors 0.10 and 0.25, estimated IPW remains better than catalogue-only at floor 0.10 but is worse at floor 0.25 at h=0.016; it also falls well behind known-probability IPW in both stress conditions.
3. The known-probability IPW arm remains close to the parent-oracle arm across floors. This isolates the gap to estimating/misspecifying inclusion probabilities rather than a generic inability of inverse-probability weighting to recover the parent density.
4. Arm E's posterior-mean density loss is close to but consistently above the MAP plug-in loss in this pilot. The posterior expected node-specific loss remains a distinct secondary estimand and must not be conflated with loss of the posterior-mean density.
5. These results are exploratory synthetic evidence only: 20 replicates per floor, one parent-population setup, and no held-out bandwidth selection. They do not establish universal estimator superiority, survey validity, or observational H0 performance.

## Next gate

Before extending to the seven-gamma quadratic-logit curvature design or H0 propagation:

- freeze the common-bandwidth training/holdout split and selection objective;
- ensure all five arms are scored on the same independent holdout seeds and physical bandwidth grid;
- retain the independent raw-record auditor and complete lifecycle ledger;
- report paired contrasts with Monte Carlo uncertainty and sensitivity to the selected bandwidth;
- only then proceed to density-to-H0 propagation, followed by separate H0-grid, quadrature, and calibration gates.

No bandwidth is selected by this document. No H0 inference is reported.
