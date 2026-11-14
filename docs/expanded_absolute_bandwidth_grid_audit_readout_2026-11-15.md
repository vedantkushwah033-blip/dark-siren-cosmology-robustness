# Expanded physical-bandwidth density pilot: independent audit readout

**Date:** 2026-11-15  
**Experiment:** 20 paired replicates × 3 catalogue-floor scenarios  
**Stage:** synthetic density recovery only  
**Status:** all three expanded-grid artifacts independently audited; exploratory bandwidth sensitivity, not a bandwidth-selection result.

## Audit provenance

The expanded-grid workflow reran the same 20 paired replicates for each of the three catalogue-floor scenarios on an 11-point physical bandwidth grid:
`0.006, 0.008, 0.010, 0.012, 0.016, 0.024, 0.032, 0.040, 0.050, 0.060, 0.080`.

All three workflow jobs completed successfully. The audit JSONs independently report, for each scenario:

- 20 planned and 20 successful replicates; zero failed replicates;
- 20/20 independent child raw-input audits passed;
- 244 manifest files verified;
- 20,460 aggregate metric rows reconstructed;
- 1,740 realized-Scott-width rows reconstructed;
- overall audit status `PASS`.

The three realized-Scott-width ranges were:
- no floor: 0.01217–0.03302;
- floor 0.10: 0.01325–0.02741;
- floor 0.25: 0.01356–0.02444.

A manifest/audit PASS establishes that the declared simulation artifacts and aggregations were reconstructed consistently. It does not establish external validity or scientific adequacy.

## Mean IAE at a common physical bandwidth

At validation size (n=1000) and (h=0.016), mean IAE across 20 paired replicates was:

| Arm | No floor | Floor 0.10 | Floor 0.25 |
|---|---:|---:|---:|
| A catalogue-only | 0.53049 | 0.37809 | 0.21349 |
| B parent oracle | 0.18378 | 0.18378 | 0.18378 |
| C known-probability IPW | 0.20391 | 0.18653 | 0.18065 |
| D MAP estimated IPW | 0.20720 | 0.29618 | 0.27331 |
| E posterior-mean density | 0.21032 | 0.30006 | 0.27605 |

Paired replicate-level contrasts at this same bandwidth (difference is first arm minus second arm; MCSE is across the 20 paired differences):

| Floor | Contrast | Mean difference | Paired MCSE |
|---|---|---:|---:|
| No floor | D − A | -0.32330 | 0.01145 |
| No floor | D − C | +0.00329 | 0.01189 |
| No floor | C − B | +0.02013 | 0.00654 |
| No floor | E − D | +0.00312 | 0.00035 |
| 0.10 | D − A | -0.08192 | 0.01470 |
| 0.10 | D − C | +0.10964 | 0.01359 |
| 0.10 | C − B | +0.00275 | 0.00401 |
| 0.10 | E − D | +0.00388 | 0.00012 |
| 0.25 | D − A | +0.05983 | 0.00907 |
| 0.25 | D − C | +0.09266 | 0.01012 |
| 0.25 | C − B | -0.00313 | 0.00256 |
| 0.25 | E − D | +0.00274 | 0.00008 |

These contrasts show that the estimated-IPW arm's relative performance depends strongly on the catalogue-selection scenario. In the no-floor scenario it substantially improves over catalogue-only and is close to known-probability IPW at (h=0.016). At floor 0.10 it still improves over catalogue-only, but remains substantially worse than known-probability IPW. At floor 0.25 it is worse than catalogue-only at this bandwidth. Known-probability IPW remains close to the parent oracle in all three scenarios.

The very small E−D MCSEs should not be read as evidence of broad posterior uncertainty calibration; they quantify only the paired density-loss difference in this finite simulation design.

## Bandwidth sensitivity: minima are not a selection result

At (n=1000), the empirical minimum mean IAE over the 11-point grid occurred at:

| Scenario | Arm A | Arm B | Arm C | Arm D | Arm E |
|---|---:|---:|---:|---:|---:|
| No floor | 0.080 (edge) | 0.016 | 0.016 | 0.016 | 0.016 |
| Floor 0.10 | 0.080 (edge) | 0.016 | 0.016 | 0.016 | 0.012 |
| Floor 0.25 | 0.012 | 0.016 | 0.016 | 0.012 | 0.012 |

The no-floor and floor-0.10 catalogue-only minima occur at the largest tested bandwidth (0.080), so those minima are not bracketed by the grid. They cannot be treated as identified optima. Several estimated-IPW curves are also fairly flat near their minima. Choosing a bandwidth by taking the smallest observed mean from these same 20 replicates would be an optimistic, data-reused procedure.

Therefore this experiment does **not** select a production bandwidth. A valid next design needs a declared training/holdout split (or independent tuning and evaluation seed sets), a bandwidth grid wide enough to bracket the objective where possible, and a frozen selection rule applied identically across arms and scenarios.

## Interpretation and limitations

1. The density target is the normalized latent parent intensity multiplied by the declared merger/host weight. This is a synthetic benchmark target.
2. Arm C uses truth-only inclusion probabilities and remains an oracle diagnostic.
3. Arms D/E fit the declared degree-one logistic inclusion model on independent validation outcomes.
4. Catalogue inclusion, host/merger weighting, and GW detection selection are separate mechanisms.
5. These are 20-replicate synthetic density results. They are not an observational (H_0) result, a real galaxy-survey validation, or an LVK/GWTC analysis.
6. Density-to-(H_0) propagation remains gated on a defensible bandwidth protocol, held-out validation, and the separate selected-likelihood/calibration audits.

## Decision

**Density-stage status: expanded-grid audit PASS; bandwidth selection NOT YET JUSTIFIED; proceed to a separate preregistered training/holdout bandwidth study, not H0 propagation.**
