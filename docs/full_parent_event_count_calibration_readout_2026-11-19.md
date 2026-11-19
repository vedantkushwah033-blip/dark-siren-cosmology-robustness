# Paired full-parent event-count calibration: preliminary readout

**Status:** 500 paired parent/event realizations completed; independent raw-seed reconstruction audit pending.

## Provenance

- Source PR: #101
- Workflow: #44434616291
- Artifact: #14487983108
- Artifact SHA-256: `306680d5547d3a2c793585334c44d609fdaa577836afa32c17a20887891d754d`
- 500 independent parent/event realizations; 1,000 posterior rows; all SUCCESS.
- Each 30-event realization contains the corresponding four-event realization as its first four detected events; tests and the workflow's paired-denominator assertions passed.
- Parent size 1,200; H0 truth 70 km s^-1 Mpc^-1; uniform H0 prior [20,300]; 561-point grid.
- Both arms use the full finite parent population and the same selected-event likelihood. There is no catalogue thinning or inclusion-probability fitting in this experiment.

## Results

| Detected events | Mean posterior-median bias | MCSE | RMSE | 68% coverage | 90% coverage | 95% coverage | Mean 68% interval width |
|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | +5.060 | 0.717 | 16.793 | 0.696 | 0.910 | 0.954 | 23.436 |
| 30 | +0.866 | 0.182 | 4.160 | 0.662 | 0.860 | 0.920 | 7.355 |

All biases and interval widths are in km s^-1 Mpc^-1.

The paired 30-event minus nested 4-event posterior-median bias difference is −4.194 (MCSE 0.699; paired 95% t interval [−5.567, −2.820]). Mean 68% interval width contracts from 23.436 to 7.355 km s^-1 Mpc^-1.

The 4-event arm has a positive mean median bias, but its coverage in this 500-replicate experiment is reasonably close to nominal. The 30-event arm has much smaller median bias and RMSE and narrower intervals; its coverage is somewhat below nominal at all three levels. These results support a strong finite-event-count effect in this synthetic full-parent setting, but do not establish joint SBC.

## Reconciliation with the earlier 100-parent oracle readout

The earlier 100-replicate catalogue-inclusion experiment's full-parent oracle showed mean median bias +7.048 (MCSE 1.809) and 68/90/95% coverage 0.54/0.83/0.91 across its 100 fixed parent/event realizations. The new 500-parent event-count experiment's 4-event arm gives +5.060 bias and 0.696/0.910/0.954 coverage.

The difference in coverage is material and should not be hidden or averaged away. A separate 20,000-draw experiment has now completed: 200 fresh event realizations per each of the same 100 archived parent populations (workflow #44435303396, artifact #14489935029). It found:
- 68% coverage 0.6535 (event-level Wilson 95% [0.6469, 0.6601]; parent-cluster t 95% [0.6477, 0.6593]);
- 90% coverage 0.8794 (parent-cluster t 95% [0.8743, 0.8845]);
- 95% coverage 0.9335 (parent-cluster t 95% [0.9298, 0.9372]);
- mean median bias +6.4388 (parent-cluster MCSE 0.1408; t 95% [6.1595, 6.7181]).

This repeated-draw result shows that the original 100 single-event-realization coverage (0.54/0.83/0.91) was unusually low; the same parent populations give 0.6535/0.8794/0.9335 when event-sample variability is averaged over 200 draws per parent. The independent 500-parent experiment gives 0.696/0.910/0.954, somewhat higher but broadly consistent with a modest finite-sample difference. These are distinct designs and must not be pooled or presented as paired estimates. The 20,000-draw result is the primary conditional-on-these-parent-populations baseline; the 500-parent result estimates unconditional performance over a new parent-population sample.

## Independent audit pending

The companion PR #104 regenerates all 500 parent populations and 30-event sequences from the archived seeds, independently reconstructs all 1,000 selected-event posteriors, and checks every reported quantile, coverage indicator, interval width, aggregate metric, paired contrast, and manifest hash. This readout should be upgraded to independently audited status only after that workflow passes.

## Scope

This is a synthetic low-redshift Euclidean fixed-truth diagnostic with a full parent catalogue. It is not an observational H0 measurement, joint SBC, a catalogue-inclusion result, LVK/GWTC analysis, or validation against real galaxy surveys.
