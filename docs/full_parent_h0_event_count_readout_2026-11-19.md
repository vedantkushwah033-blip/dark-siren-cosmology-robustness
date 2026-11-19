# Full-parent H0 calibration by detected-event count

**Source run:** `44434616291`  
**Source artifact:** `14487983108`  
**Independent audit run:** `44441848104`  
**Independent audit artifact:** `14491418583`

## Design and integrity

- Synthetic low-redshift full-parent selected-likelihood fixed-truth experiment.
- 500 independent parent/event realizations, each with 1,200 parent galaxies.
- True H0 = 70 km s^-1 Mpc^-1; uniform analysis prior [20,300].
- Compare 4 detected events against 30 detected events, where the four-event sample is a nested prefix of the same 30-event realization.
- 561-point H0 grid.
- The four-entry source manifest was independently SHA-256 checked with zero mismatches.
- Independent auditor regenerated 500 parent populations, 15,000 detected events, and 1,000 posterior curves; all 1,000 raw rows, two summary arms, and the paired contrast were reconstructed. Audit status: PASS.

## Fixed-truth results

| Detected events | Planned / successful | Mean posterior-median bias | MCSE of bias | RMSE | 68% coverage (Wilson 95%) | 90% coverage (Wilson 95%) | 95% coverage (Wilson 95%) | Mean 68% interval width |
|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| 4 | 500 / 500 | +5.060 | 0.717 | 16.793 | 0.696 (0.654–0.735) | 0.910 (0.882–0.932) | 0.954 (0.932–0.969) | 23.436 |
| 30 | 500 / 500 | +0.866 | 0.182 | 4.160 | 0.662 (0.619–0.702) | 0.860 (0.827–0.888) | 0.920 (0.893–0.941) | 7.355 |

All biases and interval widths are in km s^-1 Mpc^-1. Coverage Wilson intervals are based on 500 independent realizations.

The paired 30-event minus nested 4-event median-bias difference is -4.194 km s^-1 Mpc^-1 (MCSE 0.699; paired 95% t interval [-5.567, -2.820]). The 30-event experiment has substantially smaller median-bias RMSE and narrower posterior intervals in this synthetic setup. However, 68% coverage is 0.662 for 30 events and 0.696 for four events; both Wilson intervals include 0.68. The 30-event 90% and 95% coverage estimates are below nominal (0.86 and 0.92), with their uncertainty intervals to be considered rather than interpreting point estimates alone.

## Interpretation limits

This isolates event-count information under the full-parent model. It does **not** measure the effect of catalogue incompleteness, does not test estimated inclusion probabilities, and is not joint SBC or observational H0 inference. Because event counts are nested within the same generated realization, the paired contrast is useful for precision and bias comparison, not a comparison of independent experiments. The 561-point grid is a known numerical limitation to revisit against the 1,121-point reference before confirmatory calibration claims.
