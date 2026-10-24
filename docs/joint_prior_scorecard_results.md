# Joint prior-performance scorecard: validated results

The joint scorecard workflow completed successfully on Python 3.10, 3.11, and 3.12 (run 40075072714). It produced 18 rows for true H0 values 52, 70, 88; priors [50,90], [40,100], [30,110]; corrected and uncorrected methods; 300 replicates per truth; seed 20261024; 4,001-point grid.

## Key results

- At H0=52 with prior [50,90], 68% coverage was 0/300 for both methods. The 95% Wilson interval for coverage is 0–1.26%.
- At H0=88 with prior [50,90], 68% coverage was 13/300 (4.33%; Wilson 95% interval 2.55–7.27%) for both methods.
- Widening the prior increases credible interval widths substantially. For uncorrected inference at H0=52, mean 68% width rises from 26.15 ([50,90]) to 37.58 ([40,100]) and 45.43 ([30,110]).
- Selection correction is not uniformly better: at H0=88 with prior [40,100], 68% coverage was 81.0% uncorrected and 58.7% corrected; RMSE was 17.25 and 19.96 respectively.
- At H0=70 with prior [50,90], correction reduced bias (1.60 to 0.15) but slightly increased RMSE (5.90 to 6.11).

## Interpretation limits

These are synthetic toy-model results, not a real gravitational-wave measurement. Only three truths and three prior intervals were tested. Prior widening trades boundary sensitivity against larger posterior intervals; it does not by itself establish a preferred prior. Selection correction's effect depends on the scenario and metric. Interpret coverage alongside Wilson intervals, bias, RMSE, and interval width. Fixed-truth coverage and prior-predictive SBC answer different questions.

## Reproduction

`python scripts/run_joint_prior_scorecard.py --replicates 300 --seed 20261024 --grid-size 4001 --output results/joint-prior-scorecard.csv`

`python -m pytest tests/test_joint_prior_scorecard.py -q`
