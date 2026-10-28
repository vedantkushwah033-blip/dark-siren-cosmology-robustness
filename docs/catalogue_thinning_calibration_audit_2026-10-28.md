# Catalogue-thinning calibration audit

**Status:** deterministic synthetic software-calibration experiment; not a real-survey validation.  
**Seed:** 20261028 (default)  
**Default parent catalogue:** 40,000 fixed distance values over [0, 1200) Mpc  
**Selection law:** p(d) = 0.05 + 0.95 exp(-d / 450 Mpc)

## What this experiment checks

The runner assigns each object its declared inclusion probability and uses the production Bernoulli-thinning primitive to draw a retained/not-retained flag. It preserves object-level raw records, the seed, declared model, and a SHA-256 digest. Four prespecified distance strata are summarized using the expected retained count \(\sum_i p_i\), observed count, binomial standard deviation \(\sqrt{\sum_i p_i(1-p_i)}\), and standardized residual.

The separate audit script reloads raw records, verifies the checksum and formula, and reconstructs totals and stratum summaries without importing the production completeness function. A six-standard-deviation guard is only a gross implementation-error check; it is not a formal scientific hypothesis test.

## Reproduce

```bash
python scripts/run_catalogue_thinning_audit.py --n-parent 40000 --seed 20261028 --output-dir artifacts/catalogue_thinning
python scripts/audit_catalogue_thinning_audit.py artifacts/catalogue_thinning
```

The output contains:
- `catalogue_thinning_records.csv`: object-level distance, inclusion probability, and retention flag;
- `catalogue_thinning_summary.json`: declared configuration and summaries;
- `catalogue_thinning_audit.json`: independent reconstruction report.

## Interpretation limits

This experiment conditions on a fixed, evenly spaced parent distance grid. It isolates the Bernoulli-thinning mechanism and does not model a Poisson galaxy population, latent redshift clustering, uncertain catalogue selection, galaxy photometric-redshift errors, GW host selection, or event-detection selection. It does not test the selection-aware H0 likelihood and does not imply that the chosen probability law describes any real survey. The next scientific gate is a correctly specified catalogue likelihood, followed by paired synthetic H0 inference and calibration; the full GP-versus-homogeneous benchmark remains NO-GO until those gates pass.
