# Persisted validation-bootstrap sensitivity workflow

**Status:** implementation and Python 3.12 CI smoke step added on a dedicated research branch; the first end-to-end CI attempt exposed and is correcting a shell handoff issue.

The reusable bootstrap helpers previously existed only as library functions. This workflow makes one selected successful estimated-inclusion pilot replicate reproducible as a persisted artifact:

1. Select one replicate and validation-size arm from the pilot CSV.
2. Read only the inference-audit payload (catalogue redshifts, validation redshifts and retention labels, grids, and fixed GW distances).
3. Resample paired validation records with replacement, refit the logistic inclusion model, reconstruct the weighted host density, and propagate each successful density through the fixed toy selection-conditioned likelihood.
4. Persist numerical arrays in a compressed NPZ file, primitive inputs in JSON, and attempt/failure accounting plus interpretation limits in a summary JSON.
5. Run a separate auditor that reads the saved artifacts and pilot CSV, verifies provenance and shape/count invariants, and independently reconstructs every saved posterior summary from each persisted density draw using direct event-by-event quadrature.

## Commands

First create the ordinary non-confirmatory pilot artifacts:

    python scripts/run_estimated_inclusion_pilot.py --replicates 20 --n-parent 4000 --n-detected 30 --validation-sizes 100,300,1000 --seed 20261030 --output-dir artifacts/estimated_inclusion_smoke

Then select a replicate/arm that completed successfully:

    python scripts/run_validation_bootstrap_sensitivity.py --pilot-csv artifacts/estimated_inclusion_smoke/replicates.csv --replicate 0 --validation-size 300 --n-bootstrap 500 --seed 20261031 --output-dir artifacts/validation_bootstrap_rep0_n300

Audit the saved output:

    python scripts/audit_validation_bootstrap_sensitivity.py --artifact-dir artifacts/validation_bootstrap_rep0_n300 --pilot-csv artifacts/estimated_inclusion_smoke/replicates.csv

The chosen replicate and arm must exist and be successful; the runner fails closed otherwise. Use a separate output directory for each selected replicate/arm/seed combination.

## Audit boundary

The auditor independently replays the declared seed to reconstruct every validation-resample index vector, refits the bounded logistic model with a separate implementation, reconstructs inverse-weight diagnostics and reflected weighted KDE densities, and recomputes posterior means, medians, and 68% interval endpoints by direct quadrature. Its PASS therefore covers the persisted resampling/refit/density/posterior chain subject to numerical tolerances; it is still an internal independent implementation, not external scientific validation. The summary states this explicitly.

The bootstrap is an approximate paired-record nonparametric bootstrap with replacement. The actual synthetic validation sample is drawn without replacement from a finite parent frame, but no finite-population correction is applied. The observed catalogue and GW distances are held fixed. The logistic model family and KDE construction are held fixed. The percentile range across posterior means is a validation-sample sensitivity distribution, not an H0 posterior credible interval, not a coverage guarantee, and not a full uncertainty interval.

No scientific conclusion should be drawn from a single selected replicate. Confirmatory claims require a preregistered replicate-level design, sufficient Monte Carlo precision, failure accounting, and independent review.
