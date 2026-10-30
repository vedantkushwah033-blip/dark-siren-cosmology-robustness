# Independent bootstrap-refit audit extension

The artifact now retains the exact validation-record index vector for every attempted bootstrap draw and a success mask in attempt order. The auditor replays the random-resampling design from these stored indices, independently fits the bounded logistic model, independently reconstructs inverse-weight diagnostics and the reflected weighted KDE, and compares those results with the saved successful draws before reconstructing posterior summaries.

This strengthens the earlier audit: a PASS now covers both the validation-refit-to-density stage and posterior propagation, subject to numerical agreement tolerances. The auditor does not rely on the production fitting or KDE functions. The stored resample indices also make each draw inspectable without relying on a seed alone.

The bootstrap remains an iid-with-replacement approximation to a validation sample drawn without replacement from a finite parent frame. This does not solve the finite-population design issue; it makes the approximation auditable. A design-aware bootstrap or replicate-level resampling study remains a separate methodological task.
