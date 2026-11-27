# D-001: production GE versus hypothetical lower-tail oracle

**Purpose:** Prevent the alternative lower-tail theory experiment from being confused with the repository's production generative model.

## Audited production convention

The D-001 source audit supplied on 27 November 2026 traced the one-sided GE path to

\[
A_{\mathrm{GE}}=\{d_{\mathrm{obs}}\geq\tau\}.
\]

For \(d_{\mathrm{obs}}\mid z,H_0\sim\mathcal N(\mu,\sigma^2)\), its host-specific efficiency is

\[
P(A_{\mathrm{GE}}\mid z,H_0)
=\Phi\left(\frac{\mu-\tau}{\sigma}\right)
=1-\Phi\left(\frac{\tau-\mu}{\sigma}\right).
\]

The audit reported agreement among the GE generator, dataset efficiency, selection toy, shared-inclusion likelihood, and five-arm likelihood implementations. The reported targeted checks were 88 existing tests passed and 6 additional audit-only tests passed. Those test results belong to the prior audit report; they are not results from executing the two newly saved scripts.

## Alternative model saved in this directory

The two Python files saved alongside this note intentionally use

\[
A_{<}=\{d_{\mathrm{obs}}<\tau\},
\qquad
P(A_{<}\mid z,H_0)
=\Phi\left(\frac{\tau-\mu}{\sigma}\right).
\]

They are standalone mathematical/oracle exercises. The adapter dictionary in the regression-test file is set to None, so production-adapter tests are skipped. The test named as a generator acceptance test expresses the lower-tail oracle expectation only; it must not be wired to the existing GE generator.

## Required interpretation

- Do not change production GE formulas based on these lower-tail files.
- Do not use lower-tail toy outputs as evidence about GE production behaviour.
- Keep interval selection as a third, distinct event, \(d_{\min}<d_{\mathrm{obs}}<d_{\max}\).
- Historical experiment manifests must be checked individually before assigning their selection convention.
- Run the oracle test file and toy simulation independently before citing numerical outputs.
- Any future production adapter tests must be written for the actual audited GE convention and should be separate from this lower-tail alternative.

## Saved files

1. **test_d001_detection_convention.py** — lower-tail oracle tests with intentionally unconnected production adapters.
2. **d001_lower_tail_toy_experiment.py** — lower-tail simulation and posterior-normalisation comparison.
3. **D001_host_marginal_normalisation_theory.md** — derivation of the accepted-event density and host-marginal normalisation.
