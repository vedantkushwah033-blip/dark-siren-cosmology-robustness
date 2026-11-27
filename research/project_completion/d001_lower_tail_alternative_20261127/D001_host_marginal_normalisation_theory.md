# D-001 lower-tail alternative: mathematical theory

**Status:** Saved from the 27 November 2026 conversation. This is an independent oracle/theory exercise for a hypothetical lower-tail generator. It is **not** the repository's audited production GE convention.

## 1. Alternative generative model

Let a parent host have redshift \(z\sim g(z)\), with

\[
d_{\mathrm{obs}}\mid z,H_0\sim
\mathcal N\left(\mu(z,H_0),\sigma^2(z,H_0)\right),
\qquad
\mu(z,H_0)=\frac{cz}{H_0}.
\]

The toy experiment uses \(\sigma(z,H_0)=r\mu(z,H_0)\), with relative scale \(r=0.1\), and accepts an event when

\[
A_{<}=\{d_{\mathrm{obs}}<\tau\}.
\]

For a fixed host,

\[
P(A_{<}\mid z,H_0)
=\Phi\left(\frac{\tau-\mu(z,H_0)}{\sigma(z,H_0)}\right).
\]

The parent-population detection efficiency is the host-marginal probability

\[
\beta_{<}(H_0)
=\int g(z)
\Phi\left(\frac{\tau-\mu(z,H_0)}{\sigma(z,H_0)}\right)\,dz.
\]

## 2. Accepted-event density and the normalisation

For a parent host drawn from \(g(z)\), the joint density of an observation and its host is

\[
p(d,z\mid H_0)=g(z)\,
\mathcal N(d\mid\mu(z,H_0),\sigma^2(z,H_0)).
\]

Conditioning on acceptance gives

\[
p(d\mid A_{<},H_0)
=
\frac{
\int g(z)\mathcal N(d\mid\mu(z,H_0),\sigma^2(z,H_0))\,dz
}{
\beta_{<}(H_0)
}
\mathbf 1\{d<\tau\}.
\]

Thus, for \(N\) conditionally independent accepted events and a fixed detected count,

\[
\mathcal L(H_0)
\propto
\prod_{j=1}^{N}
\frac{
\int g(z)\mathcal N(d_j\mid\mu(z,H_0),\sigma^2(z,H_0))\,dz
}{
\beta_{<}(H_0)
}
\mathbf 1\{d_j<\tau\}.
\]

Equivalently, the log likelihood contains the term
\(-N\log\beta_{<}(H_0)\), once per dataset, not an average of the event-specific normalisations.

This result assumes the parent host distribution is \(g(z)\), the event observations are conditionally independent given the model, the detected count is conditioned upon, and the same acceptance event is used by the generator and likelihood.

## 3. Opposite-tail diagnostic

For a continuous Gaussian observation,

\[
P(d_{\mathrm{obs}}\geq\tau\mid z,H_0)
=
\Phi\left(\frac{\mu(z,H_0)-\tau}{\sigma(z,H_0)}\right)
=1-P(A_{<}\mid z,H_0).
\]

The host-marginal opposite-tail efficiency is \(1-\beta_{<}(H_0)\). Dividing the lower-tail accepted-event numerator by this opposite-tail probability is inconsistent with the lower-tail generator. Over the accepted region, its integrated mass is

\[
\frac{\beta_{<}(H_0)}{1-\beta_{<}(H_0)},
\]

rather than one.

## 4. Host-marginal versus per-host truncation

The accepted host distribution is not the original parent distribution. It is

\[
g(z\mid A_{<},H_0)
=
\frac{g(z)P(A_{<}\mid z,H_0)}{\beta_{<}(H_0)}.
\]

Given a host and an accepted observation,

\[
p(d\mid z,A_{<},H_0)
=
\frac{\mathcal N(d\mid\mu,\sigma^2)}
{P(A_{<}\mid z,H_0)}
\mathbf 1\{d<\tau\}.
\]

Marginalising over the **accepted-host distribution** recovers the host-marginal expression:

\[
\int g(z\mid A_{<},H_0)
p(d\mid z,A_{<},H_0)\,dz
=
\frac{\int g(z)\mathcal N(d\mid\mu,\sigma^2)\,dz}
{\beta_{<}(H_0)}.
\]

Therefore, dividing each host likelihood by its own acceptance probability while retaining the *parent* weights \(g(z)\) is not equivalent. The host weights must also be updated to the accepted-host distribution. The toy script includes this unadjusted per-host expression as a deliberately contrasting diagnostic, not as a valid substitute for the correctly marginalised accepted-event density.

## 5. Monotonicity in the toy model

For a fixed \(z>0\), \(\mu=cz/H_0\) decreases as \(H_0\) increases. With \(\sigma=r\mu\),

\[
\frac{\tau-\mu}{r\mu}
=\frac{\tau}{r\mu}-\frac{1}{r}
\]

increases as \(H_0\) increases. Hence the lower-tail host-specific acceptance probability, and its host-marginal average, increase with \(H_0\) (for positive \(r,z,\tau,H_0\)).

## 6. Reproducibility files

- **test_d001_detection_convention.py**: standalone oracle algebra, monotonicity, Monte Carlo, normalisation, and log-tail tests. Production adapter tests remain deliberately skipped until wired to real repository functions.
- **d001_lower_tail_toy_experiment.py**: simulation comparing host-marginal normalisation, an unadjusted per-host diagnostic, and the opposite-tail error, including an illustrative \(H_0\) posterior comparison.

These scripts were saved as separate files to avoid the original pasted message's concatenation of two Python modules.

## 7. Interpretation and limitations

This alternative lower-tail model is useful for understanding selection-conditioned likelihood normalisation and for demonstrating the consequences of a mismatched tail. It does not establish that the repository should use lower-tail acceptance. It must not be merged into production GE logic or used to relabel historical GE experiments. The toy results must be run and independently checked before numerical output is cited in a paper.
