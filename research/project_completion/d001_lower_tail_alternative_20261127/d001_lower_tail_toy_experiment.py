"""D-001 toy demonstration (NOT production code).

Alternative generator: host z~rho (Euclidean), d_obs = mu + sigma*eps,
mu=c*z/H0, sigma=rel*mu, Gaussian in distance.
ACCEPT iff d_obs < tau.

Question: which normalisation of the accepted-event density is correct,
and what does a wrong tail do to H0?
"""
import numpy as np
from scipy.special import ndtr

c = 299792.458
zmax = 0.1
G = 400
zg = (np.arange(G) + 0.5) * zmax / G
w = zg**2
w = w / w.sum()
rel = 0.1
tau = 330.0
H0 = 70.0

mu_ = lambda H: c * zg[None, :] / np.atleast_1d(H)[:, None]


def beta(H):
    m = mu_(H)
    return (w * ndtr((tau - m) / (rel * m))).sum(1)


# ---- 1) Which density does the generator imply? ----
rng = np.random.default_rng(20261129)
n = 1_500_000
idx = rng.choice(G, size=n, p=w)
m = c * zg[idx] / H0
d = m + rel * m * rng.standard_normal(n)
acc = d[d < tau]

bins = np.linspace(0, tau, 61)
h, _ = np.histogram(acc, bins)
p_emp = h / h.sum()
ctr = 0.5 * (bins[1:] + bins[:-1])
dd = bins[1] - bins[0]

mu0 = mu_(H0)[0]
sg = rel * mu0
N = (
    np.exp(-0.5 * ((ctr[:, None] - mu0[None, :]) / sg[None, :]) ** 2)
    / (sg[None, :] * np.sqrt(2 * np.pi))
)
b0 = beta(H0)[0]
Phi_i = ndtr((tau - mu0) / sg)

# (a) Host-marginal normalisation:
# p(d|accepted) = sum_z w N(d;mu_z,sigma_z) / beta, beta=P(d<tau).
fa = (N * w).sum(1) / b0 * dd

# (b) Per-host normalisation: a different host-weighting model.
fb = (N * (w / Phi_i)).sum(1) * dd

# (c) Wrong tail for this lower-tail generator: divide by P(d>=tau).
fc = (N * w).sum(1) / (1 - b0) * dd

tv = lambda f: 0.5 * np.abs(p_emp - f / f.sum()).sum()

print(
    f"acceptance fraction MC={len(acc)/n:.4f}  "
    f"analytic beta(H0)={b0:.4f}   (n draws {n})"
)
print(
    "TV distance to accepted-event histogram: "
    f"(a) host-marginal {tv(fa):.4f} | "
    f"(b) per-host {tv(fb):.4f} | "
    f"(c) wrong-tail {tv(fc):.4f}   "
    f"[noise floor ~ {np.sqrt(60/len(acc))*0.5:.4f}]"
)
print(
    f"mass of each candidate density on d<tau: "
    f"(a) {fa.sum():.4f}  (b) {fb.sum():.4f}  (c) {fc.sum():.4f}  "
    f"[c should equal beta/(1-beta)={b0/(1-b0):.4f}]"
)

# ---- 2) Effect on H0 posterior of the normalisation convention ----
Hg = np.linspace(50, 90, 401)
mu_g = mu_(Hg)
sg_g = rel * mu_g
B = beta(Hg)
R = 300
Nev = 60

res = {
    "correct (host-marginal beta^-N)": [],
    "no selection normalisation": [],
    "wrong tail ((1-beta)^-N)": [],
    "per-host Phi_z^-1 inside integral": [],
}
Phi_g = ndtr((tau - mu_g) / sg_g)

for r in range(R):
    ev = []
    while len(ev) < Nev:
        z = rng.choice(G, p=w)
        mm = c * zg[z] / H0
        x = mm + rel * mm * rng.standard_normal()
        if x < tau:
            ev.append(x)
    ev = np.array(ev)

    K = (
        np.exp(
            -0.5
            * ((ev[:, None, None] - mu_g[None]) / sg_g[None]) ** 2
        )
        / (sg_g[None] * np.sqrt(2 * np.pi))
    )  # (events, H0, z)

    num = np.log((K * w).sum(2)).sum(0)
    numb = np.log((K * (w / Phi_g)).sum(2)).sum(0)

    candidates = (
        (
            "correct (host-marginal beta^-N)",
            num - Nev * np.log(B),
        ),
        ("no selection normalisation", num),
        ("wrong tail ((1-beta)^-N)", num - Nev * np.log(1 - B)),
        ("per-host Phi_z^-1 inside integral", numb),
    )

    for key, ll in candidates:
        p = np.exp(ll - ll.max())
        p /= np.trapezoid(p, Hg)
        cdf = np.cumsum(p)
        cdf /= cdf[-1]
        res[key].append(np.interp(0.5, cdf, Hg))

print(
    f"\nposterior-median error (truth {H0}), "
    f"{R} replicates x {Nev} accepted events, mean ± MC SE and RMSE:"
)
for k, v in res.items():
    v = np.array(v) - H0
    print(
        f"  {k:36s} bias {v.mean():+7.3f} ± "
        f"{v.std(ddof=1)/np.sqrt(R):.3f}   "
        f"RMSE {np.sqrt((v**2).mean()):.3f}"
    )

print(
    "beta(H) increasing in H on grid:",
    bool(np.all(np.diff(B) > 0)),
    "| beta(50), beta(90):",
    round(B[0], 4),
    round(B[-1], 4),
)
