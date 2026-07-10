"""
C1 - Build the paper's canonical spike-wave graphs and validate the
epiphenomenality certificate on ground truth.

The paper (arXiv:2511.06602) defines W epiphenomenal to B given S iff
P(B|S; do(W)) = P(B|S), and gives a graphical certificate (Theorem 1):
W is epiphenomenal for ALL SCMs inducing graph G iff

    W  d-separated from  B  given  S   in   G_{W(S)},

where G_{W(S)} removes the arrows INTO the nodes of W that are not ancestors of
S. Here we can do what is impossible in vivo: build the SCM explicitly, run the
actual do(W) intervention, and check the measured verdict against the
certificate. W is treated as a designated autonomous node (constitution W=f(S)
is deferred to C2).

For each canonical graph we compute:
  * causal_effect = max_s | E[B | S=s; do(W=1)] - E[B | S=s; do(W=0)] |
      (ground-truth interventional test: does do(W) move B given the observed S?)
  * obs_dependence = max_s | E[B | S=s, W=1] - E[B | S=s, W=0] |
      (observational association, for contrast: correlation != causation)
  * certificate: does Theorem 1 predict epiphenomenal?
and check that (causal_effect ~ 0) <=> (certificate says epiphenomenal).

Outputs
-------
docs/figures/c1_certificates.png
result/c1/c1_data.npz
"""

import os
import sys
import numpy as np
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c1")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N = 200_000            # samples
rng = np.random.default_rng(0)


def bern(p):
    return (rng.random(N) < p).astype(int)


def noisy(x, flip=0.1):
    """pass x through with a small flip probability (soft structural equation)"""
    f = rng.random(N) < flip
    return np.where(f, 1 - x, x)


# ---------------------------------------------------------------------------
# Each graph: returns (samples dict, DAG edges, observed set S, do_W function).
# do_W(w) returns a samples dict under do(W=w).
# ---------------------------------------------------------------------------

def g_fork():
    """Fork: S -> W, S -> B  (no W->B). W epiphenomenal given S."""
    edges = [("S", "W"), ("S", "B")]

    def sample(do_w=None):
        S = bern(0.5)
        W = noisy(S) if do_w is None else np.full(N, do_w)
        B = noisy(S)
        return {"S": S, "W": W, "B": B}
    return sample, edges, ["S"]


def g_confounded():
    """Confounded: U -> W, U -> B, U latent, no W->B. Epiphenomenal under do,
    but W and B are observationally dependent."""
    edges = [("U", "W"), ("U", "B")]

    def sample(do_w=None):
        U = bern(0.5)
        W = noisy(U) if do_w is None else np.full(N, do_w)
        B = noisy(U)
        return {"U": U, "W": W, "B": B}
    return sample, edges, []          # U unobserved; S empty


def g_direct():
    """Direct cause: W -> B. Causal."""
    edges = [("W", "B")]

    def sample(do_w=None):
        W = bern(0.5) if do_w is None else np.full(N, do_w)
        B = noisy(W)
        return {"W": W, "B": B}
    return sample, edges, []


def g_mediator_unobs():
    """W -> M -> B, M unobserved. Causal (Prop 1: unblocked path)."""
    edges = [("W", "M"), ("M", "B")]

    def sample(do_w=None):
        W = bern(0.5) if do_w is None else np.full(N, do_w)
        M = noisy(W)
        B = noisy(M)
        return {"W": W, "M": M, "B": B}
    return sample, edges, []          # M unobserved


def g_mediator_obs():
    """W -> M -> B, M observed. Epiphenomenal given M (Fig 8): conditioning on
    the mediator blocks the only path."""
    edges = [("W", "M"), ("M", "B")]

    def sample(do_w=None):
        W = bern(0.5) if do_w is None else np.full(N, do_w)
        M = noisy(W)
        B = noisy(M)
        return {"W": W, "M": M, "B": B}
    return sample, edges, ["M"]       # M observed, condition on it


def g_frontdoor():
    """Front-door: C -> W, C -> B, W -> M -> B, M observed, C latent.
    Causal (Prop 2), despite conditioning on the observed mediator M."""
    edges = [("C", "W"), ("C", "B"), ("W", "M"), ("M", "B")]

    def sample(do_w=None):
        C = bern(0.5)
        W = noisy(C) if do_w is None else np.full(N, do_w)
        M = noisy(W)
        # B depends on both the mediator M and the latent confounder C
        B = noisy(((M + C) >= 1).astype(int))
        return {"C": C, "W": W, "M": M, "B": B}
    return sample, edges, ["M"]       # M observed; C latent


GRAPHS = {
    "fork": g_fork,
    "confounded": g_confounded,
    "direct": g_direct,
    "mediator (unobs)": g_mediator_unobs,
    "mediator (obs)": g_mediator_obs,
    "front-door": g_frontdoor,
}
EXPECTED = {                          # ground-truth causal verdict
    "fork": "epiphenomenal", "confounded": "epiphenomenal", "direct": "causal",
    "mediator (unobs)": "causal", "mediator (obs)": "epiphenomenal",
    "front-door": "causal",
}


def strata_means(s, cond_vars):
    """E[B | S=s] per stratum; returns dict {stratum_combo: mean}."""
    B = s["B"]
    n = len(B)
    if not cond_vars:
        return {(): B.mean()}
    keys = [s[v] for v in cond_vars]
    out = {}
    for combo in range(2 ** len(cond_vars)):
        mask = np.ones(n, bool)
        for bit, arr in enumerate(keys):
            mask &= (arr == ((combo >> bit) & 1))
        if mask.sum() > 200:
            out[combo] = B[mask].mean()
    return out


def certificate_epiphenomenal(edges, S):
    """Theorem 1: W epiphenomenal to B given S iff W ⊥d B | S in G_{W(S)}."""
    G = nx.DiGraph(edges)
    for n in ["W", "B"] + list(S):
        G.add_node(n)
    anc_S = set()
    for x in S:
        anc_S |= nx.ancestors(G, x) | {x}
    WX = [w for w in ["W"] if w not in anc_S]     # W-nodes not ancestor of S
    Gm = G.copy()
    Gm.remove_edges_from([(p, w) for w in WX for p in list(G.predecessors(w))])
    return nx.is_d_separator(Gm, {"W"}, {"B"}, set(S))


def main():
    rows = []
    for name, gfn in GRAPHS.items():
        sample, edges, S = gfn()
        obs = sample(do_w=None)
        do0 = sample(do_w=0)
        do1 = sample(do_w=1)

        # Definition 1: epiphenomenal iff P(B|S; do(W=w)) == observational P(B|S)
        # for all w. Causal effect = how far any interventional stratum mean
        # departs from the OBSERVATIONAL stratum mean (not do1 vs do0 -- that
        # misses the front-door case, where intervening breaks W-confounding
        # without varying with w).
        obs_m = strata_means(obs, S)
        d0, d1 = strata_means(do0, S), strata_means(do1, S)
        common = set(obs_m) & set(d0) & set(d1)
        causal_effect = max((max(abs(d0[c] - obs_m[c]), abs(d1[c] - obs_m[c]))
                             for c in common), default=0.0)

        # observational association for contrast (does B track W given S?)
        ow0 = strata_means({k: v[obs["W"] == 0] for k, v in obs.items()}, S)
        ow1 = strata_means({k: v[obs["W"] == 1] for k, v in obs.items()}, S)
        cc = set(ow0) & set(ow1)
        obs_dep = max((abs(ow0[c] - ow1[c]) for c in cc), default=0.0)

        cert_epi = certificate_epiphenomenal(edges, S)
        measured = "epiphenomenal" if causal_effect < 0.03 else "causal"
        cert = "epiphenomenal" if cert_epi else "causal"
        ok = (measured == cert == EXPECTED[name])
        rows.append((name, causal_effect, obs_dep, cert, measured, EXPECTED[name], ok))
        print(f"{name:18s} do-effect={causal_effect:.3f} obs-assoc={obs_dep:.3f} "
              f"| cert={cert:13s} measured={measured:13s} expected={EXPECTED[name]:13s} "
              f"{'OK' if ok else 'MISMATCH'}")

    all_ok = all(r[6] for r in rows)
    print(f"\nAll graphs agree (certificate == intervention == expected): {all_ok}")

    # figure: do-effect vs observational association, colored by certificate
    fig, ax = plt.subplots(figsize=(9, 5.5))
    names = [r[0] for r in rows]
    x = np.arange(len(names))
    ce = [r[1] for r in rows]
    od = [r[2] for r in rows]
    ax.bar(x - 0.2, ce, 0.4, label="causal effect of do(W) given S", color="crimson")
    ax.bar(x + 0.2, od, 0.4, label="observational assoc. B–W given S", color="gray")
    ax.axhline(0.03, ls="--", color="k", alpha=0.5, label="epiphenomenal threshold")
    for i, r in enumerate(rows):
        ax.text(i, max(r[1], r[2]) + 0.01, "epi" if r[3] == "epiphenomenal" else "causal",
                ha="center", fontsize=8,
                color="steelblue" if r[3] == "epiphenomenal" else "crimson")
    ax.set_xticks(x)
    ax.set_xticklabels(names, rotation=20, ha="right")
    ax.set_ylabel("effect size")
    ax.set_title("C1: Theorem-1 certificate matches ground-truth do(W) on every graph\n"
                 "(note 'confounded': strong association, zero causal effect)")
    ax.legend()
    fig.tight_layout()
    p = os.path.join(FIGDIR, "c1_certificates.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c1_data.npz"),
             names=np.array(names),
             causal_effect=np.array([r[1] for r in rows]),
             obs_dependence=np.array([r[2] for r in rows]),
             certificate=np.array([r[3] for r in rows]),
             measured=np.array([r[4] for r in rows]),
             expected=np.array([r[5] for r in rows]),
             all_ok=all_ok)
    print("wrote", os.path.join(DATADIR, "c1_data.npz"))


if __name__ == "__main__":
    main()
