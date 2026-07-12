"""Method-agnostic causal certificates and the Definition-1 epiphenomenality test.

These take SAMPLES (not a live substrate), so a user can score their own
estimator against the registry without touching the dynamics. The one
load-bearing subtlety (the C1 catch): Definition 1 is scored
interventional-vs-OBSERVATIONAL, never do(W=1)-vs-do(W=0) — the latter is blind
to the front-door case, where intervening breaks W-confounding without the
effect varying with the set value.
"""

from dataclasses import dataclass, field
import numpy as np

try:
    import networkx as nx
except Exception:  # networkx is only needed for the graphical certificate
    nx = None


@dataclass
class EpiResult:
    effect: float
    verdict: str            # 'epiphenomenal' | 'causal'
    threshold: float
    per_stratum: dict = field(default_factory=dict)


def strata_means(samples, cond, target="B", min_stratum=200):
    """E[target | cond=c] per stratum c (binary conditioning vars). {combo: mean}."""
    B = np.asarray(samples[target])
    n = len(B)
    if not cond:
        return {(): float(B.mean())}
    keys = [np.asarray(samples[v]) for v in cond]
    out = {}
    for combo in range(2 ** len(cond)):
        mask = np.ones(n, bool)
        for bit, arr in enumerate(keys):
            mask &= (arr == ((combo >> bit) & 1))
        if mask.sum() > min_stratum:
            out[combo] = float(B[mask].mean())
    return out


def observational_association(samples, cond, wave="W", target="B", min_stratum=200):
    """max_c | E[B | c, W=1] - E[B | c, W=0] | — correlation, for contrast."""
    w = np.asarray(samples[wave])
    m0 = strata_means({k: np.asarray(v)[w == 0] for k, v in samples.items()},
                      cond, target, min_stratum)
    m1 = strata_means({k: np.asarray(v)[w == 1] for k, v in samples.items()},
                      cond, target, min_stratum)
    common = set(m0) & set(m1)
    return max((abs(m0[c] - m1[c]) for c in common), default=0.0)


def epiphenomenality_test(samples_obs, samples_do, *, target="B", cond=("S",),
                          min_stratum=200, threshold=0.03):
    """Definition 1: is W epiphenomenal to `target` given `cond`?

    `samples_do` is a dict {w: samples_under_do(W=w)} — dict-shaped precisely so
    the contrast cannot degenerate to do(1)-vs-do(0). Effect = max over set-values
    w AND strata c of | E[target | c; do(W=w)] - E_obs[target | c] |. Verdict
    'epiphenomenal' iff effect < threshold.
    """
    cond = tuple(cond)
    obs_m = strata_means(samples_obs, cond, target, min_stratum)
    effect = 0.0
    per = {}
    for w, sd in samples_do.items():
        dm = strata_means(sd, cond, target, min_stratum)
        for c in set(obs_m) & set(dm):
            gap = abs(dm[c] - obs_m[c])
            per[(w, c)] = gap
            effect = max(effect, gap)
    verdict = "epiphenomenal" if effect < threshold else "causal"
    return EpiResult(effect=float(effect), verdict=verdict,
                     threshold=threshold, per_stratum=per)


def theorem1_certificate(edges, S, W="W", B="B"):
    """Theorem 1: W epiphenomenal to B given S iff W ⊥d B | S in G_{W(S)}.

    G_{W(S)} removes arrows into W-nodes that are not ancestors of S. Returns
    'epiphenomenal' | 'causal'. Requires networkx.
    """
    if nx is None:
        raise ImportError("theorem1_certificate requires networkx")
    G = nx.DiGraph(edges)
    for node in [W, B] + list(S):
        G.add_node(node)
    anc_S = set()
    for x in S:
        anc_S |= nx.ancestors(G, x) | {x}
    WX = [w for w in [W] if w not in anc_S]
    Gm = G.copy()
    Gm.remove_edges_from([(p, w) for w in WX for p in list(G.predecessors(w))])
    d_sep = nx.is_d_separator(Gm, {W}, {B}, set(S))
    return "epiphenomenal" if d_sep else "causal"
