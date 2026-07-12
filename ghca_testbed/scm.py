"""The six canonical spike-wave SCM graphs (C1), as executable ground truth.

These are the pure-SCM scenarios from `experiments/c1_graph_certificates.py`,
lifted here so the testbed can *run* them (fast, no substrate) and validate an
epiphenomenality method end-to-end against a known verdict. `W` is a designated
autonomous node here (constitution `W = f(S)` is the C2+ story).

Each builder returns `(obs, do_samples, edges, S)`:
  obs         : dict of arrays under no intervention (observational)
  do_samples  : {0: samples under do(W=0), 1: samples under do(W=1)}
  edges       : DAG edge list (for the Theorem-1 certificate)
  S           : the observed conditioning set (list of variable names)
"""

import numpy as np

# ground-truth causal verdict per graph (paper + C1)
EXPECTED_VERDICT = {
    "fork": "epiphenomenal",
    "confounded": "epiphenomenal",
    "direct": "causal",
    "mediator_unobs": "causal",
    "mediator_obs": "epiphenomenal",
    "frontdoor": "causal",
}

EDGES = {
    "fork": [("S", "W"), ("S", "B")],
    "confounded": [("U", "W"), ("U", "B")],
    "direct": [("W", "B")],
    "mediator_unobs": [("W", "M"), ("M", "B")],
    "mediator_obs": [("W", "M"), ("M", "B")],
    "frontdoor": [("C", "W"), ("C", "B"), ("W", "M"), ("M", "B")],
}

COND = {  # observed conditioning set S per graph
    "fork": ["S"], "confounded": [], "direct": [], "mediator_unobs": [],
    "mediator_obs": ["M"], "frontdoor": ["M"],
}


def _bern(rng, p, n):
    return (rng.random(n) < p).astype(int)


def _noisy(rng, x, flip=0.1):
    f = rng.random(len(x)) < flip
    return np.where(f, 1 - x, x)


def _sample(name, rng, n, do_w=None):
    """Structural equations for `name` (matches c1_graph_certificates.py)."""
    if name == "fork":
        S = _bern(rng, 0.5, n)
        W = _noisy(rng, S) if do_w is None else np.full(n, do_w)
        B = _noisy(rng, S)
        return {"S": S, "W": W, "B": B}
    if name == "confounded":
        U = _bern(rng, 0.5, n)
        W = _noisy(rng, U) if do_w is None else np.full(n, do_w)
        B = _noisy(rng, U)
        return {"U": U, "W": W, "B": B}
    if name == "direct":
        W = _bern(rng, 0.5, n) if do_w is None else np.full(n, do_w)
        B = _noisy(rng, W)
        return {"W": W, "B": B}
    if name in ("mediator_unobs", "mediator_obs"):
        W = _bern(rng, 0.5, n) if do_w is None else np.full(n, do_w)
        M = _noisy(rng, W)
        B = _noisy(rng, M)
        return {"W": W, "M": M, "B": B}
    if name == "frontdoor":
        C = _bern(rng, 0.5, n)
        W = _noisy(rng, C) if do_w is None else np.full(n, do_w)
        M = _noisy(rng, W)
        B = _noisy(rng, ((M + C) >= 1).astype(int))
        return {"C": C, "W": W, "M": M, "B": B}
    raise ValueError(f"unknown graph {name}")


def build_scm(name, n=200_000, seed=0):
    """Return (obs, do_samples, edges, S) for one canonical graph.

    A single RNG stream is consumed obs -> do(0) -> do(1) (as in C1); with the
    default n and seed this reproduces the C1 numbers.
    """
    rng = np.random.default_rng(seed)
    obs = _sample(name, rng, n, do_w=None)
    do0 = _sample(name, rng, n, do_w=0)
    do1 = _sample(name, rng, n, do_w=1)
    return obs, {0: do0, 1: do1}, EDGES[name], COND[name]


SCM_GRAPHS = list(EXPECTED_VERDICT)
