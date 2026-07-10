"""
Causal / constitution instrumentation for the C-series.

Provides, on top of the ``ghca_net.Network`` substrate, the machinery the
C-series needs (docs/causal_experiments.md):

  * an observation model: partial spikes S_obs (a subset of nodes),
  * wave variables W = f(S) as explicit deterministic coarse-grainings of the
    FULL node state (coherence, active fraction),
  * three intervention operators: do(S), do(W) with a family of realization
    policies (the fat-handedness probe), and do(theta) on generating parameters.

Nothing here is behaviour-specific; experiments supply the substrate and the
behaviour readout B.
"""

import numpy as np


# ----------------------------------------------------------------------------
# Wave variables W = f(S): deterministic coarse-grainings of the full state.
# ----------------------------------------------------------------------------

def wave_coherence(net, nodes=None):
    """Kuramoto phase-coherence of the cyclic phases over `nodes` (default all).

    W = | mean_j exp(i * 2*pi * phi_j / tau_j) |  -- a scalar in [0, 1].
    A deterministic function of the node state phi (hence of S).
    """
    idx = np.arange(net.N) if nodes is None else np.asarray(nodes)
    ang = 2 * np.pi * net.phi[idx] / net.tau[idx]
    return float(np.abs(np.exp(1j * ang).mean()))


def wave_active_fraction(net, nodes=None):
    """Fraction of `nodes` (default all) currently active (spiking)."""
    idx = np.arange(net.N) if nodes is None else np.asarray(nodes)
    am = net.active_mask()
    return float(am[idx].mean())


WAVE_FNS = {"coherence": wave_coherence, "active_fraction": wave_active_fraction}


# ----------------------------------------------------------------------------
# Observation model: partial spikes S_obs.
# ----------------------------------------------------------------------------

def make_observation_mask(N, observed_nodes, frac=0.4, seed=0):
    """Boolean mask selecting the observed subset S_obs of the node set.

    `observed_nodes` restricts the candidate pool (e.g. hidden nodes); a random
    `frac` of them is marked observed.
    """
    rng = np.random.default_rng(seed)
    pool = np.asarray(observed_nodes)
    k = max(1, int(round(frac * len(pool))))
    obs = rng.choice(pool, size=k, replace=False)
    mask = np.zeros(N, dtype=bool)
    mask[obs] = True
    return mask


def read_spikes(net, obs_mask):
    """S_obs feature vector: active-mask over the observed nodes (0/1)."""
    return net.active_mask()[obs_mask].astype(float)


# ----------------------------------------------------------------------------
# Intervention operators.
# ----------------------------------------------------------------------------

def do_S(net, nodes, values):
    """do(S): clamp `nodes` to phases `values` for the current step."""
    net.phi[np.asarray(nodes)] = np.asarray(values)


def do_theta(net, tau=None, coupling_scale=None, theta=None, nodes=None):
    """do(theta): intervene on generating parameters (modular, well-posed).

    Sets per-node timescale tau, scales couplings, or sets thresholds on `nodes`
    (default all). Returns nothing; mutates the net in place.
    """
    idx = np.arange(net.N) if nodes is None else np.asarray(nodes)
    if tau is not None:
        net.tau[idx] = tau
    if theta is not None:
        net.theta[idx] = theta
    if coupling_scale is not None:
        net.W[np.ix_(idx, idx)] *= coupling_scale


def do_W(net, target, wave="active_fraction", policy="random", nodes=None,
         rng=None):
    """do(W = target): force the constituted aggregate f(S) to `target` by
    editing the node state S, under a realization `policy`.

    Because many micro-states S satisfy f(S)=target, the choice of `policy` is a
    free parameter -- this is exactly the fat-handedness probe of C2. Currently
    supports the active-fraction wave (a clean level set: "make exactly
    target*|nodes| nodes active"); coherence targeting is left to C2 where the
    realization family is the object of study.

    policy:
      'random'    -- activate a uniformly random subset of size target*|nodes|
      'clustered' -- activate a contiguous block (spatially structured)
      'min_edit'  -- flip the fewest nodes from the current state to hit target
    """
    if wave != "active_fraction":
        raise NotImplementedError("do_W realization implemented for active_fraction")
    rng = np.random.default_rng() if rng is None else rng
    idx = np.arange(net.N) if nodes is None else np.asarray(nodes)
    n = len(idx)
    k = int(round(target * n))
    am = net.active_mask()
    active_val = 1                      # phase for "active"
    rest_val = 0
    if policy == "random":
        chosen = rng.choice(idx, size=k, replace=False)
    elif policy == "clustered":
        start = rng.integers(0, n)
        chosen = idx[(start + np.arange(k)) % n]
    elif policy == "min_edit":
        cur_active = idx[am[idx]]
        cur_rest = idx[~am[idx]]
        if k <= len(cur_active):
            chosen = rng.choice(cur_active, size=k, replace=False)
        else:
            extra = rng.choice(cur_rest, size=k - len(cur_active), replace=False)
            chosen = np.concatenate([cur_active, extra])
    else:
        raise ValueError(f"unknown policy {policy}")
    net.phi[idx] = rest_val
    net.phi[chosen] = active_val
