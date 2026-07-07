"""
C2 - Constitution: do(W) is fat-handed when W = f(S) (headline).

C1 validated the paper's certificates treating W as a manipulable node distinct
from S. Here W is the genuine CONSTITUTED aggregate: W = f(S) = active fraction
of the hidden population, an explicit deterministic coarse-graining. Because many
micro-states S satisfy f(S) = w, the intervention do(W = w) is not unique -- it
must pick a realization via some policy pi. We show the causal response
E[B | do_pi(W = w)] DEPENDS ON pi for a behaviour that reads micro-structure,
i.e. do(W) is fat-handed / ill-posed once W is constituted. A purely collective
behaviour (that reads only the aggregate) is pi-invariant, delimiting exactly
when the wave intervention is well-defined.

The sharpest statement: a single do(W=w) does not define a unique behaviour --
it admits a whole ACHIEVABLE BAND of B, spanned by choosing which micro-state S
realizes f(S)=w. At the moment of intervention (t=0) we compute that band
exactly (activate the k = w*N nodes with the smallest / largest readout weight
for the band edges; random subsets for the typical spread). A collective readout
(near-uniform weights) has a near-zero band -- the macro value pins it -- while a
labeled-line readout (structured weights) has a wide band -- do(W=w) leaves B
essentially free. We also track how the band survives propagation.

Outputs
-------
docs/figures/c2_fat_handed.png
result/c2/c2_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network, smallworld
from ghca_causal import do_W, wave_active_fraction

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c2")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_H = 120
ACT, PAS, THETA, P_S = 2, 8, 1.0, 0.01
POLICIES = ["random", "clustered", "min_edit"]
W_GRID = np.linspace(0.1, 0.7, 7)
T_GRID = [0, 1, 2, 4, 8]   # propagation steps before readout
N_TRIALS = 250


def build(seed=0):
    W = smallworld(N_H, k=6, beta=0.15, seed=seed)
    net = Network(W, act=ACT, pas=PAS, theta=THETA, p_s=P_S, seed=seed)
    rng = np.random.default_rng(seed + 1)
    w_coll = 1.0 + 0.1 * rng.standard_normal(N_H)      # collective readout
    w_lab = rng.standard_normal(N_H); w_lab -= w_lab.mean()   # labeled-line
    return net, w_coll, w_lab


def readouts_from_active(a, w_coll, w_lab):
    return float(w_coll @ a), float(w_lab @ a)


def baseline_std(net, w_coll, w_lab, steps=800):
    """Observational std of each readout, to standardize the band across behaviours."""
    net.phi[:] = 0
    net.seed_random(0.1, 0.1)
    rc, rl = [], []
    for _ in range(steps):
        net.step(None)
        c, l = readouts_from_active(net.active_mask().astype(float), w_coll, w_lab)
        rc.append(c); rl.append(l)
    return np.std(rc) + 1e-9, np.std(rl) + 1e-9


def achievable_band_t0(w_vec, w):
    """Exact band of the readout w_vec . active over all micro-states with
    active-fraction w, at the moment of intervention (t=0). Edges = activate the
    k smallest / largest weights; also a random-subset mean."""
    k = int(round(w * N_H))
    s = np.sort(w_vec)
    lo = s[:k].sum()          # k smallest weights active
    hi = s[-k:].sum()         # k largest weights active
    mean = w * w_vec.sum()    # E over random subsets of size k
    return lo, hi, mean


def prop_spread(net, w_vec, is_coll, w_coll, w_lab, t_prop, rng, n=200):
    """Std of the readout across RANDOM realizations of do(W=w), after t_prop
    steps of propagation, averaged over w -- how much do(W=w) leaves B free."""
    hidden = np.arange(N_H)
    spreads = []
    for w in W_GRID:
        vals = []
        for _ in range(n):
            net.phi[:] = 0
            do_W(net, target=w, wave="active_fraction", policy="random",
                 nodes=hidden, rng=rng)
            for _ in range(t_prop):
                net.step(None)
            a = net.active_mask().astype(float)
            c, l = readouts_from_active(a, w_coll, w_lab)
            vals.append(c if is_coll else l)
        spreads.append(np.std(vals))
    return np.mean(spreads)


def main():
    net, w_coll, w_lab = build(seed=0)
    std_c, std_l = baseline_std(net, w_coll, w_lab)
    rng = np.random.default_rng(11)

    # Panel A: exact achievable band of B under do(W=w) at t=0
    lo_c, hi_c, m_c = zip(*[achievable_band_t0(w_coll, w) for w in W_GRID])
    lo_l, hi_l, m_l = zip(*[achievable_band_t0(w_lab, w) for w in W_GRID])
    band_c = (np.array(hi_c) - np.array(lo_c)) / std_c
    band_l = (np.array(hi_l) - np.array(lo_l)) / std_l
    print(f"t=0 achievable band width / obs-std (mean over w): "
          f"collective={band_c.mean():.2f}  labeled-line={band_l.mean():.2f}")

    # Panel B: random-realization spread vs propagation time
    sp_c, sp_l = [], []
    for t in T_GRID:
        sp_c.append(prop_spread(net, w_coll, True, w_coll, w_lab, t, rng) / std_c)
        sp_l.append(prop_spread(net, w_lab, False, w_coll, w_lab, t, rng) / std_l)
        print(f"t_prop={t}: random-realization spread/std  "
              f"collective={sp_c[-1]:.3f}  labeled-line={sp_l[-1]:.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    axes[0].fill_between(W_GRID, np.array(lo_c) / std_c, np.array(hi_c) / std_c,
                         color="crimson", alpha=0.25, label="collective: achievable band")
    axes[0].plot(W_GRID, np.array(m_c) / std_c, "-", color="crimson")
    axes[0].fill_between(W_GRID, np.array(lo_l) / std_l, np.array(hi_l) / std_l,
                         color="steelblue", alpha=0.25, label="labeled-line: achievable band")
    axes[0].plot(W_GRID, np.array(m_l) / std_l, "-", color="steelblue")
    axes[0].set_xlabel("do(W=w): target active fraction")
    axes[0].set_ylabel("behaviour B under do(W=w)  (standardized)")
    axes[0].set_title("t=0: a single do(W=w) admits a whole band of B\n"
                      "(collective pinned; labeled-line left free)")
    axes[0].legend(fontsize=8)

    axes[1].plot(T_GRID, sp_c, "o-", color="crimson", label="collective B")
    axes[1].plot(T_GRID, sp_l, "s-", color="steelblue", label="labeled-line B")
    axes[1].set_xlabel("propagation steps after do(W)")
    axes[1].set_ylabel("realization spread of B (std / obs-std)")
    axes[1].set_title("How the under-determination survives dynamics\n"
                      "(this relaxing medium partly forgets the realization)")
    axes[1].legend()
    fig.suptitle("C2: do(W) is fat-handed for a constituted W=f(S) -- one macro "
                 "target, many behaviours")
    fig.tight_layout()
    pth = os.path.join(FIGDIR, "c2_fat_handed.png")
    fig.savefig(pth, dpi=110)
    print("wrote", pth)

    np.savez(os.path.join(DATADIR, "c2_data.npz"),
             w_grid=W_GRID, t_grid=np.array(T_GRID),
             lo_c=lo_c, hi_c=hi_c, m_c=m_c, lo_l=lo_l, hi_l=hi_l, m_l=m_l,
             band_c=band_c, band_l=band_l,
             sp_c=np.array(sp_c), sp_l=np.array(sp_l), std_c=std_c, std_l=std_l)
    print("wrote", os.path.join(DATADIR, "c2_data.npz"))


if __name__ == "__main__":
    main()
