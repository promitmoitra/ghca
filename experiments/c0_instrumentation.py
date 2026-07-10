"""
C0 - Instrument the C-series variables and operators (gate).

Confirms, on a recurrent GH substrate, that:
  1. W = f(S) is a deterministic coarse-graining of the FULL node state;
  2. under PARTIAL observation S_obs, the wave W carries predictive information
     about behaviour B beyond the observed spikes -- the analog of the paper's
     Eq. (1), P(B | S_obs, W) != P(B | S_obs);
  3. given the FULL state, W adds nothing (it is a function of S) -- so the
     question is only non-trivial under partial observation.

Substrate ("c-net"): a recurrent small-world GH medium of hidden nodes with
spontaneous firing (variable trajectories across trials). We define TWO
behaviours, both median-thresholded linear readouts over ALL hidden nodes:
  * B_collective -- a near-uniform readout (~ total activity): a COLLECTIVE code
    that the wave captures;
  * B_labeled    -- a structured zero-mean readout: a LABELED-LINE code that
    depends on WHICH nodes fire, not how many.
This lets C0 show the key instrumentation fact: whether W carries behavioural
information beyond partial spikes depends on whether behaviour reads a collective
or a labeled-line code -- the substrate hosts both, and which regime you are in
is what later determines the causal verdict (the paper's core message).

Predictive information is measured operationally, by decoder accuracy: a decoder
that also sees W should beat one that sees only S_obs for the COLLECTIVE code,
but not for the labeled-line code, and never beat one that sees the full S.

Outputs
-------
docs/figures/c0_predictive_info.png : decoder accuracy by feature set
result/c0/c0_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import cross_val_score
from sklearn.preprocessing import StandardScaler
from sklearn.pipeline import make_pipeline

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network, smallworld
from ghca_causal import (wave_coherence, wave_active_fraction,
                         make_observation_mask, read_spikes)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c0")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_H = 120           # hidden nodes (the full spike population S)
OBS_FRAC = 0.35     # fraction observed as S_obs
ACT, PAS = 2, 8
THETA = 1.0
P_S = 0.02          # spontaneous firing -> variable trajectories
SETTLE = 40
N_TRIALS = 1500


def build_cnet(seed=0):
    W = smallworld(N_H, k=6, beta=0.15, seed=seed)
    net = Network(W, act=ACT, pas=PAS, theta=THETA, p_s=P_S, seed=seed)
    rng = np.random.default_rng(seed + 1)
    # two readouts over ALL hidden nodes:
    w_coll = 1.0 + 0.1 * rng.standard_normal(N_H)     # near-uniform: collective
    w_lab = rng.standard_normal(N_H)                  # structured: labeled-line
    w_lab -= w_lab.mean()
    return net, w_coll, w_lab


def collect(seed=0):
    net, w_coll, w_lab = build_cnet(seed)
    rng = np.random.default_rng(seed + 7)
    S_full, Wc, Wa, r_coll, r_lab = [], [], [], [], []
    for _ in range(N_TRIALS):
        net.phi[:] = 0
        idx = rng.permutation(N_H)
        net.phi[idx[:int(0.1 * N_H)]] = rng.integers(1, net.act + 1, size=int(0.1 * N_H))
        for _ in range(SETTLE):
            net.step(None)
        a = net.active_mask().astype(float)
        S_full.append(a.copy())
        Wc.append(wave_coherence(net))
        Wa.append(wave_active_fraction(net))
        r_coll.append(float(w_coll @ a))
        r_lab.append(float(w_lab @ a))
    r_coll, r_lab = np.array(r_coll), np.array(r_lab)
    B_coll = (r_coll >= np.median(r_coll)).astype(int)   # balanced by median split
    B_lab = (r_lab >= np.median(r_lab)).astype(int)
    return np.array(S_full), np.array(Wc), np.array(Wa), B_coll, B_lab


def decode(X, y):
    """Cross-validated accuracy of a logistic decoder of y from X."""
    clf = make_pipeline(StandardScaler(), LogisticRegression(max_iter=2000, C=1.0))
    return cross_val_score(clf, X, y, cv=5).mean()


def main():
    S_full, Wc, Wa, B_coll, B_lab = collect(seed=0)
    W2 = np.c_[Wc, Wa]
    print(f"trials={len(B_coll)}  N_H={N_H}  "
          f"P(B_coll=1)={B_coll.mean():.2f} P(B_lab=1)={B_lab.mean():.2f}  "
          f"W_coh in [{Wc.min():.2f},{Wc.max():.2f}]  "
          f"W_act in [{Wa.min():.2f},{Wa.max():.2f}]")

    # determinism check: W adds ~nothing beyond the FULL state
    gain_full_coll = decode(np.c_[S_full, W2], B_coll) - decode(S_full, B_coll)
    gain_full_lab = decode(np.c_[S_full, W2], B_lab) - decode(S_full, B_lab)
    print(f"\nW=f(S) determinism check (gain beyond FULL spikes):"
          f" collective={gain_full_coll:+.3f}  labeled={gain_full_lab:+.3f}"
          f"  (both ~0 expected)")

    # sweep observation sparsity: W's gain beyond partial spikes vs |S_obs|
    fracs = [0.05, 0.10, 0.20, 0.35, 0.60, 1.00]
    rng = np.random.default_rng(0)
    gain = {"collective": [], "labeled-line": []}
    accS = {"collective": [], "labeled-line": []}
    for f in fracs:
        k = max(1, int(round(f * N_H)))
        cols = rng.choice(N_H, size=k, replace=False)
        Sobs = S_full[:, cols]
        for name, B in [("collective", B_coll), ("labeled-line", B_lab)]:
            a_s = decode(Sobs, B)
            a_sw = decode(np.c_[Sobs, W2], B)
            accS[name].append(a_s)
            gain[name].append(a_sw - a_s)
        print(f"  |S_obs|={k:3d}/{N_H}: "
              f"coll gain={gain['collective'][-1]:+.3f}  "
              f"lab gain={gain['labeled-line'][-1]:+.3f}")

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    x = [int(round(f * N_H)) for f in fracs]
    for name, c in [("collective", "crimson"), ("labeled-line", "steelblue")]:
        axes[0].plot(x, gain[name], "o-", color=c, label=name)
        axes[1].plot(x, accS[name], "o-", color=c, label=f"{name}: S_obs only")
    axes[0].axhline(0, color="k", alpha=0.5)
    axes[0].set_xlabel("|S_obs| (observed nodes)")
    axes[0].set_ylabel("W's decode gain beyond partial spikes")
    axes[0].set_title("Wave adds info beyond partial spikes\n(collective code only)")
    axes[0].legend()
    axes[1].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[1].set_xlabel("|S_obs| (observed nodes)")
    axes[1].set_ylabel("decode accuracy from S_obs")
    axes[1].set_title("Partial-spike decoding vs observation size")
    axes[1].legend()
    fig.suptitle("C0: W = f(S) is informative beyond PARTIAL spikes for a "
                 "collective code; the gain grows as observation gets sparser")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "c0_predictive_info.png")
    fig.savefig(p, dpi=110)
    print("\nwrote", p)

    np.savez(os.path.join(DATADIR, "c0_data.npz"),
             fracs=np.array(fracs), obs_sizes=np.array(x),
             gain_coll=np.array(gain["collective"]),
             gain_lab=np.array(gain["labeled-line"]),
             accS_coll=np.array(accS["collective"]),
             accS_lab=np.array(accS["labeled-line"]),
             gain_full_coll=gain_full_coll, gain_full_lab=gain_full_lab,
             Wc=Wc, Wa=Wa, B_coll=B_coll, B_lab=B_lab)
    print("wrote", os.path.join(DATADIR, "c0_data.npz"))


if __name__ == "__main__":
    main()
