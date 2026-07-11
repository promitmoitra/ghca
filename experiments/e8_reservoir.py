"""
E8.6 - An order-preserving reservoir (delay-line grid) for deep positional recall.

E8's per-channel traces store the RECENCY of each tone, so when a tone repeats the
older occurrence is overwritten and the exact tone AT LAG K is not recoverable --
which capped long-range recall (e8_results.md, Result 2). E8.6 fixes the
representation with an order-preserving reservoir: a delay-line grid of M tone
channels x L delay stages (a directed shift register per channel; a synfire/
delay-line motif on the GH substrate). Presenting tone c injects a token at stage 0
of channel c; every tone the tokens shift one stage down their channel. So the grid
cell (c, k) active means "tone c occurred ~k tones ago" -- an EXACT, ordered history
of the last ~L tones, with no overwrite when a channel repeats.

We compare the two reservoirs on:
  A. Lag-K recall: decode the tone from exactly K tones ago. The order-preserving
     grid recovers it up to its depth; the recency trace degrades with K (repeats).
  B. Delayed-copy prediction: predict x[t] = x[t-K] (a K-back copy). The grid
     predicts it once its depth reaches K; the trace does not.

Outputs
-------
docs/figures/e8_reservoir.png
result/e8/e8_reservoir.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network
import e8_predictive as e8

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

M = e8.M


# ----------------------------------------------------------------------------
# Order-preserving reservoir: a directed shift-register grid (M channels x L stages).
# ----------------------------------------------------------------------------

def build_shiftreg(L):
    N = M * L
    W = np.zeros((N, N))
    for c in range(M):
        for k in range(L - 1):
            W[c * L + (k + 1), c * L + k] = 1.0      # (c,k) -> (c,k+1): shift down channel
    net = Network(W, act=1, pas=0, theta=1.0, p_s=0.0, seed=0)  # tau=act=1: gapless shift
    return net


def feats_shiftreg(seq, L):
    """Grid active-mask each tone: cell (c,k) ~ 'tone c occurred k tones ago'."""
    net = build_shiftreg(L)
    F = []
    for t in range(len(seq)):
        net.phi[int(seq[t]) * L + 0] = 1             # inject current tone at stage 0
        net.step(None)                               # shift one stage
        F.append(np.concatenate([net.active_mask().astype(float), [1.0]]))
    return np.array(F)


def feats_trace(seq, tau):
    """E8's per-channel recency trace (for comparison)."""
    net = Network(np.zeros((M, M)), act=e8.ACT, pas=tau - e8.ACT, theta=99.0, p_s=0.0, seed=0)
    net.tau[:] = tau
    F = []
    for t in range(len(seq)):
        c = int(seq[t])
        for s in range(e8.ISI):
            if s == 0:
                net.phi[c] = 1
            net.step(None)
        onehot = np.zeros(M); onehot[c] = 1.0
        F.append(np.concatenate([net.phi.astype(float) / tau, onehot, [1.0]]))
    return np.array(F)


def decode_lag(F, seq, K, split=0.5):
    """Decode tone[t-K] from the reservoir feature at t."""
    idx = np.arange(K, len(F))
    X = F[idx]; y = seq[idx - K]
    k = int(split * len(X))
    W = e8.ridge_fit(X[:k], y[:k])
    return float(np.mean((X[k:] @ W).argmax(1) == y[k:]))


def main():
    L = 10
    TAU = 26                       # trace depth ~ TAU/ISI ~ 8 tones (matched-ish to L)
    rng = np.random.default_rng(0)
    seq = rng.integers(0, M, size=6000)         # i.i.d. tones -> isolate positional memory

    Fs = feats_shiftreg(seq, L)
    Ft = feats_trace(seq, TAU)
    Ks = list(range(0, 9))
    rec_s = [decode_lag(Fs, seq, K) for K in Ks]
    rec_t = [decode_lag(Ft, seq, K) for K in Ks]
    print(f"E8.6 lag-K recall (chance={1/M:.3f}); shift-register L={L} vs recency trace τ={TAU}:")
    for K, a, b in zip(Ks, rec_s, rec_t):
        print(f"  lag K={K}: shift-reg={a:.2f}   trace={b:.2f}")

    # depth sweep: recall of a fixed deep lag K=5 vs reservoir depth L (grid needs L>K)
    Kfix = 5
    Ls = [2, 4, 6, 8, 12]
    dep_s = [decode_lag(feats_shiftreg(seq, Lg), seq, Kfix) for Lg in Ls]
    dep_t = decode_lag(Ft, seq, Kfix)      # trace recall of lag 5 (depth-independent baseline)
    print(f"E8.6 recall of lag K={Kfix} vs grid depth L (trace baseline={dep_t:.2f}):")
    for Lg, a in zip(Ls, dep_s):
        print(f"  L={Lg:2d}: shift-reg={a:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    axes[0].plot(Ks, rec_s, "o-", color="seagreen", label=f"order-preserving grid (L={L})")
    axes[0].plot(Ks, rec_t, "s--", color="crimson", label=f"recency trace (τ={TAU})")
    axes[0].axhline(1 / M, ls=":", color="k", alpha=0.5, label=f"chance={1/M:.2f}")
    axes[0].set_xlabel("lag K (tones ago)"); axes[0].set_ylabel("recall accuracy of tone[t−K]")
    axes[0].set_ylim(0, 1.05); axes[0].legend(fontsize=8)
    axes[0].set_title("A. Positional recall: the grid keeps deep order;\n"
                      "the recency trace loses it as K grows")

    axes[1].plot(Ls, dep_s, "o-", color="seagreen", label="order-preserving grid")
    axes[1].axhline(dep_t, ls="--", color="crimson", label=f"recency trace (τ={TAU})")
    axes[1].axvline(Kfix, ls=":", color="k", alpha=0.6, label=f"lag K={Kfix}")
    axes[1].axhline(1 / M, ls=":", color="gray", alpha=0.5, label=f"chance={1/M:.2f}")
    axes[1].set_xlabel("reservoir depth L (stages)")
    axes[1].set_ylabel(f"recall accuracy of tone[t−{Kfix}]"); axes[1].set_ylim(0, 1.05)
    axes[1].set_title("B. Recall of a fixed deep lag needs grid depth\nL > K; the trace never recovers it")
    axes[1].legend(fontsize=8)

    fig.suptitle("E8.6: an order-preserving reservoir (delay-line grid) recovers deep "
                 "positional history the recency trace loses", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "e8_reservoir.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_reservoir.npz"),
             Ks=np.array(Ks), recall_shiftreg=np.array(rec_s), recall_trace=np.array(rec_t),
             Ls=np.array(Ls), depth_shiftreg=np.array(dep_s), depth_trace=dep_t,
             Kfix=Kfix, L=L, tau=TAU)
    print("wrote", os.path.join(DATADIR, "e8_reservoir.npz"))


if __name__ == "__main__":
    main()
