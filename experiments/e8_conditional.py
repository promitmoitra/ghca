"""
E8.7 - Conditional long-range prediction: positional memory AND conjunction together.

E8.3's long-range recurrence x[t] = (x[t-1] + x[t-K]) mod M was only muddily
predictable (~0.33->0.67), for two independent reasons the later experiments each
isolated:
  * E8.6: the per-channel recency trace cannot recall the EXACT tone at lag K
    (positional memory) -- the order-preserving delay-line grid can.
  * E8.5: predicting a value that is a (context-dependent) TRANSFORM of two symbols
    is not linearly separable from the concatenation -- it needs the CONJUNCTION
    (interaction) of the two, the same feature E5's conjunction cells compute.

E8.7 puts them together and shows BOTH are necessary and jointly sufficient. The
recurrence needs (a) positional recall of x[t] and x[t+1-K] and (b) their sum mod M
= a pairwise conjunction. We cross the two reservoirs {order-preserving grid,
recency trace} with {conjunction feature on/off} on the K-back mod-M recurrence:

  only  grid + conjunction  solves it (~1.0 for depth L > K); the other three fail.

and sweep the grid depth to show the clean threshold at L > K that E8.3 lacked.

Outputs
-------
docs/figures/e8_conditional.png
result/e8/e8_conditional.npz
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
import e8_predictive as e8
import e8_reservoir as er

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

M = e8.M


def reservoir_state(seq, kind, param):
    """Per-tone reservoir state (no bias column)."""
    if kind == "grid":
        return er.feats_shiftreg(seq, param)[:, :-1]        # (n, M*L)
    return er.feats_trace(seq, param)[:, :-1]               # (n, 2M): recency + current


def build_features(states, seq, conj):
    """feature[t] = [reservoir state, current-tone one-hot, (optional current x state
    conjunction), bias]. The conjunction is the current tone crossed with the whole
    reservoir -- it supplies (x[t] (x) x[t-lag]) pairs for every lag, from which a
    linear readout can pick the lag it needs and compute the mod-M sum."""
    n = len(states)
    feats = []
    for t in range(n):
        cur = np.zeros(M); cur[int(seq[t])] = 1.0
        parts = [states[t], cur]
        if conj:
            parts.append(np.outer(cur, states[t]).ravel())
        parts.append([1.0])
        feats.append(np.concatenate(parts))
    return np.array(feats)


def eval_condition(kind, param, conj, K, n=6000, seed=0, split=0.5):
    """i.i.d. carrier x; predict target[t] = (x[t] + x[t-K]) mod M. The carrier is
    i.i.d. (high entropy, no cycle to memorise), so the target genuinely requires
    (a) positional recall of x[t-K] and (b) the conjunction to form the mod-M sum --
    it cannot be fit by a linear readout on separate one-hots (additive scores cannot
    realise a modular sum) nor from a reservoir that loses the lag-K symbol."""
    rng = np.random.default_rng(seed)
    x = rng.integers(0, M, size=n)
    target = (x + np.roll(x, K)) % M                 # target[t] = (x[t] + x[t-K]) % M
    states = reservoir_state(x, kind, param)
    F = build_features(states, x, conj)
    idx = np.arange(K, len(F))                       # drop the wrapped first K
    X, y = F[idx], target[idx]
    k = int(split * len(X))
    W = e8.ridge_fit(X[:k], y[:k])
    return float(np.mean((X[k:] @ W).argmax(1) == y[k:]))


def main():
    K = 4
    L = 10
    TAU = 26
    conds = [
        ("grid + conj", "grid", L, True),
        ("grid only", "grid", L, False),
        ("trace + conj", "trace", TAU, True),
        ("trace only", "trace", TAU, False),
    ]
    accs = {}
    for label, kind, param, conj in conds:
        a = np.mean([eval_condition(kind, param, conj, K, seed=s) for s in range(3)])
        accs[label] = a
    print(f"E8.7 conditional long-range recurrence (K={K}, chance={1/M:.3f}):")
    for label, *_ in conds:
        print(f"  {label:14s} acc={accs[label]:.2f}")

    # depth sweep for grid+conj: threshold at L > K
    Ls = [2, 4, 6, 8, 12]
    depth = [np.mean([eval_condition("grid", Lg, True, K, seed=s) for s in range(3)]) for Lg in Ls]
    print(f"E8.7 grid+conj accuracy vs grid depth L (K={K}):")
    for Lg, a in zip(Ls, depth):
        print(f"  L={Lg:2d}: acc={a:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    labels = [c[0] for c in conds]
    colors = ["seagreen", "goldenrod", "indianred", "gray"]
    axes[0].bar(range(len(labels)), [accs[l] for l in labels], color=colors)
    axes[0].axhline(1 / M, ls=":", color="k", alpha=0.6, label=f"chance={1/M:.2f}")
    axes[0].set_xticks(range(len(labels))); axes[0].set_xticklabels(labels, rotation=20, fontsize=8)
    axes[0].set_ylabel("prediction accuracy"); axes[0].set_ylim(0, 1.05)
    axes[0].set_title(f"A. Only order-preserving grid + conjunction\nsolves the K={K} mod-M recurrence")
    axes[0].legend(fontsize=8)

    axes[1].plot(Ls, depth, "o-", color="seagreen")
    axes[1].axvline(K, ls=":", color="k", alpha=0.6, label=f"needed depth K={K}")
    axes[1].axhline(1 / M, ls=":", color="gray", alpha=0.5, label=f"chance={1/M:.2f}")
    axes[1].set_xlabel("grid depth L"); axes[1].set_ylabel("grid+conj accuracy"); axes[1].set_ylim(0, 1.05)
    axes[1].set_title("B. Clean threshold at L > K\n(the sharp step E8.3 lacked)")
    axes[1].legend(fontsize=8)

    fig.suptitle("E8.7: conditional long-range prediction needs BOTH positional memory "
                 "(order-preserving grid) AND the conjunction feature", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "e8_conditional.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_conditional.npz"),
             labels=np.array(labels), accs=np.array([accs[l] for l in labels]),
             Ls=np.array(Ls), depth=np.array(depth), K=K, L=L, tau=TAU)
    print("wrote", os.path.join(DATADIR, "e8_conditional.npz"))


if __name__ == "__main__":
    main()
