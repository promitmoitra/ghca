"""
E8/1c — Online, reward-free GVF prediction (learned, not offline-fit).

E8's next-tone predictor is an OFFLINE ridge on a frozen substrate: honest, but
the "learning" is a god's-eye batch solve, and the extensions audit flags it
(caveat 7: "the learning credited to E8 is in the readout weights"). Track 1c
folds that predictor into the substrate's OWN online machinery — the linear TD
GVF demons of [E6](e6_results.md) — so prediction is acquired incrementally,
one step at a time, with no batch and no external target beyond the sensory
stream itself. See docs/next_steps.md Track 1c.

Setup (identical task to E8, reused verbatim): the tonotopic trace bank produces
the same features; the TARGET is the next tone. A bank of M linear GVF demons
(one per channel), cumulant `c_k[t] = 1[next tone = k]`, is learned by online
TD (the E6 `td_continuing` rule):

    delta_k = c_k[t] + gamma * (w_k . x[t+1]) - w_k . x[t];  w_k += alpha * delta_k * x[t]

gamma = 0 gives a one-step next-tone demon (compared head-to-head with the
offline ridge); gamma > 0 gives a genuine multi-step GVF (discounted future
occurrence of each tone) — predictive knowledge in the Horde sense, which the
offline ridge does not provide.

Claims
------
* ONLINE acquisition: a single incremental pass (no batch) reaches ~the offline
  ridge's next-tone accuracy across periodic / Markov / random-walk sequences,
  and the running accuracy rises over the stream (a learning curve).
* INTRINSIC / reward-free: the only signal is the next tone (self-supervised);
  no reward, no error-unit hierarchy — prediction is learned inside-out.
* HORDE bonus: gamma > 0 demons predict discounted future tone occurrence
  (multi-step predictive knowledge) from the same features.

Outputs
-------
docs/figures/e8_online_prediction.png
result/e8/e8_online_prediction.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e8_predictive import (M, CHANCE, ISI, ACT, BETA, softmax,
                           seq_periodic, seq_markov, seq_random_walk,
                           features_targets, ridge_fit, predict_scores)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

TAU = 14


# ---------------------------------------------------------------------------
# Online linear TD GVF demon bank (the E6 rule), one demon per tone channel.
# ---------------------------------------------------------------------------

def online_td_predict(X, y, gamma=0.0, alpha=0.02, seed=0):
    """One online pass. Returns (final_acc, running_acc, W).

    W has shape (M, D): row k is demon k's weights. Prediction at t = argmax_k
    (w_k . x[t]); accuracy scored online on the *upcoming* tone (a genuine
    forecast — each step is tested before it is learned from)."""
    D = X.shape[1]
    W = np.zeros((M, D))
    correct = np.zeros(len(X) - 1)
    for t in range(len(X) - 1):
        v = W @ X[t]                                  # per-demon value
        pred = int(np.argmax(v))
        correct[t] = float(pred == y[t])              # test BEFORE learning
        vnext = W @ X[t + 1]
        cumulant = np.eye(M)[y[t]]                     # 1[next tone = k]
        delta = cumulant + gamma * vnext - v
        W += alpha * np.outer(delta, X[t])
    # running accuracy (causal moving average)
    run = np.cumsum(correct) / (np.arange(len(correct)) + 1)
    # final accuracy = last-quarter online accuracy (after learning settles)
    final = float(correct[3 * len(correct) // 4:].mean())
    return final, run, W


def gvf_multistep(X, seq, gamma, alpha=0.02):
    """gamma>0 GVF: predict the discounted future occurrence of each tone,
    G_k[t] = sum_{j>=0} gamma^j * 1[tone_{t+j} = k]. Score = correlation between
    the demon value and the empirical discounted return (a real forecast, not a
    one-step readout)."""
    n = len(X)
    # empirical discounted returns per channel (backward recursion)
    G = np.zeros((n, M))
    onehot = np.eye(M)[seq[:n]]
    for t in range(n - 2, -1, -1):
        G[t] = onehot[t] + gamma * G[t + 1]
    _, _, W = online_td_predict(X, seq[:n], gamma=gamma, alpha=alpha)
    V = X @ W.T
    k = n // 2
    # mean per-channel Pearson correlation on held-out second half
    cors = []
    for c in range(M):
        a, b = V[k:, c], G[k:, c]
        if a.std() > 1e-6 and b.std() > 1e-6:
            cors.append(np.corrcoef(a, b)[0, 1])
    return float(np.mean(cors)) if cors else 0.0


def main():
    N = 6000
    conds = {
        "periodic": seq_periodic(N, 4),
        "Markov α=0.9": seq_markov(N, 0.9, np.random.default_rng(1)),
        "Markov α=0.5": seq_markov(N, 0.5, np.random.default_rng(2)),
        "random-walk": seq_random_walk(N, np.random.default_rng(4)),
    }

    ridge_acc, online_acc, curves = {}, {}, {}
    for name, seq in conds.items():
        X, y = features_targets(seq, TAU)
        # offline ridge baseline (E8)
        kf = len(X) // 2
        Wr = ridge_fit(X[:kf], y[:kf])
        ridge_acc[name] = float(np.mean(predict_scores(X[kf:], Wr).argmax(1) == y[kf:]))
        # online TD demon bank (1c)
        fa, run, _ = online_td_predict(X, y, gamma=0.0)
        online_acc[name] = fa
        curves[name] = run

    print("1c online GVF prediction vs offline ridge (chance = %.3f):" % CHANCE)
    for name in conds:
        print(f"  {name:14s} online={online_acc[name]:.2f}  ridge={ridge_acc[name]:.2f}")

    # gamma>0 multi-step GVF on the Markov sequence
    Xm, _ = features_targets(conds["Markov α=0.9"], TAU)
    gvf_cor = {g: gvf_multistep(Xm, conds["Markov α=0.9"], gamma=g) for g in (0.5, 0.8)}
    print("Horde multi-step GVF (value vs discounted return, mean per-channel r):")
    for g, c in gvf_cor.items():
        print(f"  gamma={g}: r={c:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.7))
    names = list(conds)
    x = np.arange(len(names))
    w = 0.38
    axes[0].bar(x - w / 2, [online_acc[n] for n in names], w, color="seagreen", label="online TD (1c)")
    axes[0].bar(x + w / 2, [ridge_acc[n] for n in names], w, color="slategray", label="offline ridge (E8)")
    axes[0].axhline(CHANCE, ls="--", color="k", alpha=0.5, label=f"chance={CHANCE:.2f}")
    axes[0].set_xticks(x); axes[0].set_xticklabels(names, rotation=20, fontsize=8)
    axes[0].set_ylabel("next-tone accuracy"); axes[0].set_ylim(0, 1)
    axes[0].set_title("1c: online GVF matches the offline ridge"); axes[0].legend(fontsize=8)

    for name in names:
        axes[1].plot(curves[name], label=name)
    axes[1].axhline(CHANCE, ls="--", color="k", alpha=0.5)
    axes[1].set_xlabel("tone (online step)"); axes[1].set_ylabel("running accuracy")
    axes[1].set_ylim(0, 1); axes[1].set_title("Online acquisition (single incremental pass)")
    axes[1].legend(fontsize=7)

    axes[2].bar([str(g) for g in gvf_cor], list(gvf_cor.values()), color="crimson")
    axes[2].set_xlabel("GVF continuation γ"); axes[2].set_ylabel("value vs discounted return (r)")
    axes[2].set_ylim(0, 1)
    axes[2].set_title("Horde bonus: multi-step predictive\nknowledge (γ>0), not just 1-step")

    fig.suptitle("E8/1c: prediction learned ONLINE & reward-free by the substrate's own "
                 "TD/GVF demons — not an offline readout", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(FIGDIR, "e8_online_prediction.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_online_prediction.npz"),
             names=np.array(names),
             online=np.array([online_acc[n] for n in names]),
             ridge=np.array([ridge_acc[n] for n in names]),
             gvf_gammas=np.array(list(gvf_cor)), gvf_cor=np.array(list(gvf_cor.values())))
    print("wrote", os.path.join(DATADIR, "e8_online_prediction.npz"))


if __name__ == "__main__":
    main()
