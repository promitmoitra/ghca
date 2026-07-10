"""
E8 - Predictive dynamics on naturalistic tone sequences.

Asks whether the substrate makes time-forward PREDICTIONS of upcoming input, and
how that prediction differs from predictive coding. See
docs/predictive_dynamics_experiments.md. Motivated by the auditory-sequence
prediction literature (Chait/Bianco: neural activity integrates ~the last 7 tones
to predict the next pitch; anticipatory tuning sharpens for more predictable
sequences).

Substrate (minimal, on the GH clock):
  * Tonotopic map of M tone channels; a tone = a clamp on its channel.
  * A per-channel TRACE bank -- each trace is a GH-clocked node: presenting tone c
    sets trace c to phase 1, and the phase then advances 1..tau and wraps to 0. So
    a tone is represented for ~tau steps after it occurs, with the phase encoding
    time-since-onset. The set-with-recency of non-rested traces is a tau-deep
    history window: the effective memory depth is a do(tau) property (Line B).
  * A next-tone PREDICTIVE READOUT (a GVF-style linear demon, ridge, self-
    supervised) reads the WHOLE trace state (not a single locus -- heeding C5:
    read the medium's integrated state, not a fixed locus) and predicts the
    upcoming tone. No reward, no error-unit hierarchy: prediction is the forward
    readout of the medium's own dynamics.

Results:
  E8.1 periodic     -> high next-tone prediction (anticipation).
  E8.2 Markov       -> prediction tracks transition predictability.
  E8.3 random-walk  -> the integration window (how many past tones inform the
                       prediction) grows with tau: do(tau) sets the history depth.
  E8.4 oddball      -> surprise is a GLOBAL scalar (prediction-error magnitude)
                       spiking on deviants, and anticipatory confidence sharpens
                       with predictability -- not a per-channel error field
                       (the predictive-coding disambiguation).

Outputs
-------
docs/figures/e8_predictive.png
result/e8/e8_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

M = 8                     # number of tone channels (tonotopic map)
ISI = 3                  # steps between successive tones
ACT = 2
CHANCE = 1.0 / M


# ----------------------------------------------------------------------------
# Sequence generators.
# ----------------------------------------------------------------------------

def seq_periodic(n, period=4, rng=None):
    base = np.arange(period) % M
    return np.array([base[t % period] for t in range(n)])


def seq_markov(n, alpha, rng):
    """Row-stochastic transitions biased to nearest neighbours by `alpha`
    (alpha=1 -> deterministic step; alpha=0 -> uniform/unpredictable)."""
    seq = [int(rng.integers(M))]
    for _ in range(n - 1):
        c = seq[-1]
        p = np.ones(M) * (1 - alpha) / M
        p[(c + 1) % M] += alpha
        p /= p.sum()
        seq.append(int(rng.choice(M, p=p)))
    return np.array(seq)


def seq_random_walk(n, rng, drift=(-1, 0, 1)):
    seq = [M // 2]
    for _ in range(n - 1):
        seq.append(int(np.clip(seq[-1] + rng.choice(drift), 0, M - 1)))
    return np.array(seq)


def seq_oddball(n, rng, standard_period=4, p_dev=0.12):
    seq = seq_periodic(n, standard_period)
    dev = rng.random(n) < p_dev
    seq[dev] = rng.integers(0, M, size=int(dev.sum()))
    return seq, dev


def seq_longrange(n, K, rng):
    """A long-range recurrence x[t] = (x[t-1] + x[t-K]) mod M. Deterministic given
    the last K values; predicting x[t] needs BOTH lag-1 and lag-K, so it is only
    predictable when the history window reaches K -- the do(tau) threshold test."""
    x = list(rng.integers(0, M, size=K))
    for t in range(K, n):
        x.append((x[t - 1] + x[t - K]) % M)
    return np.array(x)


# ----------------------------------------------------------------------------
# Trace reservoir + predictive readout.
# ----------------------------------------------------------------------------

def make_traces(tau):
    net = Network(np.zeros((M, M)), act=ACT, pas=tau - ACT, theta=99.0, p_s=0.0, seed=0)
    net.tau[:] = tau
    return net


def features_targets(seq, tau):
    """Run the sequence through the trace bank; feature = (trace recency vector,
    current-tone one-hot); target = the NEXT tone."""
    net = make_traces(tau)
    feats, targ = [], []
    for t in range(len(seq) - 1):
        c = int(seq[t])
        for s in range(ISI):
            if s == 0:
                net.phi[c] = 1              # ignite trace c
            net.step(None)                  # advance all trace clocks
        recency = net.phi.astype(float) / tau        # 0 = rested, else time-since
        onehot = np.zeros(M); onehot[c] = 1.0
        feats.append(np.concatenate([recency, onehot, [1.0]]))
        targ.append(int(seq[t + 1]))
    return np.array(feats), np.array(targ)


def ridge_fit(X, y, lam=1.0):
    Y = np.eye(M)[y]
    return np.linalg.solve(X.T @ X + lam * np.eye(X.shape[1]), X.T @ Y)


def predict_scores(X, W):
    return X @ W


BETA = 6.0                       # softmax temperature (calibrates confidence/surprise)


def softmax(z, beta=1.0):
    z = beta * (z - z.max(axis=1, keepdims=True))
    e = np.exp(z)
    return e / e.sum(axis=1, keepdims=True)


def eval_sequence(seq, tau, split=0.5):
    X, y = features_targets(seq, tau)
    n = len(X); k = int(split * n)
    W = ridge_fit(X[:k], y[:k])
    sc = predict_scores(X[k:], W)
    pred = sc.argmax(1)
    acc = float(np.mean(pred == y[k:]))
    conf = float(softmax(sc, BETA).max(1).mean())      # anticipatory confidence
    return acc, conf, X, y, k, W


# ----------------------------------------------------------------------------
# E8.3 integration window: deepest past-tone lag decodable from the trace state.
# ----------------------------------------------------------------------------

def integration_window(tau, n=4000, thresh=0.30, seed=0):
    rng = np.random.default_rng(seed)
    seq = rng.integers(0, M, size=n)              # i.i.d. tones -> isolate memory
    X, _ = features_targets(seq, tau)
    depth = 0
    accs = []
    for lag in range(1, 12):
        yl = seq[:len(X)]                          # tone at t-lag aligns below
        if lag >= len(X):
            break
        Xl, yy = X[lag:], seq[np.arange(len(X))][:-lag]  # feature at t vs tone t-lag
        kk = len(Xl) // 2
        W = ridge_fit(Xl[:kk], yy[:kk])
        a = float(np.mean((Xl[kk:] @ W).argmax(1) == yy[kk:]))
        accs.append(a)
        if a > thresh:
            depth = lag
    return depth, accs


def main():
    rng = np.random.default_rng(0)
    N = 6000

    # ---- E8.1 / E8.2: prediction by sequence type ----
    conds = {
        "periodic": seq_periodic(N, 4),
        "Markov α=0.9": seq_markov(N, 0.9, np.random.default_rng(1)),
        "Markov α=0.5": seq_markov(N, 0.5, np.random.default_rng(2)),
        "Markov α=0.1": seq_markov(N, 0.1, np.random.default_rng(3)),
        "random-walk": seq_random_walk(N, np.random.default_rng(4)),
    }
    acc_by = {}; conf_by = {}
    for name, seq in conds.items():
        a, c, *_ = eval_sequence(seq, tau=14)
        acc_by[name] = a; conf_by[name] = c
    print("E8.1/8.2 next-tone prediction (chance = %.3f):" % CHANCE)
    for name in conds:
        print(f"  {name:14s} acc={acc_by[name]:.2f}  anticipatory-conf={conf_by[name]:.2f}")

    # ---- E8.3: integration window vs tau (do(tau)) + a long-range threshold ----
    K_LR = 4                              # long-range dependency depth (needs window >= K)
    taus = [4, 8, 14, 20, 26]
    windows = []
    lr_acc = []
    for tau in taus:
        windows.append(float(np.mean([integration_window(tau, seed=s)[0] for s in range(3)])))
        lr_acc.append(float(np.mean([eval_sequence(seq_longrange(N, K_LR, np.random.default_rng(7 + s)), tau=tau)[0]
                                     for s in range(3)])))
    print("E8.3 do(tau): integration window (tones) and long-range (K=%d) accuracy:" % K_LR)
    for tau, w, a in zip(taus, windows, lr_acc):
        print(f"  tau={tau:2d}: window={w} tones  longrange_acc={a:.2f}")

    # ---- E8.4: oddball surprise (global) + anticipatory sharpening ----
    seq_ob, dev = seq_oddball(N, np.random.default_rng(5))
    X, y = features_targets(seq_ob, tau=14)
    dev = dev[:len(X)]
    k = len(X) // 2
    W = ridge_fit(X[:k], y[:k])
    P = softmax(predict_scores(X[k:], W), BETA)
    yt = y[k:]; devt = dev[k:]
    surprise = 1.0 - P[np.arange(len(yt)), yt]      # global scalar: 1 - p(actual)
    s_std = float(surprise[~devt].mean())
    s_dev = float(surprise[devt].mean())
    print(f"E8.4 oddball surprise (1 - p(actual)):  standard={s_std:.2f}  deviant={s_dev:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    names = list(conds.keys())
    axes[0].bar(range(len(names)), [acc_by[n] for n in names], color="steelblue")
    axes[0].axhline(CHANCE, ls="--", color="k", alpha=0.6, label=f"chance={CHANCE:.2f}")
    axes[0].set_xticks(range(len(names))); axes[0].set_xticklabels(names, rotation=25, fontsize=8)
    axes[0].set_ylabel("next-tone prediction accuracy"); axes[0].set_ylim(0, 1)
    axes[0].set_title("E8.1/8.2: prediction tracks predictability"); axes[0].legend(fontsize=8)

    ax = axes[1]; ax.plot(taus, windows, "o-", color="seagreen", label="window (tones)")
    ax.axhline(K_LR, ls=":", color="k", alpha=0.6, label=f"needed depth K={K_LR}")
    ax.set_xlabel("trace timescale τ  (do(τ))"); ax.set_ylabel("window (tones)", color="seagreen")
    ax.tick_params(axis="y", labelcolor="seagreen"); ax.legend(fontsize=7, loc="upper left")
    ax2 = ax.twinx(); ax2.plot(taus, lr_acc, "s--", color="crimson", label=f"long-range (K={K_LR}) acc")
    ax2.axhline(CHANCE, ls="--", color="crimson", alpha=0.3)
    ax2.set_ylabel("long-range accuracy", color="crimson"); ax2.set_ylim(0, 1)
    ax2.tick_params(axis="y", labelcolor="crimson")
    ax.set_title("E8.3: history depth is a do(τ) property\n(K-back predictable only once window ≥ K)")

    axes[2].bar(["standard", "deviant"], [s_std, s_dev], color=["gray", "crimson"])
    axes[2].set_ylabel("surprise  (1 − p(actual), global)")
    axes[2].set_title(f"E8.4: global surprise spikes on deviants\n"
                      f"(anticipatory conf: periodic {conf_by['periodic']:.2f} "
                      f"vs random-walk {conf_by['random-walk']:.2f})")

    fig.suptitle("E8: time-forward prediction as forward dynamics + a do(τ) history "
                 "window — no error-unit hierarchy (predictive-coding disambiguation)",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "e8_predictive.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_data.npz"),
             names=np.array(names), acc=np.array([acc_by[n] for n in names]),
             conf=np.array([conf_by[n] for n in names]),
             taus=np.array(taus), windows=np.array(windows), lr_acc=np.array(lr_acc),
             K_longrange=K_LR, surprise_std=s_std, surprise_dev=s_dev)
    print("wrote", os.path.join(DATADIR, "e8_data.npz"))


if __name__ == "__main__":
    main()
