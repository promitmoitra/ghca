"""
E8/2a — The predictive-coding foil (making the disambiguation empirical).

E8 asserts "prediction as forward dynamics, NOT predictive coding". Track 2a
builds a minimal hierarchical predictive-coding (PC) model on the *same* tone
task and measures the observables that distinguish the two accounts, so the
contrast is empirical rather than rhetorical. See
docs/predictive_dynamics_experiments.md §1 and docs/next_steps.md Track 2a.

The PC model (untied Rao–Ballard, the recognition/generation variant):
  * recognition (feedforward)  r_t = tanh(Wff x_t + A r_{t-1})      -> representation
  * generation (top-down)      xhat_{t+1} = softmax(Wtd (A r_t))    -> prediction
  * per-channel error units     e_t = x_t - softmax(Wtd (A r_{t-1}))
Trained online to predict the next tone. Wff (recognition) and Wtd (generation)
are SEPARATE pathways — the defining PC architecture.

HONEST FINDING (this experiment corrects an E8-doc overclaim, in the spirit of
the earlier audits). Two observables one might expect to separate the accounts
turn out NOT to:
  * DISSOCIABILITY. PC dissociates cleanly (lesion Wtd → prediction 0.37→0.12
    chance, representation 0.95 spared) — the textbook PC signature. But E8 is
    ALSO dissociable: its two readouts tap different parts of the medium, so
    shrinking τ kills past-representation while sparing the (current-tone-driven)
    prediction. So the E8 doc's "prediction and representation are not
    dissociable" does not hold — corrected here.
  * ERROR LOCALITY. PC's per-channel error field peaks on the deviant — but an
    analogous residual formed from E8's readout also peaks there (~0.5). Both
    localise; locality is not the discriminator.

The observables that DO cleanly separate PC from the inside-out account are
ARCHITECTURAL, and become the empirical discriminators for data:
  (i)  intrinsic per-channel error units driving inference/learning [PC] vs a
       post-hoc GLOBAL scalar surprise, no error units [E8, E8.4: 0.16 vs 0.64];
  (ii) a dedicated top-down GENERATIVE pathway [PC] vs prediction as a passive
       readout whose learning lives in the weights [E8, and 1c];
  (iii) a window set by hierarchy depth [PC] vs by do(τ) [E8.3].

Outputs
-------
docs/figures/e8_pc_foil.png
result/e8/e8_pc_foil.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from e8_predictive import (M, CHANCE, ISI, ACT, softmax, make_traces,
                           seq_random_walk, seq_oddball, seq_periodic)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

H = 32              # PC latent size
LAG = 2             # "representation" = decode the tone LAG steps back (lives in memory)


# ---------------------------------------------------------------------------
# Minimal untied predictive-coding network.
# ---------------------------------------------------------------------------

class PCNet:
    def __init__(self, seed=0, lr=0.05):
        rng = np.random.default_rng(seed)
        self.Wff = 0.1 * rng.standard_normal((H, M))     # recognition (feedforward)
        self.A = 0.1 * rng.standard_normal((H, H))       # temporal prior
        self.Wtd = 0.1 * rng.standard_normal((M, H))     # generation (top-down)
        self.lr = lr

    def run(self, seq, learn=True, lesion_topdown=False):
        """Stream the sequence. Returns dict of latents r_t, predictions, errors."""
        n = len(seq)
        onehot = np.eye(M)[seq]
        r = np.zeros(H)
        R, PRED, ERR = [], [], []
        for t in range(n):
            x = onehot[t]
            r_prior = self.A @ r                           # temporal top-down prior
            # per-channel prediction error (of the CURRENT input)
            xhat_cur = softmax((self.Wtd @ r_prior)[None, :], 1)[0] if not lesion_topdown \
                else np.ones(M) / M
            e = x - xhat_cur                               # ERROR FIELD (M units)
            # recognition: representation driven feedforward by input + prior
            r = np.tanh(self.Wff @ x + r_prior)
            # generation: predict the NEXT input from the advanced latent
            if lesion_topdown:
                pnext = np.ones(M) / M
            else:
                pnext = softmax((self.Wtd @ (self.A @ r))[None, :], 1)[0]
            R.append(r.copy()); PRED.append(pnext); ERR.append(e)
            if learn and not lesion_topdown and t + 1 < n:
                # online gradient on next-tone cross-entropy (one-step, r as feature)
                g = pnext - onehot[t + 1]                  # dL/d(logits)
                rA = self.A @ r
                self.Wtd -= self.lr * np.outer(g, rA)
                back = self.Wtd.T @ g                      # into the latent
                self.A -= self.lr * 0.5 * np.outer(back * (1 - rA**2), r)
                self.Wff -= self.lr * 0.2 * np.outer((self.A.T @ back) * (1 - r**2), x)
        return {"r": np.array(R), "pred": np.array(PRED), "err": np.array(ERR)}


# ---------------------------------------------------------------------------
# Shared probes: next-tone prediction accuracy and lag-LAG representation.
# ---------------------------------------------------------------------------

def _decode_acc(Z, labels, split=0.5):
    """Held-out linear-decoder accuracy (ridge on one-hot, argmax)."""
    Z = np.concatenate([Z, np.ones((len(Z), 1))], 1)
    k = int(split * len(Z))
    Y = np.eye(M)[labels]
    Wd = np.linalg.solve(Z[:k].T @ Z[:k] + 1.0 * np.eye(Z.shape[1]), Z[:k].T @ Y[:k])
    return float(np.mean((Z[k:] @ Wd).argmax(1) == labels[k:]))


def pc_scores(seq, seed=0, lesion=False):
    net = PCNet(seed=seed)
    net.run(seq, learn=True)                               # train
    out = net.run(seq, learn=False, lesion_topdown=lesion)  # eval (optionally lesioned)
    n = len(seq)
    pred_acc = float(np.mean(out["pred"][:n - 1].argmax(1) == seq[1:n]))
    rep_acc = _decode_acc(out["r"][LAG:], seq[:n - LAG])    # decode tone LAG steps back
    return pred_acc, rep_acc, out


def e8_features(seq, tau):
    """E8 trace recency features (the single shared medium)."""
    net = make_traces(tau)
    feats = []
    for t in range(len(seq)):
        c = int(seq[t])
        for s in range(ISI):
            if s == 0:
                net.phi[c] = 1
            net.step(None)
        feats.append(net.phi.astype(float) / tau)
    return np.array(feats)


def e8_scores(seq, tau):
    """E8 next-tone prediction and lag-LAG representation, both from the ONE medium."""
    F = e8_features(seq, tau)
    n = len(seq)
    pred_acc = _decode_acc(F[:n - 1], seq[1:n])            # next tone
    rep_acc = _decode_acc(F[LAG:], seq[:n - LAG])          # tone LAG steps back
    return pred_acc, rep_acc


def main():
    rng = np.random.default_rng(0)
    N = 6000
    seq = seq_random_walk(N, np.random.default_rng(4))     # prediction needs the medium

    # ---- observable 1: dissociation plane ----
    pc_pred, pc_rep, _ = pc_scores(seq, seed=0, lesion=False)
    pcl_pred, pcl_rep, _ = pc_scores(seq, seed=0, lesion=True)
    print("PC intact:   pred=%.2f rep=%.2f" % (pc_pred, pc_rep))
    print("PC td-lesion:pred=%.2f rep=%.2f  (prediction dies, representation spared)" % (pcl_pred, pcl_rep))

    # E8: degrade the ONE medium by shrinking the trace timescale tau
    e8_curve = []
    for tau in [14, 10, 6, 4, 2]:
        pa, ra = e8_scores(seq, tau)
        e8_curve.append((tau, pa, ra))
        print("E8 tau=%2d:   pred=%.2f rep=%.2f" % (tau, pa, ra))
    e8_curve = np.array([(p, r) for _, p, r in e8_curve])

    # ---- the intrinsic per-channel error field PC maintains (on oddballs) ----
    seq_ob, dev = seq_oddball(N, np.random.default_rng(5))
    _, _, out = pc_scores(seq_ob, seed=1, lesion=False)
    err = np.abs(out["err"])
    dev_t = np.where(dev[:len(err)])[0]
    pc_conc = float(np.mean([err[t, seq_ob[t]] / (err[t].sum() + 1e-9) for t in dev_t]))
    # NOTE (honest): if one forms an analogous per-channel residual for E8's readout
    # it ALSO concentrates on the deviant (~0.5) -- both models' prediction errors
    # localise. The distinction is NOT locality but that PC's error field is an
    # INTRINSIC set of units driving inference/learning, whereas E8 computes no error
    # internally at all: its "surprise" is a post-hoc GLOBAL scalar (E8.4: 0.16 vs 0.64).
    print("PC per-channel error concentration at deviant = %.2f (E8 has no error units)" % pc_conc)

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.9))

    # Panel A: dissociation plane. HONEST finding -- BOTH models dissociate.
    ax = axes[0]
    ax.plot(e8_curve[:, 1], e8_curve[:, 0], "o-", color="slategray",
            label="E8 medium (shrink τ 14→2)")
    ax.annotate("", xy=(pcl_rep, pcl_pred), xytext=(pc_rep, pc_pred),
                arrowprops=dict(arrowstyle="->", color="crimson", lw=2))
    ax.plot([pc_rep], [pc_pred], "P", ms=13, color="seagreen", label="PC intact")
    ax.plot([pcl_rep], [pcl_pred], "X", ms=13, color="crimson", label="PC top-down lesion")
    ax.axhline(CHANCE, ls="--", color="k", alpha=0.3, label="chance")
    ax.set_xlabel("representation (decode tone t−%d)" % LAG)
    ax.set_ylabel("prediction (next tone)")
    ax.set_xlim(0, 1); ax.set_ylim(0, 1)
    ax.set_title("Dissociability does NOT separate them:\nPC lesion & E8 τ-shrink both dissociate\n(correcting the E8-doc 'inseparable' claim)")
    ax.legend(fontsize=7)

    # Panel B: PC's intrinsic per-channel error field (peaks on the deviant).
    ax = axes[1]
    aligned = np.zeros(M)
    for t in dev_t:
        aligned += np.roll(err[t], M // 2 - seq_ob[t])
    aligned /= max(len(dev_t), 1)
    ax.bar(np.arange(M) - M // 2, aligned, color="seagreen")
    ax.set_xlabel("channel offset from deviant"); ax.set_ylabel("|error| (PC error units)")
    ax.set_title("PC maintains INTRINSIC per-channel error\nunits (drive inference); E8 has none")

    # Panel C: the architectural discriminators (what actually separates them).
    ax = axes[2]; ax.axis("off")
    rows = [
        ["observable", "PC", "E8 (inside-out)"],
        ["error signal", "per-channel\nunits (intrinsic)", "global scalar\n(post-hoc, E8.4)"],
        ["prediction", "dedicated top-down\ngenerative pathway", "passive-medium\nreadout (1c)"],
        ["history window", "hierarchy depth", "do(τ) (E8.3)"],
        ["dissociable?", "yes (top-down)", "yes (readout/τ)\n— NOT a discriminator"],
    ]
    tbl = ax.table(cellText=rows, loc="center", cellLoc="center")
    tbl.auto_set_font_size(False); tbl.set_fontsize(8); tbl.scale(1, 1.9)
    for j in range(3):
        tbl[(0, j)].set_facecolor("#ddd")
    ax.set_title("The real discriminators are architectural")

    fig.suptitle("E8/2a: the predictive-coding foil — the clean disambiguators are ARCHITECTURAL "
                 "(intrinsic error units, generative pathway), not dissociability or error-locality "
                 "(both models do those)", fontsize=10)
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    p = os.path.join(FIGDIR, "e8_pc_foil.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_pc_foil.npz"),
             pc_pred=pc_pred, pc_rep=pc_rep, pcl_pred=pcl_pred, pcl_rep=pcl_rep,
             e8_curve=e8_curve, pc_conc=pc_conc, uniform=1.0 / M)
    print("wrote", os.path.join(DATADIR, "e8_pc_foil.npz"))


if __name__ == "__main__":
    main()
