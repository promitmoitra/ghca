"""Track 3c / P2 figure: the apparent causal-credit effect, and the frontier control
that explains it as an effective-learning-rate difference (honest null)."""

import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import ghca_stats as st  # noqa: E402
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SD = os.path.join(ROOT, "result", "stats")
FIG = os.path.join(ROOT, "docs", "figures")


def main():
    p2 = np.load(os.path.join(SD, "continual_p2.npz"))
    fr = np.load(os.path.join(SD, "continual_p2_frontier.npz"))
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    # left: the apparent headline effect (backward transfer at fixed lr)
    labels = ["correlational", "causal do(θ)"]
    keys = ["frozen_correlational_bwt", "frozen_causal_bwt"]
    means, los, his = [], [], []
    for k in keys:
        m, lo, hi = st.bootstrap_ci(p2[k])
        means.append(m); los.append(m - lo); his.append(hi - m)
    a0.bar([0, 1], means, 0.5, yerr=[los, his], capsize=5,
           color=["crimson", "steelblue"])
    a0.set_xticks([0, 1]); a0.set_xticklabels(labels)
    a0.axhline(0, color="k", lw=0.8)
    a0.set_ylabel("backward transfer (forgetting)")
    a0.set_title("Apparent effect: causal forgets less (fixed lr)")

    # right: the control — stability-plasticity frontier (retention vs new-task acq)
    corr, caus = fr["corr"], fr["caus"]     # columns: [new-task acq, retention]
    a1.plot(corr[:, 0], corr[:, 1], "o-", color="crimson", label="correlational (sweep η)")
    a1.plot(caus[:, 0], caus[:, 1], "s-", color="steelblue", label="causal (sweep lr)")
    # headline operating points
    a1.scatter([0.82], [0.18], s=140, facecolors="none", edgecolors="crimson",
               linewidths=2, zorder=5, label="corr @ default η")
    a1.scatter([0.71], [0.29], s=140, facecolors="none", edgecolors="steelblue",
               linewidths=2, zorder=5, label="causal @ lr=2")
    a1.set_xlabel("new-task acquisition (R[1,1])")
    a1.set_ylabel("retention of task 0 (R[1,0])")
    a1.set_title("Control: same frontier → effect is a learning-rate artifact")
    a1.legend(fontsize=8)

    fig.suptitle("3c/P2 — causal do(θ) credit does NOT beat correlational credit at "
                 "matched plasticity (K=2 reversal, n=30)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_p2.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
