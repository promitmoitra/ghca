"""Track 3c / P4 figure: the capacity ladder (credit fixed) resolves interference."""

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
RUNGS = ["shared", "context", "per-task"]
NICE = {"shared": "shared head\n(stimulus only)", "context": "shared head\n+ task context",
        "per-task": "per-task\nheads"}


def main():
    d = np.load(os.path.join(SD, "continual_capacity.npz"))
    chance = float(d["chance"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    xs = np.arange(len(RUNGS)); w = 0.36
    for j, (metric, col, lab) in enumerate([("avg_acc", "steelblue", "avg accuracy"),
                                            ("bwt", "crimson", "backward transfer")]):
        means, los, his = [], [], []
        for r in RUNGS:
            m, lo, hi = st.bootstrap_ci(d[f"{r}_{metric}"])
            means.append(m); los.append(m - lo); his.append(hi - m)
        a0.bar(xs + (j - 0.5) * w, means, w, yerr=[los, his], capsize=4, color=col, label=lab)
    a0.axhline(chance, ls=":", color="0.5", label=f"chance {chance:.2f}")
    a0.axhline(0, color="k", lw=0.8)
    a0.set_xticks(xs); a0.set_xticklabels([NICE[r] for r in RUNGS])
    a0.set_title("Capacity ladder (credit rule fixed = correlational)")
    a0.legend(fontsize=8)

    # accuracy matrices per rung
    for i, r in enumerate(RUNGS):
        Rm = d[f"{r}_Rmean"]
        ax = a1.inset_axes([i / 3.0 + 0.02, 0.15, 0.28, 0.6])
        im = ax.imshow(Rm, cmap="magma", vmin=0, vmax=1)
        ax.set_xticks([]); ax.set_yticks([])
        ax.set_title(r, fontsize=9)
        for u in range(2):
            for v in range(2):
                ax.text(v, u, f"{Rm[u,v]:.2f}", ha="center", va="center",
                        color="w" if Rm[u, v] < 0.6 else "k", fontsize=8)
    a1.axis("off")
    a1.set_title("accuracy matrix R[after, eval]  (rows=after training task, cols=eval task)",
                 fontsize=9)

    fig.suptitle("3c/P4 — capacity, not credit: task-context (conjunctive) and per-task "
                 "heads resolve interference (K=2, n=30)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_capacity.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
