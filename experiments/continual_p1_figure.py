"""Track 3c / P1 figure: interference matrices + summary metrics."""

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
    d = np.load(os.path.join(SD, "continual_p1.npz"))
    chance = float(d["chance"])
    fig, ax = plt.subplots(1, 3, figsize=(15, 4.6))

    for i, reg in enumerate(["frozen", "plastic"]):
        Rm = d[f"{reg}_Rmean"]
        im = ax[i].imshow(Rm, cmap="magma", vmin=0, vmax=1)
        ax[i].set_xticks(range(3)); ax[i].set_yticks(range(3))
        ax[i].set_xlabel("evaluated task"); ax[i].set_ylabel("after training task")
        ax[i].set_title(f"{reg}: accuracy matrix (n={int(d['n'])})")
        for a in range(3):
            for b in range(3):
                ax[i].text(b, a, f"{Rm[a,b]:.2f}", ha="center", va="center",
                           color="w" if Rm[a, b] < 0.6 else "k", fontsize=11)
        fig.colorbar(im, ax=ax[i], fraction=0.046)

    # summary: avg_acc and backward transfer per regime, mean +/- 95% CI
    labels = ["frozen", "plastic"]
    xs = np.arange(len(labels)); w = 0.35
    for j, (metric, col) in enumerate([("avg_acc", "steelblue"), ("bwt", "crimson")]):
        means, los, his = [], [], []
        for reg in labels:
            arr = d[f"{reg}_{metric}"]
            m, lo, hi = st.bootstrap_ci(arr)
            means.append(m); los.append(m - lo); his.append(hi - m)
        ax[2].bar(xs + (j - 0.5) * w, means, w,
                  yerr=[los, his], capsize=4, color=col,
                  label={"avg_acc": "avg accuracy", "bwt": "backward transfer"}[metric])
    ax[2].axhline(chance, ls=":", color="0.5", label=f"chance {chance:.2f}")
    ax[2].axhline(0, color="k", lw=0.8)
    ax[2].set_xticks(xs); ax[2].set_xticklabels(labels)
    ax[2].set_title("continual metrics (mean ± 95% CI)")
    ax[2].legend(fontsize=8)

    fig.suptitle("3c/P1 — catastrophic interference on one shared substrate "
                 "(3 conflicting remappings, correlational credit)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_p1.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
