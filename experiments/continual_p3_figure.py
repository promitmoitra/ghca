"""Track 3c / P3 figure: three credit rules trace one stability-plasticity frontier."""

import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SD = os.path.join(ROOT, "result", "stats")
FIG = os.path.join(ROOT, "docs", "figures")
STYLE = {"correlational": ("crimson", "o", "correlational (sweep η)"),
         "causal": ("steelblue", "s", "causal do(θ) (sweep lr)"),
         "causal_lowvar": ("seagreen", "^", "causal low-variance (sweep lr)")}


def main():
    d = np.load(os.path.join(SD, "continual_p3_frontier.npz"))
    fig, ax = plt.subplots(figsize=(8, 6))
    for name, (col, mk, lab) in STYLE.items():
        pts = d[f"{name}_pts"]
        order = np.argsort(pts[:, 0])
        ax.plot(pts[order, 0], pts[order, 1], mk + "-", color=col, label=lab, alpha=0.85)
    ax.axhline(0.5, ls=":", color="0.6", lw=0.8)
    ax.set_xlabel("new-task acquisition  R[1,1]")
    ax.set_ylabel("retention of task 0  R[1,0]")
    ax.set_title("3c/P3 — three credit rules trace ONE stability–plasticity frontier\n"
                 "(low-variance causal credit slides along it, does not move it) — K=2, n=15",
                 fontsize=11)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_p3.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
