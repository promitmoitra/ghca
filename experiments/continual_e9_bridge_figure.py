"""E9 ↔ 3c bridge figure: a learned conjunction basis resolves interference."""

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
BASES = ["frozen", "emergent", "wired"]
NICE = {"frozen": "frozen\n(no conjunction)", "emergent": "emergent\n(learned conj.)",
        "wired": "wired\n(hand conj.)"}


def main():
    d = np.load(os.path.join(SD, "continual_e9_bridge.npz"))
    chance = float(d["chance"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    xs = np.arange(len(BASES)); w = 0.36
    for j, (metric, col, lab) in enumerate([("avg_acc", "steelblue", "avg accuracy"),
                                            ("bwt", "crimson", "backward transfer")]):
        means, los, his = [], [], []
        for b in BASES:
            m, lo, hi = st.bootstrap_ci(d[f"{b}_{metric}"])
            means.append(m); los.append(m - lo); his.append(hi - m)
        a0.bar(xs + (j - 0.5) * w, means, w, yerr=[los, his], capsize=4, color=col, label=lab)
    a0.axhline(chance, ls=":", color="0.5", label=f"chance {chance:.2f}")
    a0.axhline(0, color="k", lw=0.8)
    a0.set_xticks(xs); a0.set_xticklabels([NICE[b] for b in BASES])
    a0.set_title("Continual reversal on the E9 substrate (readout = Line A, fixed)")
    a0.legend(fontsize=8)

    for i, b in enumerate(BASES):
        Rm = d[f"{b}_Rmean"]
        ax = a1.inset_axes([i / 3.0 + 0.02, 0.15, 0.28, 0.62])
        ax.imshow(Rm, cmap="magma", vmin=0, vmax=1)
        ax.set_xticks([]); ax.set_yticks([]); ax.set_title(b, fontsize=9)
        for u in range(2):
            for v in range(2):
                ax.text(v, u, f"{Rm[u,v]:.2f}", ha="center", va="center",
                        color="w" if Rm[u, v] < 0.6 else "k", fontsize=8)
    a1.axis("off")
    a1.set_title("accuracy matrix R[after, eval]", fontsize=9)

    fig.suptitle("E9 ↔ 3c — a LEARNED (stimulus × context) conjunction basis resolves "
                 "continual-learning interference (n=10)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_e9_bridge.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
