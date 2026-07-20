"""3d emergent-arm figure: a *grown* τ basis matches the hand-set one."""

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
BASES = [("homog", "crimson", "homogeneous τ"),
         ("emergent", "seagreen", "emergent τ (grown)"),
         ("graded", "slategray", "wired τ (hand-set)")]


def main():
    d = np.load(os.path.join(SD, "continual_temporal_emergent.npz"))
    T = d["T_list"]; chance = float(d["chance"])
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(16, 4.8))

    # panel 0: τ distributions (the emergence)
    bins = np.linspace(2, 40, 24)
    for b, col, lab in BASES:
        a0.hist(d[f"tau_{b}"], bins=bins, alpha=0.55, color=col, label=lab)
    a0.set_xlabel("hidden-unit refractory timescale τ")
    a0.set_ylabel("count (seed 0)")
    a0.set_title("The basis: homogeneous vs grown vs hand-set τ")
    a0.legend(fontsize=8)

    # panels 1,2: CL avg accuracy vs T, per-task then shared
    for ax, head, ttl in [(a1, "per-task", "per-task heads (does the basis REPRESENT it?)"),
                          (a2, "shared", "shared head (interference)")]:
        for b, col, lab in BASES:
            arr = d[f"{b}+{head}_avg"]
            means, los, his = [], [], []
            for row in arr:
                m, lo, hi = st.bootstrap_ci(row)
                means.append(m); los.append(m - lo); his.append(hi - m)
            ax.errorbar(T, means, yerr=[los, his], marker="o", capsize=4,
                        color=col, label=lab)
        ax.axhline(chance, ls=":", color="0.5", label=f"chance {chance:.2f}")
        ax.set_xlabel("number of sequential temporal tasks T")
        ax.set_ylabel("final average accuracy")
        ax.set_xticks(T); ax.set_title(ttl); ax.legend(fontsize=8)

    fig.suptitle("3d (emergent) — a reward-free, input-timing-driven rule GROWS a τ basis "
                 "that matches the hand-set one (n=%d)" % int(d["n"]), fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_temporal_emergent.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
