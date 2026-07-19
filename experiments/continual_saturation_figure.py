"""P5 saturation figure: a fixed conjunction basis runs out of room as T grows."""

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
MODES = [("shared", "crimson", "shared (no conjunction)"),
         ("context", "seagreen", "context (fixed conj. basis)"),
         ("per-task", "slategray", "per-task heads (capacity ∝ T)")]


def main():
    d = np.load(os.path.join(SD, "continual_saturation.npz"))
    T = d["T_list"]; chance = float(d["chance"]); n_h = int(d["n_h"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    for metric, ax, ttl, ylab in [("bwt", a0, "Backward transfer vs #tasks", "backward transfer"),
                                   ("avg", a1, "Final average accuracy vs #tasks", "avg accuracy")]:
        for mode, col, lab in MODES:
            arr = d[f"{mode}_{metric}"]        # (len(T), N)
            means, los, his = [], [], []
            for row in arr:
                m, lo, hi = st.bootstrap_ci(row)
                means.append(m); los.append(m - lo); his.append(hi - m)
            ax.errorbar(T, means, yerr=[los, his], marker="o", capsize=4,
                        color=col, label=lab)
        ax.set_xlabel("number of sequential tasks T")
        ax.set_ylabel(ylab)
        ax.set_title(ttl)
        ax.set_xticks(T)
        if metric == "bwt":
            ax.axhline(0, color="k", lw=0.8, ls="--")
        else:
            ax.axhline(chance, color="0.5", ls=":", label=f"chance {chance:.2f}")
        ax.legend(fontsize=8)

    fig.suptitle(f"P5 — a fixed (stimulus × context) conjunction basis (n_h={n_h}) "
                 f"saturates as tasks accumulate", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_saturation.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
