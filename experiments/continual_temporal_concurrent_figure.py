"""3e.3 figure: concurrent co-adaptation vs phase-split vs homogeneous."""

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
MODES = [("homogeneous", "crimson", "homogeneous τ (readout only)"),
         ("phase-split", "slategray", "phase-split (grow τ → freeze → learn)"),
         ("concurrent", "seagreen", "concurrent (grow τ + learn together)")]


def main():
    d = np.load(os.path.join(SD, "continual_temporal_concurrent.npz"))
    T = d["T_list"]
    fig, ax = plt.subplots(figsize=(7.5, 5.2))
    for mode, col, lab in MODES:
        arr = d[f"{mode}_avg"]
        means, los, his = [], [], []
        for row in arr:
            m, lo, hi = st.bootstrap_ci(row)
            means.append(m); los.append(m - lo); his.append(hi - m)
        ax.errorbar(T, means, yerr=[los, his], marker="o", capsize=4, color=col, label=lab)
    ax.axhline(0.5, ls=":", color="0.5", label="chance 0.50")
    ax.set_xlabel("number of sequential temporal tasks T")
    ax.set_ylabel("final average accuracy (per-task heads)")
    ax.set_xticks(T)
    ax.set_title("3e.3 — growing the timescale basis *concurrently* with the readout:\n"
                 "well above floor, but below phase-split (a real co-adaptation cost, n=%d)"
                 % int(d["n"]))
    ax.legend(fontsize=8)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_temporal_concurrent.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
