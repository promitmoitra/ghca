"""3d figure: timescale diversity is the continual-learning capacity axis for
temporal tasks (a homogeneous basis has no temporal code at any head capacity)."""

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
ARMS = [("graded+per-task", "seagreen", "-", "graded τ + per-task heads (temporal code, capacity ∝ T)"),
        ("graded+shared", "seagreen", "--", "graded τ + shared head (temporal code, one head)"),
        ("homog+per-task", "crimson", "-", "homogeneous τ + per-task heads (no temporal code)"),
        ("homog+shared", "crimson", "--", "homogeneous τ + shared head")]


def main():
    d = np.load(os.path.join(SD, "continual_temporal_saturation.npz"))
    T = d["T_list"]; chance = float(d["chance"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    for arm, col, ls, lab in ARMS:
        avg = d[f"{arm}_avg"]              # (len(T), N)
        means, los, his = [], [], []
        for row in avg:
            m, lo, hi = st.bootstrap_ci(row)
            means.append(m); los.append(m - lo); his.append(hi - m)
        a0.errorbar(T, means, yerr=[los, his], marker="o", ls=ls, capsize=4,
                    color=col, label=lab)
    a0.axhline(chance, color="0.5", ls=":", label=f"chance {chance:.2f}")
    a0.set_xlabel("number of sequential temporal tasks T")
    a0.set_ylabel("final average accuracy")
    a0.set_xticks(T)
    a0.set_title("Timescale diversity is the lever (only graded τ clears chance)")
    a0.legend(fontsize=7.5, loc="upper right")

    for arm, col, ls, lab in ARMS:
        bwt = d[f"{arm}_bwt"]
        means, los, his = [], [], []
        for row in bwt:
            m, lo, hi = st.bootstrap_ci(row)
            means.append(m); los.append(m - lo); his.append(hi - m)
        a1.errorbar(T, means, yerr=[los, his], marker="o", ls=ls, capsize=4,
                    color=col, label=arm)
    a1.axhline(0, color="k", lw=0.8, ls="--")
    a1.set_xlabel("number of sequential temporal tasks T")
    a1.set_ylabel("backward transfer")
    a1.set_xticks(T)
    a1.set_title("Only the shared head over the fixed graded basis forgets")
    a1.legend(fontsize=7.5)

    fig.suptitle("3d (wired) — a hand-set timescale-diverse basis buys continual-learning "
                 "capacity on temporal tasks (n=20)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_temporal_saturation.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
