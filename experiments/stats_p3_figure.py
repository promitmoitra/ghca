"""3a / P3 — E8 seed/CI figure (four headline panels at n=50)."""

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


def _pt(name):
    d = np.load(os.path.join(SD, f"{name}_stats.npz"))
    return {k[5:]: np.asarray(d[k], float) for k in d.files if k.startswith("seed_")}


def main():
    fig, ax = plt.subplots(1, 4, figsize=(17, 4.3))

    pred = _pt("e8_prediction")
    order = ["periodic", "markov0.9", "markov0.5", "markov0.1", "randomwalk"]
    st.strip_ci(ax[0], [pred[k] for k in order], ["per", "M.9", "M.5", "M.1", "RW"])
    ax[0].axhline(0.125, ls=":", color="0.6"); ax[0].set_ylim(0, 1.05)
    ax[0].set_title("E8.1/8.2 prediction acc"); ax[0].set_ylabel("accuracy")

    win = _pt("e8_window")
    taus = [4, 8, 14, 20, 26]
    means = [win[f"tau{t}"].mean() for t in taus]
    ax[1].plot(taus, means, "o-", color="seagreen")
    ax[1].set_title("E8.3 integration window vs τ")
    ax[1].set_xlabel("τ"); ax[1].set_ylabel("decodable depth (tones)")

    nes = _pt("e8_nested")
    st.strip_ci(ax[2], [nes["nested_intact"], nes["nested_ablate"], nes["within_intact"]],
                ["nest\nintact", "nest\nablate", "within\nintact"])
    ax[2].axhline(0.125, ls=":", color="0.6"); ax[2].set_ylim(0, 1.05)
    ax[2].set_title("E8.5 nested (τ_ctx)")

    con = _pt("e8_conditional")
    st.strip_ci(ax[3], [con["grid_conj"], con["grid_noconj"], con["trace_conj"]],
                ["grid\n+conj", "grid\n-conj", "trace\n+conj"])
    ax[3].axhline(0.125, ls=":", color="0.6"); ax[3].set_ylim(0, 1.05)
    ax[3].set_title("E8.7 conditional (K=4)")

    fig.suptitle("3a/P3 — E8 predictive series at n=50 (points = seeds; red = mean ± 95% CI)",
                 fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "stats_p3_e8.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
