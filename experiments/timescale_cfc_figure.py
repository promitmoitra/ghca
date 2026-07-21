"""3e.2b figure: cross-frequency coupling on the emergent hierarchy (closes 4a)."""

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
    d = np.load(os.path.join(SD, "timescale_cfc.npz"))
    nb = int(d["n_phase"])
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(12, 4.8))
    phase = (np.arange(nb) + 0.5) / nb

    # panel 0: phase–amplitude profile (fast amplitude vs slow phase), normalised
    for mode, col in [("uncoupled", "slategray"), ("coupled", "seagreen")]:
        prof = d[f"{mode}_prof"]
        norm = prof / (prof.sum(1, keepdims=True) + 1e-12)   # per-seed distribution
        m = norm.mean(0)
        sd = norm.std(0) / np.sqrt(norm.shape[0])
        a0.plot(phase, m, "o-", color=col, label=mode)
        a0.fill_between(phase, m - sd, m + sd, color=col, alpha=0.2)
    a0.axhline(1.0 / nb, ls=":", color="0.5", label="uniform (no coupling)")
    a0.set_xlabel(f"slow-population phase (period P_s={int(d['P_S'])})")
    a0.set_ylabel("fast-population amplitude (normalised)")
    a0.set_title("Phase–amplitude profile:\nfast amplitude rides on slow phase when coupled")
    a0.legend(fontsize=8)

    # panel 1: PAC modulation index, uncoupled vs coupled
    modes = ["uncoupled", "coupled"]
    means, los, his = [], [], []
    for mode in modes:
        m, lo, hi = st.bootstrap_ci(d[f"{mode}_mi"])
        means.append(m); los.append(m - lo); his.append(hi - m)
    a1.bar([0, 1], means, yerr=[los, his], capsize=5,
           color=["slategray", "seagreen"], width=0.6)
    a1.set_xticks([0, 1]); a1.set_xticklabels(modes)
    a1.set_ylabel("PAC modulation index (Tort)")
    a1.set_title("Cross-frequency coupling (n=%d)\nslow phase → fast amplitude" % int(d["n"]))

    fig.suptitle("3e.2b — theta–gamma-style cross-frequency coupling on the emergent "
                 "timescale hierarchy (closes 4a)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "timescale_cfc.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
