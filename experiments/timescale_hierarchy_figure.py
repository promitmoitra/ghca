"""3e.2 figure: an emergent fast/slow timescale hierarchy (closes 4a)."""

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
    d = np.load(os.path.join(SD, "timescale_hierarchy.npz"))
    P_F, P_S = int(d["P_F"]), int(d["P_S"])
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    bins = np.linspace(3, 34, 26)

    # panel 0: the headline histograms — new rule (emergent) bimodal vs old rule ratchet
    a = axes[0]
    a.hist(d["input+emergent+two_tau0"], bins=bins, alpha=0.6, color="seagreen",
           label="input rule (new) — emergent")
    a.hist(d["own+emergent+two_tau0"], bins=bins, alpha=0.5, color="crimson",
           label="own-fire rule (old, e10)")
    for p, lab in [(P_F, f"P_f={P_F}"), (P_S, f"P_s={P_S}")]:
        a.axvline(p, ls="--", color="0.4"); a.text(p + 0.3, a.get_ylim()[1] * 0.9, lab, fontsize=8)
    a.set_xlabel("learned τ"); a.set_ylabel("count (seed 0)")
    a.set_title("Two-rhythm drive: new rule splits at P_f, P_s;\nold rule ratchets up")
    a.legend(fontsize=8)

    # panel 1: single-rhythm controls stay unimodal
    a = axes[1]
    a.hist(d["input+emergent+fast_tau0"], bins=bins, alpha=0.6, color="royalblue",
           label="fast-only drive")
    a.hist(d["input+emergent+slow_tau0"], bins=bins, alpha=0.6, color="darkorange",
           label="slow-only drive")
    for p in (P_F, P_S):
        a.axvline(p, ls="--", color="0.4")
    a.set_xlabel("learned τ"); a.set_title("Single-rhythm controls (new rule):\nunimodal at the driven period")
    a.legend(fontsize=8)

    # panel 2: the honest discriminator — fraction of τ AT each true period.
    # (BC alone misleads: the old ratcheting rule also clears 5/9, but with NO fast
    # cluster — near_f = 0. Placement at the fundamentals is what matters.)
    a = axes[2]
    keys = ["input+emergent+two", "input+wired+two", "own+emergent+two",
            "own+wired+two", "input+emergent+fast", "input+emergent+slow"]
    labs = ["new/emerg\n(two)", "new/wired\n(two)", "old/emerg\n(two)",
            "old/wired\n(two)", "new\n(fast)", "new\n(slow)"]
    nf = [float(d[f"{k}_near_f"].mean()) for k in keys]
    ns = [float(d[f"{k}_near_s"].mean()) for k in keys]
    x = np.arange(len(keys)); w = 0.38
    a.bar(x - w / 2, nf, w, color="royalblue", label=f"near P_f={P_F}")
    a.bar(x + w / 2, ns, w, color="darkorange", label=f"near P_s={P_S}")
    a.set_xticks(x); a.set_xticklabels(labs, fontsize=7)
    a.set_ylabel("fraction of τ at the period"); a.set_ylim(0, 1)
    a.set_title("Clusters AT the drive periods (n=%d)\n(old rule: no fast cluster)" % int(d["n"]))
    a.legend(fontsize=8)

    fig.suptitle("3e.2 — an input-timing-driven τ rule GROWS a fast/slow hierarchy where the "
                 "old rule cannot (closes 4a)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "timescale_hierarchy.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
