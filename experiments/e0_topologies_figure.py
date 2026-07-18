"""Track 3b figure: self-sustaining band + period~tau across topologies."""

import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt  # noqa: E402

SD = os.path.join(ROOT, "result", "e0_topo")
FIG = os.path.join(ROOT, "docs", "figures")
TOPOS = ["lattice2d", "smallworld", "rgg"]
COL = {"lattice2d": "steelblue", "smallworld": "seagreen", "rgg": "crimson"}


def main():
    d = np.load(os.path.join(SD, "e0_topo.npz"))
    thetas, taus = d["thetas"], d["taus"]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    for t in TOPOS:
        g = d[f"band_{t}"]                       # (seeds, thetas)
        m = g.mean(0)
        a0.plot(thetas, m, "o-", color=COL[t], label=t)
        if g.shape[0] > 1:
            sem = g.std(0) / np.sqrt(g.shape[0])
            a0.fill_between(thetas, m - 1.96 * sem, m + 1.96 * sem, color=COL[t], alpha=0.15)
    a0.set_xlabel("threshold θ"); a0.set_ylabel("late active fraction A_ss (p_s=0)")
    a0.set_title("Self-sustaining band (matched mean degree ≈ 12)")
    a0.legend(fontsize=9)

    for t in TOPOS:
        g = d[f"period_{t}"]
        m = np.nanmean(g, 0)
        sl, ic, r = d[f"fit_{t}"]
        a1.plot(taus, m, "o", color=COL[t])
        a1.plot(taus, sl * taus + ic, "-", color=COL[t],
                label=f"{t}: {sl:.2f}·τ+{ic:.2f} (r={r:.3f})")
    a1.plot(taus, taus, "k--", alpha=0.4, label="period = τ")
    a1.set_xlabel("local timescale τ"); a1.set_ylabel("dominant global period")
    a1.set_title("period ~ τ in the loop-dominated regime (θ=1)")
    a1.legend(fontsize=8)

    fig.suptitle("3b — the E0 excitable story generalises off the lattice "
                 "(smallworld, rgg)", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "e0_topologies.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
