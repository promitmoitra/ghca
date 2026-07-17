"""3a / P2 — figures for the operating-point sweeps.

Reads result/stats/e3_latency_sweep.npz and e7_theta_sweep.npz (skips whichever
is absent) and writes docs/figures/stats_p2_*.png.
"""

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
os.makedirs(FIG, exist_ok=True)


def e3_latency():
    p = os.path.join(SD, "e3_latency_sweep.npz")
    if not os.path.exists(p):
        return
    d = np.load(p)
    lats, rate = d["lats"], d["joint_rate"]
    lo, hi = d["joint_lo"], d["joint_hi"]
    taus, reson = d["taus"], d["resonance"]
    fig, (a0, a1) = plt.subplots(1, 2, figsize=(13, 5))

    yerr = np.clip(np.vstack([rate - lo, hi - rate]), 0, None)
    a0.errorbar(lats, rate, yerr=yerr, fmt="o-", color="crimson", capsize=5, ms=7)
    a0.axhline(0.25, ls=":", color="0.6", label="≈chance floor")
    a0.set_xlabel("target latency (curriculum)")
    a0.set_ylabel("joint-success rate (identity ∧ timing)")
    a0.set_title("E3 composition vs target latency (n=%d, Wilson 95%%)" % int(d["n"]))
    a0.set_ylim(-0.02, 0.6)
    for L, r in zip(lats, rate):
        a0.annotate(f"τ≈{int(L)+2}", (L, r), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=8, color="0.4")
    a0.legend(fontsize=8)

    a1.plot(taus, reson, "s-", color="steelblue")
    a1.axhline(0.5, ls=":", color="0.6")
    for L in lats:
        a1.axvline(int(L) + 2, color="crimson", alpha=0.25, lw=1)
    a1.set_xlabel("fixed gate τ")
    a1.set_ylabel("identity accuracy")
    a1.set_title("resonance map (identity vs fixed τ); red = swept latencies' τ")
    a1.set_ylim(0, 1.05)

    fig.suptitle("3a/P2 — E3 composition is operating-point-contingent: "
                 "joint-success tracks the substrate resonance map", fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "stats_p2_e3_latency.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


def e7_theta():
    p = os.path.join(SD, "e7_theta_sweep.npz")
    if not os.path.exists(p):
        return
    d = np.load(p)
    th = d["thetas"]
    im, am = d["intact_mean"], d["ablate_mean"]
    lo, hi = d["intact_lo"], d["intact_hi"]
    fig, ax = plt.subplots(figsize=(7, 5))
    ax.errorbar(th, im, yerr=np.vstack([im - lo, hi - im]), fmt="o-", color="seagreen",
                capsize=5, ms=7, label="intact spiral")
    ax.plot(th, am, "s--", color="crimson", label="ablated (θ_dead=5.0)")
    ax.axhline(0.5, ls=":", color="0.6", label="chance")
    for x, y in zip(th, im):
        ax.annotate(f"≥{int(np.ceil(x))} nbrs", (x, y), textcoords="offset points",
                    xytext=(0, 10), ha="center", fontsize=8, color="0.4")
    ax.set_xlabel("spiral threshold θ_live  (excites on inp ≥ θ; inp is integer)")
    ax.set_ylabel("switching accuracy (last 4 blocks)")
    ax.set_title("E7 switching across the (integer) spiral band (n=%d)" % int(d["n"]))
    ax.set_ylim(0, 1.02)
    ax.legend(fontsize=9)
    fig.tight_layout()
    out = os.path.join(FIG, "stats_p2_e7_theta.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    e3_latency()
    e7_theta()
