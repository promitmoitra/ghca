"""3e.1 figure: the emergent τ basis re-tiles under a shifting delay distribution,
with a *graceful* (not catastrophic) representation-level interference."""

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
    d = np.load(os.path.join(SD, "continual_temporal_retile.npz"), allow_pickle=True)
    stages = [str(s) for s in d["stages"]]
    regimes = [str(r) for r in d["regimes"]]
    short, long_ = d["short"], d["long"]
    xs = np.arange(len(stages))
    fig, (a0, a1, a2) = plt.subplots(1, 3, figsize=(16, 4.8))

    # panel 0: τ migration (emergent seed-0 snapshots after each stage)
    tau = d["emergent_tau"]                       # (n_stages, n_h)
    bins = np.linspace(2, 40, 22)
    cols = ["#1b9e77", "#d95f02", "#7570b3"]
    for si, stage in enumerate(stages):
        a0.hist(tau[si], bins=bins, alpha=0.5, color=cols[si % 3],
                label=f"after {stage} (stage {si})")
    a0.axvspan(short[0], short[-1], color="steelblue", alpha=0.12)
    a0.axvspan(long_[0], long_[-1], color="crimson", alpha=0.12)
    a0.set_xlabel("hidden-unit τ  (blue=SHORT range, red=LONG range)")
    a0.set_ylabel("count (seed 0)")
    a0.set_title("τ migrates to tile the current regime")
    a0.legend(fontsize=7)

    # panels 1,2: decode vs stage for each eval regime
    rc = {"SHORT": "steelblue", "LONG": "crimson"}
    for ax, cond, ttl in [(a1, "emergent", "emergent τ — adapts (and gently forgets)"),
                          (a2, "wired-static", "wired-static τ — cannot adapt")]:
        M = d[f"{cond}_M"]                         # (N, n_stages, n_regimes)
        for ri, rn in enumerate(regimes):
            means, los, his = [], [], []
            for si in range(len(stages)):
                m, lo, hi = st.bootstrap_ci(M[:, si, ri])
                means.append(m); los.append(m - lo); his.append(hi - m)
            ax.errorbar(xs, means, yerr=[los, his], marker="o", capsize=4,
                        color=rc.get(rn, "k"), label=f"decode on {rn}")
        ax.axhline(1.0 / len(short), ls=":", color="0.5", label="chance")
        ax.set_xticks(xs); ax.set_xticklabels([f"after\n{s}" for s in stages])
        ax.set_ylim(0.15, 1.02)
        ax.set_xlabel("self-organisation stage (delay regime)")
        ax.set_ylabel("delay-decode accuracy")
        ax.set_title(ttl); ax.legend(fontsize=8)

    fig.suptitle("3e.1 — a plastic (emergent) timescale basis RE-TILES under a shifting "
                 "delay distribution, with graceful representation-level interference "
                 "(n=%d)" % int(d["n"]), fontsize=12)
    fig.tight_layout()
    out = os.path.join(FIG, "continual_temporal_retile.png")
    fig.savefig(out, dpi=110)
    print("wrote", out)


if __name__ == "__main__":
    main()
