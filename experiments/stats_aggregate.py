"""3a / P1 — aggregate result/stats/*_stats.npz into a master table + figure.

Reads the per-seed arrays saved by stats_seed_scaleup.py, emits a markdown table
(headline · n · mean · 95% CI · effect size · distribution shape) and a
multi-panel strip+CI figure. Distribution shape is read from the histogram, not
just Sarle's BC (which also flags ceiling-with-tail), so the label is honest.
"""

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
os.makedirs(FIG, exist_ok=True)

# (key, [(label, seed-array-name)], panel title)
PANELS = [
    ("e1",  [("A", "seed_A"), ("B", "seed_B")], "E1 conditioning (final acc)"),
    ("e2",  [("A", "seed_A"), ("B", "seed_B")], "E2 retention @ D=200"),
    ("e3t", [("A", "seed_A_identity"), ("B", "seed_B_identity")], "E3 timed: identity"),
    ("e3f", [("curric.", "seed_curriculum_identity"),
             ("factored", "seed_factored_identity")], "E3 factored: identity"),
    ("e5",  [("intact", "seed_switch_intact"), ("ablate", "seed_switch_ablate")], "E5 switching"),
    ("e7",  [("intact", "seed_switch_intact"), ("ablate", "seed_switch_ablate")], "E7 switching"),
    ("e9",  [("emergent", "seed_routing_emergent"),
             ("frozen", "seed_routing_frozen")], "E9 routing"),
]


def shape(arr):
    """Honest distribution label from the histogram, not BC alone."""
    a = np.asarray(arr, float); a = a[np.isfinite(a)]
    if a.size < 4 or a.std() == 0:
        return "point"
    h, _ = np.histogram(a, bins=10, range=(0, 1))
    top = h[7:].sum() / a.size            # mass near ceiling (0.7–1.0)
    bot = h[:3].sum() / a.size            # mass near floor (0.0–0.3)
    bc = st.bimodality(a)["bc"]
    if top >= 0.6 and bot <= 0.12:
        return "ceiling+tail"
    if top >= 0.25 and bot >= 0.20:
        return "bimodal"
    if bc > st.BC_UNIFORM:
        return "spread"
    return "unimodal"


def main():
    rows = []
    fig, axes = plt.subplots(2, 4, figsize=(16, 8))
    axes = axes.ravel()
    for i, (key, labs, title) in enumerate(PANELS):
        p = os.path.join(SD, f"{key}_stats.npz")
        if not os.path.exists(p):
            axes[i].set_visible(False); continue
        d = np.load(p)
        samples = [np.asarray(d[nm], float) for _, nm in labs]
        st.strip_ci(axes[i], samples, [l for l, _ in labs])
        axes[i].set_title(title); axes[i].set_ylim(-0.05, 1.05)
        axes[i].axhline(0.5, ls=":", color="0.7", lw=0.8)
        for (lab, nm) in labs:
            r = st.summarise(f"{key}:{lab}", np.asarray(d[nm], float))
            r["shape"] = shape(d[nm])
            rows.append(r)
    for j in range(len(PANELS), len(axes)):
        axes[j].set_visible(False)
    fig.suptitle("3a/P1 — headline dissociations at n=50 (points = seeds; red = mean ± 95% CI)",
                 fontsize=13)
    fig.tight_layout()
    out = os.path.join(FIG, "stats_sweeps.png")
    fig.savefig(out, dpi=110); print("wrote", out)

    # markdown table
    print("\n| headline | n | mean | 95% CI | shape |")
    print("|---|:--:|:--:|:--:|---|")
    for r in rows:
        print(f"| {r['name']} | {r['n']} | {r['mean']:.3f} "
              f"| [{r['ci_lo']:.3f}, {r['ci_hi']:.3f}] | {r['shape']} |")


if __name__ == "__main__":
    main()
