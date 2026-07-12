"""
2b bridge figure — recast the C5 model result as a falsifiable *data* signature.

Not a new experiment: reloads `result/c5/c5_data.npz` (the do(chirality)
fat-handedness run) and re-plots the core-displacement series as an explicit
prediction for spiral-wave recordings — a fixed-locus (single-ROI/electrode)
decoder of the spiral's rule collapses as the core drifts away from that locus,
while a core-tracked (topology-aware) decoder does not. See
`docs/spiral_predictions.md`, prediction P4.

Run: python3 experiments/p2b_signature_figure.py
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

d = np.load(os.path.join(ROOT, "result", "c5", "c5_data.npz"), allow_pickle=True)
labels = [str(x) for x in d["labels"]]
# displacement series: centered(0), disp-4, disp-8, disp-12
disp_idx = [labels.index(x) for x in ["centered", "disp-4", "disp-8", "disp-12"]]
disp = np.array([0, 4, 8, 12])
center = d["acc_center"][disp_idx]
tracked = d["acc_tracked"][disp_idx]
glob = d["acc_global"][disp_idx]
routing = {"center": d["routing_center"], "tracked": d["routing_tracked"]}

fig, axes = plt.subplots(1, 2, figsize=(12, 4.6))

ax = axes[0]
ax.plot(disp, center, "o-", color="crimson", label="fixed-locus reader (single ROI/electrode)")
ax.plot(disp, tracked, "s-", color="seagreen", label="core-tracked reader (topology-aware)")
ax.plot(disp, glob, "^--", color="slategray", alpha=0.8, label="global winding (whole-field)")
ax.axhline(0.5, ls=":", color="k", alpha=0.5, label="chance")
ax.set_xlabel("spiral-core displacement from the reader (lattice units)")
ax.set_ylabel("rule-decode accuracy")
ax.set_ylim(-0.03, 1.05)
ax.set_title("Predicted signature (P4): a fixed-ROI decoder of spiral\n"
             "rotation collapses as the core drifts; topology-aware does not")
ax.legend(fontsize=8, loc="center left")

ax = axes[1]
x = np.arange(2)
ax.bar(x - 0.18, [routing["center"][0], routing["tracked"][0]], 0.36,
       color=["crimson", "seagreen"], label="displaced core")
ax.set_xticks(x)
ax.set_xticklabels(["fixed-locus\nreader", "core-tracked\nreader"])
ax.axhline(0.5, ls=":", color="k", alpha=0.5, label="chance")
ax.set_ylabel("behavioural routing accuracy\n(rule applied) at a displaced core")
ax.set_ylim(0, 1.0)
ax.set_title("Behavioural consequence: fixed-locus reading of a\n"
             "displaced spiral drives the wrong rule (0.55 vs 0.78)")
ax.legend(fontsize=8)

fig.tight_layout()
out = os.path.join(FIGDIR, "spiral_predictions_signature.png")
fig.savefig(out, dpi=110)
print("wrote", out)
print("displacement series (disp, center, tracked, global):")
for i, dd in enumerate(disp):
    print(f"  d={dd:2d}: fixed={center[i]:.2f}  tracked={tracked[i]:.2f}  global={glob[i]:.2f}")
