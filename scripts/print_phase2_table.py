"""Emit Table 2 of docs/closed_loop_plasticity_results.md from the archive.

The Phase-2 table was hand-transcribed and drifted: the Trained Readout Acc
column disagreed with `phase2_single_task.npz` in every row (E1 read 88.2% for
a value that is 99.4%; E5 read 31.5% for 49.4%), and the E5 fixed-accuracy
standard deviations were ~4x too large. Nothing downstream caught it because
the numbers are prose, not asserted values.

This script is the single source of the table. Regenerate and paste after any
rerun of `experiments/closed_loop_plasticity.py`:

    python3 scripts/print_phase2_table.py

`--check` instead verifies the committed doc still matches the archive, and
exits non-zero if not (suitable for reproduce_all.py).

Note on $RIR$: the archive stores `mean_rir` as the per-seed mean of
`fixed_acc / trained_acc`, which by Jensen's inequality is NOT the ratio of the
column means. On E1 the gap is large because `trained_acc` sits at ceiling
while `fixed_acc` varies from 0.007 to 0.98 across seeds. The table reports the
stored per-seed statistic; the columns are context, not its derivation.
"""

import argparse
import os
import re
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPZ = os.path.join(ROOT, "result", "closed_loop_plasticity", "phase2_single_task.npz")
DOC = os.path.join(ROOT, "docs", "closed_loop_plasticity_results.md")

TASK_LABEL = {"E1_Sensorimotor": "**E1 Sensorimotor Routing**",
              "E5_Executive_Switch": "**E5 Executive Switch**"}
COND_LABEL = {"Control (No Plasticity)": "Control (No Plasticity)",
              "Axis W Only": "Axis $W$ Only",
              "Axis Tau Only": r"Axis $\tau$ Only",
              "Axis Theta Only": r"Axis $\theta$ Only",
              "Closed-Loop Multi-Axis": "**Closed-Loop Multi-Axis**"}


def rows():
    data = np.load(NPZ, allow_pickle=True)["data"].item()
    out = []
    for task, conds in data.items():
        first = True
        for cond, v in conds.items():
            # ddof=0, matching the archive's own stored mean_rir / std_rir
            fa, fs = v["mean_fixed_acc"] * 100, v["std_fixed_acc"] * 100
            ta = float(np.mean(v["trained_acc"])) * 100
            ts = float(np.std(v["trained_acc"])) * 100
            r, rs = v["mean_rir"], v["std_rir"]
            bold = cond == "Closed-Loop Multi-Axis"
            if bold:
                c1 = "**%.1f%% ± %.1f%%**" % (fa, fs)
                c2 = "**%.1f%% ± %.1f%%**" % (ta, ts)
                c3 = "**%.3f ± %.3f**" % (r, rs)
            else:
                c1 = "$%.1f\\%% \\pm %.1f\\%%$" % (fa, fs)
                c2 = "$%.1f\\%% \\pm %.1f\\%%$" % (ta, ts)
                c3 = "$%.3f \\pm %.3f$" % (r, rs)
            out.append("| %s | %s | %s | %s | %s |"
                       % (TASK_LABEL[task] if first else "", COND_LABEL[cond], c1, c2, c3))
            first = False
    return out


def max_jensen_gap():
    """Largest |ratio-of-column-means - mean-of-per-seed-ratios| over all rows."""
    data = np.load(NPZ, allow_pickle=True)["data"].item()
    return max(abs(float(np.mean(v["fixed_acc"])) / float(np.mean(v["trained_acc"]))
                   - float(v["mean_rir"]))
               for conds in data.values() for v in conds.values())


def check():
    """Verify every generated row appears verbatim in the committed doc, and
    that the footnote's stated agreement bound is not tighter than the data."""
    doc = open(DOC, encoding="utf-8").read()
    m = re.search(r"agree to\s+within (\d*\.\d+) in every row", doc)
    if m:
        stated, actual = float(m.group(1)), max_jensen_gap()
        if actual > stated:
            print("Footnote claims rows agree to within %.3f, but the largest gap "
                  "in the archive is %.4f." % (stated, actual))
            return 1
    norm = lambda s: re.sub(r"\s+", " ", s).strip()
    doc_norm = norm(doc)
    missing = [r for r in rows() if norm(r) not in doc_norm]
    if missing:
        print("Phase-2 table in %s is STALE. Rows not found in the doc:"
              % os.path.relpath(DOC, ROOT))
        for m in missing:
            print("  " + m)
        print("\nRegenerate with: python3 scripts/print_phase2_table.py")
        return 1
    print("Phase-2 table matches the archive (%d rows)." % len(rows()))
    return 0


if __name__ == "__main__":
    ap = argparse.ArgumentParser()
    ap.add_argument("--check", action="store_true",
                    help="verify the doc matches the archive; exit 1 if stale")
    args = ap.parse_args()
    if args.check:
        sys.exit(check())
    print("| Benchmark | Condition | Fixed Readout Acc (%) | Trained Readout Acc (%)"
          " | Readout Independence Ratio ($RIR$) |")
    print("| :--- | :--- | :---: | :---: | :---: |")
    print("\n".join(rows()))
