"""
Emit / verify Phase 3 table in docs/closed_loop_extensions_results.md from the archive.

Reads result/closed_loop_plasticity/structural_plasticity.npz and formats
Task 1 Initial Acc, Post-Test Acc, Raw Retention %, Chance-Corrected Retention %,
and Modularity Q into a markdown table.

Usage:
    python3 scripts/print_extensions_table.py         # print markdown table
    python3 scripts/print_extensions_table.py --check  # verify docs match archive
"""

import argparse
import os
import re
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
NPZ = os.path.join(ROOT, "result", "closed_loop_plasticity", "structural_plasticity.npz")
DOC = os.path.join(ROOT, "docs", "closed_loop_plasticity_extensions.md")


def generate_rows():
    if not os.path.exists(NPZ):
        raise FileNotFoundError(f"Archive file not found: {NPZ}")

    data = np.load(NPZ, allow_pickle=True)
    out = []
    
    # Standard conditions
    conditions = [
        "Multi-Axis Base (Tau, Theta, W)",
        "Multi-Axis + Axis G Rewiring (Tau, Theta, W, G)",
        "Multi-Axis + Axis G + Consolidation"
    ]

    for cond in conditions:
        if cond not in data:
            continue
        c_data = data[cond].item() if hasattr(data[cond], "item") else data[cond]
        
        # Pull or calculate metrics
        if "raw_init_accs" in c_data:
            init_arr = c_data["raw_init_accs"]
            test_arr = c_data["raw_test_accs"]
            ret_raw_arr = c_data["raw_retentions"]
            ret_cc_arr = c_data["raw_retentions_cc"]
            q_in_arr = c_data["raw_q_inits"]
            q_out_arr = c_data["raw_q_finals"]
            
            n = len(init_arr)
            init_m, init_s = np.mean(init_arr), np.std(init_arr)
            test_m, test_s = np.mean(test_arr), np.std(test_arr)
            ret_raw_m, ret_raw_s = np.mean(ret_raw_arr), np.std(ret_raw_arr) / np.sqrt(n)
            ret_cc_m, ret_cc_s = np.mean(ret_cc_arr), np.std(ret_cc_arr) / np.sqrt(n)
            q_in_m = np.mean(q_in_arr)
            q_out_m = np.mean(q_out_arr)
        else:
            # Fallback to summary tuples
            init_m, init_s = c_data["init_acc"]
            test_m, test_s = c_data["post_acc"]
            ret_raw_m, ret_raw_s = c_data["retention"][0], c_data["retention"][1] / np.sqrt(30)
            ret_cc_m, ret_cc_s = c_data.get("retention_chance_corr", (0.0, 0.0))[0], c_data.get("retention_chance_corr", (0.0, 0.0))[1] / np.sqrt(30)
            q_in_m, q_out_m = c_data["modularity"]

        c1 = f"{init_m:.3f} ± {init_s:.3f}"
        c2 = f"{test_m:.3f} ± {test_s:.3f}"
        c3 = f"{ret_raw_m:.1f}% ± {ret_raw_s:.1f}%"
        c4 = f"{ret_cc_m:.1f}% ± {ret_cc_s:.1f}%"
        c5 = f"{q_in_m:.4f} → {q_out_m:.4f}"

        row = f"| **{cond}** | {c1} | {c2} | {c3} | {c4} | {c5} |"
        out.append(row)

    return out


def print_table():
    headers = [
        "| Plasticity & Consolidation Regime | Task 1 Initial Acc | Task 1 Post-Test Acc | Raw Retention % | Chance-Corr Retention % | Modularity Q (Init → Final) |",
        "| :--- | :---: | :---: | :---: | :---: | :---: |"
    ]
    print("\n".join(headers))
    for r in generate_rows():
        print(r)


def check():
    if not os.path.exists(DOC):
        print(f"Doc file {DOC} does not exist.")
        return 1

    doc = open(DOC, encoding="utf-8").read()
    norm = lambda s: re.sub(r"\s+", " ", s).strip()
    doc_norm = norm(doc)

    missing = [r for r in generate_rows() if norm(r) not in doc_norm]
    if missing:
        print(f"Phase-3 table in {os.path.relpath(DOC, ROOT)} is STALE. Rows not found in doc:")
        for m in missing:
            print("  " + m)
        print("\nRegenerate and paste into doc, or update doc.")
        return 1
    print("Phase-3 table matches the archive.")
    return 0


if __name__ == "__main__":
    parser = argparse.ArgumentParser()
    parser.add_argument("--check", action="store_true", help="Verify doc matches archive")
    args = parser.parse_args()

    if args.check:
        sys.exit(check())
    else:
        print_table()
