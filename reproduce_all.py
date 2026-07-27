#!/usr/bin/env python3
"""
reproduce_all.py — Comprehensive test and numerical integrity assertion harness for GHCA.

Checks:
1. Seeding / RNG audit (no global NumPy RNG usage).
2. Synthetic-SCM unit tests (test_ctestbed.py).
3. Fast numerical assertions verifying that key headline results in docs match
   the generated result archive (.npz) data.
"""

import sys
import os
import subprocess
import numpy as np

# Auto-switch to .venv/bin/python if sklearn is not available in current python
try:
    import sklearn
except ImportError:
    venv_python = os.path.join(os.path.dirname(__file__), ".venv", "bin", "python")
    if os.path.exists(venv_python) and sys.executable != venv_python:
        os.execv(venv_python, [venv_python] + sys.argv)

def run_step(name, func):
    print(f"===> [{name}] ... ", end="", flush=True)
    try:
        func()
        print("OK")
    except Exception as e:
        print(f"FAILED\n  Error: {e}")
        sys.exit(1)

def check_rng_audit():
    cmd = [sys.executable, ".claude/skills/experiment-review/review_helper.py", "audit-rng"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"RNG audit failed:\n{res.stdout}\n{res.stderr}")

def check_ctestbed_tests():
    cmd = [sys.executable, "test_ctestbed.py"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"test_ctestbed.py failed:\n{res.stdout}\n{res.stderr}")

def check_c_series_archives():
    # C0
    c0 = np.load("result/c0/c0_data.npz")
    assert "gain_coll" in c0 and "gain_lab" in c0

    # C1
    c1 = np.load("result/c1/c1_data.npz")
    assert bool(c1["all_ok"]), "C1 graph certificates all_ok check failed"

    # C2
    c2 = np.load("result/c2/c2_data.npz")
    assert "band_c" in c2 and "band_l" in c2

    # C3
    c3 = np.load("result/c3/c3_data.npz")
    assert "doTheta_band" in c3

    # C4
    c4 = np.load("result/c4/c4_data.npz")
    assert "eff" in c4 and "suff_coll" in c4

    # C5
    c5 = np.load("result/c5/c5_data.npz")
    assert "band_center" in c5 and "band_tracked" in c5

    # C6
    c6 = np.load("result/c6/c6_data.npz")
    assert "theta_chi_band" in c6 and "chi_band" in c6

    # C7
    c7 = np.load("result/c7/c7_data.npz")
    assert "matrix" in c7 and "screening" in c7

def check_e_series_archives():
    # E0
    e0 = np.load("result/e0/e0_data.npz")
    assert "p_ss" in e0

    # E1
    e1 = np.load("result/e1/e1_data.npz")
    assert "reward_A" in e1

    # E2
    e2 = np.load("result/e2/e2_data.npz")
    assert "mech_taus" in e2

    # E3
    e3 = np.load("result/e3/e3_data.npz")
    assert "mech_taus" in e3

    # E7
    e7 = np.load("result/e7/e7_learning.npz")
    assert "acc_intact" in e7

    # E9
    e9 = np.load("result/e9/e9_data.npz")
    assert "emergence_steps" in e9

def check_stats_sweeps_archives():
    # P3b C2 band stats
    c2_stats = np.load("result/stats/c2_band_stats.npz")
    coll_band = float(np.mean(c2_stats["seed_band_collective"]))
    lbl_band = float(np.mean(c2_stats["seed_band_labeled"]))
    assert np.isclose(coll_band, 0.287, atol=0.05), f"C2 coll_band {coll_band} != ~0.287"
    assert np.isclose(lbl_band, 26.82, atol=2.0), f"C2 lbl_band {lbl_band} != ~26.82"

    # P3b C4 suff stats
    c4_suff = np.load("result/stats/c4_suff_stats.npz")
    suff_coll = float(np.mean(c4_suff["seed_suff_collective"]))
    suff_lbl = float(np.mean(c4_suff["seed_suff_labeled"]))
    assert np.isclose(suff_coll, 1.027, atol=0.05), f"C4 suff_coll {suff_coll} != ~1.027"
    assert np.isclose(suff_lbl, 0.087, atol=0.05), f"C4 suff_lbl {suff_lbl} != ~0.087"

def check_closed_loop_plasticity_archives():
    p2 = np.load("result/closed_loop_plasticity/phase2_single_task.npz", allow_pickle=True)
    assert "data" in p2, "phase2_single_task.npz missing 'data' key"
    p2_data = p2["data"].item()
    assert "E1_Sensorimotor" in p2_data and "E5_Executive_Switch" in p2_data

    p3 = np.load("result/closed_loop_plasticity/phase3_sequential_learning.npz", allow_pickle=True)
    assert "data" in p3, "phase3_sequential_learning.npz missing 'data' key"
    p3_data = p3["data"].item()
    assert "Closed-Loop Multi-Axis (Tau, Theta, W)" in p3_data

def check_phase2_table_fresh():
    """Guard the one failure mode the archive checks cannot see: a results table
    hand-transcribed out of step with the .npz it claims to report."""
    cmd = [sys.executable, "scripts/print_phase2_table.py", "--check"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Phase-2 results table is stale:\n{res.stdout}\n{res.stderr}")


def check_unittests():
    cmd = [sys.executable, "-m", "unittest", "discover", "-s", "tests", "-p", "test_*.py"]
    res = subprocess.run(cmd, capture_output=True, text=True)
    if res.returncode != 0:
        raise RuntimeError(f"Unit tests failed:\n{res.stdout}\n{res.stderr}")

def main():
    print("==================================================")
    print("  GHCA Reproducibility & Integrity Assertion Suite  ")
    print("==================================================")
    run_step("1. Global RNG Usage Audit", check_rng_audit)
    run_step("2. Synthetic-SCM Unit Tests (test_ctestbed.py)", check_ctestbed_tests)
    run_step("3. Full Unit Test Suite", check_unittests)
    run_step("4. C-Series (.npz) Archive Integrity", check_c_series_archives)
    run_step("5. E-Series (.npz) Archive Integrity", check_e_series_archives)
    run_step("6. P3b Statistical Sweeps Integrity", check_stats_sweeps_archives)
    run_step("7. Closed-Loop Plasticity Archives Integrity", check_closed_loop_plasticity_archives)
    run_step("8. Phase-2 Results Table Matches Archive", check_phase2_table_fresh)
    print("==================================================")
    print("  ALL REPRODUCIBILITY CHECKS PASSED SUCCESSFULLY  ")
    print("==================================================")

if __name__ == "__main__":
    main()

