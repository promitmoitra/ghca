"""
Phase 3 Sequential Learning & Anti-Forgetting Protocol.

Evaluates Task A (Identity Routing) -> Task B (Reversal Routing) -> Task A (Retention Test)
under closed-loop substrate plasticity across n=30 independent seeds.

Quantifies catastrophic forgetting mitigation without task heads or weight freezing:
    Task A Retention = Accuracy(Task A after Task B) / Accuracy(Task A initial) * 100%

Target Acceptance Criteria:
    - Task A Retention > 80% under Closed-Loop Multi-Axis Plasticity.
    - Demonstrates topological protection via intrinsic timescale tau distribution shifts.
    - Standard weight-only plasticity (Axis W) exhibits catastrophic forgetting (< 40% retention).
    - Zero global NumPy RNG usage (explicit default_rng seed threading).
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import ClosedLoopPlasticityEngine, make_closed_loop_learner

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "closed_loop_plasticity")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_SEEDS = 30
TRIALS_TASK_A = 350
TRIALS_TASK_B = 250
TRIALS_TEST_A = 100
SETTLE, CUE, WWIN = 10, 3, 6


def run_sequential_trial(net, x, target_action, rng, train=True):
    """Run a single trial in the sequential learning protocol."""
    net.reset_traces()
    for _ in range(SETTLE):
        net.step_learn(None)

    feats_before = net.features()
    V = net.value()

    for _ in range(CUE):
        net.step_learn(net.sensory_drive(x))

    sc = np.zeros(net.roles["A"])
    for _ in range(WWIN):
        net.step_learn(None)
        sc += net.motor_scores()

    if sc.sum() > 0:
        action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc))))
    else:
        action = -1

    reward = 1.0 if action == target_action else 0.0
    delta = reward - V

    if train:
        net.learn(delta)
        net.update_critic(delta, feats_before)

    return reward


def evaluate_sequential_seed(seed, axes=("tau", "theta", "w")):
    """Run Task A -> Task B -> Task A retention test for a single seed."""
    net = make_closed_loop_learner(axes=axes, seed=seed, cfg=dict(
        act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.08, eta_tau=0.02, eta_theta=0.01
    ))
    rng = np.random.default_rng(seed + 1000)

    # Phase 1: Train Task A (Identity Mapping: 0->0, 1->1)
    rewards_A1 = np.zeros(TRIALS_TASK_A)
    for t in range(TRIALS_TASK_A):
        x = int(rng.integers(net.roles["K"]))
        target = x  # Identity rule
        r = run_sequential_trial(net, x, target, rng, train=True)
        rewards_A1[t] = r

    acc_A1 = float(rewards_A1[-50:].mean())
    tau_after_A1 = net.tau.copy()

    # Phase 2: Train Task B (Reversal Mapping: 0->1, 1->0)
    rewards_B = np.zeros(TRIALS_TASK_B)
    for t in range(TRIALS_TASK_B):
        x = int(rng.integers(net.roles["K"]))
        target = 1 - x  # Reversal rule
        r = run_sequential_trial(net, x, target, rng, train=True)
        rewards_B[t] = r

    acc_B = float(rewards_B[-50:].mean())
    tau_after_B = net.tau.copy()

    # Phase 3: Test Task A Retention (Identity Mapping, plastic updates frozen)
    rewards_A2 = np.zeros(TRIALS_TEST_A)
    for t in range(TRIALS_TEST_A):
        x = int(rng.integers(net.roles["K"]))
        target = x  # Identity rule
        r = run_sequential_trial(net, x, target, rng, train=False)
        rewards_A2[t] = r

    acc_A2 = float(rewards_A2.mean())
    if acc_A1 >= 0.40:
        retention_pct = (acc_A2 / acc_A1) * 100.0
    else:
        retention_pct = 0.0

    return {
        "acc_A1": acc_A1,
        "acc_B": acc_B,
        "acc_A2": acc_A2,
        "retention_pct": retention_pct,
        "rewards_A1": rewards_A1,
        "rewards_B": rewards_B,
        "rewards_A2": rewards_A2,
        "tau_after_A1": tau_after_A1,
        "tau_after_B": tau_after_B
    }


def run_phase3_suite():
    """Execute Phase 3 sequential learning suite across n=30 seeds."""
    print("=" * 70)
    print(f"Running Phase 3 Sequential Learning & Anti-Forgetting Suite (n={N_SEEDS} seeds)")
    print("=" * 70)

    conditions = {
        "Weight-Only (Axis W)": ("w",),
        "Timescale + Weight (Tau, W)": ("tau", "w"),
        "Closed-Loop Multi-Axis (Tau, Theta, W)": ("tau", "theta", "w")
    }

    results = {}

    for cond_name, axes in conditions.items():
        print(f"\n---> Condition: {cond_name}")
        acc_A1s, acc_Bs, acc_A2s, retentions = [], [], [], []

        for s in range(N_SEEDS):
            res = evaluate_sequential_seed(s, axes=axes)
            acc_A1s.append(res["acc_A1"])
            acc_Bs.append(res["acc_B"])
            acc_A2s.append(res["acc_A2"])
            retentions.append(res["retention_pct"])

        acc_A1s = np.array(acc_A1s)
        acc_Bs = np.array(acc_Bs)
        acc_A2s = np.array(acc_A2s)
        retentions = np.array(retentions)

        mean_ret = retentions.mean()
        std_ret = retentions.std()
        ci95_ret = 1.96 * std_ret / np.sqrt(N_SEEDS)

        results[cond_name] = {
            "acc_A1": acc_A1s,
            "acc_B": acc_Bs,
            "acc_A2": acc_A2s,
            "retention": retentions,
            "mean_retention": mean_ret,
            "std_retention": std_ret,
            "ci95_retention": ci95_ret,
            "mean_acc_A1": acc_A1s.mean(),
            "mean_acc_B": acc_Bs.mean(),
            "mean_acc_A2": acc_A2s.mean()
        }

        print(f"  Task A Initial Acc : {acc_A1s.mean():.3f} +/- {acc_A1s.std():.3f}")
        print(f"  Task B Final Acc   : {acc_Bs.mean():.3f} +/- {acc_Bs.std():.3f}")
        print(f"  Task A Post Test Acc: {acc_A2s.mean():.3f} +/- {acc_A2s.std():.3f}")
        print(f"  Task A Retention % : {mean_ret:.1f}% +/- {ci95_ret:.1f}%")

    # Save results
    np.savez_compressed(
        os.path.join(DATADIR, "phase3_sequential_learning.npz"),
        data=results
    )

    # Plot retention summary comparison figure
    fig, ax = plt.subplots(figsize=(8, 5))
    names = list(conditions.keys())
    x_pos = np.arange(len(names))
    mean_rets = [results[c]["mean_retention"] for c in names]
    ci95s = [results[c]["ci95_retention"] for c in names]

    bars = ax.bar(x_pos, mean_rets, yerr=ci95s, capsize=5, color=["crimson", "darkorange", "forestgreen"], edgecolor="black", alpha=0.85)
    ax.axhline(80.0, color="navy", linestyle="--", linewidth=1.5, label="Target Retention Threshold (80%)")
    ax.set_xticks(x_pos)
    ax.set_xticklabels(names, rotation=15, ha="right")
    ax.set_ylabel("Task A Retention (%)")
    ax.set_ylim(0, 115)
    ax.set_title("Sequential Task A -> Task B -> Task A Retention")
    ax.legend(loc="lower right")
    ax.grid(axis="y", linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "closed_loop_sequential_retention.png"), dpi=200)
    plt.close()

    print("\nSaved result data to result/closed_loop_plasticity/phase3_sequential_learning.npz")
    print("Saved summary figure to docs/figures/closed_loop_sequential_retention.png")

    return results


if __name__ == "__main__":
    run_phase3_suite()
