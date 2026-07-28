"""
Phase 3 Benchmark: Homeostatic Tau-Relaxation & Offline Consolidation.

Evaluates refractoriness recycling and dynamic range protection across K=8 sequential
tasks in n=30 seeds, comparing continuous online learning vs interleaved offline tau consolidation.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import ClosedLoopPlasticityEngine, compute_graph_modularity
from ghca_learn import layered_graph

# Simulation parameters
SETTLE = 10
CUE = 20
WWIN = 25
K_TASKS = 8
TRIALS_PER_TASK = 180
TRIALS_TEST = 50


def task_mapping(task_id, x):
    """K=8 sequential task definitions (alternating mapping rules)."""
    mapping_type = task_id % 4
    if mapping_type == 0:
        return x            # Identity
    elif mapping_type == 1:
        return 1 - x        # Reversal
    elif mapping_type == 2:
        return 0            # Constant 0
    else:
        return 1            # Constant 1


def make_scaled_closed_loop_learner(n_h=100, axes=("tau", "theta", "w", "g"), seed=0):
    """Build a layered graph scaled with n_h hidden nodes."""
    n_s, n_m, K, A = 12, 12, 2, 2
    fanin_sh = min(8, n_s)
    fanin_hm = min(14, n_h)
    hh_k = min(6, max(2, n_h // 4))

    W, plastic, roles = layered_graph(
        n_s=n_s, n_h=n_h, n_m=n_m, K=K, A=A,
        fanin_sh=fanin_sh, fanin_hm=fanin_hm, hh_k=hh_k,
        w_hm=0.6, w_hh=0.25, seed=seed
    )

    N_actual = W.shape[0]

    cfg = dict(act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.08, eta_tau=0.02, eta_theta=0.01)

    engine = ClosedLoopPlasticityEngine(
        W, plastic, roles,
        axes=axes,
        act=cfg["act"], pas=cfg["pas"], theta=cfg["theta"], p_s=cfg["p_s"],
        eta_w=cfg["eta_w"], eta_tau=cfg["eta_tau"], eta_theta=cfg["eta_theta"],
        w_prune=0.05, w_sprout=0.10, max_sprouts=5,
        seed=seed
    )

    return engine, N_actual


def run_trial(net, x, target_action, rng, train=True):
    """Execute trial step sequence."""
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


def evaluate_consolidation_seed(seed, n_h=100, use_consolidation=False):
    """
    Run 8-task sequential learning protocol for a single seed.
    Track tau saturation (% nodes with tau >= 18.0), Task 1 retention, Task 8 accuracy.
    """
    axes = ("tau", "theta", "w", "g")
    net, N_actual = make_scaled_closed_loop_learner(n_h=n_h, axes=axes, seed=seed)
    rng = np.random.default_rng(seed + 8000)

    task_final_accs = np.zeros(K_TASKS)
    tau_means = np.zeros(K_TASKS)
    tau_saturations = np.zeros(K_TASKS)

    for task_idx in range(K_TASKS):
        if task_idx > 0 and use_consolidation:
            # Offline sleep / tau consolidation phase
            net.consolidate_tau(decay_rate=0.03)

        task_rewards = np.zeros(TRIALS_PER_TASK)
        for t in range(TRIALS_PER_TASK):
            x = int(rng.integers(2))
            target = task_mapping(task_idx, x)
            r = run_trial(net, x, target, rng, train=True)
            task_rewards[t] = r

        task_final_accs[task_idx] = float(task_rewards[-40:].mean())
        tau_means[task_idx] = float(net.tau.mean())
        tau_saturations[task_idx] = float(np.mean(net.tau >= 18.0) * 100.0)

    # Test Task 1 retention
    test_rewards = np.zeros(TRIALS_TEST)
    for t in range(TRIALS_TEST):
        x = int(rng.integers(2))
        target = task_mapping(0, x)
        r = run_trial(net, x, target, rng, train=False)
        test_rewards[t] = r

    acc_t1_test = float(test_rewards.mean())
    acc_t1_initial = task_final_accs[0]

    if acc_t1_initial >= 0.40:
        retention_t1_pct = (acc_t1_test / acc_t1_initial) * 100.0
    else:
        retention_t1_pct = 0.0

    return {
        "acc_t1_initial": acc_t1_initial,
        "acc_t1_test": acc_t1_test,
        "acc_t8_final": task_final_accs[-1],
        "retention_t1_pct": retention_t1_pct,
        "tau_means": tau_means,
        "tau_saturations": tau_saturations
    }


def run_consolidation_experiment(n_seeds=30, n_h=100):
    """Run Phase 3 Offline Tau Consolidation experiment across n=30 seeds."""
    print("=" * 80)
    print(f"Phase 3: Homeostatic Tau-Relaxation & Consolidation (K=8 tasks, n={n_seeds} seeds)")
    print("=" * 80)

    regimes = [
        ("No Consolidation (Continuous Online)", False),
        ("With Tau-Relaxation (Offline Sleep)", True)
    ]

    results = {}

    for regime_name, use_consol in regimes:
        print(f"\n---> Regime: {regime_name}")

        init_accs = []
        test_accs = []
        t8_accs = []
        retentions = []
        all_tau_means = []
        all_tau_sats = []

        for s in range(n_seeds):
            res = evaluate_consolidation_seed(seed=s, n_h=n_h, use_consolidation=use_consol)
            init_accs.append(res['acc_t1_initial'])
            test_accs.append(res['acc_t1_test'])
            t8_accs.append(res['acc_t8_final'])
            retentions.append(res['retention_t1_pct'])
            all_tau_means.append(res['tau_means'])
            all_tau_sats.append(res['tau_saturations'])

        init_m, init_s = np.mean(init_accs), np.std(init_accs)
        test_m, test_s = np.mean(test_accs), np.std(test_accs)
        t8_m, t8_s = np.mean(t8_accs), np.std(t8_accs)
        ret_m, ret_s = np.mean(retentions), np.std(retentions)
        sat_final_m = np.mean([sat[-1] for sat in all_tau_sats])

        print(f"    Task 1 Initial Acc      : {init_m:.3f} +/- {init_s:.3f}")
        print(f"    Task 1 Retention %      : {ret_m:.1f}% +/- {ret_s / np.sqrt(n_seeds):.1f}%")
        print(f"    Task 8 Final Acc        : {t8_m:.3f} +/- {t8_s:.3f}")
        print(f"    Final Tau Saturation %  : {sat_final_m:.1f}%")

        results[regime_name] = {
            'init_acc': (init_m, init_s),
            'test_acc': (test_m, test_s),
            't8_acc': (t8_m, t8_s),
            'retention': (ret_m, ret_s),
            'tau_means': np.mean(all_tau_means, axis=0),
            'tau_sats': np.mean(all_tau_sats, axis=0)
        }

    # Save archived data
    out_dir = "/home/dognosis/Documents/ghca/result/closed_loop_plasticity"
    os.makedirs(out_dir, exist_ok=True)
    np.savez(os.path.join(out_dir, "tau_consolidation.npz"), **results)
    print(f"\n[OK] Archived simulation data to {out_dir}/tau_consolidation.npz")

    # Generate visualization
    fig_dir = "/home/dognosis/Documents/ghca/docs/figures"
    os.makedirs(fig_dir, exist_ok=True)

    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5), dpi=300)

    tasks_x = np.arange(1, K_TASKS + 1)

    ax1.plot(tasks_x, results["No Consolidation (Continuous Online)"]['tau_sats'], 'o-',
             color='#E74C3C', linewidth=2.0, label='No Consolidation')
    ax1.plot(tasks_x, results["With Tau-Relaxation (Offline Sleep)"]['tau_sats'], 's-',
             color='#2ECC71', linewidth=2.0, label='With Tau-Relaxation')
    ax1.set_xlabel("Sequential Task Index (K=1..8)", fontsize=11, fontweight='bold')
    ax1.set_ylabel("Nodes with Tau >= 18.0 (%)", fontsize=11, fontweight='bold')
    ax1.set_title("Refractoriness Saturation Dynamics", fontsize=12, fontweight='bold')
    ax1.grid(linestyle=':', alpha=0.6)
    ax1.legend(loc='upper left', frameon=True)

    labels = ["No Consolidation", "With Tau-Relaxation"]
    ret_vals = [results["No Consolidation (Continuous Online)"]['retention'][0],
                results["With Tau-Relaxation (Offline Sleep)"]['retention'][0]]
    ret_errs = [results["No Consolidation (Continuous Online)"]['retention'][1] / np.sqrt(n_seeds),
                results["With Tau-Relaxation (Offline Sleep)"]['retention'][1] / np.sqrt(n_seeds)]

    bars = ax2.bar(labels, ret_vals, yerr=ret_errs, capsize=6, color=['#E74C3C', '#2ECC71'], alpha=0.85, edgecolor='black')
    ax2.set_ylabel("Task 1 Retention (%) after 8 Tasks", fontsize=11, fontweight='bold')
    ax2.set_title("Long-Horizon Retention (K=8)", fontsize=12, fontweight='bold')
    ax2.set_ylim(0, 110)
    ax2.axhline(80, color='crimson', linestyle='--', linewidth=1.5, label='Target Retention (80%)')
    ax2.grid(axis='y', linestyle=':', alpha=0.6)
    ax2.legend(loc='upper left', frameon=True)

    for bar, val in zip(bars, ret_vals):
        height = bar.get_height()
        ax2.annotate(f'{val:.1f}%',
                     xy=(bar.get_x() + bar.get_width() / 2, height),
                     xytext=(0, 5), textcoords="offset points",
                     ha='center', va='bottom', fontsize=11, fontweight='bold')

    plt.tight_layout()
    fig_path = os.path.join(fig_dir, "tau_consolidation.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"[OK] Saved figure to {fig_path}")


if __name__ == "__main__":
    run_consolidation_experiment(n_seeds=30, n_h=100)
