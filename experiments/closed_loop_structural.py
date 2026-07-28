"""
Phase 2 Benchmark: Axis G Structural Plasticity & Homeostatic Consolidation.

Evaluates reward-modulated rewiring (edge pruning & co-activity sprouting) and
offline tau consolidation in K=5 sequential task learning across n=30 seeds.
"""

import os
import sys
import numpy as np
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import ClosedLoopPlasticityEngine, compute_graph_modularity, compute_readout_independence_ratio
from ghca_learn import layered_graph

# Simulation trial parameters
SETTLE = 10
CUE = 10
WWIN = 5
K_TASKS = 5
TRIALS_PER_TASK = 200
TRIALS_TEST = 50


def task_mapping(task_id, x):
    """
    K=5 task mapping definitions:
    T1: Identity (0->0, 1->1)
    T2: Reversal (0->1, 1->0)
    T3: Constant (0->0, 1->0)
    T4: Shift (0->1, 1->1)
    T5: Inverted Reversal (0->1, 1->0)
    """
    if task_id == 0:
        return x
    elif task_id == 1:
        return 1 - x
    elif task_id == 2:
        return 0
    elif task_id == 3:
        return 1
    elif task_id == 4:
        return 1 - x
    return x


def make_scaled_closed_loop_learner(n_h=100, axes=("tau", "theta", "w"), seed=0):
    """Build a layered graph scaled with n_h hidden nodes and return learner with graph topology metrics."""
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
    m_edges = int(np.count_nonzero(W > 0))
    beta_1 = m_edges - N_actual + 1

    cfg = dict(act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.08, eta_tau=0.02, eta_theta=0.01)

    engine = ClosedLoopPlasticityEngine(
        W, plastic, roles,
        axes=axes,
        act=cfg["act"], pas=cfg["pas"], theta=cfg["theta"], p_s=cfg["p_s"],
        eta_w=cfg["eta_w"], eta_tau=cfg["eta_tau"], eta_theta=cfg["eta_theta"],
        w_prune=0.05, w_sprout=0.10, max_sprouts=5,
        seed=seed
    )

    return engine, N_actual, beta_1


def compute_chance_corrected_retention(acc_test, acc_init, chance=0.50):
    """
    Computes chance-corrected retention ratio: R = (acc_test - chance) / (acc_init - chance).
    Returns 0.0 if acc_init <= chance.
    """
    denom = acc_init - chance
    if denom <= 1e-6:
        return 0.0
    return float((acc_test - chance) / denom)


def run_trial(net, x, target_action, rng, train=True):
    """Execute a single trial step sequence with sensory drive and RPE update."""
    net.reset_traces()
    for _ in range(SETTLE):
        net.step_learn(None)

    feats_before = net.features()
    V = net.value()

    sc = np.zeros(net.roles["A"])
    drive = net.sensory_drive(x)

    for _ in range(CUE):
        net.step_learn(drive)
        sc += net.motor_scores()

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


def evaluate_structural_seed(seed, n_h=100, axes=("tau", "theta", "w"), use_consolidation=False):
    """
    Run 5-task sequential learning protocol for a single seed under a given plasticity/consolidation regime.
    """
    net, N_actual, beta_1 = make_scaled_closed_loop_learner(n_h=n_h, axes=axes, seed=seed)
    rng = np.random.default_rng(seed + 5000)

    Q_initial = compute_graph_modularity(net.W)

    task_final_accs = np.zeros(K_TASKS)

    # Sequential Training Phase across 5 tasks
    for task_idx in range(K_TASKS):
        if task_idx > 0 and use_consolidation:
            # Apply offline tau-consolidation between task transitions
            net.consolidate_tau(decay_rate=0.03)

        task_rewards = np.zeros(TRIALS_PER_TASK)
        for t in range(TRIALS_PER_TASK):
            x = int(rng.integers(2))
            target = task_mapping(task_idx, x)
            r = run_trial(net, x, target, rng, train=True)
            task_rewards[t] = r

        task_final_accs[task_idx] = float(task_rewards[-50:].mean())

    Q_final = compute_graph_modularity(net.W)

    # Retention Testing Phase for Task 1 (Identity) with plastic updates frozen
    test_rewards = np.zeros(TRIALS_TEST)
    for t in range(TRIALS_TEST):
        x = int(rng.integers(2))
        target = task_mapping(0, x)
        r = run_trial(net, x, target, rng, train=False)
        test_rewards[t] = r

    acc_t1_test = float(test_rewards.mean())
    acc_t1_initial = task_final_accs[0]

    retention_chance_corr = compute_chance_corrected_retention(acc_t1_test, acc_t1_initial, chance=0.50) * 100.0

    if acc_t1_initial >= 0.40:
        retention_t1_pct = (acc_t1_test / acc_t1_initial) * 100.0
    else:
        retention_t1_pct = 0.0

    return {
        "acc_t1_initial": acc_t1_initial,
        "acc_t1_test": acc_t1_test,
        "retention_t1_pct": retention_t1_pct,
        "retention_chance_corr": retention_chance_corr,
        "Q_initial": Q_initial,
        "Q_final": Q_final
    }


def run_structural_experiment(n_seeds=30, n_h=100):
    """Run Phase 2 Structural Plasticity & Consolidation Sweep across seeds."""
    print("=" * 80)
    print(f"Phase 2: Axis G Structural Plasticity & Offline Consolidation (n={n_seeds} seeds)")
    print("=" * 80)

    conditions = [
        ("Multi-Axis Base (Tau, Theta, W)", ("tau", "theta", "w"), False),
        ("Multi-Axis + Axis G Rewiring (Tau, Theta, W, G)", ("tau", "theta", "w", "g"), False),
        ("Multi-Axis + Axis G + Consolidation", ("tau", "theta", "w", "g"), True)
    ]

    results = {}

    for cond_name, axes, use_consolidation in conditions:
        print(f"\n---> Condition: {cond_name}")

        init_accs = []
        test_accs = []
        retentions = []
        retentions_cc = []
        q_inits = []
        q_finals = []

        for s in range(n_seeds):
            res = evaluate_structural_seed(seed=s, n_h=n_h, axes=axes, use_consolidation=use_consolidation)
            init_accs.append(res['acc_t1_initial'])
            test_accs.append(res['acc_t1_test'])
            retentions.append(res['retention_t1_pct'])
            retentions_cc.append(res['retention_chance_corr'])
            q_inits.append(res['Q_initial'])
            q_finals.append(res['Q_final'])

        init_m, init_s = np.mean(init_accs), np.std(init_accs)
        test_m, test_s = np.mean(test_accs), np.std(test_accs)
        ret_m, ret_s = np.mean(retentions), np.std(retentions)
        ret_cc_m, ret_cc_s = np.mean(retentions_cc), np.std(retentions_cc)
        q_in_m = np.mean(q_inits)
        q_out_m = np.mean(q_finals)

        print(f"    Task 1 Initial Acc   : {init_m:.3f} +/- {init_s:.3f}")
        print(f"    Task 1 Post Test Acc : {test_m:.3f} +/- {test_s:.3f}")
        print(f"    Task 1 Raw Ret %     : {ret_m:.1f}% +/- {ret_s / np.sqrt(n_seeds):.1f}% (std={ret_s:.1f}%)")
        print(f"    Task 1 Chance-Corr % : {ret_cc_m:.1f}% +/- {ret_cc_s / np.sqrt(n_seeds):.1f}% (std={ret_cc_s:.1f}%)")
        print(f"    Graph Modularity Q   : {q_in_m:.4f} (initial) -> {q_out_m:.4f} (final)")

        results[cond_name] = {
            'init_acc': (init_m, init_s),
            'post_acc': (test_m, test_s),
            'retention': (ret_m, ret_s),
            'retention_chance_corr': (ret_cc_m, ret_cc_s),
            'modularity': (q_in_m, q_out_m),
            'raw_init_accs': np.array(init_accs),
            'raw_test_accs': np.array(test_accs),
            'raw_retentions': np.array(retentions),
            'raw_retentions_cc': np.array(retentions_cc),
            'raw_q_inits': np.array(q_inits),
            'raw_q_finals': np.array(q_finals),
        }

    # Save archived data
    out_dir = "/home/dognosis/Documents/ghca/result/closed_loop_plasticity"
    os.makedirs(out_dir, exist_ok=True)
    np.savez(os.path.join(out_dir, "structural_plasticity.npz"), **results)
    print(f"\n[OK] Archived simulation data to {out_dir}/structural_plasticity.npz")

    # Generate visualization
    fig_dir = "/home/dognosis/Documents/ghca/docs/figures"
    os.makedirs(fig_dir, exist_ok=True)

    fig, ax = plt.subplots(figsize=(9, 5), dpi=300)
    cond_labels = [c[0] for c in conditions]
    ret_means = [results[c]['retention'][0] for c in cond_labels]
    ret_errs = [results[c]['retention'][1] / np.sqrt(n_seeds) for c in cond_labels]

    colors = ['#4A90E2', '#36B37E', '#9B51E0']
    bars = ax.bar(cond_labels, ret_means, yerr=ret_errs, capsize=6, color=colors, alpha=0.85, edgecolor='black', linewidth=1.2)

    ax.set_ylabel("Task 1 Retention (%) after 5 Tasks", fontsize=12, fontweight='bold')
    ax.set_title("Axis G Structural Plasticity & Offline Tau Consolidation (n=30 seeds)", fontsize=13, fontweight='bold', pad=12)
    ax.set_ylim(0, 100)
    ax.axhline(80, color='crimson', linestyle='--', linewidth=1.5, label='Target Retention (80%)')
    ax.grid(axis='y', linestyle=':', alpha=0.6)

    for bar, mean in zip(bars, ret_means):
        height = bar.get_height()
        ax.annotate(f'{mean:.1f}%',
                    xy=(bar.get_x() + bar.get_width() / 2, height),
                    xytext=(0, 5), textcoords="offset points",
                    ha='center', va='bottom', fontsize=11, fontweight='bold')

    ax.legend(loc='upper left', frameon=True)
    plt.tight_layout()

    fig_path = os.path.join(fig_dir, "structural_plasticity_g.png")
    plt.savefig(fig_path)
    plt.close()
    print(f"[OK] Saved figure to {fig_path}")


if __name__ == "__main__":
    run_structural_experiment(n_seeds=30, n_h=100)
