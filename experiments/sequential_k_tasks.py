"""
Phase 1: Long-Horizon Sequential Task Saturation (K=5 Tasks) & Capacity Scaling Sweep.

Evaluates 5-task sequential learning and catastrophic forgetting mitigation:
    Task 1 (Identity) -> Task 2 (Reversal) -> Task 3 (Constant 0) -> Task 4 (Constant 1) -> Task 5 (Reversal) -> Task 1 Retention Test
under closed-loop multi-axis substrate plasticity (tau, theta, W) vs baselines.

Sweeps substrate graph sizes N in {68, 98, 148, 248} (varying hidden recurrent nodes n_h in {20, 50, 100, 200})
and correlates memory retention with graph topological circuit rank beta_1 = m - N + 1.

House Rules Compliance:
    - Explicit default_rng(seed) threading everywhere.
    - Report per-seed spreads, standard deviations, and 95% confidence intervals across n=30 seeds.
    - Substrate-only dynamics; readout features used solely for metric logging.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import ClosedLoopPlasticityEngine
from ghca_learn import layered_graph

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "closed_loop_plasticity")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

N_SEEDS = 30
K_TASKS = 5
TRIALS_PER_TASK = 300
TRIALS_TEST = 100
SETTLE, CUE, WWIN = 10, 3, 6


def task_mapping(task_id, x):
    """5 sequential task rule definitions."""
    if task_id == 0:    # T1: Identity (0->0, 1->1)
        return x
    elif task_id == 1:  # T2: Reversal (0->1, 1->0)
        return 1 - x
    elif task_id == 2:  # T3: Constant 0 (0->0, 1->0)
        return 0
    elif task_id == 3:  # T4: Constant 1 (0->1, 1->1)
        return 1
    elif task_id == 4:  # T5: Reversal (0->1, 1->0)
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
    beta_1 = m_edges - N_actual + 1  # Circuit rank (m - N + 1)

    cfg = dict(act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.08, eta_tau=0.02, eta_theta=0.01)

    engine = ClosedLoopPlasticityEngine(
        W, plastic, roles,
        axes=axes,
        act=cfg["act"], pas=cfg["pas"], theta=cfg["theta"], p_s=cfg["p_s"],
        eta_w=cfg["eta_w"], eta_tau=cfg["eta_tau"], eta_theta=cfg["eta_theta"],
        seed=seed
    )

    return engine, N_actual, beta_1


def run_trial(net, x, target_action, rng, train=True):
    """Execute a single trial step sequence."""
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


def evaluate_k_task_seed(seed, n_h=100, axes=("tau", "theta", "w")):
    """
    Run 5-task sequential learning protocol for a single seed:
    Train T1 -> Train T2 -> Train T3 -> Train T4 -> Train T5 -> Test T1 Retention
    """
    net, N_actual, beta_1 = make_scaled_closed_loop_learner(n_h=n_h, axes=axes, seed=seed)
    rng = np.random.default_rng(seed + 5000)

    task_final_accs = np.zeros(K_TASKS)

    # Sequential Training Phase across 5 tasks
    for task_idx in range(K_TASKS):
        task_rewards = np.zeros(TRIALS_PER_TASK)
        for t in range(TRIALS_PER_TASK):
            x = int(rng.integers(2))
            target = task_mapping(task_idx, x)
            r = run_trial(net, x, target, rng, train=True)
            task_rewards[t] = r

        task_final_accs[task_idx] = float(task_rewards[-50:].mean())

    # Retention Testing Phase for Task 1 (Identity) with plastic updates frozen
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
        "N_actual": N_actual,
        "beta_1": beta_1,
        "task_final_accs": task_final_accs,
        "acc_t1_test": acc_t1_test,
        "retention_t1_pct": retention_t1_pct
    }


def run_sequential_k_suite():
    """Execute complete 5-task sequential saturation & capacity sweep across n=30 seeds."""
    print("=" * 80)
    print(f"Phase 1: Sequential 5-Task Saturation & Capacity Scaling Sweep (n={N_SEEDS} seeds)")
    print("=" * 80)

    n_h_list = [20, 50, 100, 200]
    conditions = {
        "Weight-Only (Axis W)": ("w",),
        "Timescale + Weight (Tau, W)": ("tau", "w"),
        "Closed-Loop Multi-Axis (Tau, Theta, W)": ("tau", "theta", "w")
    }

    results = {}

    for nh in n_h_list:
        print(f"\n==================== Hidden Recurrent Nodes n_h = {nh} ====================")
        results[nh] = {}

        for cond_name, axes in conditions.items():
            print(f"  ---> Condition: {cond_name}")
            ret_t1_list = []
            task_accs_list = []
            test_accs_list = []
            beta_1_val = 0
            N_act_val = 0

            for s in range(N_SEEDS):
                res = evaluate_k_task_seed(s, n_h=nh, axes=axes)
                ret_t1_list.append(res["retention_t1_pct"])
                task_accs_list.append(res["task_final_accs"])
                test_accs_list.append(res["acc_t1_test"])
                beta_1_val = res["beta_1"]
                N_act_val = res["N_actual"]

            ret_t1_arr = np.array(ret_t1_list)
            task_accs_arr = np.array(task_accs_list)
            test_accs_arr = np.array(test_accs_list)

            m_ret_t1, std_ret_t1 = ret_t1_arr.mean(), ret_t1_arr.std()
            ci95_t1 = 1.96 * std_ret_t1 / np.sqrt(N_SEEDS)

            results[nh][cond_name] = {
                "N_actual": N_act_val,
                "beta_1": beta_1_val,
                "retention_t1": ret_t1_arr,
                "mean_retention_t1": m_ret_t1,
                "std_retention_t1": std_ret_t1,
                "ci95_t1": ci95_t1,
                "mean_task_accs": task_accs_arr.mean(axis=0),
                "std_task_accs": task_accs_arr.std(axis=0),
                "mean_test_acc": test_accs_arr.mean(),
                "std_test_acc": test_accs_arr.std()
            }

            print(f"    Substrate N={N_act_val} (Circuit Rank beta_1 = {beta_1_val})")
            print(f"    Task 1 Initial Acc   : {task_accs_arr[:, 0].mean():.3f} +/- {task_accs_arr[:, 0].std():.3f}")
            print(f"    Task 1 Post Test Acc : {test_accs_arr.mean():.3f} +/- {test_accs_arr.std():.3f}")
            print(f"    Task 1 Retention %   : {m_ret_t1:.1f}% +/- {ci95_t1:.1f}% (std={std_ret_t1:.1f}%)")

    # Save complete dataset
    output_data_path = os.path.join(DATADIR, "sequential_k_tasks.npz")
    np.savez_compressed(output_data_path, data=results)
    print(f"\n[OK] Archived simulation data to {output_data_path}")

    # Generate Publication Capacity Scaling Figure
    fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(13, 5))

    colors = {"Weight-Only (Axis W)": "crimson",
              "Timescale + Weight (Tau, W)": "darkorange",
              "Closed-Loop Multi-Axis (Tau, Theta, W)": "forestgreen"}

    # Subplot 1: Task 1 Retention % vs Substrate Size N
    for cond_name in conditions.keys():
        n_vals = [results[nh][cond_name]["N_actual"] for nh in n_h_list]
        ret_vals = [results[nh][cond_name]["mean_retention_t1"] for nh in n_h_list]
        ci_vals = [results[nh][cond_name]["ci95_t1"] for nh in n_h_list]

        ax1.errorbar(n_vals, ret_vals, yerr=ci_vals, fmt="-o", color=colors[cond_name],
                     label=cond_name, linewidth=2, capsize=4, markersize=6)

    ax1.axhline(80.0, color="navy", linestyle="--", linewidth=1.5, label="Target Threshold (80%)")
    ax1.set_xlabel("Substrate Node Count N")
    ax1.set_ylabel("Task 1 Retention (%) after 5 Tasks")
    ax1.set_title("Capacity Scaling: Retention vs Substrate Size N")
    ax1.set_ylim(0, 115)
    ax1.grid(True, linestyle=":", alpha=0.6)
    ax1.legend(loc="lower right", fontsize=9)

    # Subplot 2: Task 1 Retention % vs Circuit Rank beta_1
    for cond_name in conditions.keys():
        beta_vals = [results[nh][cond_name]["beta_1"] for nh in n_h_list]
        ret_vals = [results[nh][cond_name]["mean_retention_t1"] for nh in n_h_list]

        ax2.plot(beta_vals, ret_vals, "-s", color=colors[cond_name],
                 label=cond_name, linewidth=2, markersize=6)

    ax2.axhline(80.0, color="navy", linestyle="--", linewidth=1.5)
    ax2.set_xlabel("Graph Circuit Rank beta_1 (m - N + 1)")
    ax2.set_ylabel("Task 1 Retention (%) after 5 Tasks")
    ax2.set_title("Topological Capacity: Retention vs Circuit Rank beta_1")
    ax2.set_ylim(0, 115)
    ax2.grid(True, linestyle=":", alpha=0.6)

    plt.tight_layout()
    fig_path = os.path.join(FIGDIR, "sequential_k_saturation.png")
    plt.savefig(fig_path, dpi=200)
    plt.close()
    print(f"[OK] Saved figure to {fig_path}")

    return results


if __name__ == "__main__":
    run_sequential_k_suite()
