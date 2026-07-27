"""
Phase 2 Single-Task Benchmarks & Readout Independence Evaluation.

Evaluates E1 sensorimotor conditioning and E5 executive contextual switching
under multi-axis closed-loop substrate plasticity across n=30 independent seeds.

Quantifies substrate credit assignment via the Readout Independence Ratio:
    RIR = Accuracy(Fixed Identity Readout) / Accuracy(Trained Linear Decoder)

Target Acceptance Criteria:
    - RIR > 0.85 across n=30 seeds on single-task benchmarks.
    - Full reporting of per-seed spreads, standard deviations, and 95% CIs.
    - Zero global NumPy RNG usage (explicit default_rng seed threading).
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

try:
    from sklearn.linear_model import LogisticRegression
    HAS_SKLEARN = True
except ImportError:
    HAS_SKLEARN = False

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_plasticity import (
    ClosedLoopPlasticityEngine,
    compute_readout_independence_ratio,
    make_closed_loop_learner
)
from ghca_learn import layered_graph
from experiments.e5_executive import build as build_e5_graph, CTX_CUE, L_RING, N_S, N_H, N_M

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "closed_loop_plasticity")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# Evaluation parameters
N_SEEDS = 30
E1_TRIALS = 300
E5_BLOCKS = 6
E5_TRIALS_PER_BLOCK = 50
E5_TOTAL_TRIALS = E5_BLOCKS * E5_TRIALS_PER_BLOCK
SETTLE, CUE, WWIN = 10, 3, 6


def fit_linear_decoder(X_train, y_train, X_test, y_test, seed=0):
    """
    Fits a linear decoder on hidden activity traces. Uses scikit-learn if available,
    otherwise falls back to a pure NumPy ridge classifier.
    """
    n_classes = len(np.unique(y_train))
    if n_classes <= 1 or X_train.shape[0] == 0:
        return 1.0

    if HAS_SKLEARN:
        clf = LogisticRegression(max_iter=200, random_state=seed)
        clf.fit(X_train, y_train)
        return float(clf.score(X_test, y_test))

    # Pure NumPy Ridge Classification Fallback
    Y_tr = np.eye(n_classes)[y_train]
    X_tr = np.hstack([X_train, np.ones((X_train.shape[0], 1))])
    X_te = np.hstack([X_test, np.ones((X_test.shape[0], 1))])
    
    alpha = 1.0
    I = np.eye(X_tr.shape[1])
    I[-1, -1] = 0.0  # do not penalize bias
    
    W = np.linalg.solve(X_tr.T @ X_tr + alpha * I, X_tr.T @ Y_tr)
    y_pred = np.argmax(X_te @ W, axis=1)
    return float(np.mean(y_pred == y_test))


# ----------------------------------------------------------------------------
# E1 Single-Task Benchmark Runner
# ----------------------------------------------------------------------------

def run_e1_trial(net, x, rng):
    """Run a single E1 trial and record hidden state features for readout comparison."""
    net.reset_traces()
    for _ in range(SETTLE):
        net.step_learn(None)
    
    feats_before = net.features()
    V = net.value()

    # Cue presentation
    for _ in range(CUE):
        net.step_learn(net.sensory_drive(x))

    sc = np.zeros(net.roles["A"])
    hidden_act_sum = np.zeros(len(net.roles["hidden"]))

    for _ in range(WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
        hidden_act_sum += net.active_mask()[net.roles["hidden"]].astype(float)

    # Substrate fixed readout action (argmax motor channel score)
    if sc.sum() > 0:
        action_fixed = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc))))
    else:
        action_fixed = -1

    reward = 1.0 if action_fixed == x else 0.0
    delta = reward - V

    net.learn(delta)
    net.update_critic(delta, feats_before)

    return reward, hidden_act_sum / float(WWIN)


def evaluate_e1_seed(seed, axes=("tau", "theta", "w")):
    """Run E1 benchmark for a single seed and return fixed readout acc, trained readout acc, and RIR."""
    net = make_closed_loop_learner(axes=axes, seed=seed)
    rng = np.random.default_rng(seed + 100)

    rewards_fixed = np.zeros(E1_TRIALS)
    hidden_features = np.zeros((E1_TRIALS, len(net.roles["hidden"])))
    labels = np.zeros(E1_TRIALS, dtype=int)

    for t in range(E1_TRIALS):
        x = int(rng.integers(net.roles["K"]))
        labels[t] = x
        r, h_feat = run_e1_trial(net, x, rng)
        rewards_fixed[t] = r
        hidden_features[t] = h_feat

    # Evaluate trained linear decoder on the second half of trials
    eval_start = E1_TRIALS // 2
    X_train, y_train = hidden_features[:eval_start], labels[:eval_start]
    X_test, y_test = hidden_features[eval_start:], labels[eval_start:]

    trained_acc = fit_linear_decoder(X_train, y_train, X_test, y_test, seed=seed)
    fixed_acc = float(rewards_fixed[eval_start:].mean())
    rir = compute_readout_independence_ratio(fixed_acc, trained_acc)

    return fixed_acc, trained_acc, rir, rewards_fixed


# ----------------------------------------------------------------------------
# E5 Executive Contextual Switch Runner
# ----------------------------------------------------------------------------

def evaluate_e5_seed(seed, axes=("tau", "theta", "w")):
    """Run E5 executive switching for a single seed using ClosedLoopPlasticityEngine."""
    W, plastic, roles, theta_init = build_e5_graph(seed=seed)
    net = ClosedLoopPlasticityEngine(
        W, plastic, roles,
        axes=axes,
        act=2, pas=8, theta=theta_init, p_s=3e-3,
        eta_w=0.04, eta_tau=0.10, eta_theta=0.02,
        seed=seed
    )

    rng = np.random.default_rng(seed + 500)
    rewards_fixed = np.zeros(E5_TOTAL_TRIALS)
    hidden_features = np.zeros((E5_TOTAL_TRIALS, len(roles["hidden"])))
    labels = np.zeros(E5_TOTAL_TRIALS, dtype=int)

    trial_idx = 0
    for block in range(E5_BLOCKS):
        rule = block % 2  # Alternate identity (rule 0) and reversal (rule 1)
        
        # Context ring ignition cue
        ctx_drive = np.zeros(net.N, dtype=bool)
        ctx_drive[roles["rings"][rule][:2]] = True
        for _ in range(CTX_CUE):
            net.step_learn(ctx_drive)

        for _ in range(E5_TRIALS_PER_BLOCK):
            x = int(rng.integers(roles["K"]))
            correct_action = x if rule == 0 else (1 - x)
            labels[trial_idx] = correct_action

            net.reset_traces()
            for _ in range(SETTLE):
                net.step_learn(None)

            feats_before = net.features()
            V = net.value()

            sens_drive = np.zeros(net.N, dtype=bool)
            sens_drive[roles["sensory"][x]] = True
            for _ in range(CUE):
                net.step_learn(sens_drive)

            sc = np.zeros(roles["A"])
            h_sum = np.zeros(len(roles["hidden"]))
            for _ in range(WWIN):
                net.step_learn(None)
                sc += net.motor_scores()
                h_sum += net.active_mask()[roles["hidden"]].astype(float)

            action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
            r = 1.0 if action == correct_action else 0.0
            delta = r - V

            net.learn(delta)
            net.update_critic(delta, feats_before)

            rewards_fixed[trial_idx] = r
            hidden_features[trial_idx] = h_sum / float(WWIN)
            trial_idx += 1

    # Evaluate second half performance
    eval_start = E5_TOTAL_TRIALS // 2
    X_train, y_train = hidden_features[:eval_start], labels[:eval_start]
    X_test, y_test = hidden_features[eval_start:], labels[eval_start:]

    trained_acc = fit_linear_decoder(X_train, y_train, X_test, y_test, seed=seed)
    fixed_acc = float(rewards_fixed[eval_start:].mean())
    rir = compute_readout_independence_ratio(fixed_acc, trained_acc)

    return fixed_acc, trained_acc, rir, rewards_fixed


# ----------------------------------------------------------------------------
# Suite Execution & Plotting
# ----------------------------------------------------------------------------

def run_phase2_suite():
    """Execute Phase 2 suite across n=30 seeds and log statistics."""
    print("=" * 70)
    print(f"Running Phase 2 Single-Task Plasticity Suite (n={N_SEEDS} seeds)")
    print("=" * 70)

    axis_conditions = {
        "Control (No Plasticity)": (),
        "Axis W Only": ("w",),
        "Axis Tau Only": ("tau",),
        "Axis Theta Only": ("theta",),
        "Closed-Loop Multi-Axis": ("tau", "theta", "w")
    }

    results = {}

    for task_name, eval_fn in [("E1_Sensorimotor", evaluate_e1_seed), ("E5_Executive_Switch", evaluate_e5_seed)]:
        results[task_name] = {}
        print(f"\n---> Task Benchmark: {task_name}")

        for cond_label, axes in axis_conditions.items():
            fixed_accs = np.zeros(N_SEEDS)
            trained_accs = np.zeros(N_SEEDS)
            rirs = np.zeros(N_SEEDS)

            for s in range(N_SEEDS):
                f_acc, t_acc, rir, _ = eval_fn(s, axes=axes)
                fixed_accs[s] = f_acc
                trained_accs[s] = t_acc
                rirs[s] = rir

            mean_rir = rirs.mean()
            std_rir = rirs.std()
            ci95_rir = 1.96 * std_rir / np.sqrt(N_SEEDS)

            results[task_name][cond_label] = {
                "fixed_acc": fixed_accs,
                "trained_acc": trained_accs,
                "rir": rirs,
                "mean_rir": mean_rir,
                "std_rir": std_rir,
                "ci95_rir": ci95_rir,
                "mean_fixed_acc": fixed_accs.mean(),
                "std_fixed_acc": fixed_accs.std()
            }

            print(f"  [{cond_label:<25}] Fixed Acc: {fixed_accs.mean():.3f} +/- {fixed_accs.std():.3f} | "
                  f"Trained Acc: {trained_accs.mean():.3f} | RIR: {mean_rir:.3f} +/- {ci95_rir:.3f}")

    # Save quantitative data
    np.savez_compressed(
        os.path.join(DATADIR, "phase2_single_task.npz"),
        data=results
    )

    # Plot RIR Comparison Figure
    fig, axes_fig = plt.subplots(1, 2, figsize=(12, 5))
    cond_labels = list(axis_conditions.keys())
    x_pos = np.arange(len(cond_labels))

    for idx, (task_name, ax) in enumerate(zip(["E1_Sensorimotor", "E5_Executive_Switch"], axes_fig)):
        mean_rirs = [results[task_name][cond]["mean_rir"] for cond in cond_labels]
        ci95s = [results[task_name][cond]["ci95_rir"] for cond in cond_labels]

        bars = ax.bar(x_pos, mean_rirs, yerr=ci95s, capsize=5, color="steelblue", edgecolor="black", alpha=0.85)
        ax.axhline(0.85, color="crimson", linestyle="--", linewidth=1.5, label="Target RIR Threshold (0.85)")
        ax.set_xticks(x_pos)
        ax.set_xticklabels(cond_labels, rotation=25, ha="right")
        ax.set_ylim(0.0, 1.15)
        ax.set_ylabel("Readout Independence Ratio (RIR)")
        ax.set_title(f"Readout Independence — {task_name}")
        ax.legend(loc="lower right")
        ax.grid(axis="y", linestyle=":", alpha=0.6)

    plt.tight_layout()
    plt.savefig(os.path.join(FIGDIR, "closed_loop_rir_comparison.png"), dpi=200)
    plt.close()

    print("\nSaved result data to result/closed_loop_plasticity/phase2_single_task.npz")
    print("Saved summary figure to docs/figures/closed_loop_rir_comparison.png")

    return results


if __name__ == "__main__":
    run_phase2_suite()
