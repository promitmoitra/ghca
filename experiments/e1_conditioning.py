"""
E1 - Stimulus->response conditioning (identity mapping, no delay).

Tests whether a strict scalar reward can carve a stimulus->action mapping into
the GH network substrate, and dissociates the two plasticity lines:
    Line A (conduction weights) should learn the identity mapping;
    Line B (timescales) should stay near chance, since identity routing needs
    spatial credit assignment, not timing.
See docs/learning_experiments.md, experiment E1.

Protocol (per trial): settle -> read critic baseline V -> present cue for
stimulus x (population code on sensory channel x) -> response window ->
action = argmax motor-channel activity -> reward r = 1[action == x] ->
broadcast delta = r - V to plasticity and critic.

Outputs
-------
docs/figures/e1_learning_curves.png : accuracy vs trial for lines A / B / AB
docs/figures/e1_summary.png         : final accuracy and critic correlation
result/e1/e1_data.npz               : per-trial rewards and baselines
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_learn import layered_graph, GHLearner

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e1")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# operating point inherited from E0 / tuning (see e1_results.md)
CFG = dict(act=3, pas=5, theta=4.0, p_s=3e-3, eta_w=0.05, eta_tau=0.12,
           w_hm=0.6, w_hh=0.25)
SETTLE, CUE, WWIN = 20, 3, 6
N_TRIALS, N_SEEDS = 500, 6


def make(line, seed):
    W, plastic, roles = layered_graph(seed=seed, w_hm=CFG["w_hm"], w_hh=CFG["w_hh"])
    return GHLearner(W, plastic, roles, line=line, act=CFG["act"], pas=CFG["pas"],
                     theta=CFG["theta"], p_s=CFG["p_s"], eta_w=CFG["eta_w"],
                     eta_tau=CFG["eta_tau"], seed=seed + 100)


def trial(net, x, rng):
    net.reset_traces()
    for _ in range(SETTLE):
        net.step_learn(None)
    feats = net.features()
    V = net.value()
    for _ in range(CUE):
        net.step_learn(net.sensory_drive(x))
    sc = np.zeros(net.roles["A"])
    for _ in range(WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    r = 1.0 if action == x else 0.0
    delta = r - V
    net.learn(delta)
    net.update_critic(delta, feats)
    return r, V


def run_condition(line, seed):
    net = make(line, seed)
    rng = np.random.default_rng(seed + 7)
    R = np.zeros(N_TRIALS)
    Vb = np.zeros(N_TRIALS)
    for t in range(N_TRIALS):
        R[t], Vb[t] = trial(net, int(rng.integers(net.roles["K"])), rng)
    return R, Vb


def moving_avg(x, w=25):
    k = np.ones(w) / w
    return np.convolve(x, k, mode="valid")


def main():
    lines = ["A", "B", "AB"]
    rewards = {L: np.zeros((N_SEEDS, N_TRIALS)) for L in lines}
    baselines = {L: np.zeros((N_SEEDS, N_TRIALS)) for L in lines}
    for L in lines:
        for s in range(N_SEEDS):
            rewards[L][s], baselines[L][s] = run_condition(L, s)
        print(f"line {L}: final acc = {rewards[L][:, -120:].mean():.3f}")

    # learning curves
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"A": "crimson", "B": "steelblue", "AB": "seagreen"}
    for L in lines:
        curves = np.array([moving_avg(rewards[L][s]) for s in range(N_SEEDS)])
        m, sem = curves.mean(0), curves.std(0) / np.sqrt(N_SEEDS)
        x = np.arange(len(m))
        ax.plot(x, m, color=colors[L], label=f"Line {L}")
        ax.fill_between(x, m - sem, m + sem, color=colors[L], alpha=0.2)
    ax.axhline(0.5, ls="--", color="k", alpha=0.5, label="chance")
    ax.set_xlabel("trial (25-trial moving average)")
    ax.set_ylabel("accuracy")
    ax.set_ylim(0, 1)
    ax.set_title("E1 - conditioning: Line A learns, Line B stays at chance")
    ax.legend()
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "e1_learning_curves.png")
    fig.savefig(p1, dpi=110)
    print("wrote", p1)

    # summary: final accuracy + critic correlation
    fig, axes = plt.subplots(1, 2, figsize=(11, 4.5))
    fin = {L: rewards[L][:, -120:].mean(1) for L in lines}
    axes[0].bar(lines, [fin[L].mean() for L in lines],
                yerr=[fin[L].std() / np.sqrt(N_SEEDS) for L in lines],
                color=[colors[L] for L in lines], capsize=4)
    axes[0].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[0].set_ylabel("final accuracy (last 120 trials)")
    axes[0].set_ylim(0, 1)
    axes[0].set_title("Final accuracy by plasticity line")

    # critic: corr(baseline V, realised reward) over last half, Line A
    corr = []
    for s in range(N_SEEDS):
        half = slice(N_TRIALS // 2, None)
        v, r = baselines["A"][s][half], rewards["A"][s][half]
        corr.append(np.corrcoef(v, r)[0, 1] if v.std() > 0 else 0.0)
    axes[1].bar(range(N_SEEDS), corr, color="crimson")
    axes[1].axhline(0, color="k")
    axes[1].set_xlabel("seed")
    axes[1].set_ylabel("corr(critic V, reward)")
    axes[1].set_title(f"Critic tracks reward (Line A), mean r={np.mean(corr):.2f}")
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "e1_summary.png")
    fig.savefig(p2, dpi=110)
    print("wrote", p2)

    np.savez(os.path.join(DATADIR, "e1_data.npz"),
             **{f"reward_{L}": rewards[L] for L in lines},
             **{f"baseline_{L}": baselines[L] for L in lines},
             critic_corr=np.array(corr))
    print("wrote", os.path.join(DATADIR, "e1_data.npz"))
    print("\nfinal accuracy: " + ", ".join(
        f"{L}={rewards[L][:, -120:].mean():.2f}" for L in lines))


if __name__ == "__main__":
    main()
