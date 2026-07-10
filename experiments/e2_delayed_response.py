"""
E2 - Delayed response (working memory).

Tests memory as a persistent reentrant loop that bridges a stimulus-free delay,
and the predicted inversion of the E1 dissociation: here Line B (timescale
plasticity) is critical, because loop persistence is set by the local timescale
tau, while Line A (conduction weights) cannot extend memory beyond the loop's
intrinsic lifetime. See docs/learning_experiments.md, experiment E2.

Substrate: K stimulus-specific directed rings (length L) form the recurrent
hidden medium; a cue ignites ring x; each motor channel reads its own ring over
a full loop period (so readout is phase-insensitive). A directed ring sustains a
rotating pulse indefinitely when tau < L and dies in ~L steps when tau >= L, so
tau is the memory-duration control variable.

Two results are produced:
  1. MECHANISM (no learning): accuracy vs delay for fixed uniform tau, showing
     memory duration is tau-controlled.
  2. LEARNING: lines A / B / AB with a strict reward at the end of the delay;
     Line B uses a shared regional timescale learned by reward-gated
     perturbation.

Outputs
-------
docs/figures/e2_mechanism.png       : accuracy vs delay for fixed tau values
docs/figures/e2_learning.png        : retention curve (accuracy vs delay) by line
result/e2/e2_data.npz               : swept arrays
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network
from ghca_learn import GHLearner

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e2")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

L = 24                 # ring length (loop transit time)
ACT = 6
IGNITE, FANIN, N_M = 6, 18, 12
CUE, WWIN = 4, L       # readout window = one loop period (phase-insensitive)
DELAYS_EVAL = [0, 20, 50, 100, 200]


def two_ring(w_hm=0.6, theta_m=1.0, jitter=0.15, seed=0):
    """Two stimulus-specific directed rings + per-ring motor readout."""
    rng = np.random.default_rng(seed)
    K = A = 2
    N = K * L + A * N_M
    ring = [np.arange(x * L, (x + 1) * L) for x in range(K)]
    m0 = K * L
    motor = [np.arange(m0 + c * N_M, m0 + (c + 1) * N_M) for c in range(A)]
    W = np.zeros((N, N))
    plastic = np.zeros((N, N), dtype=bool)
    theta = np.zeros(N)
    for x in range(K):
        for i in range(L):
            W[ring[x][i], ring[x][(i - 1) % L]] = 1.0     # directed ring
        theta[ring[x]] = 1.0
    for c in range(A):
        for m in motor[c]:
            src = rng.choice(ring[c], FANIN, replace=False)
            W[m, src] = w_hm * (1 + jitter * rng.standard_normal(FANIN))
            plastic[m, src] = True
        theta[motor[c]] = theta_m
    hidden = np.arange(K * L)
    roles = {"sensory": [ring[x][:IGNITE] for x in range(K)], "motor": motor,
             "hidden": hidden, "K": K, "A": A, "N": N, "n_m": N_M, "L": L}
    tau_mask = np.zeros(N, dtype=bool)
    tau_mask[hidden] = True
    return np.clip(W, 0, None), plastic, roles, theta, tau_mask


# --------------------------------------------------------------------------
# 1. Mechanism sweep (no learning): accuracy vs delay for fixed uniform tau
# --------------------------------------------------------------------------
def mechanism(taus=(16, 20, 22, 24, 28), delays=(0, 10, 20, 40, 80, 150), ntr=40):
    W, plastic, roles, theta, _ = two_ring(seed=0)
    ring, motor = None, roles["motor"]
    acc = np.zeros((len(taus), len(delays)))
    for i, tau in enumerate(taus):
        for j, D in enumerate(delays):
            correct = 0
            for tr in range(ntr):
                Wt, _, r2, th, _ = two_ring(seed=tr)
                net = Network(Wt, act=ACT, pas=int(tau - ACT), theta=th, p_s=0.0, seed=tr)
                x = tr % 2
                net.phi[:] = 0
                for _ in range(CUE):
                    d = np.zeros(net.N, bool); d[r2["sensory"][x]] = True
                    net.step(d)
                for _ in range(D):
                    net.step(None)
                sc = np.zeros(2)
                for _ in range(WWIN):
                    net.step(None)
                    am = net.active_mask()
                    sc += [am[r2["motor"][k]].sum() for k in range(2)]
                if sc.sum() > 0 and np.argmax(sc) == x:
                    correct += 1
            acc[i, j] = correct / ntr
        print(f"  mechanism tau={tau} done")
    return np.array(taus), np.array(delays), acc


# --------------------------------------------------------------------------
# 2. Learning: lines A / B / AB, retention curve
# --------------------------------------------------------------------------
def make(line, seed, tau0=26, eta_tau=0.4, tau_sigma=3.0):
    W, plastic, roles, theta, tmask = two_ring(seed=seed)
    net = GHLearner(W, plastic, roles, line=line, act=ACT, pas=int(tau0 - ACT),
                    theta=theta, p_s=0.0, eta_w=0.05, eta_tau=eta_tau,
                    tau_min=8, tau_max=40, tau_sigma=tau_sigma, tau_mask=tmask,
                    tau_shared=True, seed=seed + 50)
    return net, roles


def trial(net, roles, x, D, rng, learn=True):
    net.reset_traces()
    net.phi[:] = 0
    if learn:
        net.perturb_tau()
    for _ in range(CUE):
        net.step_learn(net.sensory_drive(x))
    for _ in range(D):
        net.step_learn(None)
    V = net.value()
    feats = net.features()
    sc = np.zeros(roles["A"])
    for _ in range(WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    r = 1.0 if action == x else 0.0
    if learn:
        net.learn(r - V)
        net.update_critic(r - V, feats)
    return r


def run_condition(line, seed, ntr=1800, Dlo=15, Dhi=70):
    net, roles = make(line, seed)
    rng = np.random.default_rng(seed + 3)
    for _ in range(ntr):
        trial(net, roles, int(rng.integers(2)), int(rng.integers(Dlo, Dhi)), rng, learn=True)
    net.tau[roles["hidden"]] = net.tau_scalar
    net._eps_s = 0.0
    ret = np.array([np.mean([trial(net, roles, int(rng.integers(2)), D, rng, learn=False)
                             for _ in range(60)]) for D in DELAYS_EVAL])
    return ret, net.tau_scalar


def main():
    print("E2: mechanism sweep (fixed tau, no learning) ...")
    taus, delays, acc = mechanism()

    print("E2: learning (lines A / B / AB) ...")
    lines = ["A", "B", "AB"]
    n_seeds = 5
    ret = {L_: np.zeros((n_seeds, len(DELAYS_EVAL))) for L_ in lines}
    tau_learned = {L_: np.zeros(n_seeds) for L_ in lines}
    for L_ in lines:
        for s in range(n_seeds):
            ret[L_][s], tau_learned[L_][s] = run_condition(L_, s)
        print(f"  line {L_}: mean tau={tau_learned[L_].mean():.1f} (L={L}), "
              f"retention={ret[L_].mean(0).round(2)}")

    # mechanism figure
    fig, ax = plt.subplots(figsize=(7, 5))
    for i, tau in enumerate(taus):
        ax.plot(delays, acc[i], "o-", label=f"tau={tau}")
    ax.axvline(0, color="k", alpha=0)
    ax.set_xlabel("delay D (steps)")
    ax.set_ylabel("accuracy")
    ax.set_ylim(-0.03, 1.05)
    ax.set_title(f"E2 mechanism: memory duration is tau-controlled (ring L={L})")
    ax.legend(title="fixed uniform tau")
    fig.tight_layout()
    p1 = os.path.join(FIGDIR, "e2_mechanism.png")
    fig.savefig(p1, dpi=110)
    print("wrote", p1)

    # learning retention figure
    fig, ax = plt.subplots(figsize=(7, 5))
    colors = {"A": "crimson", "B": "steelblue", "AB": "seagreen"}
    for L_ in lines:
        m = ret[L_].mean(0)
        sem = ret[L_].std(0) / np.sqrt(n_seeds)
        ax.plot(DELAYS_EVAL, m, "o-", color=colors[L_],
                label=f"Line {L_} (tau*={tau_learned[L_].mean():.0f})")
        ax.fill_between(DELAYS_EVAL, m - sem, m + sem, color=colors[L_], alpha=0.2)
    ax.axhline(0.5, ls="--", color="k", alpha=0.5, label="chance")
    ax.set_xlabel("delay D (steps)")
    ax.set_ylabel("accuracy (retention)")
    ax.set_ylim(0, 1.05)
    ax.set_title("E2 learning: Line B holds memory across delay, Line A cannot")
    ax.legend()
    fig.tight_layout()
    p2 = os.path.join(FIGDIR, "e2_learning.png")
    fig.savefig(p2, dpi=110)
    print("wrote", p2)

    np.savez(os.path.join(DATADIR, "e2_data.npz"),
             mech_taus=taus, mech_delays=delays, mech_acc=acc,
             delays_eval=np.array(DELAYS_EVAL),
             **{f"ret_{L_}": ret[L_] for L_ in lines},
             **{f"tau_{L_}": tau_learned[L_] for L_ in lines})
    print("wrote", os.path.join(DATADIR, "e2_data.npz"))


if __name__ == "__main__":
    main()
