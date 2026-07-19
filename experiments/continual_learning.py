"""Track 3c / P1 — continual learning on one shared substrate (harness + baselines).

See docs/continual_learning_plan.md. v1 uses a remapping task triple on ONE shared
E1-style substrate (three conflicting stimulus->action mappings, K=A=3) rather than
the heterogeneous E1/E2/E5 — this isolates catastrophic interference from
architecture-swapping and reuses E1's validated conditioning machinery. Tasks are
cyclic relabellings that maximally conflict on the shared readout:

    T0: action = stim            T1: action = stim+1 (mod 3)   T2: action = stim+2

P1 establishes the harness, the continual-learning metrics, and the two "what
adapts" baselines:
  - frozen : the S->H representation is fixed (reservoir); only the shared H->M
             head adapts. Reservoir theory predicts little interference.
  - plastic: S->H and H->M both adapt (the shared representation is rewritten by
             each task) — the hard case where interference is expected.

Both use the current correlational eligibility-trace credit (Line A). The causal
do(theta) credit rule is P2.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_learn import layered_graph, GHLearner  # noqa: E402
import e1_conditioning as e1  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "30"))
N_TRAIN = int(os.environ.get("STATS_NTRAIN", "500"))   # trials per task
K = A = 3
CHANCE = 1.0 / A
TASKS = [(lambda x: x), (lambda x: (x + 1) % A), (lambda x: (x + 2) % A)]
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def make(seed, regime):
    W, plastic, roles = layered_graph(K=K, A=A, seed=seed,
                                      w_hm=e1.CFG["w_hm"], w_hh=e1.CFG["w_hh"])
    if regime == "frozen":
        # freeze the S->H representation; keep only the H->M head plastic
        h0, h1 = roles["hidden"][0], roles["hidden"][-1] + 1
        plastic = plastic.copy()
        plastic[h0:h1, :] = False          # no S->H updates (rows = hidden nodes)
    net = GHLearner(W, plastic, roles, line="A", act=e1.CFG["act"], pas=e1.CFG["pas"],
                    theta=e1.CFG["theta"], p_s=e1.CFG["p_s"], eta_w=e1.CFG["eta_w"],
                    eta_tau=e1.CFG["eta_tau"], seed=seed + 100)
    return net


def trial(net, x, target, rng, learn=True):
    """One conditioning trial (mirrors e1.trial) with reward vs an explicit target."""
    net.reset_traces()
    for _ in range(e1.SETTLE):
        net.step_learn(None)
    feats = net.features()
    V = net.value()
    for _ in range(e1.CUE):
        net.step_learn(net.sensory_drive(x))
    sc = np.zeros(net.roles["A"])
    for _ in range(e1.WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(len(sc)))) if sc.sum() > 0 else -1
    r = 1.0 if action == target else 0.0
    if learn:
        delta = r - V
        net.learn(delta)
        net.update_critic(delta, feats)
    return r


def evaluate(net, task_idx, rng, n=120):
    m = TASKS[task_idx]
    hits = 0
    for _ in range(n):
        x = int(rng.integers(K))
        hits += trial(net, x, m(x), rng, learn=False)
    return hits / n


def run_sequence(seed, regime):
    """Train T0->T1->T2 on one net; return accuracy matrix R[after_task, eval_task]."""
    net = make(seed, regime)
    rng = np.random.default_rng(seed + 7)
    T = len(TASKS)
    R = np.zeros((T, T))
    just_trained = np.zeros(T)              # acc on task k right after training it
    for k in range(T):
        m = TASKS[k]
        for _ in range(N_TRAIN):
            x = int(rng.integers(K))
            trial(net, x, m(x), rng, learn=True)
        just_trained[k] = evaluate(net, k, rng)
        for j in range(T):
            R[k, j] = evaluate(net, j, rng)
    return R, just_trained


def cl_metrics(R, just_trained):
    T = R.shape[0]
    avg_acc = float(R[-1].mean())                                   # after full sequence
    # backward transfer: change on task k (k<T-1) from its just-trained level to the end
    bwt = float(np.mean([R[-1, k] - just_trained[k] for k in range(T - 1)]))
    # forward transfer: acc on task k before training it (after training k-1), vs chance
    fwt = float(np.mean([R[k - 1, k] - CHANCE for k in range(1, T)]))
    return avg_acc, bwt, fwt


def main():
    regimes = ["frozen", "plastic"]
    out = {}
    for regime in regimes:
        avg = np.zeros(N); bwt = np.zeros(N); fwt = np.zeros(N)
        Rsum = np.zeros((len(TASKS), len(TASKS)))
        for s in range(N):
            R, jt = run_sequence(s, regime)
            a, b, f = cl_metrics(R, jt)
            avg[s], bwt[s], fwt[s] = a, b, f
            Rsum += R
        out[regime] = {"avg_acc": avg, "bwt": bwt, "fwt": fwt, "Rmean": Rsum / N}
        print(f"\n=== {regime} (n={N}, chance={CHANCE:.3f}) ===", flush=True)
        for lab, arr in [("avg_acc", avg), ("backward_transfer", bwt), ("forward_transfer", fwt)]:
            print("  " + st.fmt_row(st.summarise(lab, arr)), flush=True)
        print("  Rmean[after_task, eval_task]:\n" +
              "\n".join("    " + " ".join(f"{v:.2f}" for v in row) for row in out[regime]["Rmean"]),
              flush=True)

    np.savez(os.path.join(OUT, "continual_p1.npz"), n=N, chance=CHANCE,
             **{f"{r}_{k}": out[r][k] for r in regimes for k in ("avg_acc", "bwt", "fwt")},
             **{f"{r}_Rmean": out[r]["Rmean"] for r in regimes})
    with open(os.path.join(OUT, "continual_p1.json"), "w") as f:
        json.dump({"n": N, "chance": CHANCE,
                   "rows": {f"{r}:{k}": st.summarise(k, out[r][k])
                            for r in regimes for k in ("avg_acc", "bwt", "fwt")}},
                  f, indent=2, default=float)
    print("\nwrote continual_p1.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
