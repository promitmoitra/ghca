"""Track 3c / P4 — capacity, not credit: what actually resolves the interference.

P2/P3 showed no credit rule beats the stability-plasticity frontier — the limit is
representational capacity. This is the positive complement: holding the credit rule
fixed (the correlational eligibility trace that catastrophically forgot in P1/P2),
vary *capacity* and show the interference collapses.

Because the two tasks are anti-correlated (identity vs reversal), raw dimensionality
cannot help — a stimulus-only hidden pattern forces the shared head to emit opposite
actions for the same pattern. The capacity that matters is **task-context /
conjunctive** capacity. Three rungs (all: frozen representation, shared readout
unless noted, correlational credit; K_stim=2 reversal, n=30):

  1. shared   — stimulus only, shared head (the P2 interference baseline)
  2. context  — add a task-context input channel so hidden patterns are
                (stimulus × context) conjunctive; one shared head (E9-flavoured fix)
  3. per-task — separate H->M heads per task (the trivial full-capacity bracket)
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
N_TRAIN = int(os.environ.get("STATS_NTRAIN", "500"))
K_STIM, A, T = 2, 2, 2                       # 2 stimuli, 2 actions, 2 tasks (reversal)
CHANCE = 1.0 / A
TASKS = [lambda x: x, lambda x: (x + 1) % A]  # identity, reversal
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def make(seed, context):
    K = K_STIM + (T if context else 0)       # extra channels = task-context inputs
    W, plastic, roles = layered_graph(K=K, A=A, seed=seed,
                                      w_hm=e1.CFG["w_hm"], w_hh=e1.CFG["w_hh"])
    # frozen representation: freeze S->H, only the H->M head adapts (as P2)
    h0, h1 = roles["hidden"][0], roles["hidden"][-1] + 1
    plastic = plastic.copy(); plastic[h0:h1, :] = False
    net = GHLearner(W, plastic, roles, line="A", act=e1.CFG["act"], pas=e1.CFG["pas"],
                    theta=e1.CFG["theta"], p_s=e1.CFG["p_s"], eta_w=e1.CFG["eta_w"],
                    eta_tau=e1.CFG["eta_tau"], seed=seed + 100)
    return net


def trial(net, x, task, target, rng, context, learn=True):
    net.reset_traces()
    for _ in range(e1.SETTLE):
        net.step_learn(None)
    feats = net.features(); V = net.value()
    drive = net.sensory_drive(x)
    if context:
        drive = drive | net.sensory_drive(K_STIM + task)   # co-active task-context input
    for _ in range(e1.CUE):
        net.step_learn(drive)
    sc = np.zeros(A)
    for _ in range(e1.WWIN):
        net.step_learn(None)
        sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(A))) if sc.sum() > 0 else -1
    r = 1.0 if action == target else 0.0
    if learn:
        delta = r - V
        net.learn(delta)
        net.update_critic(delta, feats)
    return r


def evaluate(net, task, rng, context, n=120):
    m = TASKS[task]
    return sum(trial(net, (x := int(rng.integers(K_STIM))), task, m(x), rng, context,
                     learn=False) for _ in range(n)) / n


def motor_plastic_idx(net):
    P = net.plastic.copy()
    mask = np.zeros_like(P)
    mask[np.concatenate(net.roles["motor"]), :] = True
    return np.where(P & mask)


def run_sequence(seed, mode):
    context = (mode == "context")
    net = make(seed, context)
    rng = np.random.default_rng(seed + 7)
    R = np.zeros((T, T)); just = np.zeros(T)
    heads = None
    if mode == "per-task":                    # per-task copies of the H->M head
        pidx = motor_plastic_idx(net)
        heads = [net.W[pidx].copy() for _ in range(T)]
    for k in range(T):
        m = TASKS[k]
        if mode == "per-task":
            net.W[pidx] = heads[k]
        for _ in range(N_TRAIN):
            x = int(rng.integers(K_STIM))
            trial(net, x, k, m(x), rng, context, learn=True)
        if mode == "per-task":
            heads[k] = net.W[pidx].copy()
        just[k] = _eval_head(net, k, rng, context, mode, heads, pidx if heads else None)
        for j in range(T):
            R[k, j] = _eval_head(net, j, rng, context, mode, heads, pidx if heads else None)
    return R, just


def _eval_head(net, j, rng, context, mode, heads, pidx):
    if mode == "per-task":
        net.W[pidx] = heads[j]                 # evaluate with task j's own head
    return evaluate(net, j, rng, context)


def cl_metrics(R, just):
    avg = float(R[-1].mean())
    bwt = float(np.mean([R[-1, k] - just[k] for k in range(T - 1)]))
    fwt = float(np.mean([R[k - 1, k] - CHANCE for k in range(1, T)]))
    return avg, bwt, fwt


def main():
    modes = ["shared", "context", "per-task"]
    out = {}
    for mode in modes:
        avg = np.zeros(N); bwt = np.zeros(N); fwt = np.zeros(N)
        Rsum = np.zeros((T, T))
        for s in range(N):
            R, jt = run_sequence(s, mode)
            a, b, f = cl_metrics(R, jt)
            avg[s], bwt[s], fwt[s] = a, b, f
            Rsum += R
        out[mode] = {"avg_acc": avg, "bwt": bwt, "fwt": fwt, "Rmean": Rsum / N}
        print(f"\n=== {mode} (n={N}, chance={CHANCE:.3f}) ===", flush=True)
        for lab, arr in [("avg_acc", avg), ("backward_transfer", bwt)]:
            print("  " + st.fmt_row(st.summarise(lab, arr)), flush=True)
        print("  Rmean:\n" + "\n".join("    " + " ".join(f"{v:.2f}" for v in row)
                                       for row in out[mode]["Rmean"]), flush=True)
    np.savez(os.path.join(OUT, "continual_capacity.npz"), n=N, chance=CHANCE,
             **{f"{k}_{m}": out[k][m] for k in out for m in ("avg_acc", "bwt", "fwt")},
             **{f"{k}_Rmean": out[k]["Rmean"] for k in out})
    with open(os.path.join(OUT, "continual_capacity.json"), "w") as f:
        json.dump({"n": N, "chance": CHANCE,
                   "rows": {f"{k}:{m}": st.summarise(m, out[k][m])
                            for k in out for m in ("avg_acc", "bwt")}}, f, indent=2, default=float)
    print("\nwrote continual_capacity.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
