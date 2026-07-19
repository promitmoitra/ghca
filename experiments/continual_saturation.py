"""Track 3c / P5 — a *fixed* conjunction basis must saturate.

The bridge (E9 ↔ 3c) showed a learned (stimulus × context) conjunction basis
resolves interference for T=2 tasks, and P4 showed the same for a frozen-random
context basis (avg 0.66) — but both bases are **fixed-size**. Any fixed basis has
finite capacity: as tasks accumulate, the K_stim × T conjunction combinations it
must tile eventually exceed what n_h hidden units can keep separable, and
interference must return. This is the honest bound on the bridge's claim, and the
mechanistic prediction of "capacity, not credit": when capacity is *fixed* the
demand eventually overruns it; when capacity *grows with demand* (per-task heads)
it never does.

Protocol (P4's substrate, generalised past T=2):
  - K_stim stimuli, A=2 actions, T sequential tasks. Each task is a distinct
    BALANCED DICHOTOMY of the stimuli (half -> action 0, half -> action 1), drawn
    once from a fixed task seed so every network seed sees the SAME nested task
    sequence (tasks[:T] for T=4 is the prefix of the T=8 sequence). Binary actions
    keep each task easily learnable (chance 0.5), so backward transfer has full
    dynamic range and the saturation signal is not masked by a many-way-readout
    floor. Conflicting labels for the same stimulus across tasks are what force
    interference at a shared head. Frozen S->H representation, correlational credit
    (the P1/P2 eligibility trace) throughout.
  - Sweep T. Three modes, identical everywhere else:
      shared   — stimulus only, one head (no conjunction; the interference floor)
      context  — + one task-context input channel per task (a fixed-size
                 (stimulus × context) basis over the SAME n_h hidden units)
      per-task — separate H->M head per task (capacity that GROWS with T; control)

Prediction: `context` defers interference at small T but its backward transfer
degrades toward the `shared` floor as T grows past the basis's tiling capacity;
`per-task` stays flat near zero for all T. The crossover locates the capacity of
the fixed conjunction basis.
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

N = int(os.environ.get("STATS_N", "20"))
N_TRAIN = int(os.environ.get("STATS_NTRAIN", "500"))
N_H = int(os.environ.get("STATS_NH", "50"))          # the FIXED conjunction basis
K_STIM = int(os.environ.get("STATS_KSTIM", "4"))     # stimuli (dichotomised per task)
A = 2                                                # binary action; chance 0.5
TASK_SEED = int(os.environ.get("STATS_TASKSEED", "12345"))
# KS=4 admits C(4,2)=6 distinct balanced dichotomies, so T is capped at 6.
T_LIST = [int(t) for t in os.environ.get("STATS_TLIST", "2,3,4,5,6").split(",")]
CHANCE = 1.0 / A
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def make_tasks(n_tasks, k_stim=K_STIM, seed=TASK_SEED):
    """`n_tasks` distinct balanced dichotomies of `k_stim` stimuli, generated from
    a fixed seed so the sequence is shared across network seeds and NESTED across T
    (the first k of any longer sequence are identical)."""
    from math import comb
    n_avail = comb(k_stim, k_stim // 2)
    if n_tasks > n_avail:
        raise ValueError(f"only {n_avail} balanced dichotomies exist for k_stim="
                         f"{k_stim}; requested {n_tasks}")
    rng = np.random.default_rng(seed)
    half = k_stim // 2
    seen, tasks = set(), []
    while len(tasks) < n_tasks:
        t = np.zeros(k_stim, dtype=int)
        t[rng.choice(k_stim, half, replace=False)] = 1
        key = tuple(t)
        if key in seen:
            continue
        seen.add(key); tasks.append(t)
    return tasks


def make(seed, mode, n_tasks):
    context = (mode == "context")
    K = K_STIM + (n_tasks if context else 0)         # extra channels = task-context inputs
    W, plastic, roles = layered_graph(K=K, A=A, n_h=N_H, seed=seed,
                                      w_hm=e1.CFG["w_hm"], w_hh=e1.CFG["w_hh"])
    # frozen representation: freeze S->H, only the H->M head adapts (as P2/P4)
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


def evaluate(net, tasks, task, rng, context, reps=20):
    tgt = tasks[task]
    hits = tot = 0
    for _ in range(reps):                         # deterministic coverage of all stimuli
        for x in range(K_STIM):
            hits += trial(net, x, task, int(tgt[x]), rng, context, learn=False)
            tot += 1
    return hits / tot


def motor_plastic_idx(net):
    P = net.plastic.copy()
    mask = np.zeros_like(P)
    mask[np.concatenate(net.roles["motor"]), :] = True
    return np.where(P & mask)


def _eval_head(net, tasks, j, rng, context, mode, heads, pidx):
    if mode == "per-task":
        net.W[pidx] = heads[j]
    return evaluate(net, tasks, j, rng, context)


def run_sequence(seed, mode, n_tasks):
    context = (mode == "context")
    net = make(seed, mode, n_tasks)
    tasks = make_tasks(n_tasks)
    rng = np.random.default_rng(seed + 7)
    R = np.zeros((n_tasks, n_tasks)); just = np.zeros(n_tasks)
    heads = pidx = None
    if mode == "per-task":
        pidx = motor_plastic_idx(net)
        heads = [net.W[pidx].copy() for _ in range(n_tasks)]
    for k in range(n_tasks):
        tgt = tasks[k]
        if mode == "per-task":
            net.W[pidx] = heads[k]
        for _ in range(N_TRAIN):
            x = int(rng.integers(K_STIM))
            trial(net, x, k, int(tgt[x]), rng, context, learn=True)
        if mode == "per-task":
            heads[k] = net.W[pidx].copy()
        just[k] = _eval_head(net, tasks, k, rng, context, mode, heads, pidx)
        for j in range(n_tasks):
            R[k, j] = _eval_head(net, tasks, j, rng, context, mode, heads, pidx)
    return R, just


def cl_metrics(R, just):
    T = R.shape[0]
    avg = float(R[-1].mean())
    bwt = float(np.mean([R[-1, k] - just[k] for k in range(T - 1)]))
    fwt = float(np.mean([R[k - 1, k] - CHANCE for k in range(1, T)]))
    return avg, bwt, fwt


def main():
    modes = ["shared", "context", "per-task"]
    out = {m: {"T": [], "avg_acc": [], "bwt": [], "Rmean": {}} for m in modes}
    for mode in modes:
        for T in T_LIST:
            avg = np.zeros(N); bwt = np.zeros(N)
            Rsum = np.zeros((T, T))
            for s in range(N):
                R, jt = run_sequence(s, mode, T)
                a, b, _ = cl_metrics(R, jt)
                avg[s], bwt[s] = a, b
                Rsum += R
            out[mode]["T"].append(T)
            out[mode]["avg_acc"].append(avg)
            out[mode]["bwt"].append(bwt)
            out[mode]["Rmean"][T] = Rsum / N
            am, alo, ahi = st.bootstrap_ci(avg)
            bm, blo, bhi = st.bootstrap_ci(bwt)
            print(f"{mode:9s} T={T}: avg={am:.3f}[{alo:.3f},{ahi:.3f}] "
                  f"bwt={bm:+.3f}[{blo:+.3f},{bhi:+.3f}]", flush=True)
    # save
    save = {"n": N, "chance": CHANCE, "n_h": N_H, "k_stim": K_STIM, "A": A,
            "T_list": np.array(T_LIST)}
    for mode in modes:
        save[f"{mode}_avg"] = np.array(out[mode]["avg_acc"])   # (len(T_LIST), N)
        save[f"{mode}_bwt"] = np.array(out[mode]["bwt"])
    np.savez(os.path.join(OUT, "continual_saturation.npz"), **save)
    summ = {"n": N, "chance": CHANCE, "n_h": N_H, "A": A, "k_stim": K_STIM,
            "rows": {}}
    for mode in modes:
        for i, T in enumerate(T_LIST):
            am, alo, ahi = st.bootstrap_ci(out[mode]["avg_acc"][i])
            bm, blo, bhi = st.bootstrap_ci(out[mode]["bwt"][i])
            summ["rows"][f"{mode}:T{T}"] = {
                "avg_acc": [am, alo, ahi], "bwt": [bm, blo, bhi]}
    with open(os.path.join(OUT, "continual_saturation.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote continual_saturation.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
