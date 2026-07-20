"""Track 3d (wired arm) — timescale diversity as a continual-learning capacity axis.

3c/P5 showed a fixed *spatial* (stimulus × context) conjunction basis has a finite
continual-learning ceiling. 3d asks whether the substrate's other representational
axis — **timescale** — buys capacity a spatial basis can't, on tasks whose
discriminative variable is *timing* rather than stimulus identity.

This is the **wired** (afforded) arm of 3d's ladder: timescale diversity is
hand-set, not grown (the *emergent* arm waits on Track 4a's input-tracked-τ rule).
It is the cheap go/no-go gate: if a hand-set timescale-diverse basis does not raise
the ceiling here, the emergent hierarchy is not worth building.

Task family (the E3 timed-response regime, stripped to its temporal core): a single
stimulus pulse, then a variable delay drawn from K_DELAY bins, then a probe. The
correct action depends on the *delay bin*, not the stimulus — so the tasks share
stimulus and (trivial) context and differ only in timing. Each task is a balanced
dichotomy over the delay bins (same nested/shared construction as P5). A purely
spatial basis is blind here; only a temporal code can separate the tasks.

Mechanism the wired basis exploits: hidden activity is not self-sustaining, so a
stimulus leaves each fired node **refractory for τ steps**. At delay Δ the set of
*rested* nodes (τ < Δ) is a thermometer code of elapsed time — but only if τ is
diverse. A probe pulse at Δ re-excites exactly the rested nodes, making that code
readable at the motor pool. With homogeneous τ the code collapses to ~1 bit
(everything flips rested at Δ=τ); with graded τ it tiles the delay axis (confirmed:
delay-decode 1.00 vs 0.53 for homogeneous).

Four arms isolate representation from capacity — {homogeneous, graded τ} × {shared
head, per-task heads}, swept over the number of sequential tasks T:
  homog+shared    — no temporal code, one head: the floor.
  graded+shared   — the fixed timescale basis + one head: works at low T, must
                    saturate as tasks accumulate (the 3d analogue of P5's context arm).
  homog+per-task  — unlimited heads but no temporal code: stays at the floor, proving
                    the homogeneous failure is *representational*, not head-conflict.
  graded+per-task — temporal code + capacity that grows with T: the flat control.

Prediction: graded ≫ homogeneous (timescale diversity is the lever); graded+shared
falls toward the floor as T grows (a fixed temporal basis is finite too); the two
per-task arms bracket it (homog low, graded flat-high).
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
N_H = int(os.environ.get("STATS_NH", "50"))
A = 2                                              # binary action; chance 0.5
DELAYS = [2, 6, 10, 16, 24, 34]                    # delay bins (steps after cue)
K_DELAY = len(DELAYS)
TAU_LO, TAU_HI = 4.0, 36.0                         # graded thermometer range
TAU_HOMOG = 8.0
CUE, PROBE, WWIN, SETTLE = e1.CUE, 2, e1.WWIN, e1.SETTLE
TASK_SEED = int(os.environ.get("STATS_TASKSEED", "12345"))
T_LIST = [int(t) for t in os.environ.get("STATS_TLIST", "2,3,4,5,6").split(",")]
CHANCE = 1.0 / A
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def make_tasks(n_tasks, k=K_DELAY, seed=TASK_SEED):
    """`n_tasks` distinct balanced dichotomies of the delay bins, from a fixed seed
    (shared across network seeds, nested across T)."""
    from math import comb
    if n_tasks > comb(k, k // 2):
        raise ValueError(f"only {comb(k, k // 2)} balanced dichotomies for k={k}")
    rng = np.random.default_rng(seed)
    half = k // 2
    seen, tasks = set(), []
    while len(tasks) < n_tasks:
        t = np.zeros(k, dtype=int)
        t[rng.choice(k, half, replace=False)] = 1
        key = tuple(t)
        if key in seen:
            continue
        seen.add(key); tasks.append(t)
    return tasks


def make(seed, basis):
    W, plastic, roles = layered_graph(K=1, A=A, n_h=N_H, seed=seed,
                                      w_hm=e1.CFG["w_hm"], w_hh=e1.CFG["w_hh"])
    h0, h1 = roles["hidden"][0], roles["hidden"][-1] + 1
    plastic = plastic.copy(); plastic[h0:h1, :] = False       # frozen representation
    net = GHLearner(W, plastic, roles, line="A", act=e1.CFG["act"], pas=e1.CFG["pas"],
                    theta=e1.CFG["theta"], p_s=e1.CFG["p_s"], eta_w=e1.CFG["eta_w"],
                    eta_tau=e1.CFG["eta_tau"], tau_max=40, seed=seed + 100)
    hid = roles["hidden"]
    if basis == "homog":
        net.tau[hid] = TAU_HOMOG
    else:                                                     # graded thermometer
        net.tau[hid] = np.linspace(TAU_LO, TAU_HI, len(hid))
    net.tau_base = net.tau.copy()
    return net


def trial(net, dbin, target, rng, learn=True):
    net.reset_traces()
    net.phi[:] = 0                                            # clean, self-contained trial
    net.t = 0
    for _ in range(SETTLE):
        net.step_learn(None)
    feats = net.features(); V = net.value()
    drive = net.sensory_drive(0)
    for _ in range(CUE):
        net.step_learn(drive)                                 # stimulus pulse
    for _ in range(DELAYS[dbin]):
        net.step_learn(None)                                  # variable delay
    sc = np.zeros(A)
    for k in range(PROBE + WWIN):
        net.step_learn(drive if k < PROBE else None)          # probe, then read
        if k >= 1:
            sc += net.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(A))) if sc.sum() > 0 else -1
    r = 1.0 if action == target else 0.0
    if learn:
        delta = r - V
        net.learn(delta)
        net.update_critic(delta, feats)
    return r


def evaluate(net, tasks, task, rng, reps=20):
    tgt = tasks[task]
    hits = tot = 0
    for _ in range(reps):
        for d in range(K_DELAY):
            hits += trial(net, d, int(tgt[d]), rng, learn=False)
            tot += 1
    return hits / tot


def motor_plastic_idx(net):
    P = net.plastic.copy()
    mask = np.zeros_like(P)
    mask[np.concatenate(net.roles["motor"]), :] = True
    return np.where(P & mask)


def _eval_head(net, tasks, j, rng, per_task, heads, pidx):
    if per_task:
        net.W[pidx] = heads[j]
    return evaluate(net, tasks, j, rng)


def run_sequence(seed, basis, head_mode, n_tasks):
    per_task = (head_mode == "per-task")
    net = make(seed, basis)
    tasks = make_tasks(n_tasks)
    rng = np.random.default_rng(seed + 7)
    R = np.zeros((n_tasks, n_tasks)); just = np.zeros(n_tasks)
    heads = pidx = None
    if per_task:
        pidx = motor_plastic_idx(net)
        heads = [net.W[pidx].copy() for _ in range(n_tasks)]
    for kk in range(n_tasks):
        tgt = tasks[kk]
        if per_task:
            net.W[pidx] = heads[kk]
        for _ in range(N_TRAIN):
            d = int(rng.integers(K_DELAY))
            trial(net, d, int(tgt[d]), rng, learn=True)
        if per_task:
            heads[kk] = net.W[pidx].copy()
        just[kk] = _eval_head(net, tasks, kk, rng, per_task, heads, pidx)
        for j in range(n_tasks):
            R[kk, j] = _eval_head(net, tasks, j, rng, per_task, heads, pidx)
    return R, just


def cl_metrics(R, just):
    T = R.shape[0]
    avg = float(R[-1].mean())
    bwt = float(np.mean([R[-1, k] - just[k] for k in range(T - 1)])) if T > 1 else 0.0
    return avg, bwt


ARMS = [("homog", "shared"), ("graded", "shared"),
        ("homog", "per-task"), ("graded", "per-task")]


def main():
    out = {f"{b}+{h}": {"avg": [], "bwt": []} for b, h in ARMS}
    for basis, head in ARMS:
        arm = f"{basis}+{head}"
        for T in T_LIST:
            avg = np.zeros(N); bwt = np.zeros(N)
            for s in range(N):
                R, jt = run_sequence(s, basis, head, T)
                avg[s], bwt[s] = cl_metrics(R, jt)
            out[arm]["avg"].append(avg); out[arm]["bwt"].append(bwt)
            am, alo, ahi = st.bootstrap_ci(avg)
            bm, blo, bhi = st.bootstrap_ci(bwt)
            print(f"{arm:16s} T={T}: avg={am:.3f}[{alo:.3f},{ahi:.3f}] "
                  f"bwt={bm:+.3f}[{blo:+.3f},{bhi:+.3f}]", flush=True)
    save = {"n": N, "chance": CHANCE, "n_h": N_H, "k_delay": K_DELAY,
            "delays": np.array(DELAYS), "T_list": np.array(T_LIST)}
    for arm in out:
        save[f"{arm}_avg"] = np.array(out[arm]["avg"])
        save[f"{arm}_bwt"] = np.array(out[arm]["bwt"])
    np.savez(os.path.join(OUT, "continual_temporal_saturation.npz"), **save)
    summ = {"n": N, "chance": CHANCE, "n_h": N_H, "delays": DELAYS, "rows": {}}
    for arm in out:
        for i, T in enumerate(T_LIST):
            am, alo, ahi = st.bootstrap_ci(out[arm]["avg"][i])
            bm, blo, bhi = st.bootstrap_ci(out[arm]["bwt"][i])
            summ["rows"][f"{arm}:T{T}"] = {"avg": [am, alo, ahi], "bwt": [bm, blo, bhi]}
    with open(os.path.join(OUT, "continual_temporal_saturation.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote continual_temporal_saturation.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
