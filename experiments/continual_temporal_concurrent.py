"""Track 3e.3 — concurrent co-adaptation: grow τ *while* the readout learns.

3d-emergent and E9 both use a **phase split**: self-organise the representation
(reward-free), freeze it, then let reward carve the readout. The oldest deferred item
across the programme (E9 caveats, 1c) is the *concurrent* case — representation and
readout plastic at the same time, the strongest form of "one homogeneous machine
learning end-to-end". This tests it on the temporal axis: does growing the timescale
basis *concurrently* with the reward readout reach the phase-split performance, or do
the two plastic loops interfere (early readout learning on a not-yet-tiled basis
poisoning the result)?

Three arms, per-task heads, on the 3d delay-keyed tasks (n=20):
  * **homogeneous** — τ fixed at the homogeneous value, readout only (the floor).
  * **phase-split** — τ grown reward-free first (3d-emergent), frozen, then readout
    learns (the reference).
  * **concurrent** — τ starts homogeneous and is updated by the input-timing rule *on
    every reward trial*, while the reward readout learns simultaneously.

The τ rule's teaching signal (the probe delay) is present on every task trial, so no
separate self-organisation phase is needed — the same trial drives both the reward
readout (Line A) and the competitive τ update.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import continual_temporal_saturation as B      # make, tasks, trial constants  # noqa: E402
import continual_temporal_emergent as B2        # grow_emergent_tau (phase-split arm)  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "20"))
N_TRAIN = int(os.environ.get("STATS_NTRAIN", "500"))
T_LIST = [int(t) for t in os.environ.get("STATS_TLIST", "2,3,4,5,6").split(",")]
ETA_TAU, WIN_FRAC, CONSC_C, CONSC_RATE = 0.20, 0.25, 6.0, 0.05
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def concurrent_trial(net, dbin, target, rng, p_win, plastic_tau):
    """One reward trial that also (optionally) applies the competitive τ update."""
    hid = net.roles["hidden"]
    net.reset_traces()
    net.phi[:] = 0
    net.t = 0
    for _ in range(B.SETTLE):
        net.step_learn(None)
    feats = net.features(); V = net.value()
    drive = net.sensory_drive(0)
    for _ in range(B.CUE):
        net.step_learn(drive)
    for _ in range(B.DELAYS[dbin]):
        net.step_learn(None)
    sc = np.zeros(B.A); acc = np.zeros(len(hid))
    for k in range(B.PROBE + B.WWIN):
        net.step_learn(drive if k < B.PROBE else None)
        if k >= 1:
            sc += net.motor_scores()
        acc += net.active_mask()[hid].astype(float)      # hidden response for the τ rule
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(B.A))) if sc.sum() > 0 else -1
    r = 1.0 if action == target else 0.0
    delta = r - V
    net.learn(delta); net.update_critic(delta, feats)     # Line A reward learning
    if plastic_tau:                                        # concurrent input-timing τ update
        resp = acc / (B.PROBE + B.WWIN)
        elig = resp > 0
        score = resp - CONSC_C * (p_win - WIN_FRAC)
        score[~elig] = -1e9
        n_win = max(1, int(WIN_FRAC * len(hid)))
        win = np.argsort(-score)[:n_win]; win = win[elig[win]]
        won = np.zeros(len(hid), bool); won[win] = True
        p_win += CONSC_RATE * (won.astype(float) - p_win)
        d = B.DELAYS[dbin]
        net.tau[hid[win]] = np.clip(net.tau[hid[win]] + ETA_TAU * (d - net.tau[hid[win]]), 2, 40)
        net.tau_base = net.tau.copy()
    return r


_GROWN = {}                                                # cache: seed -> grown τ (phase-split)


def run_sequence(seed, mode, n_tasks):
    net = B.make(seed, "homog")                            # homogeneous τ start
    hid = net.roles["hidden"]
    if mode == "phase-split":
        if seed not in _GROWN:
            _GROWN[seed] = B2.grow_emergent_tau(seed)      # grow reward-free (cached per seed)
        net.tau[hid] = _GROWN[seed]                        # then freeze
        net.tau_base = net.tau.copy()
    plastic = (mode == "concurrent")
    tasks = B.make_tasks(n_tasks)
    rng = np.random.default_rng(seed + 7)
    p_win = np.full(len(hid), WIN_FRAC)
    pidx = B.motor_plastic_idx(net)                        # per-task heads
    heads = [net.W[pidx].copy() for _ in range(n_tasks)]
    R = np.zeros((n_tasks, n_tasks)); just = np.zeros(n_tasks)
    for k in range(n_tasks):
        tgt = tasks[k]
        net.W[pidx] = heads[k]
        for _ in range(N_TRAIN):
            d = int(rng.integers(B.K_DELAY))
            concurrent_trial(net, d, int(tgt[d]), rng, p_win, plastic)
        heads[k] = net.W[pidx].copy()
        just[k] = _eval(net, tasks, k, rng, heads, pidx)
        for j in range(n_tasks):
            R[k, j] = _eval(net, tasks, j, rng, heads, pidx)
    return R, just


def _eval(net, tasks, j, rng, heads, pidx):
    net.W[pidx] = heads[j]
    return B.evaluate(net, tasks, j, rng)


def main():
    modes = ["homogeneous", "phase-split", "concurrent"]
    out = {m: {"avg": []} for m in modes}
    for mode in modes:
        for T in T_LIST:
            avg = np.zeros(N)
            for s in range(N):
                R, jt = run_sequence(s, mode, T)
                avg[s], _ = B.cl_metrics(R, jt)
            out[mode]["avg"].append(avg)
            m, lo, hi = st.bootstrap_ci(avg)
            print(f"{mode:12s} T={T}: avg={m:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
    save = {"n": N, "T_list": np.array(T_LIST)}
    for mode in modes:
        save[f"{mode}_avg"] = np.array(out[mode]["avg"])
    np.savez(os.path.join(OUT, "continual_temporal_concurrent.npz"), **save)
    summ = {"n": N, "rows": {}}
    for mode in modes:
        for i, T in enumerate(T_LIST):
            m, lo, hi = st.bootstrap_ci(out[mode]["avg"][i])
            summ["rows"][f"{mode}:T{T}"] = {"avg": [m, lo, hi]}
    with open(os.path.join(OUT, "continual_temporal_concurrent.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote continual_temporal_concurrent.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
