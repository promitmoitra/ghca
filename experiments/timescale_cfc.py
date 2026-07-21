"""Track 3e.2b / fully closes 4a — cross-frequency coupling on the emergent hierarchy.

3e.2 grew a fast/slow timescale hierarchy but its pool had no inter-population pathway,
so the theta–gamma-style **cross-frequency coupling** half of 4a was deferred. This
adds the pathway and tests for phase–amplitude coupling (PAC): does the slow
population's phase modulate the fast population's amplitude?

Mechanism (standard theta–gamma nesting): the slow population modulates the **fast
nodes' excitability** — when the slow population is active it lowers the fast nodes'
firing threshold, so the fast rhythm is expressed only around the slow-active phase.
The two populations and their timescales are *learned* (3e.2's self-organised
hierarchy, reused verbatim); the excitability pathway linking them is *structural*
(added here, not learned) — the honest scope of this result.

Periods are non-commensurate (P_f=5, P_s=23) so the baseline PAC has no
period-alignment artifact. Two conditions:
  * **uncoupled** — fast nodes supra-driven (fire every fast cycle), threshold fixed:
    fast amplitude is uniform across slow phase → PAC ≈ 0 (the control).
  * **coupled** — fast source at the excitability-gated level, threshold lowered by
    slow-population activity: fast fires preferentially at the slow-active phase → PAC.

PAC metric: Tort modulation index (normalised KL of the mean-fast-amplitude-per-
slow-phase-bin distribution from uniform; 0 = no coupling, 1 = all amplitude in one
phase bin).
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ghca_net import Network  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "20"))
NPOOL = 120
GROW_STEPS = int(os.environ.get("STATS_GROW", "10000"))
RUN_STEPS = int(os.environ.get("STATS_RUN", "4000"))
P_F, P_S = 5, 23                          # non-commensurate fast / slow periods
ACT, TMIN, TMAX = 2, 3, 34
ETA_TAU, ETA_W, CONSC = 0.15, 0.03, 0.6
TH_HI, DEPTH, S_REF = 1.35, 0.5, 0.05     # excitability modulation of fast thresholds
N_PHASE = 9
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def grow_hierarchy(seed):
    """The 3e.2 emergent fast/slow τ hierarchy (input-timing rule + conscience)."""
    rng = np.random.default_rng(seed)
    M = NPOOL + 2
    fast, slow = NPOOL, NPOOL + 1
    W = np.zeros((M, M))
    W[:NPOOL, fast] = 1.2 * (1 + 0.2 * rng.standard_normal(NPOOL))
    W[:NPOOL, slow] = 1.2 * (1 + 0.2 * rng.standard_normal(NPOOL))
    net = Network(W, act=ACT, pas=8, theta=np.full(M, 1.0), p_s=0.0, seed=seed)
    net.tau = net.tau.astype(float)
    net.tau[:NPOOL] = rng.uniform(TMIN, TMAX, NPOOL)
    net.tau[fast], net.tau[slow] = P_F, P_S
    last_in = np.full(NPOOL, -1.0)
    for t in range(GROW_STEPS):
        d = np.zeros(M, bool)
        if t % P_F == 0:
            d[fast] = True
        if t % P_S == 0:
            d[slow] = True
        net.step(d)
        ff, sf = bool(d[fast]), bool(d[slow])
        if not (ff or sf):
            continue
        wf, ws = net.W[:NPOOL, fast], net.W[:NPOOL, slow]
        pf = float((wf > ws).mean())
        domf = (wf - CONSC * (pf - 0.5)) >= (ws - CONSC * ((1 - pf) - 0.5))
        for i in range(NPOOL):
            if domf[i]:
                net.W[i, fast] = min(net.W[i, fast] + ETA_W, 2.5)
                net.W[i, slow] = max(net.W[i, slow] - ETA_W, 0.2)
            else:
                net.W[i, slow] = min(net.W[i, slow] + ETA_W, 2.5)
                net.W[i, fast] = max(net.W[i, fast] - ETA_W, 0.2)
            fd = (domf[i] and ff) or ((not domf[i]) and sf)
            if fd:
                if last_in[i] >= 0 and (t - last_in[i]) <= P_S * 1.5:
                    net.tau[i] = np.clip(net.tau[i] + ETA_TAU * ((t - last_in[i]) - net.tau[i]),
                                         TMIN, TMAX)
                last_in[i] = t
    dom = net.W[:NPOOL, fast] > net.W[:NPOOL, slow]
    return net, dom


def tort_mi(phase, amp, nb=N_PHASE):
    bins = np.clip((phase * nb).astype(int), 0, nb - 1)
    m = np.array([amp[bins == b].mean() if (bins == b).any() else 0.0 for b in range(nb)])
    if m.sum() == 0:
        return 0.0, m
    p = m / m.sum()
    p = np.clip(p, 1e-12, None)
    return float(np.sum(p * np.log(p * nb)) / np.log(nb)), m


def run_cfc(net, dom, coupled, seed):
    M = net.N
    fast, slow = NPOOL, NPOOL + 1
    fid, sid = np.where(dom)[0], np.where(~dom)[0]
    net.phi[:] = 0
    net.p_s = 0.0
    if not coupled:
        net.W[fid, fast] = 1.5           # supra-threshold: fires every fast cycle
        net.theta[fid] = 1.0
    else:
        net.W[fid, fast] = 1.0           # gated level; expressed only when threshold drops
    aF = np.zeros(RUN_STEPS)
    s_prev = 0.0
    for t in range(RUN_STEPS):
        if coupled:
            net.theta[fid] = TH_HI - DEPTH * np.clip(s_prev / S_REF, 0, 1)
        d = np.zeros(M, bool)
        if t % P_F == 0:
            d[fast] = True
        if t % P_S == 0:
            d[slow] = True
        net.step(d)
        am = net.active_mask()
        aF[t] = am[fid].mean()
        s_prev = am[sid].mean()
    warm = 800
    sl = slice(warm, RUN_STEPS)
    phase = ((np.arange(RUN_STEPS) % P_S) / P_S)[sl]
    env = np.convolve(aF, np.ones(P_F) / P_F, mode="same")[sl]
    mi, profile = tort_mi(phase, env)
    return mi, profile, float(aF[sl].mean())


def main():
    res = {"uncoupled": {"mi": np.zeros(N), "prof": np.zeros((N, N_PHASE)), "amp": np.zeros(N)},
           "coupled": {"mi": np.zeros(N), "prof": np.zeros((N, N_PHASE)), "amp": np.zeros(N)}}
    for s in range(N):
        net, dom = grow_hierarchy(s)
        base = (net.W.copy(), net.theta.copy(), net.tau.copy())
        for mode in ("uncoupled", "coupled"):
            net.W[:], net.theta[:], net.tau[:] = base[0], base[1], base[2]
            mi, prof, amp = run_cfc(net, dom, mode == "coupled", s)
            res[mode]["mi"][s] = mi
            res[mode]["prof"][s] = prof
            res[mode]["amp"][s] = amp
        print(f"seed {s+1}/{N}: uncoupled MI={res['uncoupled']['mi'][s]:.3f} "
              f"coupled MI={res['coupled']['mi'][s]:.3f}", flush=True)
    for mode in ("uncoupled", "coupled"):
        m, lo, hi = st.bootstrap_ci(res[mode]["mi"])
        print(f"=== {mode}: PAC modulation index = {m:.3f} [{lo:.3f}, {hi:.3f}] "
              f"(mean fast amp {res[mode]['amp'].mean():.3f}) ===", flush=True)
    np.savez(os.path.join(OUT, "timescale_cfc.npz"), n=N, P_F=P_F, P_S=P_S, n_phase=N_PHASE,
             **{f"{mode}_{k}": res[mode][k] for mode in res for k in ("mi", "prof", "amp")})
    summ = {"n": N, "P_F": P_F, "P_S": P_S, "rows": {}}
    for mode in res:
        m, lo, hi = st.bootstrap_ci(res[mode]["mi"])
        summ["rows"][mode] = {"pac_mi": [m, lo, hi], "fast_amp": float(res[mode]["amp"].mean())}
    with open(os.path.join(OUT, "timescale_cfc.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote timescale_cfc.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
