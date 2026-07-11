"""
E10 (WIP) — diagnostics for an emergent timescale hierarchy.

STATUS: exploratory / negative result. This script is NOT a finished experiment;
it is the evidence behind the finding written up in `docs/e10_notes.md`, kept so
the next iteration does not repeat the dead ends.

Goal (Track 4a of docs/next_steps.md): with per-node timescale tau plastic
(Line B), do learned tau distributions self-organise into a fast/slow HIERARCHY
(a bimodal tau distribution, clusters at the two drive periods) under a
two-rhythm drive?

Finding: the existing Line B resonance rule (tau += eta*(interfire - tau)) CANNOT
build such a hierarchy. It only ratchets tau UPWARD toward MULTIPLES of a drive
period (once tau exceeds the fundamental, the node skips pulses, so its observed
inter-fire interval is a larger multiple, which pulls tau up further). So no
small-tau "fast" population can form from a reasonable init. Three diagnostics
below, most decisive last.

Run: python3 experiments/e10_diagnostics.py
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network

P_F, P_S = 6, 24            # fast / slow drive periods
ACT = 2
TMIN, TMAX = 3, 34


# ---------------------------------------------------------------------------
# Diagnostic 1 — naive resonance on a recurrent pool with a two-rhythm drive.
# Expectation (hoped): tau splits into clusters at P_f and P_s.
# Result: no split at the drive periods; two-rhythm == single-rhythm (the slow
# rhythm is invisible because the fast drive fires every node first).
# ---------------------------------------------------------------------------

def diag1_naive(two_rhythm=True, N=200, steps=4000, eta=0.15, seed=0):
    rng = np.random.default_rng(seed)
    M = N + 2
    fast_src, slow_src = N, N + 1
    W = np.zeros((M, M))
    for i in range(N):
        W[i, fast_src] = 1.2
        if two_rhythm:
            W[i, slow_src] = 1.2
        W[i, rng.choice(N, 4, replace=False)] = 0.3        # sparse recurrence
    net = Network(W, act=ACT, pas=8, theta=np.full(M, 1.0), p_s=0.0, seed=seed)
    net.tau = net.tau.astype(float)
    net.tau[:N] = rng.uniform(TMIN, TMAX, N)
    net.tau[fast_src], net.tau[slow_src] = P_F, P_S
    last = np.full(M, -1.0)
    for t in range(steps):
        d = np.zeros(M, bool)
        if t % P_F == 0:
            d[fast_src] = True
        if two_rhythm and t % P_S == 0:
            d[slow_src] = True
        prev = net.phi.copy()
        net.step(d)
        for i in np.where((net.phi == 1) & (prev == 0))[0]:
            if i < N and last[i] >= 0:
                net.tau[i] = np.clip(net.tau[i] + eta * ((t - last[i]) - net.tau[i]), TMIN, TMAX)
            last[i] = t
    tau = net.tau[:N]
    return dict(mean=tau.mean(), sd=tau.std(),
                near_f=(np.abs(tau - P_F) <= 2).mean(),
                near_s=(np.abs(tau - P_S) <= 2).mean())


# ---------------------------------------------------------------------------
# Diagnostic 3 (the decisive one) — IDEAL case: isolated nodes (no recurrence),
# each hand-fed a SINGLE rhythm. If resonance cannot lock tau to a period even
# here, the rule itself is the problem, not the coupling.
# Result: fast group (fed P_f=6) -> tau ~ 24-25 (a multiple, not 6);
#         slow group (fed P_s=24) -> tau rails at TMAX. Upward ratchet confirmed.
# ---------------------------------------------------------------------------

def diag3_ideal(N=100, steps=5000, eta=0.15, seed=0):
    rng = np.random.default_rng(seed)
    M = N + 2
    fast_src, slow_src = N, N + 1
    W = np.zeros((M, M))
    W[:N // 2, fast_src] = 1.5          # first half: fast only
    W[N // 2:N, slow_src] = 1.5         # second half: slow only
    net = Network(W, act=ACT, pas=8, theta=np.full(M, 1.0), p_s=0.0, seed=seed)
    net.tau = net.tau.astype(float)
    net.tau[:N] = rng.uniform(TMIN, TMAX, N)
    net.tau[fast_src], net.tau[slow_src] = P_F, P_S
    last = np.full(M, -1.0)
    for t in range(steps):
        d = np.zeros(M, bool)
        if t % P_F == 0:
            d[fast_src] = True
        if t % P_S == 0:
            d[slow_src] = True
        prev = net.phi.copy()
        net.step(d)
        for i in np.where((net.phi == 1) & (prev == 0))[0]:
            if i < N and last[i] >= 0:
                net.tau[i] = np.clip(net.tau[i] + eta * ((t - last[i]) - net.tau[i]), TMIN, TMAX)
            last[i] = t
    fast_tau, slow_tau = net.tau[:N // 2], net.tau[N // 2:N]
    return dict(fast_mean=fast_tau.mean(), fast_sd=fast_tau.std(),
                slow_mean=slow_tau.mean(), slow_sd=slow_tau.std())


def main():
    print("=== E10 diagnostics: can Line B build a fast/slow tau hierarchy? ===\n")

    print("Diagnostic 1 — naive resonance, recurrent pool, two-rhythm drive")
    for two in (True, False):
        r = diag1_naive(two_rhythm=two)
        print(f"  {'two-rhythm ' if two else 'fast-only  '}: tau mean={r['mean']:.1f} "
              f"sd={r['sd']:.1f}  near P_f={r['near_f']:.2f}  near P_s={r['near_s']:.2f}")
    print("  -> two-rhythm == fast-only, and no cluster at P_f: the slow rhythm is\n"
          "     invisible and tau does not lock to the drive periods.\n")

    print(f"Diagnostic 3 (decisive) — IDEAL: isolated nodes, each fed ONE rhythm")
    r = diag3_ideal()
    print(f"  fast group (fed P_f={P_F}): tau mean={r['fast_mean']:.1f} sd={r['fast_sd']:.1f}  "
          f"(target {P_F})")
    print(f"  slow group (fed P_s={P_S}): tau mean={r['slow_mean']:.1f} sd={r['slow_sd']:.1f}  "
          f"(target {P_S})")
    print("  -> even in the ideal case tau ratchets UP toward multiples / the ceiling;\n"
          "     it cannot settle at a fundamental period, so no fast population forms.")
    print("\nConclusion: the resonance rule is structurally wrong for 4a. A NEW,\n"
          "bidirectional tau rule is needed (see docs/e10_notes.md).")


if __name__ == "__main__":
    main()
