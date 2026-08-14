"""Blackboard session on the anchor law: the debt floor is exactly -1 via
merged cells, the merge/unmerge events are classified, the regime-carrying
atomic fact is isolated ("a merged receptive cell never sees an active
v-neighbour"), and the live debt field is far more rigid than the anchor law
needs -- at most one merged cell, at least two zeros, and no (-1, >=1) edge.

The oracle-refereed attempt at proving debt_anchor_gradient.py's ANCHOR LAW.
Every step below is a claim I formulated at the board and the oracle
(exhaustive enumeration) checked. All hold at every tau_a >= tau_p 2x2 cell
tested; the atomic fact FAILS at tau_p > tau_a, as it must (it carries the
regime).

THE LADDER (live clock-shift pairs, u = orbit(c+1), v = orbit(c)):

  B1 (floor). Delta_i ranges in [-1, ceil(S/2)]: the negative side is
     EXACTLY -1, never -2. Two-sided but asymmetric -- u (the shifted copy,
     one tick "older") can out-dwell v by at most one waiting event, ever.
  B2 (merge identity). Delta_i = -1  <=>  u_i = v_i. The ledger hitting -1
     IS local merging: u_i = v_i + 1 + Delta_i (the decomposition identity),
     so Delta_i = -1 iff the copies agree at i. 0 violations.
  B3 (event classification). The 0 -> -1 transition happens ONLY via
     (u_i dwells, v_i wraps S -> 0): the young copy waits one beat while the
     old copy comes off refractory -- and they are equal from then on until
     unmerged. The -1 -> 0 transition (unmerge) is v_i dwelling at a merged
     receptive cell while u_i fires. Merge and unmerge counts balance
     exactly on live pairs at every dwell-free-attractor cell (96/96 at
     (2,1) ... 1536/1536 at (4,4)); at (1,1) -- the dwelling-attractor
     exception yet again -- merges lead by exactly 16 (1016/1000): the
     attractor itself dwells, so merge/unmerge events recur forever and
     some are mid-flight at any finite horizon.
  B4 (the atomic fact, REGIME-CARRYING). At a merged receptive cell
     (Delta_i = -1, u_i = v_i = 0) at tau_a >= tau_p, NO v-neighbour is
     active -- the strong form holds: every v-neighbour is receptive or
     passive (measured distributions: only R and P bands, e.g. (3,2):
     480 R / 320 P1 / 320 P2). Hence v cannot fire at a merged cell before
     u does, hence -1 -> -2 is impossible: B4 => B1's floor. At
     tau_p > tau_a the fact FAILS (64 events at (1,2), 160 at (2,3)) --
     this is where the regime enters the ledger.
  B5 (rigidity of live shapes). The live debt field's sorted multiset only
     ever takes shapes with AT LEAST TWO zeros, AT MOST ONE -1, and never
     an edge between a -1 cell and a >= 1 cell; all-merged never occurs and
     no-zero never occurs ((0,0,0,0), (-1,0,0,0), (0,0,k,l) with k,l >= 0
     exhaust what is seen). The anchor law holds with multiplicity two.

WHAT REMAINS (the honest gap): B4 is verified, not proven -- the hand
argument needs WHY a merged receptive cell's v-neighbourhood is quiet
exactly when tau_a >= tau_p. Sketch worth pursuing: at the merge moment
v_i = S, so the wave last left i a full excursion ago; a v-neighbour active
now would have to have fired within the last tau_a steps, i.e. the wave
re-approached i while u_i (one tick ahead) was still receptive-dwelling; at
tau_a >= tau_p the young copy's receptive window is short enough that the
returning wavefront reaches u and v within the same beat. Also open: B5's
two-zero rigidity (stronger than needed; likely the loop parity of the 2x2).

House Rules Compliance: exhaustive, deterministic, no RNG; asserts B1-B5 and
the tau_p > tau_a failure of B4; aborts on regression.
Output: result/topology/anchor_law_blackboard.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
EDGES2 = [(0, 1), (0, 2), (1, 3), (2, 3)]
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
WRONG = [(1, 2), (2, 3)]
B4_FAIL_EXPECT = {(1, 2): 64, (2, 3): 160}


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def alive(cfg, ta, tp, cache):
    if cfg in cache:
        return cache[cfg]
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    f = any(any(x) for x in [c for c, t in seen.items() if t >= seen[cur]])
    for c in seen:
        cache[c] = f
    return f


def audit(ta, tp, T=None):
    S = ta + tp
    B = S + 1
    if T is None:
        T = 6 * (ta + 2 * tp + 1) + 40
    cache = {}
    minD = maxD = 0
    b2v = b3v_merge = b3v_never = 0
    merges = unmerges = 0
    b4_strong_viol = 0
    nbr_bands = defaultdict(int)
    all_merged = no_zero = edge_jump = 0
    for c in product(range(B), repeat=4):
        if not alive(c, ta, tp, cache):
            continue
        u, v = tuple((x + 1) % B for x in c), c
        D = [0] * 4
        for t in range(T):
            minD = min(minD, min(D))
            maxD = max(maxD, max(D))
            for i in range(4):
                if (D[i] == -1) != (u[i] == v[i]):
                    b2v += 1
                if D[i] == -1 and u[i] == 0:
                    for j in ADJ2[i]:
                        band = ("R" if v[j] == 0 else
                                "A" if v[j] <= ta else "P")
                        nbr_bands[band] += 1
                        if 1 <= v[j] <= ta:
                            b4_strong_viol += 1
            if all(x == -1 for x in D):
                all_merged += 1
            if 0 not in D:
                no_zero += 1
            for a, b in EDGES2:
                if (D[a] == -1 and D[b] >= 1) or (D[b] == -1 and D[a] >= 1):
                    edge_jump += 1
            nu, nv = step2(u, ta, tp), step2(v, ta, tp)
            for i in range(4):
                du_ = (u[i] == 0 and nu[i] == 0)
                dv_ = (v[i] == 0 and nv[i] == 0)
                dD = (1 if dv_ else 0) - (1 if du_ else 0)
                if D[i] == 0 and dD == -1:
                    merges += 1
                    if not (du_ and v[i] == S):
                        b3v_merge += 1
                if D[i] == -1 and dD == -1:
                    b3v_never += 1
                if D[i] == -1 and dD == +1:
                    unmerges += 1
                D[i] += dD
            u, v = nu, nv
    return dict(ta=ta, tp=tp, min_debt=minD, max_debt=maxD, b2_viol=b2v,
                merge_events=merges, unmerge_events=unmerges,
                b3_bad_merge=b3v_merge, b3_below_floor=b3v_never,
                b4_strong_viol=b4_strong_viol,
                nbr_R=nbr_bands["R"], nbr_A=nbr_bands["A"],
                nbr_P=nbr_bands["P"], all_merged=all_merged,
                no_zero=no_zero, edge_jump=edge_jump)


def main():
    print("=== ladder B1-B5 at tau_a >= tau_p ===")
    recs = []
    for ta, tp in CELLS:
        r = audit(ta, tp)
        recs.append(r)
        print(f"({ta},{tp}) floor={r['min_debt']} b2v={r['b2_viol']} "
              f"merge/unmerge={r['merge_events']}/{r['unmerge_events']} "
              f"b3v={r['b3_bad_merge']}+{r['b3_below_floor']} "
              f"B4 nbrs R/A/P={r['nbr_R']}/{r['nbr_A']}/{r['nbr_P']} "
              f"allm/nozero/edgejump={r['all_merged']}/{r['no_zero']}/"
              f"{r['edge_jump']}", flush=True)
        assert r["min_debt"] == -1, f"B1 floor changed at ({ta},{tp})"
        assert r["b2_viol"] == 0, f"B2 merge identity broke at ({ta},{tp})"
        assert r["b3_bad_merge"] == 0 and r["b3_below_floor"] == 0, \
            f"B3 event classification broke at ({ta},{tp})"
        if (ta, tp) == (1, 1):
            assert r["merge_events"] - r["unmerge_events"] == 16, \
                "(1,1) in-flight merge count changed"
        else:
            assert r["merge_events"] == r["unmerge_events"], \
                f"merge/unmerge balance broke at ({ta},{tp})"
        assert r["b4_strong_viol"] == 0 and r["nbr_A"] == 0, \
            f"B4 atomic fact broke at ({ta},{tp})"
        assert r["all_merged"] == 0 and r["no_zero"] == 0 \
            and r["edge_jump"] == 0, f"B5 rigidity broke at ({ta},{tp})"

    print("\n=== B4 must FAIL at tau_p > tau_a ===")
    wrecs = []
    for ta, tp in WRONG:
        r = audit(ta, tp, T=120)
        wrecs.append(r)
        print(f"({ta},{tp}) merged-receptive cells with active v-neighbour: "
              f"{r['b4_strong_viol']}", flush=True)
        assert r["b4_strong_viol"] == B4_FAIL_EXPECT[(ta, tp)], \
            f"B4 failure count changed at ({ta},{tp})"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "anchor_law_blackboard.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g_{k}": np.array([r[k] for r in recs]) for k in recs[0]},
             **{f"w_{k}": np.array([r[k] for r in wrecs]) for k in wrecs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | floor | merge=unmerge | B4 active-nbr events | "
          "all-merged | no-zero | (-1,≥1) edges |")
    print("| :---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in recs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['min_debt']} | "
              f"{r['merge_events']} | {r['b4_strong_viol']} | "
              f"{r['all_merged']} | {r['no_zero']} | {r['edge_jump']} |")
    for r in wrecs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['min_debt']} | — | "
              f"{r['b4_strong_viol']} | — | — | — |")


if __name__ == "__main__":
    main()
