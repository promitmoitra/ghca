"""The regime law reduces to clock-shift merging: at tau_a >= tau_p every orbit
merges with its global-clock-shift, and that single fact is the whole boundary.

Follow-up to docs/spectrum_sufficiency_proof.md (P6: the open problem is
signature-measurability of the live set). Two reductions, both exhaustive:

R1. DECOMPOSITION. Signature-measurability of fate = invariance under
    (A) global clock-shift (add 1 to every phase, mod S+1) and
    (B) rearrangement (same gap multiset, different cyclic arrangement).
    Measured separately:
      (B) holds at EVERY cell in BOTH regimes (0 violations) -- rearrangement
          never changes fate, regime-independently.
      (A) holds exactly at tau_a >= tau_p and fails at every tau_p > tau_a
          cell tested (112/640/384/2080 violations at (1,2)/(2,3)/(2,4)/(3,4)).
    So the ENTIRE regime law is the clock-shift law: fate is invariant under
    global phase shift iff the active band covers the passive band. Absolute
    refractory depth -- the exact information the Section-10 ladder found
    missing -- is precisely what a clock-shift moves.

R2. THE MECHANISM IS ORBIT MERGING. At every tau_a >= tau_p cell, EVERY
    configuration c reaches the SAME attractor set as c+1 (100% of configs;
    e.g. 6561/6561 at (4,4)). At tau_p > tau_a, orbits either merge or differ
    in fate -- there are ZERO different-attractor-same-fate cases anywhere.
    Since fate is a tail property, merging implies (A); the dichotomy makes
    merging the exact mechanism, not merely a sufficient one.

THE REMAINING HAND-PROOF OBLIGATION, localised by a step-commutation analysis
(recorded for the whiteboard): comparing step(c+1) with step(c)+1 cell by cell,
the two agree EXCEPT
  - at a cell that DWELLS in c (receptive, stays receptive): step(c)+1 shows 1,
    step(c+1) shows 2 -- the shifted copy runs one tick ahead there; and
  - at a cell at state S in c (wraps to receptive): step(c+1) must fire-or-
    dwell where step(c)+1 already shows 1.
So the shifted trajectory is the original with dwell/wrap events displaced one
step in time, and the merge theorem to prove is that at tau_a >= tau_p these
displacements heal (the orbits re-synchronise) instead of compounding. The
(2,3) counterexample pair shows compounding: (0,0,1,3) lives while
(1,1,2,4) = (0,0,1,3)+1 dies.

House Rules Compliance: exhaustive, deterministic, no RNG; asserts R1 and R2
and aborts on regression.
Output: result/topology/clock_shift_merge.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
         (1, 2), (2, 3), (2, 4), (3, 4)]


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def fates(ta, tp):
    S = ta + tp
    fate = {}
    for c in product(range(S + 1), repeat=4):
        if c in fate:
            continue
        path, cur = [], c
        while cur not in fate and cur not in path:
            path.append(cur)
            cur = step2(cur, ta, tp)
        val = fate[cur] if cur in fate else any(
            any(x) for x in path[path.index(cur):])
        for p in path:
            fate[p] = val
    return fate


def attractor(cfg, ta, tp, cache):
    if cfg in cache:
        return cache[cfg]
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    att = frozenset(c for c, t in seen.items() if t >= seen[cur])
    for c in seen:
        cache[c] = att
    return cache[cfg]


def gaps_arr(cfg, B):
    ph = [cfg[i] for i in CYC]
    return tuple((ph[(k + 1) % 4] - ph[k]) % B for k in range(4))


def audit(ta, tp):
    S = ta + tp
    B = S + 1
    fate = fates(ta, tp)

    a_viol = sum(1 for c in fate for d in range(1, B)
                 if fate[tuple((x + d) % B for x in c)] != fate[c])

    by_ms = defaultdict(lambda: defaultdict(set))
    for c, f in fate.items():
        arr = gaps_arr(c, B)
        canon = min(arr[r:] + arr[:r] for r in range(4))
        by_ms[tuple(sorted(arr))][canon].add(f)
    b_viol = 0
    for ms, groups in by_ms.items():
        pures = [next(iter(f)) for f in groups.values() if len(f) == 1]
        if len(pures) > 1 and len(set(pures)) > 1:
            b_viol += 1

    cache = {}
    same_att = diff_same = diff_fate = 0
    for c in fate:
        c1 = tuple((x + 1) % B for x in c)
        if attractor(c, ta, tp, cache) == attractor(c1, ta, tp, cache):
            same_att += 1
        elif fate[c] == fate[c1]:
            diff_same += 1
        else:
            diff_fate += 1
    return dict(ta=ta, tp=tp, n=len(fate), a_viol=a_viol, b_viol=b_viol,
                same_att=same_att, diff_att_same_fate=diff_same,
                diff_fate=diff_fate)


def main():
    print(f"{'(ta,tp)':>8} {'regime':>7} {'N':>6} | {'A-viol':>7} {'B-viol':>7} | "
          f"{'merge':>6} {'d-same':>7} {'d-fate':>7}")
    recs = []
    for ta, tp in CELLS:
        r = audit(ta, tp)
        recs.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"{str((ta,tp)):>8} {regime:>7} {r['n']:6d} | {r['a_viol']:7d} "
              f"{r['b_viol']:7d} | {r['same_att']:6d} "
              f"{r['diff_att_same_fate']:7d} {r['diff_fate']:7d}", flush=True)

    for r in recs:
        cell = (r["ta"], r["tp"])
        assert r["b_viol"] == 0, f"rearrangement invariance broke at {cell}"
        assert r["diff_att_same_fate"] == 0, f"merge dichotomy broke at {cell}"
        if r["ta"] >= r["tp"]:
            assert r["a_viol"] == 0, f"clock-shift invariance broke at {cell}"
            assert r["same_att"] == r["n"], f"full merge broke at {cell}"
        else:
            assert r["a_viol"] > 0 and r["diff_fate"] > 0, \
                f"expected clock-shift failure at {cell}"

    # the counterexample pair cited in the docstring
    fate23 = fates(2, 3)
    assert fate23[(0, 0, 1, 3)] and not fate23[(1, 1, 2, 4)], \
        "(2,3) counterexample pair changed"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "clock_shift_merge.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")
    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | clock-shift viol. | rearr. viol. | "
          "orbits merging with c+1 | split-fate |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: |")
    for r in recs:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | {r['a_viol']} | "
              f"{r['b_viol']} | {r['same_att']}/{r['n']} | {r['diff_fate']} |")


if __name__ == "__main__":
    main()
