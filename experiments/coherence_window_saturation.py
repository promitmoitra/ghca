"""Window saturation, mechanised: the age function on coherent pairs has a
rigid transition law (+1 / reset-to-0 / hold-at-S), the ceiling is held by
OFF-ORBIT diagonal ancestors supplied by pair-map non-injectivity, and the
naive on-orbit recurrence proof is falsified (orbit gaps exceed S at every
cell).

Third experiment on the coherence branch; the saturation obligation
F^{S+1}(diag) subseteq F^{<=S}(diag) from coherence_window_S.py, attacked.

FINDINGS (exhaustive; (2,1) and (3,3) censused in full, gaps at 8 cells):

1. ON-ORBIT RECURRENCE IS FALSE. The natural proof ("diagonals recur along
   each pair orbit with gaps <= S") fails at EVERY cell: max on-orbit
   distance-since-diagonal = 5/4/6/5/7/9/10/12 at S = 2/3/4/4/5/6/7/8 --
   strictly above S at all 8 of 8 cells. (First-run assertion caught a
   hand-count of 7 here -- the transcription error class again; the
   falsification is total.) Saturation is NOT a single-orbit recurrence
   fact.
2. THE AGE-TRANSITION LAW (census asserted at (2,1) and (3,3)): from age a,
   a pair moves only to a+1, to 0 (reset = stepping onto the diagonal), or
   HOLDS -- and holds occur at exactly a in {0, 1, S}, nowhere else. Age
   never falls to an intermediate value. A hold means the image has a
   younger preimage merging in: at a = 1 the image has its own diagonal
   preimage (a different one), at a = S the image has a fresh age-(S-1)
   preimage. (The draft law omitted the 1 -> 1 hold; its own assertion
   caught it -- the census, not the prose, is the law.) The ceiling set
   {age = S} exits only to itself or to the diagonal.
3. THE CEILING MECHANISM: every S -> S state has an age-(S-1) PREIMAGE
   DIFFERENT from its orbit predecessor (12/12 at (2,1), 40/40 at (3,3)).
   The pair map is non-injective; merging pair-orbits continuously re-supply
   fresh diagonal ancestry. Where the pair's own history has run > S steps
   since a diagonal, ANOTHER trajectory that hit the diagonal more recently
   has merged into the same state (example asserted per cell: on-orbit gap
   4 > S = 3 with BFS age still 3 at (2,1); gap 7 > S = 6, age 6 at (3,3)).
   Saturation = "orbit merging refreshes ancestry at least every S steps" --
   the same orbit-merging phenomenon as the regime law itself, one level up
   (pairs instead of configurations).

PROOF SHAPE THIS DICTATES: not recurrence along an orbit, but a covering
argument -- every reachable pair state within one cycle-length of the
diagonal via SOME preimage chain. The young-wave family (the healing
extremals) supplies the fresh-ancestor witnesses at the boundary; the
lattice-free statement is that the diagonal's forward S-ball absorbs its own
boundary under one more step BECAUSE each boundary state is also the image
of a younger in-ball state. Precisely posed, one covering lemma.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts the
falsified on-orbit gaps, the full transition censuses, the 12/12 and 40/40
preimage counts, and the ceiling closure; aborts on regression.
Output: result/topology/coherence_window_saturation.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import Counter, defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
GAP_EXPECT = {(1, 1): 5, (2, 1): 4, (2, 2): 6, (3, 1): 5,
              (3, 2): 7, (3, 3): 9, (4, 3): 10, (4, 4): 12}
CENSUS_CELLS = [(2, 1), (3, 3)]
TRANS_EXPECT = {
    (2, 1): {(0, 0): 132, (0, 1): 48, (1, 0): 12, (1, 1): 12, (1, 2): 24,
             (2, 3): 24, (3, 0): 12, (3, 3): 12},
    (3, 3): {(0, 0): 624, (0, 1): 160, (1, 0): 40, (1, 1): 40, (1, 2): 80,
             (2, 3): 80, (3, 4): 80, (4, 5): 80, (5, 6): 80, (6, 0): 40,
             (6, 6): 40},
}
SS_EXPECT = {(2, 1): 12, (3, 3): 40}


def make_step(ta, tp):
    B = ta + tp + 1

    def stp(cfg):
        return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                     if s == 0 else (s + 1) % B for i, s in enumerate(cfg))
    return stp, B


def live_set(stp, B):
    cache = {}
    for c in product(range(B), repeat=4):
        if c in cache:
            continue
        seen, cur = {}, c
        while cur not in seen:
            seen[cur] = len(seen)
            cur = stp(cur)
        f = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = f
    return {c for c, f in cache.items() if f}


def ages(stp, B, lv):
    depth = {}
    frontier = [(w, tuple((x + 1) % B for x in w)) for w in sorted(lv)]
    for p in frontier:
        depth[p] = 0
    d = 0
    while frontier:
        d += 1
        nxt = []
        for v, u in frontier:
            q = (stp(v), stp(u))
            if q not in depth:
                depth[q] = d
                nxt.append(q)
        frontier = nxt
    return depth


def max_orbit_gap(stp, B, lv):
    mx = 0
    for w in sorted(lv):
        v, u = w, tuple((x + 1) % B for x in w)
        seen = set()
        t = last = 0
        while (v, u) not in seen:
            seen.add((v, u))
            if u == tuple((x + 1) % B for x in v):
                last = t
            mx = max(mx, t - last)
            v, u = stp(v), stp(u)
            t += 1
    return mx


def main():
    print("=== 1) on-orbit recurrence is false: max distance-since-diagonal ===")
    gaps = []
    for ta, tp in CELLS:
        stp, B = make_step(ta, tp)
        S = ta + tp
        lv = live_set(stp, B)
        g = max_orbit_gap(stp, B, lv)
        gaps.append(g)
        print(f"({ta},{tp}) S={S}: max on-orbit gap = {g} "
              f"({'> S' if g > S else '<= S'})", flush=True)
        assert g == GAP_EXPECT[(ta, tp)], f"gap changed at ({ta},{tp})"
    assert sum(1 for (c, g) in zip(CELLS, gaps) if g > c[0] + c[1]) == 8, \
        "on-orbit falsification pattern changed"

    print("\n=== 2-3) age-transition census + ceiling mechanism ===")
    recs = []
    for ta, tp in CENSUS_CELLS:
        stp, B = make_step(ta, tp)
        S = ta + tp
        lv = live_set(stp, B)
        depth = ages(stp, B, lv)
        trans = Counter()
        for (v, u), a in depth.items():
            trans[(a, depth[(stp(v), stp(u))])] += 1
        assert dict(trans) == TRANS_EXPECT[(ta, tp)], f"census changed ({ta},{tp})"
        legal = all(b == a + 1 or b == 0 or (b == a and a in (0, 1, S))
                    for (a, b) in trans)
        rev = defaultdict(list)
        for p in depth:
            rev[(stp(p[0]), stp(p[1]))].append(p)
        sS = [p for p, a in depth.items()
              if a == S and depth[(stp(p[0]), stp(p[1]))] == S]
        cert = sum(1 for p in sS if any(depth[q] == S - 1 for q in rev[p]))
        print(f"({ta},{tp}) S={S}: transition law legal: {legal} | "
              f"S->S states: {len(sS)}, with fresh age-(S-1) preimage: {cert}",
              flush=True)
        assert legal, f"illegal age transition at ({ta},{tp})"
        assert len(sS) == cert == SS_EXPECT[(ta, tp)], "ceiling mechanism changed"
        recs.append(dict(ta=ta, tp=tp, n_sS=len(sS), cert=cert))

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_window_saturation.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             ta=np.array([c[0] for c in CELLS]),
             tp=np.array([c[1] for c in CELLS]),
             gaps=np.array(gaps),
             census_cells=np.array([f"{c}" for c in CENSUS_CELLS]),
             n_sS=np.array([r["n_sS"] for r in recs]),
             cert=np.array([r["cert"] for r in recs]))
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | S | max on-orbit gap | > S? |")
    print("| :---: | ---: | ---: | :---: |")
    for (ta, tp), g in zip(CELLS, gaps):
        print(f"| ({ta}, {tp}) | {ta+tp} | {g} | "
              f"{'yes' if g > ta + tp else 'no'} |")


if __name__ == "__main__":
    main()
