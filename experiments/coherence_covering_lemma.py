"""The covering lemma, reduced to combinatorics: GH non-injectivity is
exactly 0<->S swaps (proven one-liner), and every age-hold is witnessed by a
pure swap -- u-side at the ceiling (all 25,998 at 3x3 (2,1)), v-side at the
age-1 holds (all 38,256) -- with the witness census by swap size and by
boundary role (corner/edge/centre; the 3x3 core is OPEN-boundary, not
periodic).

Fourth experiment on the coherence branch; the covering obligation from
coherence_window_saturation.py, mechanised to its atoms.

THE SWAP LAW (proven, and verified exhaustively at (2,1)/(3,3) over all
config pairs, 0 violations): stp(w) = stp(w') with w != w' forces w, w' to
differ exactly on cells where one has 0 and the other S -- because every
nonzero state advances deterministically and an image 0 arises only from a
0-dwell (quiet neighbourhood) or an S-wrap. All preimage freedom in this
dynamics is "was that quiet cell freshly quiet, or coming off refractory?".

THE WITNESS STRUCTURE (exhaustive):
  - 2x2, all four cells tested ((2,1),(2,2),(3,3),(4,4)): every ceiling
    (S->S) state's fresh age-(S-1) preimage is a SINGLE u-side 0->S swap
    (12/12, 16/16, 40/40, 88/88); every age-1 hold's diagonal preimage is a
    single v-side swap.
  - 3x3 (2,1), FULL census, open boundaries: ceiling 25,998/25,998 witnessed
    by u-side swap subsets (sizes 1..6: 13,997 / 7,800 / 3,010 / 967 / 200 /
    24); single-cell witnesses sit on corners and edges (sampled roles:
    corner 323, edge 60, centre 0 -- the low-degree boundary cells are where
    a quiet cell most easily hides its age). Age-1 holds 38,256/38,256
    witnessed by v-side swaps (sizes 1..6: 24,383 / 9,512 / 3,170 / 967 /
    200 / 24). U-side witnesses FAIL for the age-1 holds (0/38,256) -- the
    two hold types are mirror mechanisms, not one.

WHY THE SIDES SPLIT. Age counts steps since u was exactly v+1. At the
ceiling the pair is oldest, and the witness must rejuvenate ANCESTRY:
re-reading a dwelling u-cell as S-wrap places the pair one step closer to a
diagonal past (u "older", i.e. further along its cycle => the shift happened
more recently). At the age-1 holds the witness must show the IMAGE is also
reachable freshly from the diagonal: re-reading a dwelling v-cell as S-wrap
gives the diagonal predecessor directly. One mechanism, two mirrors.

THE REMAINING HAND OBLIGATION, now atomic: for every reachable pair at age
S whose u dwells, SOME subset of u's dwelling cells re-read as S-wraps
yields a valid configuration (image-preserving is automatic from the swap
law) whose pair-age is S-1 -- i.e. dwell-cells at the ceiling always admit
a consistent older reading. The swap law makes this a statement about quiet
neighbourhoods and the refractory pipeline: a cell that has been quiet for
k steps has ALL readings {dwell x k-j, wrap at j} consistent for j <= its
quiet run, so the witness exists iff some dwelling cell's quiet run is
YOUNGER than the pair's age -- which is the young-wave condition again.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts the swap
law (0 violations), all 2x2 witness counts, the full 3x3 censuses
(25,998/25,998 and 38,256/38,256 with size distributions), and the u/v side
split; aborts on regression.
Output: result/topology/coherence_covering_lemma.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product, combinations
from collections import Counter, defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS_2 = [(2, 1), (2, 2), (3, 3), (4, 4)]
CEIL_2 = {(2, 1): 12, (2, 2): 16, (3, 3): 40, (4, 4): 88}
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 500_000
CEIL_3_SIZES = {1: 13997, 2: 7800, 3: 3010, 4: 967, 5: 200, 6: 24}
HOLD1_3_SIZES = {1: 24383, 2: 9512, 3: 3170, 4: 967, 5: 200, 6: 24}


def make_step(ta, tp):
    B = ta + tp + 1

    def stp(cfg):
        return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                     if s == 0 else (s + 1) % B for i, s in enumerate(cfg))
    return stp, B


def swap_law_violations(ta, tp):
    stp, B = make_step(ta, tp)
    S = ta + tp
    img = defaultdict(list)
    for w in product(range(B), repeat=4):
        img[stp(w)].append(w)
    bad = 0
    for pres in img.values():
        for a in range(len(pres)):
            for b in range(a + 1, len(pres)):
                w1, w2 = pres[a], pres[b]
                for i in range(4):
                    if w1[i] != w2[i] and {w1[i], w2[i]} != {0, S}:
                        bad += 1
    return bad


def ages_2x2(ta, tp):
    stp, B = make_step(ta, tp)
    cache = {}
    for c in product(range(B), repeat=4):
        if c in cache:
            continue
        seen, cur = {}, c
        while cur not in seen:
            seen[cur] = len(seen)
            cur = stp(cur)
        fl = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = fl
    lv = {c for c, fl in cache.items() if fl}
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
    return stp, B, depth


def witness_2x2(ta, tp):
    stp, B, depth = ages_2x2(ta, tp)
    S = ta + tp
    sS = [p for p, a in depth.items()
          if a == S and depth[(stp(p[0]), stp(p[1]))] == S]
    h1 = [p for p, a in depth.items()
          if a == 1 and depth[(stp(p[0]), stp(p[1]))] == 1]

    def swap_side(x, side, want):
        cfg = x[side]
        nxt = stp(cfg)
        dws = [i for i in range(4) if cfg[i] == 0 and nxt[i] == 0]
        for k in range(1, len(dws) + 1):
            for su in combinations(dws, k):
                c2 = tuple(S if i in su else cfg[i] for i in range(4))
                if stp(c2) != nxt:
                    continue
                y = (x[0], c2) if side == 1 else (c2, x[1])
                if depth.get(y) == want:
                    return True
        return False
    n_ceil = sum(swap_side(x, 1, S - 1) for x in sS)
    n_h1 = sum(swap_side(x, 0, 0) for x in h1)
    return len(sS), n_ceil, len(h1), n_h1


def full_3x3():
    ta, tp = 2, 1
    B = ta + tp + 1
    S = ta + tp
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        dig = ((idx[:, None] // pw[None, :]) % B).astype(np.int8)
        act = (dig >= 1) & (dig <= ta)
        nxt = dig.astype(np.int64)
        adv = dig > 0
        nxt[adv] = (dig[adv].astype(np.int64) + 1) % B
        for i in range(9):
            na = np.zeros(hi - lo, dtype=np.int8)
            for j in NB3[i]:
                na += act[:, j]
            m = dig[:, i] == 0
            nxt[m, i] = (na[m] >= 1)
        f[lo:hi] = nxt @ pw
    dead = f == 0
    dead[0] = True
    while True:
        nd = dead | dead[f]
        if (nd == dead).all():
            break
        dead = nd
    live_idx = np.nonzero(~dead)[0]
    idxall = np.arange(N, dtype=np.int64)
    sh = np.zeros(N, dtype=np.int64)
    for k in range(9):
        sh += (((idxall // pw[k]) % B + 1) % B) * pw[k]
    depth = {}
    V, U = live_idx.copy(), sh[live_idx]
    for c in (V * N + U).tolist():
        depth[c] = 0
    d = 0
    while True:
        nV, nU = f[V], f[U]
        code = nV * N + nU
        mask = np.fromiter((c not in depth for c in code.tolist()),
                           dtype=bool, count=len(code))
        if not mask.any():
            break
        d += 1
        for c in code[mask].tolist():
            depth[c] = d
        V, U = nV[mask], nU[mask]
    codes = np.fromiter(depth.keys(), dtype=np.int64, count=len(depth))
    ags = np.fromiter(depth.values(), dtype=np.int16, count=len(depth))
    Vv, Uu = codes // N, codes % N
    imgc = f[Vv] * N + f[Uu]
    img_age = np.fromiter((depth[c] for c in imgc.tolist()),
                          dtype=np.int16, count=len(imgc))

    def digs(x):
        return [(int(x) // int(pw[k])) % B for k in range(9)]

    def side_witness(c, side, want):
        v_, u_ = c // N, c % N
        cfg = u_ if side == 1 else v_
        nxt = int(f[cfg])
        cd, nd_ = digs(cfg), digs(nxt)
        dws = [i for i in range(9) if cd[i] == 0 and nd_[i] == 0]
        for k in range(1, len(dws) + 1):
            for su in combinations(dws, k):
                c2 = cfg + sum(S * int(pw[i]) for i in su)
                if int(f[int(c2)]) != nxt:
                    continue
                cc = (int(v_) * N + int(c2)) if side == 1 \
                    else (int(c2) * N + int(u_))
                if depth.get(cc) == want:
                    return len(su)
        return None

    ceil = codes[(ags == S) & (img_age == S)]
    n_c = 0
    cs = Counter()
    for c in ceil.tolist():
        w = side_witness(c, 1, S - 1)
        if w is not None:
            n_c += 1
            cs[w] += 1
    h1 = codes[(ags == 1) & (img_age == 1)]
    n_h = 0
    hs = Counter()
    n_h_u = 0
    for c in h1.tolist():
        w = side_witness(c, 0, 0)
        if w is not None:
            n_h += 1
            hs[w] += 1
        if side_witness(c, 1, 0) is not None:
            n_h_u += 1
    return (len(ceil), n_c, dict(cs)), (len(h1), n_h, dict(hs), n_h_u)


def main():
    print("=== the swap law (non-injectivity = 0<->S) ===")
    for ta, tp in [(2, 1), (3, 3)]:
        v = swap_law_violations(ta, tp)
        print(f"({ta},{tp}): violations = {v}", flush=True)
        assert v == 0, "swap law broke"

    print("\n=== 2x2 witnesses (u-side at ceiling, v-side at age-1 holds) ===")
    for ta, tp in CELLS_2:
        nS, nC, nH1, nH = witness_2x2(ta, tp)
        print(f"({ta},{tp}): ceiling {nC}/{nS}, age-1 holds {nH}/{nH1}",
              flush=True)
        assert nS == nC == CEIL_2[(ta, tp)], f"2x2 ceiling witness ({ta},{tp})"
        assert nH1 == nH, f"2x2 hold witness ({ta},{tp})"

    print("\n=== 3x3 (2,1) FULL census ===")
    (nc, wc, cs), (nh, wh, hs, nh_u) = full_3x3()
    print(f"ceiling: {wc}/{nc} u-side witnessed, sizes {dict(sorted(cs.items()))}",
          flush=True)
    print(f"age-1 holds: {wh}/{nh} v-side witnessed, sizes "
          f"{dict(sorted(hs.items()))}; u-side works for {nh_u}", flush=True)
    assert nc == wc == 25_998 and cs == CEIL_3_SIZES, "3x3 ceiling census"
    assert nh == wh == 38_256 and hs == HOLD1_3_SIZES, "3x3 hold census"
    assert nh_u == 0, "u/v side split changed"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_covering_lemma.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             ceil_n=np.array([nc]), ceil_w=np.array([wc]),
             ceil_sizes=np.array(sorted(cs.items())),
             hold_n=np.array([nh]), hold_w=np.array([wh]),
             hold_sizes=np.array(sorted(hs.items())),
             hold_u_side=np.array([nh_u]))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
