"""Scoping the coherence invariant: Theorem 4's certificate LIFTS to 3x3
(anchor law proven at 3x3 (2,1), 483,446 pair states); the invariant cannot
be clock-shift-quotiented (absolute phase is load-bearing); and <=1-merged is
the sixth 2x2-smallness artifact (up to 6 merged cells at 3x3).

Answers "are we talking about plaquettes again, and how big do we need to
go?" for the open generalisation target of anchor_law_certificate.py.

FINDINGS:

1. THEOREM 4 AT 3x3 (2,1) -- the certificate method scales. From all 255,484
   persistent v, iterate the pair map on (v, Delta) with u := v+1+Delta,
   Delta_0 = 0 (no assumption on u's fate). Forward closure converges in 4
   steps at 483,446 states (|R|/live = 1.89 vs 1.53 at 2x2); on ALL of R:
     - ANCHOR LAW: min_i |Delta_i| = 0 at every state (with merge already
       exhaustive from damage_relaxation_3x3, the regime law at 3x3 (2,1)
       now has the same certificate status as 2x2);
     - floor Delta >= -1 universal; ceiling = 3 = S (matches the committed
       confinement constant);
     - >= 2 zeros at every state (two-zero rigidity SURVIVES the lattice);
     - <= 1 merged cell FAILS: up to 6 merged cells at once. The sixth
       2x2-smallness artifact; B5 must be restated lattice-free as
       "floor + two zeros", dropping the merged-count clause.
   Encoding lesson (bug caught in-session): the Delta digit encoding must be
   sized by the LATTICE constant S, not the 2x2 constant ceil(S/2) -- the
   first pass encoded base-4, Delta = 3 overflowed into the neighbouring
   digit, and 7 corrupted codes masqueraded as anchor violations. A scalar
   replay of a flagged trajectory exposed the corruption; base-8 with
   overflow guard fixed it. Asserted below via the overflow guard.

2. NO CLOCK-SHIFT QUOTIENT (2x2). Membership in R is NOT a function of
   (cyclic gap arrangement of v, Delta): 60 of 8,100 such keys are mixed --
   e.g. arrangement (0,2,2,0), Delta = (-1,0,0,0): v = (0,0,0,2) is IN R,
   its clock-shifts (1,1,1,3), (2,2,2,0), (3,3,3,1) are OUT. Adding u's
   arrangement does not help (still 60 mixed). So the coherence invariant
   MUST carry absolute clock phase -- consistent with the DEC-ladder result
   (sufficiency threshold at absolute state information) and with the
   compression barrier. (Note: arrangement + Delta + ONE absolute phase
   trivially decides -- it reconstructs v bijectively -- so that is a
   reparametrisation, not a compression; recorded to prevent rediscovery.)
   Spectrum (sorted) + Delta-multiset leaves 11 of 237 keys mixed.

3. SCOPE VERDICT. Plaquettes: yes for the RELATIVE part -- the arrangement
   generalises to the plaquette gap field, already validated as the right
   local object at 3x3 -- but the invariant needs (gap field, ledger, one
   absolute phase) jointly; pure relative-phase formulations are excluded by
   (2). How big: formulate at 2x2 (276 states, hand-inspectable), validate
   at 3x3 (done here for (2,1); other cells are a loop away), 4x4 only by
   sampling (B^16 configs; the successor array alone is ~34 GB -- out of
   exhaustive reach).

House Rules Compliance: exhaustive/deterministic, no RNG; asserts the 3x3
closure size, convergence, anchor, floor/ceiling, two-zero rigidity, the
6-merged maximum, the encoding overflow guard, and the 60/8100 mixed-key
count; aborts on regression.
Output: result/topology/coherence_invariant_scope.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

TA, TP = 2, 1
S = TA + TP
B = S + 1
C2 = -(-S // 2)          # 2x2 confinement constant
ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 500_000
OFF, BASE = 3, 8         # Delta digit encoding, sized for |Delta| <= S with guard


def step2(cfg):
    return tuple((1 if sum(1 <= cfg[j] <= TA for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % B for i, s in enumerate(cfg))


def live2():
    cache = {}
    for c in product(range(B), repeat=4):
        if c in cache:
            continue
        seen, cur = {}, c
        while cur not in seen:
            seen[cur] = len(seen)
            cur = step2(cur)
        f = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = f
    return {c for c, f in cache.items() if f}


def build_R2():
    live = live2()
    R = set()
    frontier = [(v, (0, 0, 0, 0)) for v in sorted(live)]
    while frontier:
        st = frontier.pop()
        if st in R:
            continue
        R.add(st)
        v, D = st
        u = tuple((v[i] + 1 + D[i]) % B for i in range(4))
        nu, nv = step2(u), step2(v)
        D2 = tuple(D[i] + (1 if (v[i] == 0 and nv[i] == 0) else 0)
                   - (1 if (u[i] == 0 and nu[i] == 0) else 0) for i in range(4))
        frontier.append((nv, D2))
    return live, R


def mixed_keys(live, R):
    def arr(v):
        ph = [v[i] for i in CYC]
        return tuple((ph[(k + 1) % 4] - ph[k]) % B for k in range(4))
    keys = defaultdict(set)
    for v in sorted(live):
        for D in product(range(-1, C2 + 1), repeat=4):
            u = tuple((v[i] + 1 + D[i]) % B for i in range(4))
            if u not in live:
                continue
            keys[(arr(v), tuple(D[i] for i in CYC))].add((v, D) in R)
    mixed = sum(1 for s in keys.values() if len(s) > 1)
    return mixed, len(keys)


def successor_array_3x3():
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        dig = ((idx[:, None] // pw[None, :]) % B).astype(np.int8)
        act = (dig >= 1) & (dig <= TA)
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
    return f, pw


def build_R3():
    N = B ** 9
    f, pw = successor_array_3x3()
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
    enc = (BASE ** np.arange(9)).astype(np.int64)
    V = live_idx.copy()
    U = sh[live_idx]
    D = np.zeros((len(V), 9), dtype=np.int8)
    seen = set()
    conv_t = -1
    for t in range(40):
        assert D.min() >= -OFF and D.max() < BASE - OFF, \
            "Delta encoding overflow -- widen BASE (the base-4 lesson)"
        code = V * (int(BASE ** 9)) + ((D.astype(np.int64) + OFF) @ enc)
        n0 = len(seen)
        seen.update(code.tolist())
        if len(seen) == n0 and t > 0:
            conv_t = t
            break
        fV, fU = f[V], f[U]
        for i in range(9):
            dv = (((V // pw[i]) % B) == 0) & (((fV // pw[i]) % B) == 0)
            du = (((U // pw[i]) % B) == 0) & (((fU // pw[i]) % B) == 0)
            D[:, i] += dv.astype(np.int8) - du.astype(np.int8)
        V, U = fV, fU
    codes = np.fromiter(seen, dtype=np.int64, count=len(seen))
    dd = np.empty((len(codes), 9), dtype=np.int8)
    rem = codes % int(BASE ** 9)
    for k in range(9):
        dd[:, k] = (rem % BASE) - OFF
        rem //= BASE
    return len(live_idx), len(seen), conv_t, dd


def main():
    print("=== 2x2: no clock-shift quotient ===")
    live, R2 = build_R2()
    assert len(R2) == 276, "2x2 certificate size changed"
    mixed, total = mixed_keys(live, R2)
    print(f"(arrangement, Delta) keys: {total}, mixed: {mixed}", flush=True)
    assert (mixed, total) == (60, 8100), "mixed-key census changed"

    print("\n=== 3x3 (2,1): Theorem 4 lifted ===")
    n_live, n_R, conv_t, dd = build_R3()
    anchor = bool((np.abs(dd).min(axis=1) == 0).all())
    floor_ok = int(dd.min()) == -1
    ceil_val = int(dd.max())
    min_zeros = int((dd == 0).sum(axis=1).min())
    max_merged = int((dd == -1).sum(axis=1).max())
    print(f"live={n_live:,} |R|={n_R:,} converged t={conv_t} anchor={anchor} "
          f"floor={int(dd.min())} ceil={ceil_val} min-zeros={min_zeros} "
          f"max-merged={max_merged}", flush=True)
    assert n_live == 255_484 and n_R == 483_446 and conv_t == 4
    assert anchor, "ANCHOR LAW failed at 3x3 (2,1)"
    assert floor_ok and ceil_val == S, "floor/ceiling changed"
    assert min_zeros == 2, "two-zero rigidity changed"
    assert max_merged == 6, "sixth-artifact census changed (<=1-merged was 2x2-only)"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_invariant_scope.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, mixed=np.array([mixed]), total=np.array([total]),
             n_live=np.array([n_live]), n_R=np.array([n_R]),
             conv_t=np.array([conv_t]), anchor=np.array([anchor]),
             floor=np.array([int(dd.min())]), ceil=np.array([ceil_val]),
             min_zeros=np.array([min_zeros]), max_merged=np.array([max_merged]))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
