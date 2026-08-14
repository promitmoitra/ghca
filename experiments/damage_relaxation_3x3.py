"""Damage relaxation at 3x3: the qualitative transition survives scaling; the
2x2 closed form does not -- and the right healing criterion is MERGE, not
return-to-uniform.

Scaling test of damage_relaxation.py (2x2). Same formalism: u = orbit of the
clock-shifted copy, v = orbit of c, damage d = u - v (mod S+1) per cell;
d uniform = the zero mode.

FINDINGS (exhaustive over all B^9 configurations per cell):

1. THE TRANSITION SURVIVES. Fate is clock-shift invariant (split-fate = 0) at
   EVERY tau_a >= tau_p cell tested -- (1,1), (2,1), (2,2), (3,1), (3,2), up
   to 10,077,696 configs -- and fails at (1,2) (77,376 split pairs). The
   regime law scales.

2. THE 2x2 CLOSED FORM DOES NOT LIFT. Max live relaxation-to-uniform at 3x3:
   (2,1): 9, (2,2): 12, (3,1): 10, (3,2): 14 -- vs tau_a + 2*tau_p + 1 =
   5, 7, 6, 8. The values are not even affine in (tau_a, tau_p) (the
   (2,1)->(3,1) increment is +1 while (2,2)->(3,2) is +2). The relaxation
   time is TOPOLOGY-DEPENDENT; tau_a + 2*tau_p + 1 is exact on the 2x2 core
   only. (The fourth 2x2-smallness lesson of this program, after automaton
   soundness, the entropy decoy, and the perimeter-signature failure.)

3. THE RIGHT CRITERION IS MERGE. At (1,1), 14,052 pairs NEVER return to
   uniform, yet split-fate = 0 and an explicit attractor-id computation shows
   EVERY orbit merges with its clock-shift. Cause: (1,1) attractors dwell
   (the 2x2 Lemma-R exception generalises), and on a dwelling attractor the
   two copies never differ by a constant. Return-to-uniform is a sufficient
   witness of merging where attractors are dwell-free; merging itself is the
   regime-law mechanism. Likewise at (1,2) never-return (105,184) exceeds
   split-fate (77,376): the excess is dwelling-attractor pairs that merge
   without uniformising, so the clean 2x2 identity never-return == split-fate
   is also a dwell-free artifact.

CONSEQUENCE FOR THE PROOF: the target theorem is
    "at tau_a >= tau_p every orbit merges with its clock-shift, within a
     topology-dependent bounded time"
with the 2x2 formula demoted to the exact constant for L=4 loops. Any
confinement argument must not lean on the 2x2 pipeline arithmetic.

House Rules Compliance: exhaustive, deterministic, no RNG; asserts findings
1-3 and aborts on regression.
Output: result/topology/damage_relaxation_3x3.npz.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 500_000
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (1, 2)]
EXPECT_MX = {(2, 1): 9, (2, 2): 12, (3, 1): 10, (3, 2): 14}


def successor_array(ta, tp):
    B = ta + tp + 1
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
    return f


def persistent_mask(f):
    dead = f == 0
    dead[0] = True
    while True:
        nd = dead | dead[f]
        if (nd == dead).all():
            return ~nd
        dead = nd


def shift_perm(B):
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    idx = np.arange(N, dtype=np.int64)
    out = np.zeros(N, dtype=np.int64)
    for k in range(9):
        out += (((idx // pw[k]) % B + 1) % B) * pw[k]
    return out


def audit(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = successor_array(ta, tp)
    pm = persistent_mask(f)
    sh = shift_perm(B)
    split = int((pm != pm[sh]).sum())
    cap = 10 * (ta + 2 * tp + 1) + 60
    mx_live = 0
    never = 0
    n_break = 0
    CH = 1_000_000
    for lo in range(0, N, CH):
        hi = min(lo + CH, N)
        U = sh[lo:hi].copy()
        V = np.arange(lo, hi, dtype=np.int64)
        broken = np.zeros(hi - lo, dtype=bool)
        ret_t = np.full(hi - lo, -1, dtype=np.int32)
        unres = np.ones(hi - lo, dtype=bool)
        for t in range(cap):
            w = np.nonzero(unres)[0]
            if len(w) == 0:
                break
            uw, vw = U[w], V[w]
            d0 = (uw % B - vw % B) % B
            uni = np.ones(len(w), dtype=bool)
            for k in range(1, 9):
                uni &= (((uw // pw[k]) % B - (vw // pw[k]) % B) % B == d0)
            nb = ~uni & ~broken[w]
            broken[w[nb]] = True
            ret = uni & broken[w]
            ret_t[w[ret]] = t
            unres[w[ret]] = False
            live_w = w[unres[w]]
            U[live_w], V[live_w] = f[U[live_w]], f[V[live_w]]
        n_break += int(broken.sum())
        nv = broken & (ret_t < 0)
        never += int(nv.sum())
        seg = pm[lo:hi] & (ret_t >= 0)
        if seg.any():
            mx_live = max(mx_live, int(ret_t[seg].max()))
    return dict(ta=ta, tp=tp, N=N, n_break=n_break, never=never, split=split,
                mx_live=mx_live, formula_2x2=ta + 2 * tp + 1)


def merge_total_11():
    """(1,1): every orbit merges with its clock-shift (attractor-id check)."""
    ta, tp = 1, 1
    B = 3
    N = B ** 9
    f = successor_array(ta, tp)
    sh = shift_perm(B)
    g = f.copy()
    for _ in range(60):
        g = g[g]
    ids = np.full(N, -1, dtype=np.int64)
    for c in range(N):
        x = g[c]
        if ids[x] >= 0:
            continue
        cyc = [x]
        y = f[x]
        while y != x:
            cyc.append(y)
            y = f[y]
        m = min(cyc)
        for z in cyc:
            ids[z] = m
    aid = ids[g]
    return bool((aid == aid[sh]).all())


def main():
    print(f"{'(ta,tp)':>8} {'regime':>7} {'N':>11} | {'scatter':>10} "
          f"{'never-ret':>10} {'split':>8} | {'max live relax':>14} "
          f"{'2x2 formula':>11}")
    recs = []
    for ta, tp in CELLS:
        t0 = time.time()
        r = audit(ta, tp)
        recs.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"{str((ta,tp)):>8} {regime:>7} {r['N']:11,d} | {r['n_break']:10,d} "
              f"{r['never']:10,d} {r['split']:8,d} | {r['mx_live']:14d} "
              f"{r['formula_2x2']:11d}  [{time.time()-t0:.0f}s]", flush=True)

    by = {(r["ta"], r["tp"]): r for r in recs}
    # 1: the transition
    for cell, r in by.items():
        if cell[0] >= cell[1]:
            assert r["split"] == 0, f"shift-invariance of fate broke at {cell}"
    assert by[(1, 2)]["split"] == 77376, "(1,2) split count changed"
    # 2: the non-lift
    for cell, want in EXPECT_MX.items():
        r = by[cell]
        assert r["mx_live"] == want, f"3x3 relaxation time changed at {cell}"
        assert r["mx_live"] != r["formula_2x2"], \
            f"2x2 formula unexpectedly lifted at {cell}"
        assert r["never"] == 0, f"dwell-free never-return broke at {cell}"
    # 3: merge vs uniform at (1,1)
    assert by[(1, 1)]["never"] == 14052, "(1,1) never-return count changed"
    assert merge_total_11(), "(1,1) merge totality broke"
    assert by[(1, 2)]["never"] > by[(1, 2)]["split"], \
        "(1,2) never-return no longer exceeds split-fate"
    print("\nall assertions pass (incl. (1,1) merge totality despite "
          "14,052 never-uniform pairs)")

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "damage_relaxation_3x3.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"wrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | N | scatter | never-return | split-fate | "
          "max live relax | 2×2 formula |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in recs:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | {r['N']:,} | "
              f"{r['n_break']:,} | {r['never']:,} | {r['split']:,} | "
              f"{r['mx_live']} | {r['formula_2x2']} |")


if __name__ == "__main__":
    main()
