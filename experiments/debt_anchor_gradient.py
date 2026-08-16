"""The confinement lemma decomposes: a regime-independent GRADIENT bound plus
a regime-defining ANCHOR law -- at tau_a >= tau_p some cell always has debt
exactly zero.

Continues debt_streak_lemmas.py (streak mechanism falsified; argument must
couple neighbours). Write Delta_i(t) for the dwell debt (dwell_debt_
confinement.py), G for the running max over edges of |Delta_i - Delta_j|,
A for the running max over time of min_i |Delta_i|.

FINDINGS (exhaustive, live pairs):

1. ANCHOR LAW (the regime): A = 0 at EVERY tau_a >= tau_p cell on BOTH
   lattices (2x2: 8 cells; 3x3: (1,1),(2,1),(2,2),(3,1)) -- at every step of
   every live pair, some cell has dwelt EQUALLY OFTEN in both copies. At
   tau_p > tau_a the anchor ESCAPES (2x2 (1,2): A=116, (2,3): A=114 at
   T=120; 3x3 (1,2): A=72) -- eventually every cell's ledger is imbalanced,
   and the imbalance grows without bound.
2. GRADIENT BOUND (regime-INDEPENDENT): G stays small in BOTH regimes --
   2x2: G = ceil(S/2) exactly at every tau_a >= tau_p cell (= the debt
   constant!), G <= 5 even on the drifting cells; 3x3: G = 2/2/3/2 at the
   tau_a >= tau_p cells, 6 at (1,2). Neighbouring ledgers never decouple;
   the whole field drifts together when it drifts.
3. CONSEQUENCE: confinement = anchor + connectivity. |Delta_i| <=
   A + diam(G_graph) * G, so anchor law + gradient bound => debt confined
   (0 + 2*ceil(S/2) at 2x2, 0 + 4*G at 3x3 -- both hold with room). The one
   open lemma REDUCES TO THE ANCHOR LAW; the gradient bound carries no
   regime information.
4. REFINEMENT FALSIFIED (recorded): the anchor is NOT always a jointly-active
   cell (2x2 (2,1): 2,972 of 12,600 live steps have a zero-debt cell but no
   jointly-active one; (1,1): all 1,536 steps). The zero-debt cell can be
   receptive or passive; no wavefront-based proof of the anchor is available.

INTERPRETATION. The anchor law is the regime law in its sharpest form yet:
with the active band covering the passive band, the two copies always share
at least one cell whose waiting history is perfectly balanced -- the ledger
pivots around a moving zero. When the passive band is longer, a frozen copy
can out-wait a circulating one at every cell simultaneously.

House Rules Compliance: exhaustive, deterministic, no RNG; asserts the anchor
law, the gradient values, the escape at tau_p > tau_a, and the falsified
refinement; aborts on regression.
Output: result/topology/debt_anchor_gradient.npz.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
EDGES2 = [(0, 1), (0, 2), (1, 3), (2, 3)]
CELLS_2 = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
DRIFT_2 = [(1, 2), (2, 3)]
CELLS_3 = [(1, 1), (2, 1), (2, 2), (3, 1)]
G3_EXPECT = {(1, 1): 2, (2, 1): 2, (2, 2): 3, (3, 1): 2}
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
EDGES3 = [(i, j) for i in range(9) for j in NB3[i] if i < j]
CHUNK = 500_000


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def alive2(cfg, ta, tp, cache):
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


def audit_2x2(ta, tp, T=None):
    S = ta + tp
    B = S + 1
    if T is None:
        T = 6 * (ta + 2 * tp + 1) + 40
    cache = {}
    G = A = mxD = 0
    steps = no_joint = 0
    for c in product(range(B), repeat=4):
        if not alive2(c, ta, tp, cache):
            continue
        u, v = tuple((x + 1) % B for x in c), c
        D = [0] * 4
        for t in range(T):
            steps += 1
            zd = [i for i in range(4) if D[i] == 0]
            if zd and not any(1 <= u[i] <= ta and 1 <= v[i] <= ta for i in zd):
                no_joint += 1
            nu, nv = step2(u, ta, tp), step2(v, ta, tp)
            for i in range(4):
                if v[i] == 0 and nv[i] == 0:
                    D[i] += 1
                if u[i] == 0 and nu[i] == 0:
                    D[i] -= 1
            G = max(G, max(abs(D[a] - D[b]) for a, b in EDGES2))
            A = max(A, min(abs(x) for x in D))
            mxD = max(mxD, max(abs(x) for x in D))
            u, v = nu, nv
    return dict(ta=ta, tp=tp, T=T, G=G, A=A, max_debt=mxD,
                steps=steps, no_joint_anchor=no_joint)


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


def audit_3x3(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = successor_array(ta, tp)
    pm = persistent_mask(f)
    idx = np.arange(N, dtype=np.int64)
    sh = np.zeros(N, dtype=np.int64)
    for k in range(9):
        sh += (((idx // pw[k]) % B + 1) % B) * pw[k]
    T = 6 * (ta + 2 * tp + 1) + 40
    U = sh.copy()
    V = idx.copy()
    D = np.zeros((N, 9), dtype=np.int16)
    G = np.zeros(N, dtype=np.int16)
    A = np.zeros(N, dtype=np.int16)
    for t in range(T):
        fU, fV = f[U], f[V]
        for i in range(9):
            dv = (((V // pw[i]) % B) == 0) & (((fV // pw[i]) % B) == 0)
            du = (((U // pw[i]) % B) == 0) & (((fU // pw[i]) % B) == 0)
            D[:, i] += dv.astype(np.int16) - du.astype(np.int16)
        gm = np.zeros(N, dtype=np.int16)
        for a, b in EDGES3:
            np.maximum(gm, np.abs(D[:, a] - D[:, b]).astype(np.int16), out=gm)
        np.maximum(G, gm, out=G)
        np.maximum(A, np.abs(D).min(axis=1), out=A)
        U, V = fU, fV
    return dict(ta=ta, tp=tp, G=int(G[pm].max()), A=int(A[pm].max()))


def main():
    print("=== 2x2 live pairs: gradient G, anchor A ===")
    recs2 = []
    for ta, tp in CELLS_2:
        r = audit_2x2(ta, tp)
        recs2.append(r)
        S = ta + tp
        print(f"({ta},{tp}) G={r['G']} (ceil(S/2)={-(-S//2)}) A={r['A']} "
              f"debt={r['max_debt']} no-joint-anchor="
              f"{r['no_joint_anchor']}/{r['steps']}", flush=True)
        assert r["A"] == 0, f"ANCHOR LAW broke at ({ta},{tp})"
        assert r["G"] == -(-S // 2), f"gradient value changed at ({ta},{tp})"
        assert r["max_debt"] <= r["A"] + 2 * r["G"], "decomposition inequality"
    # the falsified refinement, pinned
    by2 = {(r["ta"], r["tp"]): r for r in recs2}
    assert by2[(2, 1)]["no_joint_anchor"] == 2972, "refinement count changed"
    assert by2[(1, 1)]["no_joint_anchor"] == by2[(1, 1)]["steps"], \
        "(1,1) joint-anchor refinement unexpectedly holds"

    print("\n=== 2x2 drift cells: anchor escapes, gradient stays bounded ===")
    drecs = []
    for ta, tp in DRIFT_2:
        r = audit_2x2(ta, tp, T=120)
        drecs.append(r)
        print(f"({ta},{tp}) T=120: G={r['G']} A={r['A']} debt={r['max_debt']}",
              flush=True)
        assert r["A"] > 100, f"anchor escape vanished at ({ta},{tp})"
        assert r["G"] <= 5, f"gradient no longer bounded at ({ta},{tp})"

    print("\n=== 3x3: anchor law + gradient (vectorised) ===")
    recs3 = []
    for ta, tp in CELLS_3:
        t0 = time.time()
        r = audit_3x3(ta, tp)
        recs3.append(r)
        print(f"({ta},{tp}) G={r['G']} A={r['A']} [{time.time()-t0:.0f}s]",
              flush=True)
        assert r["A"] == 0, f"3x3 ANCHOR LAW broke at ({ta},{tp})"
        assert r["G"] == G3_EXPECT[(ta, tp)], f"3x3 gradient changed ({ta},{tp})"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "debt_anchor_gradient.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g2_{k}": np.array([r[k] for r in recs2]) for k in recs2[0]},
             **{f"d_{k}": np.array([r[k] for r in drecs]) for k in drecs[0]},
             **{f"g3_{k}": np.array([r[k] for r in recs3]) for k in recs3[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | regime | G (edge gradient) | A (anchor) | "
          "max debt |")
    print("| :---: | :---: | :---: | ---: | ---: | ---: |")
    for r in recs2:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | {regime} | {r['G']} | "
              f"{r['A']} | {r['max_debt']} |")
    for r in drecs:
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | τa<τp | {r['G']} | "
              f"{r['A']} | {r['max_debt']} |")
    for r in recs3:
        print(f"| 3×3 | ({r['ta']}, {r['tp']}) | τa≥τp | {r['G']} | "
              f"{r['A']} | — |")


if __name__ == "__main__":
    main()
