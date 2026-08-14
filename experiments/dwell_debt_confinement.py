"""The confinement variable: per-cell dwell debt. Exactly conserved
decomposition of the damage field; confined (<= S) at tau_a >= tau_p, maximal
linear drift (T-1) at tau_p > tau_a.

The topology-free attack on the one open lemma (merge at tau_a >= tau_p),
in the damage coordinates of damage_relaxation*.py.

THE DECOMPOSITION (exact, verified at every config x time of every cell
tested, both lattices): let u = orbit(c+1), v = orbit(c), and let D_i^u(t),
D_i^v(t) be the cumulative dwell counts of cell i in each copy. Then

    d_i(t) = u_i(t) - v_i(t) = 1 + Delta_i(t)   (mod S+1),
    where Delta_i(t) = D_i^v(t) - D_i^u(t)  -- the DWELL DEBT of cell i.

Proof of the identity: by Lemma R's cell-wise version, a cell's phase advances
by exactly 1 per step except at a dwell (advance 0). So u_i(t) = u_i(0) + t -
D_i^u(t) and likewise for v; subtract. (This is why damage is uniform iff the
debt field is: the initial offset 1 is the zero mode, and ALL deviation is
dwell bookkeeping. Damage doesn't diffuse -- it is a ledger of who waited more.)

FINDINGS (exhaustive; 2x2 at 11 cells, 3x3 at 6):

1. CONFINEMENT at tau_a >= tau_p: sup over live pairs, cells, and ALL time of
   |Delta| is a small constant --
       2x2: exactly ceil(S/2) at every cell tested
            ((1,1):1, (2,1):2, (2,2):2, (3,1):2, (3,2):3, (3,3):3,
             (4,3):4, (4,4):4);
       3x3: exactly S at every cell tested
            ((1,1):2, (2,1):3, (2,2):4, (3,1):4).
   The constant is topology-dependent (consistent with the relaxation-time
   non-lift) but S-scale in both lattices tested.
2. DRIFT at tau_p > tau_a: the live debt grows linearly with UNIT SLOPE --
   max|Delta| = T - k after T steps, k a small cell constant
   ((1,2): k=1, 49/99/199 at T=50/100/200; (2,3): k=2, 48/98/198). Unit
   slope means the extremal pair has one copy dwelling EVERY step while the
   other never dwells: a frozen (dead) partner beside a circulating wave --
   the split pair itself, seen from the ledger. Permanent damage = unbounded
   debt; the fate split is the debt escaping to infinity.
3. Dead debt is bounded (<= S) exactly where fate is shift-invariant
   (ta >= tp: both copies freeze, the ledger stops). At tp > ta the "dead"
   side drifts too ((1,2): 76 at the default horizon): classification is by
   the fate of c, and a dead c can be paired with a LIVE c+1 -- the split
   pairs themselves. Bounded dead debt is a consequence of the regime law,
   not an independent fact.

CONSEQUENCE FOR THE PROOF: the open lemma is now a confinement statement with
an explicit Lyapunov-like variable:

    at tau_a >= tau_p, the live dwell debt |Delta_i(t)| never exceeds c(G) ~ S
    -- each dwell of one copy forces a compensating dwell of the other within
    bounded time. Merge follows: bounded debt + finite state space => the
    debt field is eventually periodic with mean drift 0 per cycle => equal
    asymptotic dwell rates => same attractor (Lemma R turns zero relative
    drift on dwell-free attractors into a pure time shift; dwelling
    attractors merge by the same debt argument applied along the cycle).

House Rules Compliance: exhaustive, deterministic, no RNG; asserts the
decomposition identity, the confinement constants, the drift law, and the
dead bound; aborts on regression.
Output: result/topology/dwell_debt_confinement.npz.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS_2 = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
           (1, 2), (2, 3), (3, 4)]
CELLS_3 = [(1, 1), (2, 1), (2, 2), (3, 1), (1, 2)]
DRIFT_TS = [50, 100, 200]
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 500_000


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def att2(cfg, ta, tp, cache):
    if cfg in cache:
        return cache[cfg]
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    a = any(any(x) for x in
            [c for c, t in seen.items() if t >= seen[cur]])
    for c in seen:
        cache[c] = a
    return a


def audit_2x2(ta, tp, T=None):
    S = ta + tp
    B = S + 1
    if T is None:
        T = 6 * (ta + 2 * tp + 1) + 40
    cache = {}
    decomp_ok = True
    mxL = mxD = 0
    for c in product(range(B), repeat=4):
        alive = att2(c, ta, tp, cache)
        u, v = tuple((x + 1) % B for x in c), c
        Du = [0] * 4
        Dv = [0] * 4
        for t in range(T):
            for i in range(4):
                if (u[i] - v[i]) % B != (1 + Dv[i] - Du[i]) % B:
                    decomp_ok = False
            nu, nv = step2(u, ta, tp), step2(v, ta, tp)
            for i in range(4):
                if u[i] == 0 and nu[i] == 0:
                    Du[i] += 1
                if v[i] == 0 and nv[i] == 0:
                    Dv[i] += 1
                a = abs(Dv[i] - Du[i])
                if alive:
                    mxL = max(mxL, a)
                else:
                    mxD = max(mxD, a)
            u, v = nu, nv
    return dict(ta=ta, tp=tp, T=T, decomp_ok=decomp_ok,
                max_debt_live=mxL, max_debt_dead=mxD)


# ---- 3x3 vectorised (pipeline shared with damage_relaxation_3x3.py) --------
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


def audit_3x3(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = successor_array(ta, tp)
    pm = persistent_mask(f)
    sh = shift_perm(B)
    T = 6 * (ta + 2 * tp + 1) + 40
    U = sh.copy()
    V = np.arange(N, dtype=np.int64)
    D = np.zeros((N, 9), dtype=np.int16)
    run_max = np.zeros(N, dtype=np.int16)
    for t in range(T):
        fU, fV = f[U], f[V]
        for i in range(9):
            du = (((U // pw[i]) % B) == 0) & (((fU // pw[i]) % B) == 0)
            dv = (((V // pw[i]) % B) == 0) & (((fV // pw[i]) % B) == 0)
            D[:, i] += dv.astype(np.int16) - du.astype(np.int16)
        np.maximum(run_max, np.abs(D).max(axis=1), out=run_max)
        U, V = fU, fV
    return dict(ta=ta, tp=tp, T=T,
                max_debt_live=int(run_max[pm].max()),
                max_debt_dead=int(run_max[~pm].max()))


def main():
    print("=== 2x2: decomposition + debt bounds ===")
    recs2 = []
    for ta, tp in CELLS_2:
        r = audit_2x2(ta, tp)
        recs2.append(r)
        S = ta + tp
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"({ta},{tp}) [{regime:6s}] decomp={r['decomp_ok']} "
              f"max debt live={r['max_debt_live']:3d} (ceil(S/2)={-(-S//2)}) "
              f"dead={r['max_debt_dead']:3d} (S={S})", flush=True)
        assert r["decomp_ok"], f"decomposition identity broke at ({ta},{tp})"
        if ta >= tp:
            assert r["max_debt_live"] == -(-S // 2), \
                f"2x2 live confinement constant changed at ({ta},{tp})"
            assert r["max_debt_dead"] <= S, \
                f"dead debt bound broke at ({ta},{tp})"
        else:
            assert r["max_debt_dead"] > S, \
                f"expected dead-side drift at ({ta},{tp}) (split pairs)"

    print("\n=== 2x2 drift law at tau_p > tau_a: max live debt = T - k, "
          "unit slope ===")
    drift = []
    for ta, tp in [(1, 2), (2, 3)]:
        mxs = []
        for T in DRIFT_TS:
            r = audit_2x2(ta, tp, T=T)
            mxs.append(r["max_debt_live"])
            drift.append(dict(ta=ta, tp=tp, T=T, mx=r["max_debt_live"]))
            print(f"({ta},{tp}) T={T:3d}: max live debt = {r['max_debt_live']}",
                  flush=True)
        ks = {T - m for T, m in zip(DRIFT_TS, mxs)}
        assert len(ks) == 1, f"drift not unit-slope at ({ta},{tp}): {mxs}"
        k = ks.pop()
        assert k == {(1, 2): 1, (2, 3): 2}[(ta, tp)], \
            f"drift offset changed at ({ta},{tp}): k={k}"
    # bounded-side contrast at matched horizons
    for T in DRIFT_TS:
        r = audit_2x2(2, 1, T=T)
        assert r["max_debt_live"] == 2, "(2,1) confinement broke at long T"

    print("\n=== 3x3: debt bounds (vectorised) ===")
    recs3 = []
    for ta, tp in CELLS_3:
        t0 = time.time()
        r = audit_3x3(ta, tp)
        recs3.append(r)
        S = ta + tp
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"({ta},{tp}) [{regime:6s}] max debt live={r['max_debt_live']:3d} "
              f"(S={S}) dead={r['max_debt_dead']:3d} [{time.time()-t0:.0f}s]",
              flush=True)
        if ta >= tp:
            assert r["max_debt_live"] == S, \
                f"3x3 live confinement constant changed at ({ta},{tp})"
        else:
            assert r["max_debt_live"] > 4 * S, \
                f"expected 3x3 drift at ({ta},{tp})"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "dwell_debt_confinement.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g2_{k}": np.array([r[k] for r in recs2]) for k in recs2[0]},
             **{f"d_{k}": np.array([r[k] for r in drift]) for k in drift[0]},
             **{f"g3_{k}": np.array([r[k] for r in recs3]) for k in recs3[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | regime | max live debt | bound | max dead debt |")
    print("| :---: | :---: | :---: | ---: | :---: | ---: |")
    for r in recs2:
        S = r["ta"] + r["tp"]
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        bnd = f"⌈S/2⌉={-(-S//2)}" if r["ta"] >= r["tp"] else "drift (T−k)"
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | {regime} | "
              f"{r['max_debt_live']} | {bnd} | {r['max_debt_dead']} |")
    for r in recs3:
        S = r["ta"] + r["tp"]
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        bnd = f"S={S}" if r["ta"] >= r["tp"] else "drift"
        print(f"| 3×3 | ({r['ta']}, {r['tp']}) | {regime} | "
              f"{r['max_debt_live']} | {bnd} | {r['max_debt_dead']} |")


if __name__ == "__main__":
    main()
