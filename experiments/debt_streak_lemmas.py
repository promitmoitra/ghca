"""The whiteboard attempt on the confinement lemma: a clean 2x2 skeleton
(debt constant = max dwell streak = ceil(S/2); matching window = S), and its
FALSIFICATION at 3x3 -- confinement is not explained by short streaks.

Attempt at the one open lemma (dwell_debt_confinement.py: live debt bounded
at tau_a >= tau_p). The candidate skeleton had two local ingredients:

  C1 (streak bound)     live trajectories have dwell streaks <= ceil(S/2);
  C2 (matching window)  in a live pair, every dwell of one copy at cell i has
                        a dwell of the other copy at i within S steps;

and the observation that would have made them sufficient:

  2x2 IDENTITY (exhaustive, 8 cells): max live dwell streak = ceil(S/2)
  = the measured debt-confinement constant, exactly, at every tau_a >= tau_p
  cell. Matching window = S everywhere except (1,1) (= S+1, the dwelling-
  attractor cell). Debt looked like "at most one unmatched streak in flight".

FALSIFICATION (exhaustive, 3x3): max live dwell streak at (1,1)/(2,1)/(2,2)/
(3,1) is 5/6/8/6 -- far ABOVE the measured debt constants 2/3/4/4 (= S).
Long quiet spells exist (a corner cell can wait many steps while the wave
circulates elsewhere) yet the ledger stays S-confined. So confinement is NOT
short streaks; it is EPISODE OVERLAP: during a long dwell episode of v at
cell i, the partner u must be dwelling at i almost simultaneously (bounded
debt forces >= ell - 2*const of any ell-long v-streak to be matched by
u-dwells -- but that direction is a consequence of confinement, not a proof
of it). The fifth 2x2-smallness artifact of this program, caught before the
proof was built on it.

WHAT SURVIVES FOR THE PROOF: any valid argument must couple the dwell
EPISODES of the two copies locally (a dwell at i happens iff i's closed
neighbourhood is inactive; the partner sees the same neighbourhood shifted by
1 + the neighbours' debts, so episode boundaries differ by neighbour debts --
an induction on the damage field along graph edges), rather than bounding
streak lengths. The 2x2 identity stands as the exact L=4 answer and as the
base case.

House Rules Compliance: exhaustive, deterministic, no RNG; asserts the 2x2
identity, the window values, AND the 3x3 falsification (so the negative
result cannot rot either); aborts on regression.
Output: result/topology/debt_streak_lemmas.npz.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS_2 = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
CELLS_3 = [(1, 1), (2, 1), (2, 2), (3, 1)]
DEBT_3 = {(1, 1): 2, (2, 1): 3, (2, 2): 4, (3, 1): 4}     # from dwell_debt npz
STREAK_3 = {(1, 1): 5, (2, 1): 6, (2, 2): 8, (3, 1): 6}
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
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


def audit_2x2(ta, tp):
    S = ta + tp
    B = S + 1
    cache = {}
    T = 6 * (ta + 2 * tp + 1) + 40
    mx_streak = w_max = mxD = 0
    for c in product(range(B), repeat=4):
        if not alive2(c, ta, tp, cache):
            continue
        u, v = tuple((x + 1) % B for x in c), c
        dw_u = [[] for _ in range(4)]
        dw_v = [[] for _ in range(4)]
        Du = [0] * 4
        Dv = [0] * 4
        streak = [0] * 4
        for t in range(T):
            nu, nv = step2(u, ta, tp), step2(v, ta, tp)
            for i in range(4):
                if v[i] == 0 and nv[i] == 0:
                    dw_v[i].append(t)
                    Dv[i] += 1
                    streak[i] += 1
                else:
                    streak[i] = 0
                mx_streak = max(mx_streak, streak[i])
                if u[i] == 0 and nu[i] == 0:
                    dw_u[i].append(t)
                    Du[i] += 1
                mxD = max(mxD, abs(Dv[i] - Du[i]))
            u, v = nu, nv
        for i in range(4):
            for t in dw_v[i]:
                if dw_u[i]:
                    w_max = max(w_max, min(abs(t - s) for s in dw_u[i]))
            for t in dw_u[i]:
                if dw_v[i]:
                    w_max = max(w_max, min(abs(t - s) for s in dw_v[i]))
    return dict(ta=ta, tp=tp, max_streak=mx_streak, window=w_max,
                max_debt=mxD)


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


def audit_3x3_streak(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = successor_array(ta, tp)
    pm = persistent_mask(f)
    T = 6 * (ta + 2 * tp + 1) + 40
    V = np.arange(N, dtype=np.int64)
    streak = np.zeros((N, 9), dtype=np.int16)
    mx = np.zeros(N, dtype=np.int16)
    for t in range(T):
        fV = f[V]
        for i in range(9):
            dv = (((V // pw[i]) % B) == 0) & (((fV // pw[i]) % B) == 0)
            streak[:, i] = np.where(dv, streak[:, i] + 1, 0)
        np.maximum(mx, streak.max(axis=1), out=mx)
        V = fV
    return int(mx[pm].max())


def main():
    print("=== 2x2: streak / window / debt (live, tau_a >= tau_p) ===")
    recs2 = []
    for ta, tp in CELLS_2:
        r = audit_2x2(ta, tp)
        recs2.append(r)
        S = ta + tp
        print(f"({ta},{tp}) streak={r['max_streak']} (ceil(S/2)={-(-S//2)}) "
              f"window={r['window']} (S={S}) debt={r['max_debt']}", flush=True)
        assert r["max_streak"] == -(-S // 2), f"2x2 streak law broke ({ta},{tp})"
        assert r["max_debt"] == -(-S // 2), f"2x2 debt law broke ({ta},{tp})"
        want_w = S + 1 if (ta, tp) == (1, 1) else S
        assert r["window"] == want_w, f"2x2 window law broke ({ta},{tp})"

    print("\n=== 3x3: the falsification (streak >> debt constant) ===")
    recs3 = []
    for ta, tp in CELLS_3:
        t0 = time.time()
        s = audit_3x3_streak(ta, tp)
        recs3.append(dict(ta=ta, tp=tp, max_streak=s, debt_const=DEBT_3[(ta, tp)]))
        print(f"({ta},{tp}) streak={s} vs debt constant={DEBT_3[(ta, tp)]} "
              f"[{time.time()-t0:.0f}s]", flush=True)
        assert s == STREAK_3[(ta, tp)], f"3x3 streak changed at ({ta},{tp})"
        assert s > DEBT_3[(ta, tp)], \
            f"3x3 falsification vanished at ({ta},{tp}) -- revisit the skeleton"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "debt_streak_lemmas.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g2_{k}": np.array([r[k] for r in recs2]) for k in recs2[0]},
             **{f"g3_{k}": np.array([r[k] for r in recs3]) for k in recs3[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | max live streak | ⌈S/2⌉ | window | debt const |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: |")
    for r in recs2:
        S = r["ta"] + r["tp"]
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | {r['max_streak']} | "
              f"{-(-S//2)} | {r['window']} | {r['max_debt']} |")
    for r in recs3:
        print(f"| 3×3 | ({r['ta']}, {r['tp']}) | {r['max_streak']} | — | — | "
              f"{r['debt_const']} |")


if __name__ == "__main__":
    main()
