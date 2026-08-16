"""THE COHERENCE INVARIANT, FOUND: a pair is coherent iff it lies within S
steps downstream of an exact clock-shift (diagonal) state -- the memory window
is exactly S = tau_a + tau_p, the full cycle length, at every 2x2 cell tested
and at 3x3 (2,1).

Second experiment on the coherence branch; resolves the target left by
coherence_formulate_2x2.py (an intrinsic predicate separating the 276
certified pairs from the 192 lockstep impostors).

THE CHARACTERISATION. Let diag = {(w, w+1) : w persistent} and F the pair map
(component-wise step). Then

    R  =  union over k = 0..S of F^k(diag),

with S sharp: BFS from diag reaches all of R at depth exactly S at every cell
((1,1): 2, (2,1): 3, (2,2): 4, (3,1): 4, (3,2): 5, (3,3): 6, (4,3): 7,
(4,4): 8 -- candidates S+1, ta+2tp+1, 2S-2 all falsified by (2,1) or (3,1)).
Equivalently, INTRINSICALLY: (v, u) is coherent iff its backward window of
length S contains a diagonal state -- a finite-memory (sofic-style) condition
with window = the cycle length. The 192 impostors are precisely the
lockstep-bound pairs whose backward window is diagonal-free (they flow INTO
R -- the traced example enters R at t=3 -- but were never themselves an exact
clock-shift).

WHY S IS THE RIGHT CONSTANT -- **SUPERSEDED, see below.** The rationale in
this paragraph was falsified by coherence_window_saturation.py (commit
5ba2c5f): on-orbit recurrence gaps are 5/4/6/5/7/9/10/12, strictly ABOVE S at
all 8 cells, so a diagonal state does NOT recur on every lockstep cycle. The
window is still exactly S (that reproduces); the explanation below is wrong.
The current candidate mechanism is the quiet-run decay in
coherence_covering_lemma.py, which is an OPEN lemma, not an account.

Superseded text: A diagonal state recurs on every lockstep cycle
(the pair returns to u = v+1 within the attractor period); S is one full
excursion+refractory span, so "within S of a diagonal" is exactly "the pair's
current disagreement pattern is one the clock-shift can generate inside a
single cycle". This also retro-explains the certificate convergence depths
(2x2 (2,1): frontier exhausted at t = 4 = S+1, i.e. new states through
depth S; 3x3 (2,1): same, 483,446 states, depth 3 = S -- asserted here).

STATUS. This is the formulation step done: an intrinsic, bounded-memory
predicate, exact at eight 2x2 cells and validated at 3x3 (2,1). Open for the
lattice-free proof: show F^{S+1}(diag) is contained in union F^{<=S}(diag)
(the window saturates at S) on any graph at tau_a >= tau_p -- a statement
with the same shape as the healing bound, likely the same young-wave witness
family at the boundary.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts the window
= S at all eight 2x2 cells, |R| at each, the separation (276/0 at (2,1)),
and the 3x3 depth; aborts on regression.
Output: result/topology/coherence_window_S.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
# Every cell above has tau_a >= tau_p. The window = S law is a REGIME law and
# does not hold outside it -- asserted below as a negative control so the
# scope cannot rot. (1,2) is a coincidence at 2x2 and dies at 3x3 (window 10
# vs S 3), so it is excluded from the negative set.
CELLS_NEG = [(2, 3), (2, 4), (3, 4), (3, 5), (4, 5)]
WINDOW_NEG = {(2, 3): 9, (2, 4): 12, (3, 4): 15, (3, 5): 17, (4, 5): 20}
R_EXPECT = {(1, 1): 48, (2, 1): 276, (2, 2): 360, (3, 1): 730,
            (3, 2): 1194, (3, 3): 1344, (4, 3): 3400, (4, 4): 3600}
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 500_000


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


def bfs_from_diag(stp, B, lv):
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


def separation_check_21(depth):
    """At (2,1): window-S membership must accept all of R and reject all 192
    lockstep impostors."""
    stp, B = make_step(2, 1)
    lv = live_set(stp, B)

    def lockstep(v, u, K=60):
        for _ in range(K):
            if u == tuple((x + 1) % B for x in v):
                return True
            v, u = stp(v), stp(u)
        return False
    P = {(v, u) for v in sorted(lv) for u in sorted(lv) if lockstep(v, u)}
    Rp = set(depth)
    extras = P - Rp
    return len(Rp), len(extras & Rp), len(extras)


def depth_3x3_21():
    ta, tp = 2, 1
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
    # BFS on pair codes (v, u), both as int64
    seen = set()
    V, U = live_idx.copy(), sh[live_idx]
    code = V * N + U
    seen.update(code.tolist())
    d = 0
    n_total = len(seen)
    while True:
        V, U = f[V], f[U]
        code = V * N + U
        mask = np.fromiter((c not in seen for c in code.tolist()),
                           dtype=bool, count=len(code))
        if not mask.any():
            break
        d += 1
        seen.update(code[mask].tolist())
        V, U = V[mask], U[mask]
    return len(live_idx), len(seen), d


def main():
    print("=== window w (BFS depth from the diagonal) vs S, all 2x2 cells ===")
    ws, sizes = [], []
    for ta, tp in CELLS:
        S = ta + tp
        stp, B = make_step(ta, tp)
        lv = live_set(stp, B)
        depth = bfs_from_diag(stp, B, lv)
        w = max(depth.values())
        ws.append(w)
        sizes.append(len(depth))
        print(f"({ta},{tp}) S={S}: |R|={len(depth)} window={w}", flush=True)
        assert w == S, f"window != S at ({ta},{tp})"
        assert len(depth) == R_EXPECT[(ta, tp)], f"|R| changed at ({ta},{tp})"
        if (ta, tp) == (2, 1):
            nR, fp, nex = separation_check_21(depth)
            print(f"   separation at (2,1): |R|={nR}, impostors={nex}, "
                  f"impostors inside window-S set: {fp}", flush=True)
            assert (nR, fp, nex) == (276, 0, 192), "separation changed"

    print("\n=== negative control: tau_p > tau_a, where window != S ===")
    ws_neg = []
    for ta, tp in CELLS_NEG:
        S = ta + tp
        stp, B = make_step(ta, tp)
        depth = bfs_from_diag(stp, B, live_set(stp, B))
        w = max(depth.values())
        ws_neg.append(w)
        print(f"({ta},{tp}) S={S}: window={w}  -- window != S, as expected",
              flush=True)
        assert w == WINDOW_NEG[(ta, tp)], f"neg-control window moved at ({ta},{tp})"
        assert w != S, f"window == S at tau_p > tau_a ({ta},{tp}) -- the regime law would be WIDER than claimed; re-derive before celebrating"

    print("\n=== 3x3 (2,1): window at scale ===")
    n_live, n_R, d33 = depth_3x3_21()
    print(f"live={n_live:,} |R|={n_R:,} window={d33} (S=3)", flush=True)
    assert (n_live, n_R, d33) == (255_484, 483_446, 3), "3x3 window changed"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_window_S.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             ta_neg=np.array([c[0] for c in CELLS_NEG]),
             tp_neg=np.array([c[1] for c in CELLS_NEG]),
             window_neg=np.array(ws_neg),
             ta=np.array([c[0] for c in CELLS]),
             tp=np.array([c[1] for c in CELLS]),
             window=np.array(ws), n_R=np.array(sizes),
             d33=np.array([d33]), n_R_33=np.array([n_R]))
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | S | |R| | window |")
    print("| :---: | :---: | ---: | ---: | ---: |")
    for (ta, tp), w, n in zip(CELLS, ws, sizes):
        print(f"| 2×2 | ({ta}, {tp}) | {ta+tp} | {n:,} | {w} |")
    print(f"| 3×3 | (2, 1) | 3 | {n_R:,} | {d33} |")


if __name__ == "__main__":
    main()
