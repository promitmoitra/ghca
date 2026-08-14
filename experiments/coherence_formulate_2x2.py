"""Formulating the coherence invariant at 2x2, step 1: orbit co-membership is
falsified in both directions; +1-LOCKSTEP is proven necessary (all of R, within
4 steps) but not sufficient (468 lockstep pairs vs 276 coherent); the residual
192-pair gap is the invariant's exact target.

First experiment on the branch opened by coherence_invariant_scope.py
(workflow: formulate at 2x2, validate at 3x3, sample at 4x4). R = the
certified 276-state pair closure of Theorem 4 at (2,1), viewed as a set of
ordered pairs (v, u) -- Delta is determined by (v, u) on R (0 ambiguities).

HYPOTHESES TESTED (oracle-refereed):

  H1 u in orbit(v) + 1 (uniform translate of v's own forward orbit):
     192/276 -- FALSIFIED forward; converse also fails badly (816/1008
     orbit-translate pairs are NOT in R). Orbit co-membership with a
     clock offset is neither necessary nor sufficient.
  H2 u in orbit(v+1): 204/276 -- FALSIFIED.
  H3 u in orbit(v) + m for any m: 204/276 -- FALSIFIED.
  H5 +1-LOCKSTEP: exists k with step^k(u) = step^k(v) + 1 (uniformly).
     Holds for ALL 276 pairs of R with max k = 4 -- NECESSARY. Among all
     live x live pairs, exactly 468 reach lockstep: R is a strict subset
     (192 extras). NOT SUFFICIENT. The extras include ledger shapes that
     also occur inside R (e.g. one cell at +2), so no shape predicate
     separates them: the missing information is which lockstep pairs are
     reachable from an exact clock-shift start (Delta_0 = 0).

ON-ATTRACTOR LEDGER CENSUS (the persistent residue): 144 of R's states have
v on its attractor; only 84 have Delta = (0,0,0,0). The rest carry permanent
balanced residues (census asserted below) -- the transient bakes a nonzero
but anchored ledger into the eternal lockstep. So the invariant is NOT
"eventually Delta = 0"; it is "eventually lockstep, with an anchored residue
from a Delta_0 = 0 start".

STATE OF THE FORMULATION after step 1: coherent = lockstep-reachable AND
initialised at exact clock-shift. The open compression: an intrinsic
predicate on (v, u) alone -- no reference to the start -- separating the 276
from the 468. The 192 extras are the complete counterexample set to every
future candidate; they are archived.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts all
hypothesis counts, the lockstep max-k, the 468/276/192 split, and the
on-attractor census; aborts on regression.
Output: result/topology/coherence_formulate_2x2.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import Counter

TA, TP = 2, 1
S = TA + TP
B = S + 1
ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
LOCK_K = 60


def step2(cfg):
    return tuple((1 if sum(1 <= cfg[j] <= TA for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % B for i, s in enumerate(cfg))


def live_set():
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


def build_R(live):
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
    return R


def orbit(v):
    o, cur, seen = [], v, set()
    while cur not in seen:
        seen.add(cur)
        o.append(cur)
        cur = step2(cur)
    return o


def lockstep_k(v, u, K=LOCK_K):
    for k in range(K):
        if u == tuple((x + 1) % B for x in v):
            return k
        v, u = step2(v), step2(u)
    return None


def on_attractor(v):
    seen, cur = set(), v
    while cur not in seen:
        seen.add(cur)
        cur = step2(cur)
    cyc, x = {cur}, step2(cur)
    while x != cur:
        cyc.add(x)
        x = step2(x)
    return v in cyc


def main():
    live = live_set()
    R = build_R(live)
    assert len(R) == 276, "certificate size changed"
    pairs = {}
    for v, D in R:
        u = tuple((v[i] + 1 + D[i]) % B for i in range(4))
        assert pairs.get((v, u), D) == D, "Delta ambiguous on R"
        pairs[(v, u)] = D
    assert len(pairs) == 276

    print("=== orbit-translate hypotheses ===")
    orbs = {v: orbit(v) for v in live}
    h1 = sum(any(tuple((x + 1) % B for x in c) == u for c in orbs[v])
             for (v, u) in pairs)
    h2 = sum(u in orbs[tuple((x + 1) % B for x in v)] for (v, u) in pairs)
    h3 = sum(any(tuple((x + m) % B for x in c) == u
                 for c in orbs[v] for m in range(B)) for (v, u) in pairs)
    conv1 = sum((v, tuple((x + 1) % B for x in c)) not in pairs
                for v in sorted(live) for c in orbs[v])
    n_conv = sum(len(orbs[v]) for v in live)
    print(f"H1 {h1}/276  H2 {h2}/276  H3 {h3}/276  "
          f"H1-converse extras {conv1}/{n_conv}", flush=True)
    assert (h1, h2, h3) == (192, 204, 204), "hypothesis census changed"
    assert (conv1, n_conv) == (816, 1008), "converse census changed"

    print("\n=== +1-lockstep ===")
    ks = [lockstep_k(v, u) for (v, u) in pairs]
    assert all(k is not None for k in ks), "an R pair fails lockstep"
    print(f"all 276 reach lockstep; max k = {max(ks)}", flush=True)
    assert max(ks) == 4, "lockstep horizon changed"
    P = {(v, u) for v in sorted(live) for u in sorted(live)
         if lockstep_k(v, u) is not None}
    extras = sorted(P - set(pairs))
    print(f"lockstep pairs among live x live: {len(P)}; extras: {len(extras)}",
          flush=True)
    assert len(P) == 468 and len(extras) == 192, "lockstep split changed"
    assert set(pairs) <= P

    print("\n=== on-attractor ledger census ===")
    att = [(v, D) for (v, D) in R if on_attractor(v)]
    cen = Counter(D for v, D in att)
    print(f"on-attractor states: {len(att)}; Delta=(0,0,0,0): "
          f"{cen[(0, 0, 0, 0)]}", flush=True)
    assert len(att) == 144 and cen[(0, 0, 0, 0)] == 84, "attractor census changed"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_formulate_2x2.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             h_counts=np.array([h1, h2, h3]),
             conv=np.array([conv1, n_conv]),
             lock=np.array([len(P), len(pairs), len(extras), max(ks)]),
             att_census_keys=np.array([str(k) for k in sorted(cen)]),
             att_census_vals=np.array([cen[k] for k in sorted(cen)]),
             extras_v=np.array([v for v, u in extras]),
             extras_u=np.array([u for v, u in extras]))
    print(f"\nwrote {out}")


if __name__ == "__main__":
    main()
