"""The persistent set's hull in DIRECT phase-space coordinates: small, saturating,
necessary-but-not-sufficient.

Companion to persistent_set_structure.py, prompted by the observation that the
natural phase space of the system is the product of chains Omega = {0..S}^4 --
one axis per cell, a configuration is a point, the dynamics move the point.
That embedding has dimension 4 regardless of S (the one-hot embedding grows as
4S), so EXACT facet enumeration is tractable and the polytope question can be
answered directly instead of via separability proxies.

Findings, all exhaustive on the 2x2 core:

1. THE HULL IS SMALL AND ITS COMPLEXITY SATURATES. conv(persistent set) in R^4:

       (ta,tp)   pers pts   vertices   facets   dead pts INSIDE hull
       (1,1)         24        24        108        13
       (2,2)        200        40        202       229
       (3,3)        784        48        248      1131
       (4,4)       2160        48        242      3355

   Vertex count stalls near 48 and facet count near 240-250 while the point
   count grows by two orders of magnitude -- an O(1)-complexity outer
   description. Membership is a valid NECESSARY condition for persistence.
   (Facet counts are as reported by Qhull's triangulated output and can shift
   by a few between Qhull runs on different point orderings; VERTEX counts and
   the dead-inside counts are order-independent, and the latter are re-verified
   by exact LP. The committed archive is the authority.)

2. IT IS NOT SUFFICIENT, AND THE GAP GROWS. Dead configurations strictly inside
   the hull outnumber the persistent ones from (2,2) onward (3355 vs 2160 at
   (4,4)). The persistent set is therefore NOT the set of integer points of any
   polytope in these coordinates -- convexity is the wrong closure, consistent
   with the one-hot non-separability of persistent_set_structure.py seen from
   the other side. The sharp description remains the combinatorial gap
   signature.

3. AT (1,1) THE ATTRACTOR SET IS AN EXACT HYPERPLANE SLICE. Every attractor
   configuration satisfies x0+x1+x2+x3 = 3: the rotating wave conserves total
   phase, the attractor set is 3-dimensional inside R^4, and all 8 attractor
   points are vertices of the persistent hull. This is the cleanest polytope
   statement in the system -- and it is SPECIAL, not generic:

4. ABOVE (1,1) ATTRACTORS LEAVE THE HULL BOUNDARY. At (2,2), (3,3), (4,4) zero
   attractor points are hull vertices, attractor phase-sums span a whole range
   rather than one value, and the attractor set is full-dimensional. Attractors
   live in the interior of phase space.

5. NO LINEAR LYAPUNOV FUNCTION. Transient trajectories are not monotone in
   |phase_sum - attractor mean| (e.g. 112/1440 monotone at (4,4)): approach to
   the attractor is genuinely non-convex in these coordinates.

House Rules Compliance:
    - Deterministic and exhaustive; no RNG, nothing to seed, rerun bit-identical.
    - Substrate/analysis boundary: dynamics validated bit-exact against
      ghca_net.Network by the companion experiment; hulls, LPs and invariants
      here are analysis constructs. No claim the substrate represents them.
    - Qhull runs in floating point; all hull memberships are re-checked with an
      exact-feasibility LP (HiGHS), and vertex/facet counts are reported as
      computed. Coordinates are small integers, far from Qhull's precision
      limits, and (1,1)'s conservation law is verified in integer arithmetic.
Output: result/topology/phase_space_polytope.npz + a printed doc table.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from scipy.spatial import ConvexHull
from scipy.optimize import linprog

ADJ = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
THETA = 1
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 2), (4, 4)]


def step(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ[i]) >= THETA else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def orbit(cfg, ta, tp):
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step(cur, ta, tp)
    start = seen[cur]
    inv = {t: c for c, t in seen.items()}
    return [inv[t] for t in range(start)], [inv[t] for t in range(start, len(seen))]


def in_hull(x, P):
    """Exact-feasibility LP: is x a convex combination of rows of P?"""
    n = len(P)
    r = linprog(np.zeros(n), A_eq=np.vstack([P.T, np.ones(n)]),
                b_eq=np.append(x, 1.0), bounds=[(0, None)] * n, method="highs")
    return r.status == 0


def affine_dim(pts):
    X = np.asarray(pts, float)
    X = X - X.mean(0)
    return int((np.linalg.svd(X, compute_uv=False) > 1e-9).sum())


def main():
    recs = []
    print(f"{'(ta,tp)':>8} {'pers':>5} {'hV':>4} {'hF':>4} {'dead_in':>8} "
          f"{'attr':>5} {'adim':>5} {'attr_on_hull':>12} {'sums':>18} {'mono':>10}")
    for ta, tp in CELLS:
        S = ta + tp
        cfgs = list(product(range(S + 1), repeat=4))
        pers, dead, attr_pts, transients = [], [], set(), []
        for c in cfgs:
            tr, at = orbit(c, ta, tp)
            alive = any(any(x) for x in at)
            (pers if alive else dead).append(c)
            if alive:
                attr_pts.update(at)
                if tr:
                    transients.append((tr, at))

        P = np.array(pers, float)
        hull = ConvexHull(P)
        hv, hf = len(hull.vertices), len(hull.equations)
        Pv = P[hull.vertices]
        dead_in = sum(in_hull(np.array(d, float), Pv) for d in dead)

        A = sorted(attr_pts)
        adim = affine_dim(A)
        hull_verts = {tuple(P[v]) for v in hull.vertices}
        on_hull = sum(tuple(map(float, a)) in hull_verts for a in A)
        sums = sorted({sum(a) for a in A})

        # (1,1) conservation law in INTEGER arithmetic
        conserved = all(sum(a) == 3 for a in A) if (ta, tp) == (1, 1) else None

        mono = 0
        for tr, at in transients:
            target = np.mean([sum(a) for a in at])
            d = [abs(sum(c) - target) for c in tr] + [abs(sum(at[0]) - target)]
            mono += all(d[i + 1] <= d[i] + 1e-9 for i in range(len(d) - 1))

        recs.append(dict(ta=ta, tp=tp, n_pers=len(pers), n_dead=len(dead),
                         hull_v=hv, hull_f=hf, dead_inside=dead_in,
                         n_attr=len(A), attr_dim=adim, attr_on_hull=on_hull,
                         sum_min=min(sums), sum_max=max(sums), n_sums=len(sums),
                         mono=mono, n_transients=len(transients),
                         conserved_11=bool(conserved) if conserved is not None else False))
        r = recs[-1]
        srange = f"{r['sum_min']}..{r['sum_max']} ({r['n_sums']})"
        print(f"{str((ta,tp)):>8} {r['n_pers']:5d} {hv:4d} {hf:4d} {dead_in:8d} "
              f"{len(A):5d} {adim:5d} {on_hull:12d} {srange:>18} "
              f"{mono:4d}/{len(transients):<5d}", flush=True)

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "phase_space_polytope.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| (tau_a, tau_p) | persistent | hull vertices | facets | dead inside hull | "
          "attractor pts | attr dim | attr pts on hull | monotone transients |")
    print("| :---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: |")
    for r in recs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['n_pers']} | {r['hull_v']} | "
              f"{r['hull_f']} | {r['dead_inside']} | {r['n_attr']} | {r['attr_dim']} | "
              f"{r['attr_on_hull']} | {r['mono']}/{r['n_transients']} |")

    c11 = [r for r in recs if (r['ta'], r['tp']) == (1, 1)][0]
    print(f"\n(1,1) attractor conservation law x.1 == 3 (integer check): {c11['conserved_11']}")
    big = [r for r in recs if r['ta'] == r['tp'] and r['ta'] >= 2]
    print(f"attractor pts on hull for (ta=tp>=2): "
          f"{[r['attr_on_hull'] for r in big]} (all interior)")
    print(f"hull vertices across cells: {[r['hull_v'] for r in recs]} (saturating)")
    print(f"dead-inside exceeds persistent from (2,2) on: "
          f"{[(r['dead_inside'] > r['n_pers']) for r in recs if r['ta'] == r['tp']]}")


if __name__ == "__main__":
    main()
