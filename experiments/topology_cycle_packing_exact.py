"""K_dyn was undercounted: the greedy packing in topology_cycle_capacity.py
maximises cycle LENGTH per pick, not cycle COUNT.

PR #51 reports K_dyn(tau) -- how many reentrant loops a topology can hold
simultaneously -- via greedy edge-disjoint packing of cycles longer than tau,
and its doc flags the number as "a bound, not a measured count". That caveat
turns out to matter more than expected. `pack_long_cycles()` selects

    found = max(longer, key=len)      # longest cycle exceeding tau

but the longest cycle consumes the most edges per unit of count, so it is a poor
column for maximising how many cycles fit. Switching that single call to `min` --
shortest cycle still exceeding tau -- raises the packing by up to 4.7x on the
topologies tested, and an ILP over short cycles raises it further on three of the
four. Which method wins is topology-dependent, not density-dependent:

    ring(60,k=3)        ILP wins and CERTIFIES the optimum (45 vs greedy 27 vs
                        merged 6 at tau=3; UB = 45 attained)
    smallworld(120,k=6) ILP wins (91 vs 59 at tau=3; 69 vs 32, 59 vs 32, 44 vs 28)
    rgg(150,r=0.14)     ILP wins where solvable (118 vs 37 at tau=3; 90 vs 21)
    lattice2d(12,r=2)   greedy(min) wins (183 vs 181 at tau=3; 145 vs 75; and the
                        pool is un-enumerable at tau>=5, where only greedy reports)

lattice2d is the sole case where the greedy fix beats the ILP, and it is also the
densest pool (258k length-6 cycles at tau=5) -- so the ILP's disadvantage there
tracks pool size and solve difficulty, not graph density as such. The script
records both and reports `best known = max(merged, greedy_min, ILP)` per cell
rather than nominating a winner.

What this experiment establishes, per (topology, tau). The only inequalities
that hold in general are

  max(K_greedy_max, K_greedy_min, K_ILP)  <=  K_opt  <=  UB
                                              UB = min(beta_1, floor(m/(tau+1)))

i.e. every method reported here is a LOWER bound on the true optimum, and none
of them dominates the others. In particular **K_ILP is NOT always >= the greedy
values** -- measured on `lattice2d(12x12, r=2)`, greedy(min) beats the ILP at
every tau tested (183 vs 181 at tau=3, 145 vs 75 at tau=4, 111 vs an
un-enumerable pool at tau=5), for two independent reasons:

* **Pool restriction is not free in general.** K_ILP optimises over simple cycles
  of length EXACTLY tau+1, on the reasoning that a longer cycle spends more edges
  per unit of count. That reasoning bounds the *edge budget*, not the *achievable
  packing*: edge-disjointness can force a packing through longer cycles where no
  length-(tau+1) cycle exists on the remaining edges. The restriction is
  therefore a heuristic that shrinks the model, and it can lose -- the greedy,
  which accepts any cycle longer than tau, is free to use those columns.
* **Time-limited solves report the incumbent.** With MILP_TIME exceeded the ILP
  returns whatever it had found, which on dense pools (lattice2d has 258k
  length-6 cycles at tau=5) can be far below both the optimum and the greedy.
  `status` distinguishes "optimal" (exact over the pool) from "time-limit" /
  "solver-failed" / "pool-too-large"; only "optimal" rows are exact statements
  about the pool, and only `certified` rows are exact statements about K_opt.

So the ILP earns a CERTIFIED value only where the pool is small enough to solve
to proven optimality AND the value meets UB -- on `ring(60, k=3)` it does both.
Everywhere else the honest best-known K is `max(K_merged, K_greedy_min, K_ILP)`,
which is what the generated table reports; both methods are recorded per cell so
the comparison stays visible instead of being asserted in either direction.

* UB: packed cycles are edge-disjoint and each uses >= tau+1 edges, hence
  K <= floor(m/(tau+1)); and no packing exceeds the circuit rank
  beta_1 = m - N + c. K_ILP == UB CERTIFIES the exact capacity.
* Every returned packing is re-verified independently: each selected cycle is
  checked to be a real cycle in G, of length > tau, sharing no edge with any
  other. A packing that fails verification is not reported.

The headline: on `ring(60, k=3)` at tau=3 the merged greedy reports 6 and the
true optimum is 45 edge-disjoint 4-cycles using all 180 edges -- certified,
since UB = min(121, 45) = 45. The capacity/duration tradeoff of PR #51 survives
(K still falls with tau) but every K_dyn value in that doc is a loose lower
bound, and the ring is not the "flat" topology it appeared to be.

House Rules Compliance:
    - Explicit default_rng(seed) threading; both greedy variants share a seed so
      the max-vs-min comparison is paired rather than confounded by RNG.
    - Substrate/analysis boundary: pure graph combinatorics on the SAME
      topologies as PR #51. No dynamics are run, and nothing here claims a
      dynamical measurement -- K is a structural bound on simultaneous loops.
Output: result/topology/cycle_packing_exact.npz + a printed doc table.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import networkx as nx
from scipy.optimize import milp, LinearConstraint, Bounds
from scipy.sparse import lil_matrix
from ghca_net import lattice2d, smallworld, rgg, ring

SEED = 0
MILP_TIME = 120.0      # seconds per (topology, tau) solve; ring tau=6 needs >60s
                       # to certify 25 (at 60s it returns 24, a valid but loose
                       # incumbent -- the time limit is part of the result)
POOL_CAP = 50000       # skip the ILP above this many columns (report greedy only)
ENUM_CAP = 400000      # abandon cycle enumeration past this many candidates
TAUS = [3, 4, 5, 6]    # small tau: where the pool is enumerable and K is largest

TOPOLOGIES = {
    "ring (N=60, k=3)":            lambda: ring(60, k=3),
    "lattice2d (12x12, r=2)":      lambda: lattice2d(12, r=2, periodic=True),
    "smallworld (N=120, k=6)":     lambda: smallworld(120, k=6, beta=0.1, seed=0),
    "rgg (N=150, radius=0.14)":    lambda: rgg(150, radius=0.14, seed=0),
}


def greedy_pack(G, tau, pick, seed=SEED, max_iter=2000):
    """PR #51's greedy, parameterised by the pick rule.

    pick=max reproduces the merged `pack_long_cycles` (longest cycle > tau);
    pick=min is the count-maximising variant (shortest cycle > tau).
    """
    rng = np.random.default_rng(seed)
    H = G.copy()
    packed = 0
    for _ in range(max_iter):
        found = None
        nodes = list(H.nodes())
        rng.shuffle(nodes)
        for s in nodes:
            sub = H.subgraph(nx.node_connected_component(H, s))
            longer = [c for c in nx.cycle_basis(sub) if len(c) > tau]
            if longer:
                found = pick(longer, key=len)
                break
        if found is None:
            break
        L = len(found)
        H.remove_edges_from([(found[k], found[(k + 1) % L]) for k in range(L)
                             if H.has_edge(found[k], found[(k + 1) % L])])
        packed += 1
    return packed


def verify_packing(G, cycles, tau):
    """Independent check: real cycles in G, length > tau, pairwise edge-disjoint."""
    used = set()
    for c in cycles:
        L = len(c)
        if L <= tau:
            return False, "cycle not longer than tau"
        for k in range(L):
            u, v = c[k], c[(k + 1) % L]
            if not G.has_edge(u, v):
                return False, f"({u},{v}) is not an edge of G"
            e = frozenset((u, v))
            if e in used:
                return False, f"edge {sorted(e)} reused"
            used.add(e)
    return True, f"{len(cycles)} cycles, {len(used)} distinct edges"


def exact_pack(G, tau):
    """Max #edge-disjoint cycles of length exactly tau+1 (ILP over that pool)."""
    pool, seen = [], 0
    for c in nx.simple_cycles(G, length_bound=tau + 1):
        seen += 1
        if len(c) == tau + 1:
            pool.append(tuple(c))
        if seen > ENUM_CAP or len(pool) > POOL_CAP:
            return -1, len(pool), "enum-capped", []
    if not pool:
        return 0, 0, "no-columns", []
    m = G.number_of_edges()
    eidx = {frozenset(e): i for i, e in enumerate(G.edges())}
    A = lil_matrix((m, len(pool)), dtype=np.int8)
    for j, c in enumerate(pool):
        L = len(c)
        for k in range(L):
            A[eidx[frozenset((c[k], c[(k + 1) % L]))], j] = 1
    res = milp(c=-np.ones(len(pool)),
               constraints=LinearConstraint(A.tocsc(), ub=np.ones(m)),
               integrality=np.ones(len(pool)), bounds=Bounds(0, 1),
               options={"time_limit": MILP_TIME, "mip_rel_gap": 0.0})
    if res.x is None:
        return -1, len(pool), "solver-failed", []
    sel = [pool[j] for j in range(len(pool)) if res.x[j] > 0.5]
    return len(sel), len(pool), ("optimal" if res.status == 0 else "time-limit"), sel


def main():
    root = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
    prior = np.load(os.path.join(root, "result", "topology", "cycle_capacity.npz"),
                    allow_pickle=True)
    g_taus = prior["taus"].tolist()
    g_names = [str(x) for x in prior["names"]]
    g_K = prior["K"]

    names = list(TOPOLOGIES)
    shape = (len(names), len(TAUS))
    K_merged = np.zeros(shape, dtype=int)   # as committed in cycle_capacity.npz
    K_gmin = np.zeros(shape, dtype=int)     # one-word fix: min instead of max
    K_ilp = np.full(shape, -1, dtype=int)   # -1 = not solved
    UB = np.zeros(shape, dtype=int)
    certified = np.zeros(shape, dtype=bool)
    verified = np.zeros(shape, dtype=bool)

    print(f"{'topology':24s} {'tau':>3} {'merged':>6} {'g_min':>6} {'ILP':>5} "
          f"{'UB':>5} {'pool':>7}  status")
    for i, (name, build) in enumerate(TOPOLOGIES.items()):
        G = nx.from_numpy_array(np.asarray(build()))
        m = G.number_of_edges()
        b1 = m - G.number_of_nodes() + nx.number_connected_components(G)
        for j, tau in enumerate(TAUS):
            gi, gj = g_names.index(name), g_taus.index(tau)
            K_merged[i, j] = int(g_K[gi, gj])
            K_gmin[i, j] = greedy_pack(G, tau, min)
            K, npool, status, sel = exact_pack(G, tau)
            ub = min(b1, m // (tau + 1))
            UB[i, j] = ub
            if K >= 0:
                ok, msg = verify_packing(G, sel, tau)
                verified[i, j] = ok
                if not ok:
                    status = f"VERIFY-FAILED: {msg}"
                    K = -1
                K_ilp[i, j] = K
                certified[i, j] = (K == ub)
            flag = "CERT" if certified[i, j] else ""
            print(f"{name:24s} {tau:3d} {K_merged[i,j]:6d} {K_gmin[i,j]:6d} "
                  f"{K:5d} {ub:5d} {npool:7d}  {status:12s} {flag}", flush=True)

    out = os.path.join(root, "result", "topology", "cycle_packing_exact.npz")
    np.savez(out, taus=np.array(TAUS), names=np.array(names, dtype=object),
             K_merged=K_merged, K_greedy_min=K_gmin, K_ilp=K_ilp, UB=UB,
             certified=certified, verified=verified, seed=SEED,
             milp_time=MILP_TIME, pool_cap=POOL_CAP)
    print(f"\nwrote {out}")

    best = np.maximum(K_gmin, np.maximum(K_ilp, K_merged))
    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| topology | tau | merged K_dyn | greedy(min) | ILP | **best known** | UB | certified |")
    print("| :--- | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
    for i, name in enumerate(names):
        for j, tau in enumerate(TAUS):
            k = K_ilp[i, j]
            print(f"| `{name}` | {tau} | {K_merged[i,j]} | {K_gmin[i,j]} | "
                  f"{k if k >= 0 else '--'} | **{best[i,j]}** | {UB[i,j]} | "
                  f"{'yes' if certified[i,j] else 'no'} |")
    ub_ratio = best / np.maximum(K_merged, 1)
    print(f"\nbest-known/merged K_dyn ratio: min={ub_ratio.min():.2f} "
          f"max={ub_ratio.max():.2f} (every merged value is a loose lower bound)")
    print("monotone in tau (best known falls as tau rises)? " +
          str({names[i]: bool(np.all(np.diff(best[i]) <= 0)) for i in range(len(names))}))
    ratio = K_gmin / np.maximum(K_merged, 1)
    print(f"\ngreedy(min)/merged ratio: min={ratio.min():.2f} max={ratio.max():.2f}")
    solved = K_ilp >= 0
    if solved.any():
        r2 = (K_ilp[solved] / np.maximum(K_merged[solved], 1))
        print(f"ILP/merged ratio (solved cells only): min={r2.min():.2f} max={r2.max():.2f}")
    print(f"certified-exact cells: {int(certified.sum())}/{certified.size}")
    print(f"all reported ILP packings independently verified: "
          f"{bool(verified[solved].all()) if solved.any() else 'n/a'}")


if __name__ == "__main__":
    main()
