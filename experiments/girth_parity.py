"""Does the dwelling-attractor threshold track GIRTH? Prediction falsified at one
cell, then repaired: girth is exactly right ONCE Lemma E's obstruction is
factored out.

The open claim of coherent_core.py was "every LIVE attractor is dwell-free when
S+1 >= 4", certified only on square lattices, where 4 is the girth. Its honest
scope named the discriminator: a lattice of different girth AND parity. This is
that experiment.

Setting: arbitrary graphs, theta = 1, states 0 (receptive), 1..ta (active),
ta+1..S (passive), S = ta+tp, B = S+1. Nine graphs chosen to separate girth from
parity, every one computed (girth, bipartiteness and full cycle spectrum are
derived here, not asserted):

    C3..C8 rings   girth = m, bipartite iff m even -- the clean instrument,
                   since a ring's girth is its length and nothing else varies
    square 3x3     girth 4, bipartite      (the baseline coherent_core used)
    triangular 3x3 girth 3, NON-bipartite, cycle spectrum {3..9}
    honeycomb      girth 6, bipartite, cycle spectrum {6, 10} (two fused hexagons)

---------------------------------------------------------------------------
PREDICTION, recorded before computing (this is the house rule that makes a
falsification legible):

    H_g:  a live attractor dwells  <=>  B < girth.

RESULT: H_g is FALSIFIED -- 79/80 (graph, cell) rows, one failure:
square 3x3 at (ta,tp) = (1,3), where B = 5 >= girth 4 yet 28 of the 28 live
attractors dwell. The failure is diagnostic rather than fatal: at (1,3) the
coherent set is EMPTY for the arithmetic reason of Lemma E (ta = 1 forces the
lag sum to equal the cycle length, and no multiple of B = 5 lies in {4, 6, 8}),
so by Theorem C no dwell-free attractor can exist at all, girth notwithstanding.

REPAIRED HYPOTHESIS, certified 80/80:

    H_g':  a dwelling live attractor exists
           <=>  live attractors exist  AND  (B < girth  OR  C is empty).

which decomposes into one triviality and two substantive halves:

  A (trivial, immediate from Theorem C): C empty and some live attractor exists
    => every live attractor dwells.                        0 violations.
  B (SUBSTANTIVE): C non-empty and B >= girth => NO live attractor dwells.
    0 violations over 39 rows in scope.
  C (SUBSTANTIVE): C non-empty and B <  girth => SOME live attractor dwells.
    0 violations over 21 rows in scope.

So the girth reading is RESCUED, correctly conditioned. coherent_core.md's
verdict that "the obstruction is arithmetic, not metric" was half the story:
Lemma E's arithmetic decides whether dwell-free attractors can exist AT ALL,
and given that they can, GIRTH decides whether dwelling ones coexist. Two
different obstructions, previously conflated because the square lattice at the
tested cells never separated them.

---------------------------------------------------------------------------
THE DISCRIMINATING COMPARISON. Girth is doing real work, and one pair shows it
with everything else held fixed -- same cell (2,1), same B = 4, both bipartite:

    square 3x3  girth 4:  B >= girth  ->  12,049 live,     0 dwelling
    honeycomb   girth 6:  B <  girth  ->  20,401 live,    10 dwelling

Same dynamics, same parameters, different girth, opposite answer. And the
triangular lattice (girth 3, so B >= girth ALWAYS since S >= 2) has ZERO
dwelling attractors at every one of its six affordable cells -- the anomaly
that made (1,1) special on the square lattice simply does not exist there.

LEMMA E EXTENDED. coherent_core.py proved Lemma E for any graph but only
checked it on square lattices, both bipartite. It is re-checked here on
NON-bipartite graphs (C3, C5, C7, triangular) and on girths 3, 5, 6, 7, 8:
0 necessity violations over 80 rows, and the converse held on all 80 as well
(still reported, not claimed).

House Rules Compliance:
  - Exhaustive and deterministic throughout; no RNG, nothing to seed.
  - Graph properties (girth, bipartiteness, cycle spectrum) are COMPUTED, so a
    mis-specified graph cannot silently pass as its intended shape.
  - Substrate / analysis boundary: the step map, its cycles and dwell events are
    substrate; C, girth and the cycle spectrum are analysis constructs imposed
    on the enumeration.
  - The falsified prediction H_g is asserted AS falsified (exactly one failure,
    at the named cell), so a future change that silently "fixes" it trips.
Output: result/topology/girth_parity.npz + printed doc tables.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from collections import deque

CELLS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 1), (1, 3),
         (3, 2), (2, 3), (4, 1), (1, 4), (3, 3)]
HGP = "H_g'"          # label; avoids a quote clash inside f-strings
STATE_CAP = 4_000_000          # exhaustive census budget, in configurations
SETTLE = 128                   # iterations to land on a cycle; asserted below
H_G_EXPECTED_FAILURES = [("square 3x3", 1, 3)]


# --------------------------------------------------------------------------
# graphs, and their measured properties
# --------------------------------------------------------------------------
def ring(m):
    return tuple(tuple(sorted(((i - 1) % m, (i + 1) % m))) for i in range(m))


def box(L, offsets):
    """L x L patch with the given coordinate offsets, open boundary."""
    adj = []
    for i in range(L):
        for j in range(L):
            a = []
            for di, dj in offsets:
                x, y = i + di, j + dj
                if 0 <= x < L and 0 <= y < L:
                    a.append(x * L + y)
            adj.append(tuple(sorted(a)))
    return tuple(adj)


def from_edges(n, edges):
    adj = [[] for _ in range(n)]
    for u, v in edges:
        adj[u].append(v)
        adj[v].append(u)
    return tuple(tuple(sorted(a)) for a in adj)


SQUARE_OFF = ((-1, 0), (1, 0), (0, -1), (0, 1))
TRI_OFF = SQUARE_OFF + ((1, -1), (-1, 1))
# two hexagons fused along the edge (2,3)
HONEY_E = [(0, 1), (1, 2), (2, 3), (3, 4), (4, 5), (5, 0),
           (3, 6), (6, 7), (7, 8), (8, 9), (9, 2)]

GRAPHS = {f"C{m} ring": ring(m) for m in range(3, 9)}
GRAPHS["square 3x3"] = box(3, SQUARE_OFF)
GRAPHS["triangular 3x3"] = box(3, TRI_OFF)
GRAPHS["honeycomb 2-hex"] = from_edges(10, HONEY_E)


def girth_and_parity(adj):
    """BFS girth and 2-colourability. Both measured, never assumed."""
    n = len(adj)
    g = 10 ** 9
    for s in range(n):
        dist = [-1] * n
        par = [-1] * n
        dist[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if dist[v] == -1:
                    dist[v], par[v] = dist[u] + 1, u
                    q.append(v)
                elif v != par[u]:
                    g = min(g, dist[u] + dist[v] + 1)
    col = [-1] * n
    bip = True
    for s in range(n):
        if col[s] != -1:
            continue
        col[s] = 0
        q = deque([s])
        while q:
            u = q.popleft()
            for v in adj[u]:
                if col[v] == -1:
                    col[v] = 1 - col[u]
                    q.append(v)
                elif col[v] == col[u]:
                    bip = False
    return g, bip


def cycle_lengths(adj):
    lengths = set()

    def walk(s, u, vis, d):
        for v in adj[u]:
            if v == s and d >= 3:
                lengths.add(d)
            elif v > s and v not in vis:
                walk(s, v, vis | {v}, d + 1)

    for s in range(len(adj)):
        walk(s, s, {s}, 1)
    return lengths


def emptiness_allows_C(ta, tp, lengths):
    """LEMMA E (proven in coherent_core.py, any graph): C non-empty requires a
    cycle length l and k >= 1 with l <= k*B <= l*ta."""
    B = ta + tp + 1
    for l in lengths:
        k = -(-l // B)
        if k >= 1 and k * B <= l * ta:
            return True
    return False


# --------------------------------------------------------------------------
def pad_adj(adj):
    n = len(adj)
    d = max(len(a) for a in adj)
    nbi = np.zeros((n, d), np.int64)
    nbm = np.zeros((n, d), bool)
    for i, a in enumerate(adj):
        for k, j in enumerate(a):
            nbi[i, k], nbm[i, k] = j, True
    return nbi, nbm


def census(adj, ta, tp, cap=STATE_CAP):
    """Exhaustive attractor census on an arbitrary graph. None if unaffordable."""
    n, B = len(adj), ta + tp + 1
    if B ** n > cap:
        return None
    nbi, nbm = pad_adj(adj)
    N = B ** n
    D = np.empty((N, n), np.int8)
    idx = np.arange(N, dtype=np.int64)
    for k in range(n - 1, -1, -1):
        D[:, k] = idx % B
        idx //= B

    act = (D >= 1) & (D <= ta)
    fires = (act[:, nbi] & nbm[None, :, :]).any(axis=2)
    nxt = np.where(D == 0, fires.astype(D.dtype), (D + 1) % B)
    pw = (B ** np.arange(n - 1, -1, -1)).astype(np.int64)
    succ = nxt.astype(np.int64) @ pw

    land = np.arange(N, dtype=np.int64)
    for _ in range(SETTLE):
        land = succ[land]
    cyc = np.unique(land)
    assert np.array_equal(np.unique(succ[cyc]), cyc), "SETTLE insufficient"

    seen = np.zeros(N, bool)
    n_live = n_df = n_dead = 0
    periods = set()
    for c0 in cyc:
        if seen[c0]:
            continue
        orb = [int(c0)]
        seen[c0] = True
        x = succ[c0]
        while x != c0:
            orb.append(int(x))
            seen[x] = True
            x = succ[x]
        blk = D[orb]
        alive = bool(((blk >= 1) & (blk <= ta)).any())
        dwellfree = not bool(((blk == 0) & (nxt[orb] == 0)).any())
        if alive:
            n_live += 1
            periods.add(len(orb))
            if dwellfree:
                n_df += 1
        else:
            n_dead += 1
    assert n_dead == 1, f"Theorem Z broke: {n_dead} dead attractors"

    lag = (D[:, nbi].astype(np.int16) - D[:, :, None]) % B
    good = (lag >= 1) & (lag <= ta) & nbm[None, :, :]
    nC = int(good.any(axis=2).all(axis=1).sum())
    # Theorem C, re-verified on every graph here (not just square lattices)
    assert nC % B == 0 and nC // B == n_df, \
        f"Theorem C broke: |C|={nC}, B={B}, dwell-free={n_df}"

    return dict(nC=nC, live=n_live, dwellfree=n_df, dwelling=n_live - n_df,
                B=B, live_periods=sorted(periods))


# --------------------------------------------------------------------------
def main():
    print("=== the graphs (girth, parity and cycle spectrum all MEASURED) ===",
          flush=True)
    print(f"{'graph':>16} {'n':>3} {'girth':>6} {'bipartite':>10} "
          f"{'cycle lengths':>22}", flush=True)
    props = {}
    for name, adj in GRAPHS.items():
        g, bip = girth_and_parity(adj)
        lens = cycle_lengths(adj)
        props[name] = (g, bip, lens)
        print(f"{name:>16} {len(adj):>3} {g:>6} {str(bip):>10} "
              f"{str(sorted(lens)):>22}", flush=True)
    assert props["triangular 3x3"][0] == 3 and not props["triangular 3x3"][1], \
        "triangular patch is not girth-3 non-bipartite"
    assert props["honeycomb 2-hex"][0] == 6 and props["honeycomb 2-hex"][1], \
        "honeycomb patch is not girth-6 bipartite"
    assert props["square 3x3"][0] == 4, "square patch is not girth-4"

    print("\n=== census (exhaustive) ===", flush=True)
    print("H_g PREDICTION, recorded before computing: a live attractor dwells "
          "<=> B < girth", flush=True)
    print(f"{'graph':>16} {'g':>2} {'cell':>7} {'B':>3} {'|C|':>9} {'live':>8} "
          f"{'dwell-free':>10} {'DWELLING':>9} {'H_g':>5} {HGP:>5}", flush=True)
    rows = []
    for name, adj in GRAPHS.items():
        g, bip, lens = props[name]
        for ta, tp in CELLS:
            r = census(adj, ta, tp)
            if r is None:
                continue
            B = r["B"]
            hg_ok = ((B < g) == (r["dwelling"] > 0))
            hgp = (r["live"] > 0) and ((B < g) or (r["nC"] == 0))
            hgp_ok = (hgp == (r["dwelling"] > 0))
            rows.append(dict(graph=name, girth=g, bip=int(bip), ta=ta, tp=tp,
                             B=B, nC=r["nC"], live=r["live"],
                             dwellfree=r["dwellfree"], dwelling=r["dwelling"],
                             lemE=int(emptiness_allows_C(ta, tp, lens)),
                             hg_ok=int(hg_ok), hgp_ok=int(hgp_ok)))
            print(f"{name:>16} {g:>2} {'(' + str(ta) + ',' + str(tp) + ')':>7} "
                  f"{B:>3} {r['nC']:>9,} {r['live']:>8,} {r['dwellfree']:>10,} "
                  f"{r['dwelling']:>9,} {'OK' if hg_ok else 'FAIL':>5} "
                  f"{'OK' if hgp_ok else 'FAIL':>5}", flush=True)

    # ---- H_g is falsified, at exactly one named row
    hg_fail = [(r["graph"], r["ta"], r["tp"]) for r in rows if not r["hg_ok"]]
    print(f"\nH_g : {len(rows) - len(hg_fail)}/{len(rows)} consistent -- "
          f"FALSIFIED at {hg_fail}", flush=True)
    assert hg_fail == H_G_EXPECTED_FAILURES, \
        f"H_g failure set changed: {hg_fail} != {H_G_EXPECTED_FAILURES}"

    # ---- H_g' and its decomposition
    hgp_fail = [(r["graph"], r["ta"], r["tp"]) for r in rows if not r["hgp_ok"]]
    print(f"{HGP}: {len(rows) - len(hgp_fail)}/{len(rows)} consistent -- "
          f"failures: {hgp_fail or 'NONE'}", flush=True)
    assert not hgp_fail, f"H_g' broke at {hgp_fail}"

    scope_b = [r for r in rows if r["nC"] > 0 and r["B"] >= r["girth"]]
    scope_c = [r for r in rows if r["nC"] > 0 and r["B"] < r["girth"]]
    viol_a = [r for r in rows if r["nC"] == 0 and r["live"] > 0 and r["dwellfree"]]
    viol_b = [r for r in scope_b if r["dwelling"] > 0]
    viol_c = [r for r in scope_c if r["dwelling"] == 0]
    print(f"  A (trivial)     C empty & live>0 => all live dwell : "
          f"{len(viol_a)} violations", flush=True)
    print(f"  B (substantive) C non-empty & B >= girth => none dwell : "
          f"{len(viol_b)} violations over {len(scope_b)} rows", flush=True)
    print(f"  C (substantive) C non-empty & B <  girth => some dwell : "
          f"{len(viol_c)} violations over {len(scope_c)} rows", flush=True)
    assert not viol_a and not viol_b and not viol_c, "H_g' decomposition broke"
    assert len(scope_b) >= 30 and len(scope_c) >= 15, "sub-claims under-exercised"

    # ---- Lemma E, now on non-bipartite and higher-girth graphs
    lemE_viol = [r for r in rows if r["nC"] > 0 and not r["lemE"]]
    conv = all((r["nC"] > 0) == bool(r["lemE"]) for r in rows)
    nonbip = [r for r in rows if not r["bip"]]
    print(f"\nLEMMA E necessity: {len(lemE_viol)} violations over {len(rows)} "
          f"rows ({len(nonbip)} of them NON-bipartite); converse also held: {conv}",
          flush=True)
    assert not lemE_viol, "Lemma E necessity violated"
    assert len(nonbip) >= 20, "non-bipartite coverage too thin to claim extension"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "girth_parity.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{f"r_{k}": np.array([r[k] for r in rows]) for k in rows[0]},
             g_name=np.array(list(GRAPHS)),
             g_n=np.array([len(a) for a in GRAPHS.values()]),
             g_girth=np.array([props[k][0] for k in GRAPHS]),
             g_bip=np.array([int(props[k][1]) for k in GRAPHS]))
    print(f"\nwrote {out}", flush=True)

    print("\n--- results-doc table 1: the graphs (generated) ---")
    print("| graph | n | girth | bipartite | cycle lengths |")
    print("| :--- | ---: | ---: | :---: | :--- |")
    for name in GRAPHS:
        g, bip, lens = props[name]
        print(f"| {name} | {len(GRAPHS[name])} | {g} | {'yes' if bip else 'no'} "
              f"| {{{', '.join(map(str, sorted(lens)))}}} |")

    print("\n--- results-doc table 2: the discriminating comparison (generated) ---")
    print("| graph | girth | cell | B | B vs girth | live | dwelling |")
    print("| :--- | ---: | :---: | ---: | :---: | ---: | ---: |")
    for r in rows:
        if (r["ta"], r["tp"]) == (2, 1) and r["graph"] in (
                "square 3x3", "honeycomb 2-hex", "triangular 3x3"):
            rel = "≥" if r["B"] >= r["girth"] else "<"
            print(f"| {r['graph']} | {r['girth']} | ({r['ta']}, {r['tp']}) | "
                  f"{r['B']} | B {rel} girth | {r['live']:,} | {r['dwelling']:,} |")

    print("\n--- results-doc table 3: full census (generated) ---")
    print("| graph | girth | (τa, τp) | B | \\|C\\| | live | dwell-free | "
          "dwelling | H_g | H_g' |")
    print("| :--- | ---: | :---: | ---: | ---: | ---: | ---: | ---: | :---: | :---: |")
    for r in rows:
        print(f"| {r['graph']} | {r['girth']} | ({r['ta']}, {r['tp']}) | {r['B']} "
              f"| {r['nC']:,} | {r['live']:,} | {r['dwellfree']:,} | "
              f"{r['dwelling']:,} | {'✓' if r['hg_ok'] else '**✗**'} | "
              f"{'✓' if r['hgp_ok'] else '✗'} |")


if __name__ == "__main__":
    main()
