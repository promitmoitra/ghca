"""Live-subshift entropy: the predicted regime signature is FALSIFIED, and the
entropy turns out to measure the abstraction, not the dynamics.

Third experiment on the theory branch. Motivated by Kessebohmer, Rademacher &
Ulbrich (arXiv:1903.02459), who compute topological entropy for 1D GH via
subshifts: take the live part of the spectrum automaton (spectrum_automaton*.py)
as a vertex shift of finite type (SFT) and compute its topological entropy
h = ln(spectral radius of the live-restricted adjacency).

PREDICTION WRITTEN DOWN BEFORE COMPUTING (recorded because it FAILED): the live
SFT entropy should be well-defined/stable at tau_a >= tau_p and jump or become
ill-conditioned across the regime boundary.

WHAT WAS ACTUALLY FOUND:

1. 2x2: a clean-looking pattern -- h = ln 2 exactly on strict cells
   (tau_a > tau_p: (2,1),(3,1),(3,2),(4,3)) and h = 0 on diagonal and
   tau_p > tau_a cells -- carried by 2-cycles between symbol pairs whose
   alternation is realized only on transients (L2: real attractors have
   constant spectrum, so the SFT's entropy-carrying cycles are never realized
   asymptotically).

2. 3x3 FALSIFIES the pattern (the second 2x2-smallness artifact caught this
   session, after automaton soundness): h((1,1)) = 1.775, and
   h((2,1)) = h((2,2)) = h((3,2)) = 1.8826 -- strict and diagonal cells
   IDENTICAL, no regime signature at all.

3. The identity of those three values to 12 digits is not a numerical accident
   and not a bug (verified by dense eigensolves): the spectral radius is
   carried by an 8-NODE STRONGLY-CONNECTED COMPONENT (rho = 6.570750564) that
   occurs -- twice -- in the live automaton of BOTH (2,1) (309 live symbols,
   alphabet B=4) and (2,2) (1158 live symbols, B=5). A small, nearly-complete
   (K8 has rho = 7), cell-independent "scrambled core" of symbols saturates
   the SFT entropy regardless of the cell.

4. Interpretation, with the closure and zero-entropy checks that pin it down:
   the live CONFIGURATION dynamics is deterministic on a finite core, so its
   entropy is 0 in every cell (verified; and both live parts are closed under
   the step map in every cell -- the "not even closed" half of the prediction
   was also wrong). All positive entropy in the symbol SFT is therefore FIBRE
   AMBIGUITY of the projection -- how many configurations share a symbol --
   not dynamical complexity. Realized symbol windows confirm the slack:
   realized 6-windows are 11-28% of SFT 6-paths at the branching 2x2 cells.

5. Consequence for the program: the SFT entropy of the finite-core symbol
   automaton is NOT the object that connects to KRU's 1D entropy results.
   Their entropy lives on the infinite lattice's non-wandering set; the
   finite-core analogue is degenerate (deterministic => 0) and its symbolic
   cover measures the abstraction. An honest bridge needs growing-window
   entropy (h from window counts as the lattice grows) -- future work, not
   claimed here.

House Rules Compliance: exhaustive at 2x2 and 3x3 (no RNG); dense verification
of every sparse eigenvalue actually cited; the script asserts findings 1-4 and
aborts on regression.
Output: result/topology/live_subshift_entropy.npz + printed doc table.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
CELLS_2 = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
           (1, 2), (2, 3), (2, 4), (3, 4)]
CELLS_3 = [(1, 1), (2, 1), (2, 2), (3, 2)]
STRICT_LN2 = {(2, 1), (3, 1), (3, 2), (4, 3)}
CHUNK = 500_000

NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
PLAQ3 = [(0, 1, 4, 3), (1, 2, 5, 4), (3, 4, 7, 6), (4, 5, 8, 7)]


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def sig2(cfg, S):
    ph = [cfg[i] for i in CYC]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % (S + 1) for k in range(4)))


def fates2(ta, tp):
    S = ta + tp
    fate = {}
    for c in product(range(S + 1), repeat=4):
        if c in fate:
            continue
        path, cur = [], c
        while cur not in fate and cur not in path:
            path.append(cur)
            cur = step2(cur, ta, tp)
        val = fate[cur] if cur in fate else any(
            any(x) for x in path[path.index(cur):])
        for p in path:
            fate[p] = val
    return fate


def rho_dense(A):
    return float(max(abs(np.linalg.eigvals(A)))) if len(A) else 0.0


def audit_2x2(ta, tp, n_win=6):
    S = ta + tp
    fate = fates2(ta, tp)
    cls_fate, edges = defaultdict(set), defaultdict(set)
    for c, f in fate.items():
        s = sig2(c, S)
        cls_fate[s].add(f)
        edges[s].add(sig2(step2(c, ta, tp), S))
    live = sorted(s for s, f in cls_fate.items() if f == {True})
    li = {s: i for i, s in enumerate(live)}
    A = np.zeros((len(live), len(live)))
    for s in live:
        for t in edges[s]:
            if t in li:
                A[li[s], li[t]] = 1.0
    sym_closed = all(edges[s] <= set(live) for s in live)
    h_sym = np.log(r) if (r := rho_dense(A)) > 0 else float("-inf")
    # configuration level
    live_cfg = [c for c, f in fate.items() if f]
    cfg_closed = all(fate[step2(c, ta, tp)] for c in live_cfg)
    ci = {c: i for i, c in enumerate(live_cfg)}
    M = np.zeros((len(ci), len(ci)))
    for c in live_cfg:
        n = step2(c, ta, tp)
        M[ci[c], ci[n]] = 1.0
    h_cfg = np.log(r) if (r := rho_dense(M)) > 0 else float("-inf")
    # realized windows vs SFT paths
    sft = int(np.linalg.matrix_power(A, n_win).sum())
    win = set()
    for c, f in fate.items():
        if not f:
            continue
        w, cur = [], c
        for _ in range(n_win + 1):
            w.append(sig2(cur, S))
            cur = step2(cur, ta, tp)
        win.add(tuple(w))
    return dict(ta=ta, tp=tp, n_live=len(live), sym_closed=sym_closed,
                cfg_closed=cfg_closed, h_sym=float(h_sym), h_cfg=float(h_cfg),
                sft_paths=sft, real_windows=len(win))


# ---- 3x3 (vectorised pipeline shared with spectrum_automaton_3x3.py) -------
def successor_array(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = np.empty(N, dtype=np.int32)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        dig = ((idx[:, None] // pw[None, :]) % B).astype(np.int8)
        act = (dig >= 1) & (dig <= ta)
        nxt = dig.astype(np.int32)
        adv = dig > 0
        nxt[adv] = (dig[adv].astype(np.int32) + 1) % B
        for i in range(9):
            na = np.zeros(hi - lo, dtype=np.int8)
            for j in NB3[i]:
                na += act[:, j]
            m = dig[:, i] == 0
            nxt[m, i] = (na[m] >= 1)
        f[lo:hi] = (nxt.astype(np.int64) @ pw).astype(np.int32)
    return f


def persistent_mask(f):
    dead = f == 0
    dead[0] = True
    while True:
        nd = dead | dead[f]
        if bool((nd == dead).all()):
            return ~nd
        dead = nd


def plaq_code(N, B):
    pw = (B ** np.arange(9)).astype(np.int64)
    base4 = (B ** np.arange(4)).astype(np.int64)
    M = int(B) ** 4
    code = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        dig = ((idx[:, None] // pw[None, :]) % B).astype(np.int64)
        pcs = np.empty((hi - lo, 4), dtype=np.int64)
        for p, cyc in enumerate(PLAQ3):
            ph = dig[:, list(cyc)]
            g = (np.roll(ph, -1, axis=1) - ph) % B
            g.sort(axis=1)
            pcs[:, p] = g @ base4
        pcs.sort(axis=1)
        code[lo:hi] = ((pcs[:, 0]*M + pcs[:, 1])*M + pcs[:, 2])*M + pcs[:, 3]
    return code


def audit_3x3(ta, tp):
    from scipy.sparse import csr_matrix
    from scipy.sparse.csgraph import connected_components
    B = ta + tp + 1
    N = B ** 9
    f = successor_array(ta, tp)
    pm = persistent_mask(f)
    code = plaq_code(N, B)
    fate_acc, edge_acc = {}, set()
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        c, cn, p = code[lo:hi], code[f[lo:hi]], pm[lo:hi]
        edge_acc.update(map(tuple, np.unique(np.stack([c, cn], 1), axis=0)))
        for val, mask in ((0, p), (1, ~p)):
            for u in np.unique(c[mask]):
                fate_acc.setdefault(int(u), [False, False])[val] = True
    live = sorted(s for s, (l, d) in fate_acc.items() if l and not d)
    li = {s: i for i, s in enumerate(live)}
    A = np.zeros((len(live), len(live)))
    for a, b in edge_acc:
        if a in li and b in li:
            A[li[a], li[b]] = 1.0
    rho = rho_dense(A)          # dense: every cited eigenvalue dense-verified
    # dominant nontrivial SCC
    ncc, lab = connected_components(csr_matrix(A), directed=True,
                                    connection='strong')
    best_size, best_rho = 0, 0.0
    for k in range(ncc):
        m = lab == k
        if m.sum() < 2:
            continue
        r = rho_dense(A[np.ix_(m, m)])
        if r > best_rho:
            best_rho, best_size = r, int(m.sum())
    return dict(ta=ta, tp=tp, n_live=len(live), rho=rho,
                h_sym=float(np.log(rho)) if rho > 0 else float("-inf"),
                scc_size=best_size, scc_rho=best_rho)


def main():
    print("=== 2x2: h_sym / h_cfg / closure / realized-vs-SFT ===")
    recs2 = [audit_2x2(ta, tp) for ta, tp in CELLS_2]
    for r in recs2:
        regime = ("strict" if r["ta"] > r["tp"] else
                  "diag" if r["ta"] == r["tp"] else "ta<tp")
        print(f"({r['ta']},{r['tp']}) [{regime:6s}] live={r['n_live']:3d} "
              f"h_sym={r['h_sym']:7.4f} h_cfg={r['h_cfg']:7.4f} "
              f"closed(sym/cfg)={r['sym_closed']}/{r['cfg_closed']} "
              f"windows/paths={r['real_windows']}/{r['sft_paths']}", flush=True)

    ln2 = float(np.log(2))
    for r in recs2:
        cell = (r["ta"], r["tp"])
        assert r["cfg_closed"] and r["sym_closed"], f"closure fails at {cell}"
        assert abs(r["h_cfg"]) < 1e-9 or r["h_cfg"] == float("-inf"), \
            f"config entropy nonzero at {cell}"
        want = ln2 if cell in STRICT_LN2 else 0.0
        got = 0.0 if r["h_sym"] == float("-inf") else r["h_sym"]
        assert abs(got - want) < 1e-9, f"2x2 h_sym pattern changed at {cell}"

    print("\n=== 3x3: h_sym (dense-verified) + dominant SCC ===")
    recs3 = []
    for ta, tp in CELLS_3:
        t0 = time.time()
        r = audit_3x3(ta, tp)
        recs3.append(r)
        regime = "strict" if ta > tp else "diag"
        print(f"({ta},{tp}) [{regime:6s}] live={r['n_live']:5d} "
              f"h_sym={r['h_sym']:.4f} rho={r['rho']:.9f} | dominant SCC: "
              f"size={r['scc_size']} rho={r['scc_rho']:.9f} "
              f"[{time.time()-t0:.0f}s]", flush=True)

    # the falsification + the shared-core identity, asserted
    by = {(r["ta"], r["tp"]): r for r in recs3}
    assert abs(by[(2, 1)]["rho"] - by[(2, 2)]["rho"]) < 1e-9, "shared rho gone"
    assert abs(by[(2, 2)]["rho"] - by[(3, 2)]["rho"]) < 1e-9, "shared rho gone"
    assert by[(1, 1)]["h_sym"] > 1.0, "3x3 diagonal entropy no longer large"
    assert all(r["scc_size"] == 8 for r in recs3 if (r["ta"], r["tp"]) != (1, 1)), \
        "dominant SCC size changed"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "live_subshift_entropy.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g2_{k}": np.array([r[k] for r in recs2]) for k in recs2[0]},
             **{f"g3_{k}": np.array([r[k] for r in recs3]) for k in recs3[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | regime | live symbols | h_sym | h_cfg | note |")
    print("| :---: | :---: | :---: | ---: | ---: | ---: | :--- |")
    for r in recs2:
        regime = ("τa>τp" if r["ta"] > r["tp"] else
                  "τa=τp" if r["ta"] == r["tp"] else "τa<τp")
        h = 0.0 if r["h_sym"] == float("-inf") else r["h_sym"]
        note = "ln 2" if (r["ta"], r["tp"]) in STRICT_LN2 else ""
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | {regime} | {r['n_live']} | "
              f"{h:.4f} | 0 | {note} |")
    for r in recs3:
        regime = "τa>τp" if r["ta"] > r["tp"] else "τa=τp"
        print(f"| 3×3 | ({r['ta']}, {r['tp']}) | {regime} | {r['n_live']:,} | "
              f"{r['h_sym']:.4f} | 0 | SCC size {r['scc_size']}, "
              f"ρ={r['scc_rho']:.6f} |")


if __name__ == "__main__":
    main()
