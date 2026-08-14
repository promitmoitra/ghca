"""The spectrum automaton at 3x3 (exhaustive, 40M configs) and 4x4 (random
init): soundness is itself regime-dependent -- and 2x2 was too small to see it.

Scaling test of spectrum_automaton.py. Symbols are plaquette signatures
(persistent_set_3x3.py); the automaton has an edge s -> s' iff some
configuration with signature s steps to one with signature s'; z is the
all-zero configuration's symbol. Exhaustive at 3x3 (up to (3,3): 40,353,607
configs, 19,251 symbols, 277,036 edges); sampled at 4x4 (random initial
conditions, seeded, trajectories to attractor; 9 plaquettes, (S+1)^16 configs
-- far beyond enumeration, so 4x4 edges/fates are SAMPLED lower bounds).

FINDINGS:

1. AT tau_a >= tau_p THE AUTOMATON REMAINS FATE-EXACT AT 3x3, exhaustively:
   completeness (every dead class reaches z) and soundness (NO live class has
   any path to z) hold with zero leaks at (1,1), (2,1), (2,2), (3,2), (3,3) --
   the last over all 40M configurations. Persistence at tau_a >= tau_p is
   finite-state z-reachability over <= 19,251 symbols, one lattice size up.

2. SOUNDNESS FAILS AT tau_p > tau_a ON 3x3 -- THE 2x2 PICTURE WAS INCOMPLETE.
   At 2x2, soundness held in both regimes (12/12 cells). At 3x3 it BREAKS
   exactly in the failing regime: (1,2) has 14 pure-live classes with automaton
   paths to z; (2,3) has 288. So the fate-exactness of the symbolic abstraction
   is not a generic smallness artifact: it is a property of the tau_a >= tau_p
   regime, and the regime boundary now shows up in the automaton at TWO levels
   -- mixed classes appear AND soundness fails, both exactly at tau_p > tau_a.
   (2x2's universal soundness was the small-system accident, not the rule.)

3. 4x4 RANDOM-INIT AGREES (sampled, seeded): zero sampled-path leaks at every
   tau_a >= tau_p cell tested ((1,1) through (3,3), 2000-4000 runs each), and
   leaks appear at tau_p > tau_a ((1,2): 361, (2,3): 7). Attractor symbols of
   surviving waves never reach z in any cell. Caveats: sampled edges are a
   subgraph (missed edges can only UNDERcount leaks, so the tau_a >= tau_p
   zeros are evidence, not proof); "pure-live" at tau_p > tau_a may be
   under-sampled mixed classes; trajectories that fail to close within t_max
   are conservatively scored live.

Proof target after this: "class reaches z => all members die", now to be
proven at tau_a >= tau_p only -- 3x3 shows it is genuinely false in the other
regime, so no regime-blind argument can work. That is a real constraint on the
proof shape (it must USE the active-band-covers-passive-band structure).

House Rules Compliance: 3x3 exhaustive (no RNG); 4x4 sampling seeded
default_rng(1000 + 10*ta + tp) per cell; per-cell counts reported; the script
asserts findings 1-3 and aborts on regression.
Output: result/topology/spectrum_automaton_3x3.npz + printed doc table.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from collections import defaultdict

CHUNK = 500_000
CELLS_3 = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 2), (2, 3), (3, 3)]
CELLS_4 = [(1, 1, 4000), (2, 1, 4000), (2, 2, 3000), (3, 2, 2500),
           (3, 3, 2000), (1, 2, 4000), (2, 3, 2500)]
T_MAX = 3000

NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j%3 - i%3) == 1)
       for i in range(9)]
PLAQ3 = [(0, 1, 4, 3), (1, 2, 5, 4), (3, 4, 7, 6), (4, 5, 8, 7)]


# ---------------- 3x3 exhaustive (vectorised; see persistent_set_3x3.py) ----
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
            gaps = (np.roll(ph, -1, axis=1) - ph) % B
            gaps.sort(axis=1)
            pcs[:, p] = gaps @ base4
        pcs.sort(axis=1)
        code[lo:hi] = ((pcs[:, 0]*M + pcs[:, 1])*M + pcs[:, 2])*M + pcs[:, 3]
    return code


def audit_3x3(ta, tp):
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
    z = int(code[0])
    if fate_acc[z] != [False, True]:
        raise RuntimeError(f"z-class not uniformly dead at 3x3 ({ta},{tp})")
    radj = defaultdict(set)
    for a, b in edge_acc:
        radj[b].add(a)
    reach, stack = {z}, [z]
    while stack:
        for a in radj[stack.pop()]:
            if a not in reach:
                reach.add(a)
                stack.append(a)
    live = [s for s, (l, d) in fate_acc.items() if l and not d]
    dead = [s for s, (l, d) in fate_acc.items() if d and not l]
    mixed = [s for s, (l, d) in fate_acc.items() if l and d]
    return dict(ta=ta, tp=tp, N=N, n_cls=len(fate_acc), n_live=len(live),
                n_dead=len(dead), n_mixed=len(mixed), n_edges=len(edge_acc),
                complete=all(s in reach for s in dead),
                leaks=sum(1 for s in live if s in reach),
                mixed_rz=sum(1 for s in mixed if s in reach))


# ---------------- 4x4 sampled --------------------------------------------
L4 = 4
NB4 = [tuple(j for j in range(16) if abs(j//L4 - i//L4) + abs(j % L4 - i % L4) == 1)
       for i in range(16)]
PLAQ4 = [tuple([r*L4 + c, r*L4 + c + 1, (r+1)*L4 + c + 1, (r+1)*L4 + c])
         for r in range(L4 - 1) for c in range(L4 - 1)]


def step4(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if any(1 <= cfg[j] <= ta for j in NB4[i]) else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def sig4(cfg, B):
    out = []
    for cyc in PLAQ4:
        ph = [cfg[i] for i in cyc]
        out.append(tuple(sorted((ph[(k+1) % 4] - ph[k]) % B for k in range(4))))
    return tuple(sorted(out))


def sample_4x4(ta, tp, n_samples):
    B = ta + tp + 1
    rng = np.random.default_rng(1000 + 10 * ta + tp)
    edges, fate_of = defaultdict(set), defaultdict(set)
    att_syms_live, n_live, n_unclosed = set(), 0, 0
    for _ in range(n_samples):
        cfg = tuple(rng.integers(0, B, 16).tolist())
        seen, traj = {}, []
        while cfg not in seen:
            seen[cfg] = len(seen)
            traj.append(cfg)
            cfg = step4(cfg, ta, tp)
            if len(traj) > T_MAX:
                break
        if cfg in seen:
            att = traj[seen[cfg]:]
            alive = any(any(x) for x in att)
        else:
            att, alive = traj[-5:], True   # conservative: score unclosed as live
            n_unclosed += 1
        syms = [sig4(c, B) for c in traj]
        for a, b in zip(syms, syms[1:]):
            edges[a].add(b)
        for s in syms:
            fate_of[s].add(alive)
        if alive:
            n_live += 1
            att_syms_live.update(sig4(c, B) for c in att)
    zsym = sig4(tuple([0]*16), B)
    radj = defaultdict(set)
    for a, outs in edges.items():
        for b in outs:
            radj[b].add(a)
    reach, stack = {zsym}, [zsym]
    while stack:
        for a in radj[stack.pop()]:
            if a not in reach:
                reach.add(a)
                stack.append(a)
    pure_live = [s for s, f in fate_of.items() if f == {True}]
    return dict(ta=ta, tp=tp, n_samples=n_samples, n_live=n_live,
                n_unclosed=n_unclosed, n_sigs=len(fate_of),
                n_pure_live=len(pure_live),
                n_mixed=sum(1 for f in fate_of.values() if len(f) == 2),
                leaks=sum(1 for s in pure_live if s in reach),
                att_leaks=sum(1 for s in att_syms_live
                              if s in reach and fate_of[s] == {True}))


def main():
    print("=== 3x3 exhaustive ===")
    recs3 = []
    for ta, tp in CELLS_3:
        t0 = time.time()
        r = audit_3x3(ta, tp)
        recs3.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"({ta},{tp}) [{regime:6s}] N={r['N']:>10,} cls={r['n_cls']:6d} "
              f"live={r['n_live']:6d} mixed={r['n_mixed']:5d} "
              f"edges={r['n_edges']:7d} complete={r['complete']} "
              f"leaks={r['leaks']:3d} mix->z={r['mixed_rz']}/{r['n_mixed']} "
              f"[{time.time()-t0:.0f}s]", flush=True)

    assert all(r["complete"] for r in recs3), "3x3 completeness violated"
    for r in recs3:
        if r["ta"] >= r["tp"]:
            assert r["leaks"] == 0 and r["n_mixed"] == 0, \
                f"soundness/purity violated at tau_a>=tau_p ({r['ta']},{r['tp']})"
        else:
            assert r["n_mixed"] > 0, f"expected mixed classes at ({r['ta']},{r['tp']})"
        assert r["mixed_rz"] == r["n_mixed"], "a mixed class fails to reach z"
    assert any(r["leaks"] > 0 for r in recs3 if r["ta"] < r["tp"]), \
        "expected soundness failures at tau_p > tau_a on 3x3"

    print("\n=== 4x4 random-init (sampled; leaks are lower bounds) ===")
    recs4 = []
    for ta, tp, n in CELLS_4:
        t0 = time.time()
        r = sample_4x4(ta, tp, n)
        recs4.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"({ta},{tp}) [{regime:6s}] n={r['n_samples']} live={r['n_live']} "
              f"unclosed={r['n_unclosed']} sigs={r['n_sigs']:6d} "
              f"mixed={r['n_mixed']:4d} leaks={r['leaks']:4d} "
              f"att-leaks={r['att_leaks']} [{time.time()-t0:.0f}s]", flush=True)

    assert all(r["leaks"] == 0 for r in recs4 if r["ta"] >= r["tp"]), \
        "sampled leak at tau_a >= tau_p on 4x4"
    assert all(r["att_leaks"] == 0 for r in recs4), "attractor-symbol leak"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "spectrum_automaton_3x3.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"g3_{k}": np.array([r[k] for r in recs3]) for k in recs3[0]},
             **{f"g4_{k}": np.array([r[k] for r in recs4]) for k in recs4[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | regime | classes | mixed | leaks (live→z) | verdict |")
    print("| :---: | :---: | :---: | ---: | ---: | ---: | :--- |")
    for r in recs3:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        v = ("fate-exact" if r["leaks"] == 0 and r["n_mixed"] == 0
             else "mixed + unsound")
        print(f"| 3×3 (exhaustive) | ({r['ta']}, {r['tp']}) | {regime} | "
              f"{r['n_cls']:,} | {r['n_mixed']:,} | {r['leaks']} | {v} |")
    for r in recs4:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        v = "no sampled leaks" if r["leaks"] == 0 else "sampled leaks"
        print(f"| 4×4 (sampled, n={r['n_samples']}) | ({r['ta']}, {r['tp']}) | "
              f"{regime} | {r['n_sigs']:,} | {r['n_mixed']:,} | {r['leaks']} | {v} |")


if __name__ == "__main__":
    main()
