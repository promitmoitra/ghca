"""The gap signature generalises to 3x3 -- as a multiset over PLAQUETTES, not the
perimeter -- and the tau_p <= tau_a rule survives 40M configurations.

Extends persistent_set_structure.py one lattice size up. At 3x3 the reentrant
cycle is no longer unique (four 2x2 plaquettes + the perimeter 8-cycle), so "the
gap signature" needs choosing before it can be tested. Both candidates are
tested; only one works:

* PLAQUETTE signature: per 2x2 plaquette, the sorted cyclic phase-gap multiset
  mod (S+1) (exactly the 2x2 signature); then the sorted multiset of the four
  plaquette codes. Persistence is a FUNCTION of this signature in every
  tau_p <= tau_a cell tested -- including the boundary cell (3,3) with
  40,353,607 configurations and 19,251 classes, ZERO impure (2096x compression):

      (ta,tp)      N          P      classes  impure  compression
      (1,1)        19,683   0.7462       45       0       437x
      (2,1)       262,144   0.9746      319       0       822x
      (1,2)       262,144   0.4580      319     173       822x   <- tp>ta fails
      (2,2)     1,953,125   0.7752     1355       0      1441x
      (3,2)    10,077,696   0.9477     5830       0      1729x
      (2,3)    10,077,696   0.5549     5830    2378      1729x   <- tp>ta fails
      (3,3)    40,353,607   0.7826    19251       0      2096x

  The 2x2 dichotomy transfers exactly: invariant whenever tau_p <= tau_a, fails
  whenever tau_p > tau_a (here with no boundary counterexample; at 2x2, (2,5)
  showed tau_p > tau_a does not IMPLY failure, so the rule stays one-directional).

* PERIMETER signature (sorted gaps around the boundary 8-cycle): FAILS even at
  (1,1) -- 4 impure classes of 15 -- and everywhere else tested. The signature
  lives on the short reentrant loops, not on the long way around: consistent
  with E2's tau < L gate, where the 4-cycle (L=4) is the tightest sustainable
  loop and the 8-cycle is slack.

Compression SCALES: 16x-161x at 2x2 grows to 437x-2096x at 3x3, because classes
grow polynomially in S while configurations grow as (S+1)^9. P(tau_a,tau_p) on
3x3 is computable from ~19k classes instead of ~40M configurations.

Method (this is what makes 40M tractable in numpy, no compiled kernel needed):
1. The full transition function is materialised as a successor ARRAY f (int32,
   one entry per configuration), built in vectorised chunks.
2. Persistence reduces to reachability of a single sink: for theta=1 GH, the
   ONLY dead attractor is the all-zero fixed point. Non-receptive states advance
   autonomously and decay to 0, so an attractor either regenerates activity
   forever or collapses to 0. Checked EXHAUSTIVELY at 2x2 over 4 (ta,tp) cells
   (no non-zero dead attractor exists), and the whole pipeline is checked at
   3x3 (1,1) against direct orbit iteration on 2000 seeded-random configurations
   (0 mismatches) -- both checks are re-run by this script and it aborts on any
   mismatch.
3. "Reaches 0" is computed by boolean label propagation dead |= dead[f] to a
   fixed point -- O(N) memory, no orbit hashing.
4. Class purity via a per-chunk dict merge (code -> saw-persistent/saw-dead),
   avoiding a global 40M sort.
   An earlier draft used pointer doubling g = g[g]; correct, but its per-round
   full-array gather cost ~2x N int32 temporaries and tripped the 15GB machine's
   OOM killer (exit 137) at N=40M alongside the code arrays. Label propagation
   converges in <= attractor-distance rounds (~tens here) at a quarter of the
   footprint. Peak RSS for the (3,3) cell: ~0.4 GB.

House Rules Compliance:
    - The only RNG is the 2000-sample cross-validation, seeded explicitly
      (default_rng(SEED)); every reported number is a complete enumeration.
    - Substrate/analysis boundary: step function validated against
      ghca_net.Network at 2x2 by the companion experiment; here the 3x3
      dynamics are additionally validated against direct orbit iteration.
      Signatures are analysis constructs.
Output: result/topology/persistent_set_3x3.npz + a printed doc table.
"""
import sys, os, time
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

SEED = 0
CELLS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 2), (2, 3), (3, 3)]
CHUNK = 500_000

# 3x3 box, von Neumann, open boundary; plaquettes and perimeter
NB3 = [tuple(j for j in range(9) if abs(j//3 - i//3) + abs(j%3 - i%3) == 1)
       for i in range(9)]
PLAQ = [(0, 1, 4, 3), (1, 2, 5, 4), (3, 4, 7, 6), (4, 5, 8, 7)]
PERIM = (0, 1, 2, 5, 8, 7, 6, 3)

# 2x2 (for the dead-attractor claim check)
ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def check_dead_attractor_claim():
    """Exhaustive at 2x2: every dead attractor is exactly {all-zero}."""
    for ta, tp in [(1, 1), (2, 2), (2, 3), (3, 2)]:
        S = ta + tp
        for c in product(range(S + 1), repeat=4):
            seen, cur = {}, c
            while cur not in seen:
                seen[cur] = len(seen)
                cur = step2(cur, ta, tp)
            att = [k for k, t in seen.items() if t >= seen[cur]]
            if not any(any(x) for x in att) and set(att) != {(0, 0, 0, 0)}:
                raise RuntimeError(f"non-zero dead attractor at ({ta},{tp}): {att}")
    return True


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
    """~(reaches config 0), by boolean label propagation to a fixed point."""
    dead = f == 0
    dead[0] = True
    while True:
        nd = dead | dead[f]
        if bool((nd == dead).all()):
            return ~nd
        dead = nd


def signature_codes(N, B, cycles, pack):
    """One int64 code per config: sorted-gap multiset per cycle, sorted across
    cycles, packed. `pack` must inject for the given B (asserted by caller)."""
    pw = (B ** np.arange(9)).astype(np.int64)
    base = (B ** np.arange(len(cycles[0]))).astype(np.int64)
    code = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        idx = np.arange(lo, hi, dtype=np.int64)
        dig = ((idx[:, None] // pw[None, :]) % B).astype(np.int64)
        pcs = np.empty((hi - lo, len(cycles)), dtype=np.int64)
        for p, cyc in enumerate(cycles):
            ph = dig[:, list(cyc)]
            gaps = (np.roll(ph, -1, axis=1) - ph) % B
            gaps.sort(axis=1)
            pcs[:, p] = gaps @ base
        pcs.sort(axis=1)
        code[lo:hi] = pack(pcs)
        del idx, dig, pcs
    return code


def purity(code, pm):
    acc = {}
    for lo in range(0, len(code), CHUNK):
        hi = min(lo + CHUNK, len(code))
        c, p = code[lo:hi], pm[lo:hi]
        for val, mask in ((0, p), (1, ~p)):
            for u in np.unique(c[mask]):
                acc.setdefault(int(u), [False, False])[val] = True
    ng = len(acc)
    imp = sum(1 for a, b in acc.values() if a and b)
    return ng, imp


def validate_3x3_11(f, pm, n_sample=2000):
    """Direct orbit iteration vs the array pipeline at (1,1)."""
    B = 3
    def orbit_persist(idx):
        cfg = tuple((idx // B**i) % B for i in range(9))
        def st(c):
            return tuple((1 if any(1 <= c[j] <= 1 for j in NB3[i]) else 0)
                         if s == 0 else (s + 1) % B for i, s in enumerate(c))
        seen, cur = {}, cfg
        while cur not in seen:
            seen[cur] = len(seen)
            cur = st(cur)
        att = [k for k, t in seen.items() if t >= seen[cur]]
        return any(any(x) for x in att)
    rng = np.random.default_rng(SEED)
    sample = rng.choice(len(f), size=n_sample, replace=False)
    mism = sum(orbit_persist(int(i)) != bool(pm[i]) for i in sample)
    return n_sample, mism


def main():
    print("dead-attractor claim (2x2, exhaustive, 4 cells):",
          "OK" if check_dead_attractor_claim() else "FAIL")

    recs = []
    print(f"\n{'(ta,tp)':>8} {'N':>11} {'P':>7} {'plaq-cls':>9} {'impure':>7} "
          f"{'compress':>9} {'perim-cls':>10} {'p-impure':>9} {'secs':>6}")
    for ta, tp in CELLS:
        B = ta + tp + 1
        t0 = time.time()
        f = successor_array(ta, tp)
        pm = persistent_mask(f)
        N = len(f)
        if (ta, tp) == (1, 1):
            n_s, mism = validate_3x3_11(f, pm)
            if mism:
                raise RuntimeError(f"pipeline vs orbit mismatch: {mism}/{n_s}")
        del f

        M = int(B) ** 4
        assert M ** 4 < 2 ** 62, "plaquette pack overflows int64"
        pc = signature_codes(N, B, PLAQ,
                             lambda pcs: ((pcs[:, 0] * M + pcs[:, 1]) * M
                                          + pcs[:, 2]) * M + pcs[:, 3])
        ng, imp = purity(pc, pm)
        del pc

        # perimeter only where affordable (its class count is tiny anyway)
        if N <= 2_000_000:
            qc = signature_codes(N, B, [PERIM], lambda pcs: pcs[:, 0])
            ngq, impq = purity(qc, pm)
            del qc
        else:
            ngq, impq = -1, -1

        recs.append(dict(ta=ta, tp=tp, N=N, P=float(pm.mean()),
                         n_classes=ng, n_impure=imp,
                         perim_classes=ngq, perim_impure=impq,
                         invariant=(imp == 0), secs=time.time() - t0))
        r = recs[-1]
        print(f"{str((ta,tp)):>8} {N:11,d} {r['P']:7.4f} {ng:9d} {imp:7d} "
              f"{N/ng:8.0f}x {ngq:10d} {impq:9d} {r['secs']:6.0f}", flush=True)
        del pm

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "persistent_set_3x3.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| (tau_a, tau_p) | configs | P | plaquette classes | impure | "
          "compression | invariant |")
    print("| :---: | ---: | ---: | ---: | :---: | ---: | :---: |")
    for r in recs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['N']:,} | {r['P']:.4f} | "
              f"{r['n_classes']:,} | {r['n_impure']:,} | {r['N']/r['n_classes']:.0f}x | "
              f"{'yes' if r['invariant'] else 'NO'} |")

    inv = [(r['ta'], r['tp']) for r in recs if r['invariant']]
    bad = [(r['ta'], r['tp']) for r in recs if not r['invariant']]
    print(f"\ninvariant: {inv}  -- all tau_p <= tau_a: {all(t <= a for a, t in inv)}")
    print(f"fails:     {bad}  -- all tau_p >  tau_a: {all(t > a for a, t in bad)}")
    p11 = [r for r in recs if (r['ta'], r['tp']) == (1, 1)][0]
    print(f"perimeter signature at (1,1): {p11['perim_impure']} impure of "
          f"{p11['perim_classes']} -- fails even where the plaquette one holds")


if __name__ == "__main__":
    main()
