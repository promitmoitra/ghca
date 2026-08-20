"""Witness selection: the covering witness is UNIQUE at every ceiling hold, the
selection question is VACUOUS at 2x2 (every dwelling subset witnesses), and no
cell-local rule decides membership at 3x3 -- the compression barrier again, one
level down.

Follows covering_lemma_quiet_runs_scaling.py, which falsified quiet-run age as
the covering mechanism (old cells witness) and left the question: WHICH dwelling
subsets admit a consistent older reading? Structural, not temporal.

FINDINGS:

1. THE 2x2 CORE IS DEGENERATE FOR THIS QUESTION, not merely small. At every
   2x2 cell tested ((2,1), (2,2), (3,2), (3,3), (4,4)), EVERY dwelling subset
   of every ceiling hold witnesses: 12/24/44/64/152 witnessing, ZERO
   non-witnessing. The discriminating question does not exist there. This is
   the structural reason three separate hypotheses (the quiet-run closed form,
   H1, "young cells supply the witness") all read as universal laws at 2x2 and
   all died at 3x3 -- the core cannot distinguish them. Methodological note for
   the handoff: for witness-SELECTION questions, 2x2 is not a weak test, it is
   a vacuous one.

2. AT 3x3 THE QUESTION IS REAL AND THE ANSWER IS UNIQUE. Over all 25,998
   ceiling holds at (2,1), enumerating all dwelling subsets gives 25,998
   witnessing and 90,014 non-witnessing -- and exactly ONE witnessing subset
   per hold (census {1: 25998}, no hold has two). The covering witness is not
   merely existent, it is unique. P(witness) falls with subset size
   (0.282 at k=1 down to 0.135 at k=5), and the witness is the FULL dwelling
   set in 22,622 of 25,998 holds (87.0%); P(wit | full set) = 0.870 vs
   P(wit | proper subset) = 0.038.

3. NO CELL-LOCAL RULE DECIDES MEMBERSHIP. Six predicates were tested per
   dwelling cell against the unique witness, all fail:
     v_i active                       3,576 mismatches
     v_i != 0                         6,292
     v_i does not dwell               3,380
     u2_i == v_i + 1 after the swap  18,690
     v_i == S (wraps)                25,998
     ledger + receptive tie-break     1,972   <- best, 92.4% correct
   The ledger d_i = (u_i - v_i) mod B is strongly predictive but not
   sufficient: d=1 (the unshifted baseline) is NEVER in the witness
   (0 / 3,768), d in {2,3} is ALWAYS in it (39,379 / 39,379), and d=0 splits
   4,260 in / 2,172 out. "v_i receptive" resolves the d=0 split at the
   cell level (4,260 agreements) but the combined rule still misses 1,972
   holds, so the residual is not decidable from cell i's own data.

   This is the COMPRESSION BARRIER one level down. anchor_law_certificate.py
   showed seven ledger-local invariants leak on the pair set; here the ledger
   plus one local state bit gets 92.4% of witness selection and no further. The
   same conclusion follows: the object that decides is the wave configuration,
   not any per-cell bookkeeping.

WHAT THIS LEAVES. Uniqueness is a new exact handle: a proof of the covering
lemma no longer needs to search subsets, only to exhibit THE witness. The
ledger gives two thirds of it for free (d=1 out, d in {2,3} in), so the
remaining obligation is sharply narrowed to the d=0 cells -- those where u and v
sit at the same phase -- and to why 13.0% of holds exclude some dwelling cell.
That is a statement about coincident-phase cells, not about waiting times.

House Rules Compliance: exhaustive/deterministic at both lattices, no RNG;
asserts the 2x2 degeneracy (zero non-witnessing subsets), the 3x3 counts, the
uniqueness census, the ledger split, and the mismatch count of every failed
predicate so no future edit can promote one to a law; aborts on regression.
Caveats: 3x3 results are the (2,1) cell only -- other 3x3 cells and 4x4 are
untested for selection, and the 87% full-set rate is not claimed to generalise.
Output: result/topology/covering_witness_selection.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product, combinations
from collections import Counter

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
NB3 = [tuple(j for j in range(9) if abs(j // 3 - i // 3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
CHUNK = 1_000_000
CELLS_2 = [(2, 1), (2, 2), (3, 2), (3, 3), (4, 4)]
WIT_2_EXPECT = {(2, 1): 12, (2, 2): 24, (3, 2): 44, (3, 3): 64, (4, 4): 152}
N_CEIL_3 = 25_998
N_NON_3 = 90_014
FULL_3 = 22_622
LEDGER_EXPECT = {0: (4260, 2172), 1: (0, 3768), 2: (20534, 0), 3: (18845, 0)}
RULE_MISMATCH = {"v_i active": 3576, "v_i != 0": 6292, "v_i does not dwell": 3380,
                 "u2_i == v_i+1": 18690, "v_i == S": 25998,
                 "ledger + receptive": 1972}


def make_step2(ta, tp):
    B = ta + tp + 1

    def stp(cfg):
        return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                     if s == 0 else (s + 1) % B for i, s in enumerate(cfg))
    return stp, B


def live_2x2(stp, B):
    cache = {}
    for c in product(range(B), repeat=4):
        if c in cache:
            continue
        seen, cur = {}, c
        while cur not in seen:
            seen[cur] = len(seen)
            cur = stp(cur)
        fl = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = fl
    return {c for c, fl in cache.items() if fl}


def depths_2x2(stp, B, lv):
    dep = {}
    fr = [(w, tuple((x + 1) % B for x in w)) for w in sorted(lv)]
    for p in fr:
        dep[p] = 0
    d = 0
    while fr:
        d += 1
        nx = []
        for v, u in fr:
            q = (stp(v), stp(u))
            if q not in dep:
                dep[q] = d
                nx.append(q)
        fr = nx
    return dep


def degeneracy_2x2(ta, tp):
    """Label every dwelling subset at every ceiling hold. 2x2."""
    stp, B = make_step2(ta, tp)
    S = ta + tp
    lv = live_2x2(stp, B)
    dep = depths_2x2(stp, B, lv)
    wit = non = 0
    for w in sorted(lv):
        v, u = w, tuple((x + 1) % B for x in w)
        seen = set()
        while (v, u) not in seen:
            seen.add((v, u))
            nv, nu = stp(v), stp(u)
            if dep.get((v, u)) == S and dep.get((nv, nu)) == S:
                dws = [i for i in range(4) if u[i] == 0 and nu[i] == 0]
                for k in range(1, len(dws) + 1):
                    for su in combinations(dws, k):
                        u2 = tuple(S if i in su else u[i] for i in range(4))
                        if stp(u2) == nu and dep.get((v, u2)) == S - 1:
                            wit += 1
                        else:
                            non += 1
            v, u = nv, nu
    return wit, non


def build_3x3(ta=2, tp=1):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        dig = ((np.arange(lo, hi, dtype=np.int64)[:, None] // pw[None, :]) % B).astype(np.int8)
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
    live = np.nonzero(~dead)[0]
    idxall = np.arange(N, dtype=np.int64)
    sh = np.zeros(N, dtype=np.int64)
    for k in range(9):
        sh += (((idxall // pw[k]) % B + 1) % B) * pw[k]
    depth = {}
    V, U = live.copy(), sh[live]
    for c in (V * N + U).tolist():
        depth[c] = 0
    d = 0
    while True:
        nV, nU = f[V], f[U]
        code = nV * N + nU
        mask = np.fromiter((c not in depth for c in code.tolist()),
                           dtype=bool, count=len(code))
        if not mask.any():
            break
        d += 1
        for c in code[mask].tolist():
            depth[c] = d
        V, U = nV[mask], nU[mask]
    codes = np.fromiter(depth.keys(), dtype=np.int64, count=len(depth))
    ags = np.fromiter(depth.values(), dtype=np.int16, count=len(depth))
    o = np.argsort(codes)
    return f, live, sh, B, N, pw, codes[o], ags[o]


def main():
    print("=== 1) 2x2 is DEGENERATE: every dwelling subset witnesses ===")
    d2 = []
    for ta, tp in CELLS_2:
        wit, non = degeneracy_2x2(ta, tp)
        d2.append(dict(ta=ta, tp=tp, wit=wit, non=non))
        print(f"({ta},{tp}): witnessing={wit:4d}  NON-witnessing={non}", flush=True)
        assert wit == WIT_2_EXPECT[(ta, tp)], f"2x2 witness count ({ta},{tp})"
        assert non == 0, \
            f"2x2 must be degenerate at ({ta},{tp}) -- a non-witnessing subset appeared"

    print("\n=== 2) 3x3 (2,1): the question is real, the witness is UNIQUE ===")
    ta, tp = 2, 1
    S = ta + tp
    f, live, sh, B, N, pw, codes_s, ages_s = build_3x3(ta, tp)

    def age_of(c):
        pos = np.clip(np.searchsorted(codes_s, c), 0, len(codes_s) - 1)
        hit = codes_s[pos] == c
        out = np.full(len(c), -1, dtype=np.int16)
        out[hit] = ages_s[pos[hit]]
        return out

    def digs(x):
        return [(int(x) // int(pw[k])) % B for k in range(9)]

    holds = []
    Vv, Uu = live.copy(), sh[live]
    for _ in range(4 * S + 20):
        nV, nU = f[Vv], f[Uu]
        ceil = (age_of(Vv * N + Uu) == S) & (age_of(nV * N + nU) == S)
        for j in np.nonzero(ceil)[0]:
            holds.append((int(Vv[j]), int(Uu[j])))
        Vv, Uu = nV, nU
    holds = list(dict.fromkeys(holds))

    uniq = Counter()
    n_wit = n_non = n_full = 0
    by_size = Counter()
    ledger = Counter()
    wit_of = {}
    for v_, u_ in holds:
        ud, nud = digs(u_), digs(int(f[u_]))
        dws = [i for i in range(9) if ud[i] == 0 and nud[i] == 0]
        if not dws:
            continue
        subs, cds = [], []
        for k in range(1, len(dws) + 1):
            for su in combinations(dws, k):
                subs.append(su)
                cds.append(v_ * N + (u_ + sum(S * int(pw[i]) for i in su)))
        a = age_of(np.array(cds, dtype=np.int64))
        wins = []
        for su, ai in zip(subs, a):
            if ai == S - 1:
                wins.append(su)
                n_wit += 1
                by_size[(len(su), True)] += 1
            else:
                n_non += 1
                by_size[(len(su), False)] += 1
        uniq[len(wins)] += 1
        if len(wins) == 1:
            wit_of[(v_, u_)] = (wins[0], tuple(dws))
            if set(wins[0]) == set(dws):
                n_full += 1
            vd = digs(v_)
            for i in dws:
                ledger[((ud[i] - vd[i]) % B, i in wins[0])] += 1

    print(f"ceiling holds={len(holds):,}  witnessing subsets={n_wit:,}  "
          f"NON-witnessing={n_non:,}", flush=True)
    print(f"witnesses per hold: {dict(sorted(uniq.items()))}  -> UNIQUE")
    print(f"witness == full dwelling set: {n_full:,}/{len(wit_of):,} "
          f"({100*n_full/len(wit_of):.1f}%)")
    assert len(holds) == N_CEIL_3, "3x3 ceiling-hold count changed"
    assert n_wit == N_CEIL_3 and n_non == N_NON_3, "3x3 subset labels changed"
    assert dict(uniq) == {1: N_CEIL_3}, "witness UNIQUENESS broke"
    assert n_full == FULL_3, "full-set rate changed"

    print("\nledger d_i = (u_i - v_i) mod B vs witness membership:")
    for dv in sorted({k[0] for k in ledger}):
        a, b = ledger[(dv, True)], ledger[(dv, False)]
        print(f"  d={dv}: in={a:6,}  out={b:6,}  P(in)={a/max(a+b,1):.3f}")
        assert (a, b) == LEDGER_EXPECT[dv], f"ledger split changed at d={dv}"

    print("\n=== 3) no cell-local rule decides membership ===")
    preds = {
        "v_i active": lambda vd, nvd, ud, u2d, i: 1 <= vd[i] <= ta,
        "v_i != 0": lambda vd, nvd, ud, u2d, i: vd[i] != 0,
        "v_i does not dwell": lambda vd, nvd, ud, u2d, i: not (vd[i] == 0 and nvd[i] == 0),
        "u2_i == v_i+1": lambda vd, nvd, ud, u2d, i: (u2d[i] - vd[i]) % B == 1,
        "v_i == S": lambda vd, nvd, ud, u2d, i: vd[i] == S,
        "ledger + receptive": lambda vd, nvd, ud, u2d, i:
            ((ud[i] - vd[i]) % B in (2, 3)) or ((ud[i] - vd[i]) % B == 0 and vd[i] == 0),
    }
    rule_res = []
    for name, fn in preds.items():
        bad = 0
        for (v_, u_), (w, dws) in wit_of.items():
            vd, nvd, ud = digs(v_), digs(int(f[v_])), digs(u_)
            u2d = digs(u_ + sum(S * int(pw[i]) for i in w))
            if tuple(i for i in dws if fn(vd, nvd, ud, u2d, i)) != tuple(sorted(w)):
                bad += 1
        rule_res.append(dict(name=name, mismatch=bad))
        print(f"  {name:>20}: mismatches {bad:6,} / {len(wit_of):,} "
              f"({100*(1-bad/len(wit_of)):.1f}% correct)", flush=True)
        assert bad == RULE_MISMATCH[name], f"predicate '{name}' changed"
        assert bad > 0, f"predicate '{name}' unexpectedly became exact -- re-examine"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "covering_witness_selection.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"d2_{k}": np.array([r[k] for r in d2]) for k in d2[0]},
             n_holds=np.array([len(holds)]), n_wit=np.array([n_wit]),
             n_non=np.array([n_non]), n_full=np.array([n_full]),
             uniq_keys=np.array(sorted(uniq)),
             uniq_vals=np.array([uniq[k] for k in sorted(uniq)]),
             size_keys=np.array([k for k, _ in sorted(by_size)]),
             size_lab=np.array([int(l) for _, l in sorted(by_size)]),
             size_vals=np.array([by_size[k] for k in sorted(by_size)]),
             ledger_d=np.array([d for d in sorted({k[0] for k in ledger})]),
             ledger_in=np.array([ledger[(d, True)] for d in sorted({k[0] for k in ledger})]),
             ledger_out=np.array([ledger[(d, False)] for d in sorted({k[0] for k in ledger})]),
             rule_names=np.array([r["name"] for r in rule_res]),
             rule_mismatch=np.array([r["mismatch"] for r in rule_res]))
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | witnessing subsets | non-witnessing | question exists? |")
    print("| :---: | :---: | ---: | ---: | :---: |")
    for r in d2:
        print(f"| 2×2 | ({r['ta']}, {r['tp']}) | {r['wit']} | {r['non']} | "
              f"{'yes' if r['non'] else '**no — degenerate**'} |")
    print(f"| 3×3 | (2, 1) | {n_wit:,} | {n_non:,} | yes |")


if __name__ == "__main__":
    main()
