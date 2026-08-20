"""The covering lemma's existence step, closed empirically: live quiet runs are
bounded by min(ceil(S/2)+1, tau_p+2) -- a SINGLE-trajectory fact, no pairs --
and the young cell is always the witness.

Direction 2 of docs/research_brainstorming_coherence_horizons.md, rescoped.
That doc proposes a two-week pilot to test H1 ("every ceiling hold has a
u-side cell with quiet run Q < S"). H1 runs exhaustively in seconds, and the
reason it holds turns out to be stronger and simpler than H1 itself.

DEFINITIONS. Q(i,t) = number of consecutive steps cell i has been receptive
(state 0) up to time t. "Ceiling hold" = a pair (v,u) at age S whose image is
also at age S (age = BFS depth from the clock-shift diagonal, per
coherence_window_S.py). "Witness" = a subset of u's dwelling cells that,
re-read 0 -> S (the swap law), yields a valid pair of age S-1.

FINDINGS (exhaustive over all live configurations; 2x2 core):

1. H1 HOLDS WITH ZERO VIOLATIONS at every cell tested -- 32/32 ceiling holds
   at (1,1) through 152/152 at (4,4). Confirmed, and the Direction-2 pilot is
   therefore a one-day confirmation rather than a two-week experiment.

2. THE YOUNG CELL IS THE WITNESS. At every ceiling hold, a witness exists
   using ONLY cells with Q < S (32/32 ... 152/152, matching H1 exactly). So
   the covering lemma's existence step does not merely coincide with young
   cells; young cells supply it.

3. THE REAL REASON, AND IT NEEDS NO PAIRS: on a live orbit the quiet run is
   bounded outright,
       max Q = min(ceil(S/2)+1, tau_p+2),
   exact at all 16 tau_a >= tau_p cells tested from (1,1) to (5,5). Since that
   maximum is < S at every cell except the two smallest ((1,1): 2 = S, (2,1):
   3 = S, where it merely ties), EVERY dwelling cell on a live orbit is young
   almost everywhere, and H1's existential is nearly vacuous once the bound is
   known. The bound is a statement about ONE trajectory -- no pair, no ledger,
   no BFS -- which is exactly the reduction a hand proof wants.

   Reading the two terms: ceil(S/2)+1 binds when the passive band is long (a
   cell cannot out-wait half a cycle plus one, because the wave returns);
   tau_p+2 binds when the passive band is short (the pipeline flushes and
   re-excites within tau_p+2). The crossover is at tau_p ~ ceil(S/2)-1, and
   the min form is exact on both sides.

4. NOT REGIME-CARRYING -- stated so nobody builds on it as if it were. In the
   wrong regime (tau_p > tau_a) the quiet run is ALSO short (maxQ = 2, 3, 0, 3
   at (1,2), (2,3), (1,3), (2,4)), with min(...) a valid but loose upper
   bound. Bounded quiet runs are thus necessary for the covering lemma but
   cannot by themselves explain why coherence is signature-measurable only at
   tau_a >= tau_p. The regime content lives elsewhere (the anchor law).

SUPERSEDED IN PART -- see covering_lemma_quiet_runs_scaling.py: the closed
form below is FALSE beyond 2x2 (maxQ = 7 at 3x3 (2,1), vs 3 predicted),
and H1 itself is falsified at 3x3 (8 violations of 58,588 ceiling holds).
The covering lemma survives there, but its witnesses use OLD cells
(8/8), so "young cells supply the witness" is a 2x2 fact only.
The hand obligation stated below is WITHDRAWN.

WHAT THIS LEAVES FOR THE HAND PROOF. Prove max Q <= min(ceil(S/2)+1, tau_p+2)
on a live orbit of any graph -- a local statement about the refractory
pipeline and wave return time, with no reference to pairs, ages, or the
diagonal. Then the covering lemma's existence step follows wherever that bound
is < S, and the remaining work is the swap validity check, not the search.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts H1 counts,
young-cell witnessing, the closed form on all 16 covered cells, the two
off-form cells' actual values, and the wrong-regime control; aborts on
regression. Caveats are adjacent to every claim above; all results are 2x2
core -- the bound is untested on larger lattices, which is the next step.
Output: result/topology/covering_lemma_quiet_runs.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product, combinations

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
H1_CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]
H1_EXPECT = {(1, 1): 32, (2, 1): 12, (2, 2): 24, (3, 1): 16,
             (3, 2): 44, (3, 3): 64, (4, 3): 120, (4, 4): 152}
Q_CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2),
           (4, 3), (4, 4), (5, 1), (5, 2), (5, 3), (5, 4), (5, 5), (6, 1)]
Q_EXPECT = {(1, 1): 2, (2, 1): 3, (2, 2): 3, (3, 1): 3, (3, 2): 4, (3, 3): 4,
            (4, 1): 3, (4, 2): 4, (4, 3): 5, (4, 4): 5, (5, 1): 3, (5, 2): 4,
            (5, 3): 5, (5, 4): 6, (5, 5): 6, (6, 1): 3}
WRONG_CELLS = [(1, 2), (2, 3), (1, 3), (2, 4)]
WRONG_EXPECT = {(1, 2): 2, (2, 3): 3, (1, 3): 0, (2, 4): 3}


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
        fl = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = fl
    return {c for c, fl in cache.items() if fl}


def pair_depths(stp, B, lv):
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


def max_quiet_run(ta, tp):
    """Max consecutive receptive steps on any LIVE orbit. One trajectory."""
    stp, B = make_step(ta, tp)
    mx = 0
    for w in sorted(live_set(stp, B)):
        c = w
        Q = [0] * 4
        seen = set()
        while c not in seen:
            seen.add(c)
            Q = [Q[i] + 1 if c[i] == 0 else 0 for i in range(4)]
            mx = max(mx, max(Q))
            c = stp(c)
    return mx


def audit_h1(ta, tp):
    """H1 and young-cell witnessing at every ceiling hold."""
    stp, B = make_step(ta, tp)
    S = ta + tp
    lv = live_set(stp, B)
    dep = pair_depths(stp, B, lv)
    n_ceil = ok = viol = young_wit = 0
    for w in sorted(lv):
        v, u = w, tuple((x + 1) % B for x in w)
        Q = [0] * 4
        seen = set()
        while (v, u) not in seen:
            seen.add((v, u))
            nv, nu = stp(v), stp(u)
            if dep.get((v, u)) == S and dep.get((nv, nu)) == S:
                n_ceil += 1
                dws = [i for i in range(4) if u[i] == 0 and nu[i] == 0]
                young = [i for i in dws if Q[i] < S]
                if young:
                    ok += 1
                else:
                    viol += 1
                hit = False
                for k in range(1, len(young) + 1):
                    for su in combinations(young, k):
                        u2 = tuple(S if i in su else u[i] for i in range(4))
                        if stp(u2) == nu and dep.get((v, u2)) == S - 1:
                            hit = True
                            break
                    if hit:
                        break
                young_wit += hit
            Q = [Q[i] + 1 if u[i] == 0 else 0 for i in range(4)]
            v, u = nv, nu
    return dict(ta=ta, tp=tp, n_ceil=n_ceil, h1_ok=ok, h1_viol=viol,
                young_witness=young_wit)


def main():
    print("=== H1 + young-cell witnessing (exhaustive over live configs) ===")
    h1 = []
    for ta, tp in H1_CELLS:
        r = audit_h1(ta, tp)
        h1.append(r)
        print(f"({ta},{tp}) S={ta+tp}: ceiling holds={r['n_ceil']:4d}  "
              f"H1 ok={r['h1_ok']:4d}  violations={r['h1_viol']}  "
              f"young-cell witness={r['young_witness']:4d}", flush=True)
        assert r["n_ceil"] == H1_EXPECT[(ta, tp)], f"ceiling count ({ta},{tp})"
        assert r["h1_viol"] == 0, f"H1 FALSIFIED at ({ta},{tp})"
        assert r["young_witness"] == r["n_ceil"], \
            f"young cells failed to witness at ({ta},{tp})"

    print("\n=== the underlying bound: max quiet run on a LIVE orbit ===")
    print(f"{'cell':>7} {'S':>3} {'maxQ':>5} {'min(ceil(S/2)+1, tp+2)':>23} {'<S?':>5}")
    qs = []
    for ta, tp in Q_CELLS:
        S = ta + tp
        mx = max_quiet_run(ta, tp)
        pred = min(-(-S // 2) + 1, tp + 2)
        qs.append(dict(ta=ta, tp=tp, maxQ=mx, pred=pred))
        print(f"{str((ta,tp)):>7} {S:>3} {mx:>5} {pred:>23} {str(mx < S):>5}",
              flush=True)
        assert mx == Q_EXPECT[(ta, tp)], f"quiet run changed at ({ta},{tp})"
        assert mx == pred, f"closed form broke at ({ta},{tp}): {mx} vs {pred}"

    print("\n=== control: wrong regime (tau_p > tau_a) -- bound is loose, "
          "not regime-carrying ===")
    wr = []
    for ta, tp in WRONG_CELLS:
        S = ta + tp
        mx = max_quiet_run(ta, tp)
        pred = min(-(-S // 2) + 1, tp + 2)
        wr.append(dict(ta=ta, tp=tp, maxQ=mx, pred=pred))
        print(f"({ta},{tp}) S={S}: maxQ={mx}  min-form={pred}  "
              f"exact={mx == pred}  still < S: {mx < S}", flush=True)
        assert mx == WRONG_EXPECT[(ta, tp)], f"wrong-regime Q ({ta},{tp})"
        assert mx <= pred, "min-form must remain an upper bound"
        assert mx < S, "quiet runs are short in BOTH regimes"
    assert any(r["maxQ"] != r["pred"] for r in wr), \
        "the closed form must be INEXACT off-regime (else it is regime-blind)"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "covering_lemma_quiet_runs.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"h1_{k}": np.array([r[k] for r in h1]) for k in h1[0]},
             **{f"q_{k}": np.array([r[k] for r in qs]) for k in qs[0]},
             **{f"w_{k}": np.array([r[k] for r in wr]) for k in wr[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | S | max live quiet run | min(⌈S/2⌉+1, τp+2) | < S |")
    print("| :---: | ---: | ---: | ---: | :---: |")
    for r in qs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['ta']+r['tp']} | {r['maxQ']} | "
              f"{r['pred']} | {'yes' if r['maxQ'] < r['ta']+r['tp'] else 'ties'} |")


if __name__ == "__main__":
    main()
