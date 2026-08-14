"""The spectrum automaton: fate is z-reachability at the symbol level, and the
regime boundary is exactly the appearance of mixed classes.

Second experiment on the theory branch (after spectrum_mechanism_hunt.py, which
falsified bridging and lumpability). Construction, in the style of symbolic
dynamics: project each configuration to its gradient spectrum (the symbol), and
let the SPECTRUM AUTOMATON have one node per symbol with a nondeterministic
edge s -> s' iff some configuration with spectrum s steps to one with spectrum
s'. This is the finite symbolic factor of the dynamics induced by the (a
priori non-Markov) spectrum partition -- lumpability already failed, so the
automaton over-approximates real trajectories. The all-equal-phases symbol
z = (0,0,0,0) is the candidate dead sink.

FINDINGS, exhaustive at 12 (tau_a, tau_p) cells, both regimes:

1. THE ABSTRACTION IS FATE-EXACT ON PURE CLASSES -- EVERYWHERE. In every cell
   tested, in BOTH regimes:
     - every uniformly-dead class has an automaton path to z    (completeness)
     - NO uniformly-live class has an automaton path to z       (soundness)
   The second is the remarkable half: the automaton over-approximates (edges
   exist that no single trajectory realises in context), yet it never
   fabricates a path from a live class to the dead sink. Zero leaks in
   12/12 cells.

2. THEREFORE, AT tau_a >= tau_p, FATE IS DECIDED BY FINITE-STATE REACHABILITY.
   There, every class is pure (persistent_set_structure.md Section 3), so:

       a configuration persists  <=>  its spectrum has no automaton path to z.

   Persistence -- defined over infinite time on the full configuration space --
   is computed by graph reachability on <= 55 symbols. This is the
   sofic-style characterisation the symbolic-dynamics framing promises: the
   persistent set is (spectrum-)recognised by a finite automaton.

3. THE REGIME BOUNDARY IS EXACTLY THE APPEARANCE OF MIXED CLASSES, AND EVERY
   MIXED CLASS REACHES z. At tau_p > tau_a the mixed classes (2 at (1,2), 4 at
   (2,3), 8 at (3,4)) ALL have automaton paths to z -- consistent with their
   dead members -- while also containing live members. So the automaton never
   lies about pure classes; insufficiency at tau_p > tau_a is not an
   abstraction artifact but the genuine bifurcation of same-symbol initial
   conditions.

4. PROOF TARGET, REFORMULATED. To prove spectrum-sufficiency at tau_a >= tau_p
   it now suffices to prove: "if a class has an automaton path to z, all its
   members die." The contrapositive of soundness plus completeness then gives
   purity. This is an automaton-side statement about a finite object per
   (tau_a, tau_p) -- a large structural step down from quantifying over all
   trajectories.

Context: this is the 2D-core analogue of the symbolic-dynamics program of
Kessebohmer, Rademacher & Ulbrich (arXiv:1903.02459), who analyse 1D GH via
its non-wandering set and subshifts; their topological-entropy results were
compared against this repo's capacity curves earlier (gh_entropy comparison).
Here the symbol alphabet is not raw states but the gradient spectrum, and the
factor map is exact for fate precisely in the tau_a >= tau_p regime.

House Rules Compliance: exhaustive, deterministic, no RNG; all constructs
analytic (substrate/analysis boundary per the parent doc).
Output: result/topology/spectrum_automaton.npz + printed doc table.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
         (1, 2), (2, 3), (2, 4), (3, 4)]


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def sig(cfg, S):
    ph = [cfg[i] for i in CYC]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % (S + 1) for k in range(4)))


def fates(ta, tp):
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


def audit_cell(ta, tp):
    S = ta + tp
    fate = fates(ta, tp)
    edges, cls_fate = defaultdict(set), defaultdict(set)
    for c, f in fate.items():
        s = sig(c, S)
        edges[s].add(sig(step2(c, ta, tp), S))
        cls_fate[s].add(f)
    z = (0, 0, 0, 0)
    if cls_fate[z] != {False}:
        raise RuntimeError(f"z-class not uniformly dead at ({ta},{tp})")
    reach_z, changed = {z}, True
    while changed:
        changed = False
        for s, outs in edges.items():
            if s not in reach_z and outs & reach_z:
                reach_z.add(s)
                changed = True
    pure_live = [s for s, f in cls_fate.items() if f == {True}]
    pure_dead = [s for s, f in cls_fate.items() if f == {False}]
    mixed = [s for s, f in cls_fate.items() if len(f) == 2]
    return dict(ta=ta, tp=tp, n_classes=len(edges),
                n_live=len(pure_live), n_dead=len(pure_dead), n_mixed=len(mixed),
                complete=all(s in reach_z for s in pure_dead),
                sound=all(s not in reach_z for s in pure_live),
                mixed_reach_z=sum(1 for s in mixed if s in reach_z))


def main():
    recs = []
    print(f"{'(ta,tp)':>8} {'regime':>7} {'cls':>4} {'live':>5} {'dead':>5} "
          f"{'mixed':>6} {'complete':>9} {'sound':>6} {'mixed->z':>9}")
    for ta, tp in CELLS:
        r = audit_cell(ta, tp)
        recs.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"{str((ta,tp)):>8} {regime:>7} {r['n_classes']:4d} {r['n_live']:5d} "
              f"{r['n_dead']:5d} {r['n_mixed']:6d} {str(r['complete']):>9} "
              f"{str(r['sound']):>6} {r['mixed_reach_z']:4d}/{r['n_mixed']:<4d}",
              flush=True)

    assert all(r["complete"] for r in recs), "completeness violated"
    assert all(r["sound"] for r in recs), "soundness violated"
    assert all(r["mixed_reach_z"] == r["n_mixed"] for r in recs), \
        "a mixed class fails to reach z"
    assert all((r["n_mixed"] == 0) == (r["ta"] >= r["tp"]) for r in recs), \
        "mixed classes do not track the regime boundary on these cells"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "spectrum_automaton.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | classes | pure live | pure dead | mixed | "
          "dead⟹reach z | live⟹no path to z | mixed⟹reach z |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: | :---: | :---: | :---: |")
    for r in recs:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | {r['n_classes']} | "
              f"{r['n_live']} | {r['n_dead']} | {r['n_mixed']} | yes | yes | "
              f"{r['mixed_reach_z']}/{r['n_mixed']} |")
    print("\nfate == z-reachability on pure classes in 12/12 cells; "
          "mixed classes exist exactly at tau_p > tau_a on these cells")


if __name__ == "__main__":
    main()
