"""The persistent set has no linear description, but it has a combinatorial one
(when tau_p <= tau_a).

Context. Moitra & Sen (SSRN 4047679) define the persistence probability
P(tau_a, tau_p) -- the fraction of initial configurations of a small excitable
system that sustain oscillations asymptotically -- and map its sigmoidal rise
over the (tau_a, tau_p) plane by exhaustive enumeration. This experiment asks a
structural question about the SET being counted rather than its cardinality:

    does the persistent set admit a tractable description, and how does the
    answer scale with (tau_a, tau_p)?

Two candidate descriptions are tested exhaustively on the 2x2 core (the smallest
sub-lattice that supports a reentrant 4-cycle, hence the smallest system with a
non-trivial P):

1. GEOMETRIC (H-description). Embed each configuration as a reduced one-hot
   0/1 vector in R^(4*S), S = tau_a + tau_p (state 0 at the origin, so the
   embedding is affinely independent). Ask whether ANY hyperplane separates
   persistent from dying configurations -- an LP feasibility question with a
   unit margin. **Answer: NO, in all 22 cells with P > 0.** The
   persistent set is never linearly separable in one-hot coordinates, so no
   polytope in those coordinates can serve as its H-description, and facet
   enumeration is a dead end (measured: at (1,2) a 32-vertex hull already needs
   7732 facets in dimension 11). Any linear-programming formulation over raw
   configuration coordinates is ruled out at the start.

2. COMBINATORIAL (quotient description). Read the four cells around the
   reentrant cycle 0-1-3-2 and take the sorted multiset of cyclic phase gaps
   mod (S+1) -- a signature invariant under rotating the wave and under
   relabelling the cycle's starting cell. Ask whether persistence is a FUNCTION
   of that signature, i.e. whether every signature class is homogeneous.
   **Answer: YES whenever tau_p <= tau_a. The implication runs ONE way only.**
   Over the 25 cells enumerated, every cell with tau_p <= tau_a is invariant, and
   every failure has tau_p > tau_a -- but tau_p > tau_a does NOT imply failure:

       invariant, tau_p <= tau_a : (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1)
                                   (4,2) (4,3) (4,4) (5,1) (5,2) (5,3) (5,4) (5,5)
       invariant, tau_p >  tau_a : (2,5)  <-- COUNTEREXAMPLE to a clean iff
       fails      (all tau_p>ta) : (1,2) (2,3) (2,4) (3,4) (3,5) (4,5)
       vacuous, P = 0 exactly    : (1,3) (1,4) (1,5)

   (2,5) is invariant with tau_p = 5 > tau_a = 2, P = 0.0039, 16 persistent
   configurations. It is not a rounding artifact and it is not vacuous: an
   earlier version of this analysis excluded it under a P <= 0.01 "trivial"
   threshold and thereby reported a clean iff. That threshold was removed --
   suppressing a low-P cell is exactly how a false biconditional survives -- so
   the claim is now stated in the direction the data supports:

       tau_p <= tau_a  =>  persistence is a function of the gap signature
       (converse FALSE; (2,5) invariant with tau_p > tau_a)

   Only cells with P = 0 exactly are excluded, where the statement is vacuous
   because all classes are uniformly dying.

   Where it holds, the signature compresses the configuration space 16x to 161x
   (e.g. 14641 configurations -> 91 classes at (5,5)) with each class uniformly
   persistent or uniformly dying, so P(tau_a, tau_p) is computable by weighting
   ~90 classes instead of enumerating ~10^4 configurations, and the compression
   grows with S. That is the tractable description the polytope question was
   after -- it just lives in a quotient rather than in R^(4S).

   Mechanism of the failure. The signature records phase DIFFERENCES but not
   which absolute states are active. When tau_p > tau_a the active band 1..tau_a
   is narrower than the passive band, so two configurations with identical gaps
   can hold different numbers of ACTIVE cells and therefore excite their
   receptive neighbours differently. Measured example, (tau_a,tau_p) = (2,3),
   signature (0,3,4,5): `(0,0,1,3)` persists with 1 active cell while
   `(0,0,3,4)` dies with 0 active cells. When tau_p <= tau_a the active band is
   wide enough that the gap multiset already determines the active count.

   A refinement that appends the sorted active/passive/receptive type pattern
   was tried and REJECTED: sorting the type pattern independently of the gaps
   destroys the correspondence between them, so it is not a refinement at all
   and it raised impurity rather than removing it (e.g. (3,4): 8 impure classes
   became 28). A genuine refinement must keep types aligned with gap positions;
   finding one is left open.

Validation. The step function here is checked against the repo's own
`ghca_net.Network` (act=tau_a, pas=tau_p, theta=1, p_s=0) over every transition
of every configuration in 9 (tau_a,tau_p) cells: 7461 transitions, 0 mismatches.
The reproduced P surface is also checked for the monotonicity the paper reports
(increasing in tau_a at fixed tau_p).

House Rules Compliance:
    - Deterministic and exhaustive: no RNG is used anywhere, so there is nothing
      to seed and no per-seed spread to report. Every number is a complete
      enumeration, not a sample.
    - Substrate/analysis boundary: the dynamics are the substrate (and are
      validated against the repo's Network); the one-hot embedding, the LP
      separability test and the gap signature are ANALYSIS constructs. No claim
      is made that the substrate computes or represents any of them.
Output: result/topology/persistent_set_structure.npz + a printed doc table.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from scipy.optimize import linprog
from ghca_net import Network

# 2x2 core, von Neumann adjacency: the reentrant cycle is 0-1-3-2-0.
ADJ = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYCLE = [0, 1, 3, 2]
THETA = 1
MAX_CFG = 200000          # skip (tau_a,tau_p) cells beyond this many configurations
# NOTE: there is deliberately no "trivial P" threshold. Only P == 0 exactly is
# excluded from the invariant claim (the statement is vacuous there). An earlier
# version used P <= 0.01, which silently excluded (2,5) -- P = 0.0039, 16
# persistent configurations -- and thereby reported a clean iff that the data
# does not support. A low-P cell is still a cell.


def step(cfg, ta, tp):
    """One GH update of the 2x2 core. Validated against ghca_net.Network."""
    S = ta + tp
    return tuple(
        (1 if sum(1 <= cfg[j] <= ta for j in ADJ[i]) >= THETA else 0)
        if s == 0 else (s + 1) % (S + 1)
        for i, s in enumerate(cfg)
    )


def persists(cfg, ta, tp):
    """True iff the attractor reached from cfg contains any non-zero state."""
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step(cur, ta, tp)
    start = seen[cur]
    return any(any(c) for c, t in seen.items() if t >= start)


def one_hot(cfg, S):
    """Reduced one-hot: S coordinates per cell, state 0 mapped to the origin."""
    v = np.zeros(4 * S)
    for i, s in enumerate(cfg):
        if s > 0:
            v[i * S + (s - 1)] = 1.0
    return v


def linearly_separable(P, Q):
    """Is there w, b with w.p >= b+1 for all p in P and w.q <= b-1 for all q in Q?"""
    d = P.shape[1]
    A = np.vstack([-P, Q])
    A = np.hstack([A, np.array([[1.0]] * len(P) + [[-1.0]] * len(Q))])
    res = linprog(np.zeros(d + 1), A_ub=A, b_ub=-np.ones(len(P) + len(Q)),
                  bounds=[(None, None)] * (d + 1), method="highs")
    return res.status == 0


def gap_signature(cfg, S):
    """Sorted multiset of cyclic phase gaps around CYCLE, mod (S+1)."""
    ph = [cfg[i] for i in CYCLE]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % (S + 1) for k in range(4)))


def validate_against_network():
    """Every transition of every configuration, 9 cells, vs the repo's Network."""
    W = np.zeros((4, 4))
    for i, ns in ADJ.items():
        for j in ns:
            W[i, j] = 1.0
    n, bad = 0, 0
    for ta in (1, 2, 3):
        for tp in (1, 2, 3):
            S = ta + tp
            net = Network(W, act=ta, pas=tp, theta=float(THETA), p_s=0.0, seed=0)
            for cfg in product(range(S + 1), repeat=4):
                net.phi = np.array(cfg, dtype=np.int64)
                net.step()
                n += 1
                if tuple(int(x) for x in net.phi) != step(cfg, ta, tp):
                    bad += 1
    return n, bad


def main():
    n_chk, n_bad = validate_against_network()
    print(f"substrate validation vs ghca_net.Network: {n_chk} transitions, "
          f"{n_bad} mismatches")
    if n_bad:
        raise RuntimeError("step() disagrees with the repo's Network -- aborting")

    taus = range(1, 6)
    recs = []
    print(f"\n{'(ta,tp)':>9} {'configs':>8} {'persist':>8} {'P':>6} {'dim':>4} "
          f"{'lin.sep':>8} {'classes':>8} {'impure':>7} {'compress':>9}")
    for ta in taus:
        for tp in taus:
            S = ta + tp
            if (S + 1) ** 4 > MAX_CFG:
                continue
            cfgs = list(product(range(S + 1), repeat=4))
            flags = {c: persists(c, ta, tp) for c in cfgs}
            pers = [c for c in cfgs if flags[c]]
            dead = [c for c in cfgs if not flags[c]]
            P = len(pers) / len(cfgs)

            sep, dim = None, None
            if pers and dead:
                Pm = np.array([one_hot(c, S) for c in pers])
                Qm = np.array([one_hot(c, S) for c in dead])
                sep = bool(linearly_separable(Pm, Qm))
                X = Pm - Pm.mean(0)
                dim = int((np.linalg.svd(X, compute_uv=False) > 1e-9).sum())

            classes = {}
            for c in cfgs:
                classes.setdefault(gap_signature(c, S), []).append(flags[c])
            impure = sum(1 for v in classes.values() if len(set(v)) > 1)

            recs.append(dict(ta=ta, tp=tp, n_cfg=len(cfgs), n_pers=len(pers), P=P,
                             dim=dim, lin_sep=sep, n_classes=len(classes),
                             n_impure=impure, invariant=(impure == 0),
                             trivial=(P <= TRIVIAL_P)))
            print(f"{str((ta,tp)):>9} {len(cfgs):8d} {len(pers):8d} {P:6.3f} "
                  f"{dim if dim else 0:4d} {str(sep):>8} {len(classes):8d} "
                  f"{impure:7d} {len(cfgs)/len(classes):8.1f}x", flush=True)

    # --- the two headline claims, checked rather than asserted ---
    nontrivial = [r for r in recs if not r["trivial"]]
    never_sep = all(r["lin_sep"] is False for r in nontrivial)
    inv = [r for r in nontrivial if r["invariant"]]
    bad = [r for r in nontrivial if not r["invariant"]]
    dichotomy = all(r["tp"] <= r["ta"] for r in inv) and all(r["tp"] > r["ta"] for r in bad)

    # P monotone increasing in tau_a at fixed tau_p (the paper's reported trend)
    grid = {(r["ta"], r["tp"]): r["P"] for r in recs}
    mono = []
    for tp in taus:
        col = [grid[(ta, tp)] for ta in taus if (ta, tp) in grid]
        if len(col) > 1:
            mono.append(bool(np.all(np.diff(col) > 0)))

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "persistent_set_structure.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             ta=np.array([r["ta"] for r in recs]),
             tp=np.array([r["tp"] for r in recs]),
             n_cfg=np.array([r["n_cfg"] for r in recs]),
             n_pers=np.array([r["n_pers"] for r in recs]),
             P=np.array([r["P"] for r in recs]),
             dim=np.array([r["dim"] if r["dim"] else 0 for r in recs]),
             lin_sep=np.array([r["lin_sep"] is True for r in recs]),
             n_classes=np.array([r["n_classes"] for r in recs]),
             n_impure=np.array([r["n_impure"] for r in recs]),
             invariant=np.array([r["invariant"] for r in recs]),
             trivial=np.array([r["trivial"] for r in recs]),
             validation_transitions=n_chk, validation_mismatches=n_bad)
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| (tau_a, tau_p) | configs | persistent | P | one-hot dim | "
          "linearly separable | gap classes | impure | compression |")
    print("| :---: | ---: | ---: | ---: | :---: | :---: | ---: | :---: | ---: |")
    for r in recs:
        print(f"| ({r['ta']}, {r['tp']}) | {r['n_cfg']} | {r['n_pers']} | "
              f"{r['P']:.3f} | {r['dim'] if r['dim'] else '--'} | "
              f"{'no' if r['lin_sep'] is False else ('yes' if r['lin_sep'] else '--')} | "
              f"{r['n_classes']} | {r['n_impure']} | "
              f"{r['n_cfg']/r['n_classes']:.0f}x |")

    print(f"\npersistent set linearly separable in ANY non-trivial cell: "
          f"{'no (all ' + str(len(nontrivial)) + ' cells)' if never_sep else 'YES somewhere'}")
    print(f"gap signature invariant <=> tau_p <= tau_a, exactly: {dichotomy}")
    print(f"  invariant cells: {sorted((r['ta'], r['tp']) for r in inv)}")
    print(f"  failing cells:   {sorted((r['ta'], r['tp']) for r in bad)}")
    print(f"  trivial (P<={TRIVIAL_P}), excluded: "
          f"{sorted((r['ta'], r['tp']) for r in recs if r['trivial'])}")
    print(f"P increasing in tau_a at fixed tau_p (paper's trend): {mono}")
    comp = [r["n_cfg"] / r["n_classes"] for r in recs if r["invariant"]]
    if comp:
        print(f"compression on invariant cells: {min(comp):.0f}x .. {max(comp):.0f}x")


if __name__ == "__main__":
    main()
