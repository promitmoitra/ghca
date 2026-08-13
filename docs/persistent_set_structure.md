# The persistent set has no linear description, but it has a combinatorial one

**Status:** new result, external-facing. Structural companion to Moitra & Sen
(SSRN 4047679), which defines the persistence probability `P(tau_a, tau_p)` — the
fraction of initial configurations that sustain oscillations asymptotically — and
maps its sigmoidal rise by exhaustive enumeration. That paper counts the
persistent set; this asks what **shape** it has, and how the answer scales.

Experiment: [`experiments/persistent_set_structure.py`](../experiments/persistent_set_structure.py) ·
Archive: `result/topology/persistent_set_structure.npz`

## 1. The question

Two candidate descriptions of the persistent set, tested exhaustively on the 2x2
core — the smallest sub-lattice supporting a reentrant 4-cycle, hence the smallest
system with a non-trivial `P`:

1. **Geometric.** Embed each configuration as a reduced one-hot 0/1 vector in
   `R^(4S)`, `S = tau_a + tau_p`. Is there a hyperplane separating persistent from
   dying configurations? If yes, an H-description exists and linear-programming
   methods apply.
2. **Combinatorial.** Read the four cells around the reentrant cycle `0-1-3-2` and
   take the sorted multiset of cyclic phase gaps mod `S+1` — invariant under
   rotating the wave and under relabelling the starting cell. Is persistence a
   *function* of that signature?

## 2. Geometric: no, never

**The persistent set is not linearly separable from the dying set in any of the
22 non-trivial cells tested.** Not one. So no polytope in one-hot coordinates can
serve as its H-description, and no LP formulation over raw configuration
coordinates can decide persistence.

Facet enumeration confirms this is hopeless rather than merely hard: at
`(tau_a, tau_p) = (1,2)` a 32-vertex hull already requires **7,732 facets** in
dimension 11, and the one-hot dimension grows as `4(tau_a + tau_p)` — reaching 40
at `(5,5)`. The convex hull of the persistent set is a generic 0/1 polytope with
essentially no exploitable linear structure.

## 3. Combinatorial: yes, whenever `tau_p <= tau_a`

Persistence **is** a function of the gap signature whenever `tau_p <= tau_a`. The
implication runs **one way only** — see the counterexample below.

| (τa, τp) | configs | persistent | P | one-hot dim | lin. separable | gap classes | impure | compression |
| :---: | ---: | ---: | ---: | :---: | :---: | ---: | :---: | ---: |
| (1, 1) | 81 | 24 | 0.296 | 8 | no | 5 | 0 | 16x |
| (1, 2) | 256 | 32 | 0.125 | 11 | no | 10 | 2 | 26x |
| (2, 1) | 256 | 180 | 0.703 | 12 | no | 10 | 0 | 26x |
| (2, 2) | 625 | 200 | 0.320 | 16 | no | 14 | 0 | 45x |
| (2, 3) | 1296 | 216 | 0.167 | 20 | no | 22 | 4 | 59x |
| (2, 5) | 4096 | 16 | 0.004 | 11 | no | 43 | 0 | 95x |
| (3, 3) | 2401 | 784 | 0.327 | 24 | no | 30 | 0 | 80x |
| (3, 4) | 4096 | 800 | 0.195 | 28 | no | 43 | 8 | 95x |
| (4, 4) | 6561 | 2160 | 0.329 | 32 | no | 55 | 0 | 119x |
| (4, 5) | 10000 | 2160 | 0.216 | 36 | no | 73 | 12 | 137x |
| (5, 5) | 14641 | 4840 | 0.331 | 40 | no | 91 | 0 | 161x |

(Abridged; the archive and the experiment's generated table carry all 25 cells.)

**Scaling — the answer to "how does it scale".** Where the invariant holds, the
signature compresses the configuration space **16x to 161x**, and the compression
*grows* with `S`: 14,641 configurations collapse to 91 homogeneous classes at
`(5,5)`. `P(tau_a, tau_p)` is then computable by weighting ~90 classes instead of
enumerating ~10^4 configurations. Since the configuration count grows as
`(S+1)^(L^2)` while the class count grows far more slowly, this is the lever for
pushing the paper's enumeration to larger cores.

## 4. Where it fails, and why

The dichotomy over 25 enumerated cells:

```
invariant, tau_p <= tau_a : (1,1) (2,1) (2,2) (3,1) (3,2) (3,3) (4,1) (4,2)
                            (4,3) (4,4) (5,1) (5,2) (5,3) (5,4) (5,5)
invariant, tau_p >  tau_a : (2,5)          <-- counterexample to a clean iff
fails      (all tau_p>ta) : (1,2) (2,3) (2,4) (3,4) (3,5) (4,5)
vacuous, P = 0 exactly    : (1,3) (1,4) (1,5)
```

So the supported claim is one-directional:

> `tau_p <= tau_a` **=>** persistence is a function of the gap signature.
> The converse is **false**: `(2,5)` is invariant with `tau_p = 5 > tau_a = 2`.

**Mechanism.** The signature records phase *differences* but not which absolute
states are active. When `tau_p > tau_a` the active band `1..tau_a` is narrower
than the passive band, so two configurations with identical gaps can hold
different numbers of *active* cells and excite their receptive neighbours
differently. Measured at `(2,3)`, signature `(0,3,4,5)`: `(0,0,1,3)` persists with
1 active cell, `(0,0,3,4)` dies with 0. When `tau_p <= tau_a` the active band is
wide enough that the gap multiset already pins the active count.

## 5. Two negative results worth recording

- **A rejected refinement.** Appending the sorted active/passive/receptive type
  pattern to the signature was tried and **fails**: sorting the type pattern
  independently of the gaps destroys their correspondence, so it is not a
  refinement at all and it *raised* impurity (at `(3,4)`: 8 impure classes became
  28). A genuine refinement must keep types aligned with gap positions. Finding
  one is open.
- **A threshold that hid a counterexample.** An earlier pass used a
  `P <= 0.01` "trivial" cutoff and thereby excluded `(2,5)` — P = 0.0039, but 16
  genuinely persistent configurations — and so reported a clean iff the data does
  not support. Only `P == 0` exactly is now excluded, where the statement is
  vacuous. A low-`P` cell is still a cell.

## 6. Validation

- **Dynamics.** The step function is checked against the repo's own
  `ghca_net.Network` (`act=tau_a, pas=tau_p, theta=1, p_s=0`) over **every**
  transition of **every** configuration in 9 cells: **7,461 transitions, 0
  mismatches**. The experiment aborts if this ever fails.
- **Paper's trend.** The reproduced `P` surface is monotone **increasing in
  `tau_a`** at every fixed `tau_p` — the sigmoidal rise the paper reports.
- **No RNG.** Every number is a complete enumeration, so there is nothing to seed
  and no per-seed spread to report. Deterministic by construction.

## 7. Caveats

- **The 2x2 core only.** Whether the gap signature generalises to 3x3 (where the
  reentrant cycle is no longer unique and multiple loops can coexist) is untested,
  and it is the obvious next question. 3x3 needs `(S+1)^9` configurations — 40M at
  `S+1=7` — so it wants a compiled or GPU kernel, or the symmetry quotient applied
  *before* enumeration rather than after.
- **`theta = 1`, von Neumann adjacency.** Both fixed. The paper's core analysis
  uses the same, but the invariant is untested against `theta > 1`.
- **The negative is coordinate-specific.** "Not linearly separable" is a statement
  about the one-hot embedding. A different embedding might separate; the result
  rules out the natural deep-learning-style flattening, not all embeddings.
- **No dynamics claim.** The signature is an *analysis* construct. Nothing here
  says the substrate computes, represents, or is sensitive to it.

## 8. Addendum: the polytope in *direct* phase-space coordinates

The natural phase space of the system is the product of chains
`Omega = {0..S}^4` — one axis per cell, a configuration is a point, the
dynamics move the point. Unlike the one-hot embedding (dimension `4S`), this
space has dimension 4 regardless of `S`, so **exact facet enumeration is
tractable** and the polytope question can be answered directly.

Experiment: [`experiments/phase_space_polytope.py`](../experiments/phase_space_polytope.py) ·
Archive: `result/topology/phase_space_polytope.npz`

| (τa, τp) | persistent | hull vertices | facets | dead inside hull | attractor pts | attr dim | attr pts on hull | monotone transients |
| :---: | ---: | ---: | ---: | ---: | ---: | :---: | ---: | ---: |
| (1, 1) | 24 | 24 | 108 | 13 | 8 | 3 | 8 | 8/16 |
| (2, 1) | 180 | 48 | 220 | 44 | 84 | 4 | 12 | 8/96 |
| (2, 2) | 200 | 40 | 202 | 229 | 40 | 4 | 0 | 24/160 |
| (3, 1) | 530 | 64 | 304 | 63 | 330 | 4 | 16 | 12/200 |
| (3, 2) | 786 | 48 | 213 | 378 | 378 | 4 | 12 | 48/408 |
| (3, 3) | 784 | 48 | 248 | 1131 | 224 | 4 | 0 | 64/560 |
| (4, 2) | 1806 | 64 | 284 | 463 | 1078 | 4 | 16 | 68/728 |
| (4, 4) | 2160 | 48 | 242 | 3355 | 720 | 4 | 0 | 112/1440 |

Four facts, all exhaustive:

1. **The hull is small and its complexity saturates.** Vertex count stalls near
   48 and facets near 240–250 while the persistent count grows two orders of
   magnitude — an `O(1)`-complexity **outer** description. Membership is a valid
   *necessary* condition for persistence. (Facet counts are Qhull's triangulated
   output and can shift by a few with point ordering; vertex and dead-inside
   counts are order-independent, the latter re-verified by exact LP.)
2. **It is not sufficient, and the gap grows.** Dead configurations strictly
   inside the hull outnumber persistent ones from (2,2) onward (3355 vs 2160 at
   (4,4)). The persistent set is **not** the integer-point set of any polytope in
   these coordinates either — convexity is the wrong closure, consistent with §2
   seen from the other side. The sharp description remains the gap signature (§3).
3. **At (1,1) the attractor set is an exact hyperplane slice**: every attractor
   configuration satisfies `x0+x1+x2+x3 = 3` (verified in integer arithmetic) —
   the rotating wave conserves total phase, the attractor set is 3-dimensional,
   and all 8 attractor points are vertices of the persistent hull. This is the
   cleanest polytope statement in the system, and it is special: at
   (2,2), (3,3), (4,4) **zero** attractor points are hull vertices and attractor
   phase-sums span a range — attractors live in the *interior* of phase space.
4. **No linear Lyapunov function**: transients are not monotone in
   `|phase_sum − attractor mean|` (112/1440 at (4,4)), so approach to the
   attractor is genuinely non-convex in these coordinates.

Caveat: same scope as the rest of this doc — 2x2 core, `theta = 1`, von Neumann.
The saturating-hull observation is an empirical pattern over the 8 cells tested,
not a theorem.
