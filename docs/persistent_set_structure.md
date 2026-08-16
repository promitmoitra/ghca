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

- **A rejected refinement — argument corrected in §10.** Appending the sorted
  active/passive/receptive type pattern was tried and found insufficient. The
  *original* argument here claimed it "is not a refinement at all" because
  impure **class counts** rose (at `(3,4)`: 8 → 28). That argument was wrong:
  appending data to a signature always refines the partition, and an impure
  class *count* can legitimately rise under refinement — one impure class
  splits into several, some still impure. The metric that cannot rise is the
  number of **configurations in impure classes**, and by that correct metric
  the augmentation *was* a genuine refinement that merely remained
  insufficient. §10 redoes this properly with a full refinement ladder and
  locates the exact sufficiency threshold.
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

## 9. Addendum: the signature generalises to 3x3 — 40M configurations

At 3x3 the reentrant cycle is not unique (four 2x2 plaquettes + the perimeter
8-cycle), so the signature needs choosing. The **plaquette signature** — the 2x2
gap multiset per plaquette, then the sorted multiset of the four plaquette codes
— is the one that works.

Experiment: [`experiments/persistent_set_3x3.py`](../experiments/persistent_set_3x3.py) ·
Archive: `result/topology/persistent_set_3x3.npz`

| (τa, τp) | configs | P | plaquette classes | impure | compression | invariant |
| :---: | ---: | ---: | ---: | :---: | ---: | :---: |
| (1, 1) | 19,683 | 0.7462 | 45 | 0 | 437x | yes |
| (2, 1) | 262,144 | 0.9746 | 319 | 0 | 822x | yes |
| (1, 2) | 262,144 | 0.4580 | 319 | 173 | 822x | **NO** |
| (2, 2) | 1,953,125 | 0.7752 | 1,355 | 0 | 1441x | yes |
| (3, 2) | 10,077,696 | 0.9477 | 5,830 | 0 | 1729x | yes |
| (2, 3) | 10,077,696 | 0.5549 | 5,830 | 2,378 | 1729x | **NO** |
| (3, 3) | 40,353,607 | 0.7826 | 19,251 | 0 | 2096x | yes |

- **The τp ≤ τa rule transfers exactly**, including the boundary cell (3,3):
  **40,353,607 configurations, 19,251 classes, zero impure** — every class
  uniformly persistent or uniformly dying. Failures are precisely the τp > τa
  cells, as at 2x2. (The rule stays one-directional: 2x2's (2,5) already showed
  τp > τa does not *imply* failure.)
- **The perimeter signature fails even at (1,1)** (4 impure of 15 classes): the
  invariant lives on the shortest reentrant loops, not the long way around —
  consistent with the τ < L sustain gate, where the 4-cycle is the tight loop.
- **Compression scales**: 16–161x at 2x2 → 437–2096x at 3x3. Classes grow
  polynomially while configurations grow as (S+1)^9, so P becomes computable
  from ~19k classes instead of ~40M configurations — no GPU needed after all:
  the whole sweep runs in ~2 minutes and <0.5 GB via a successor array, boolean
  label propagation to the all-zero sink (valid because θ=1 GH has no non-zero
  dead attractor — checked exhaustively at 2x2, and the pipeline is checked
  against direct orbit iteration at 3x3 (1,1), both inside the script).
  **That validity condition is now a theorem** for every graph and every
  `(tau_a, tau_p)` — Theorem Z of [`coherent_core.md`](coherent_core.md) — so
  the propagation no longer rests on a 2x2 check.

Scope: open-boundary 3x3 box, θ = 1, von Neumann, cells up to S+1 = 7. The
(4,3)/(3,4) pair (S+1 = 8, 134M configs) is the next affordable check.

## 10. Addendum: the discrete-calculus ladder — where does decidability begin?

The gap field is a discrete-exterior-calculus object: it is the **discrete
gradient** of cell phase along the reentrant loop; its loop-sum is the
**discrete curl**, equal to `(S+1) × w` where `w` is the GGH winding number
(the invariant of PR #54); and edge-differences give a **discrete divergence**
(a phase Laplacian). Each rung is tested as a persistence invariant.

Experiment: [`experiments/persistent_set_dec.py`](../experiments/persistent_set_dec.py) ·
Archive: `result/topology/persistent_set_dec.npz`

**The curl is well-defined but insufficient.** The telescoping identity
`Σ gaps ≡ 0 (mod S+1)` holds with zero violations in every cell tested, so `w`
is always defined. `w = 0` is uniformly dead everywhere (an irrotational phase
field cannot circulate) — but intermediate windings are **mixed** in every cell
tested, including τp ≤ τa ones. The winding classifies which loops are alive on
a *running* trajectory; it under-determines which initial conditions get there.
The §3 signature is strictly finer than the curl.

**Relative phase, in its entirety, is insufficient at τp > τa.** Div-based
refinements and even the full canonical gap *arrangement* (rotation-normalised,
subsuming every invariant derivable from the gap field alone) leave impure
classes at every failing cell:

**(1, 2)**, N = 256:

| invariant | classes | impure classes | configs in impure |
| :--- | ---: | ---: | ---: |
| grad spectrum (the §3 signature) | 10 | 2 | 96 |
| + sorted div spectrum | 20 | 4 | 64 |
| grad arrangement (canonical rotation) | 20 | 4 | 64 |
| arrangement + active count | 49 | 2 | 16 |
| arrangement + absolute state multiset | 67 | 0 | 0 |
| canonical phases (loop-rotation orbit) | 70 | 0 | 0 |

**(2, 3)**, N = 1296:

| invariant | classes | impure classes | configs in impure |
| :--- | ---: | ---: | ---: |
| grad spectrum (the §3 signature) | 22 | 4 | 432 |
| + sorted div spectrum | 56 | 10 | 264 |
| grad arrangement (canonical rotation) | 58 | 10 | 240 |
| arrangement + active count | 176 | 8 | 96 |
| arrangement + absolute state multiset | 330 | 0 | 0 |
| canonical phases (loop-rotation orbit) | 336 | 0 | 0 |

**(2, 4)**, N = 2401:

| invariant | classes | impure classes | configs in impure |
| :--- | ---: | ---: | ---: |
| grad spectrum (the §3 signature) | 30 | 4 | 336 |
| + sorted div spectrum | 85 | 4 | 112 |
| grad arrangement (canonical rotation) | 88 | 4 | 112 |
| arrangement + active count | 279 | 2 | 32 |
| arrangement + absolute state multiset | 616 | 0 | 0 |
| canonical phases (loop-rotation orbit) | 616 | 0 | 0 |

**(3, 4)**, N = 4096:

| invariant | classes | impure classes | configs in impure |
| :--- | ---: | ---: | ---: |
| grad spectrum (the §3 signature) | 43 | 8 | 1152 |
| + sorted div spectrum | 130 | 18 | 608 |
| grad arrangement (canonical rotation) | 134 | 18 | 576 |
| arrangement + active count | 439 | 20 | 304 |
| arrangement + absolute state multiset | 1030 | 0 | 0 |
| canonical phases (loop-rotation orbit) | 1044 | 0 | 0 |

**The sufficiency threshold is `arrangement + absolute state multiset`** —
exact (zero configs in impure classes) at all four failing cells. Knowing *how
many* cells are active is not enough; the decider needs *which* absolute states
are present, i.e. how deep into the refractory band each passive cell sits.
And that threshold invariant is barely coarser than the loop-rotation symmetry
orbit itself (330 vs 336 classes at (2,3); 1030 vs 1044 at (3,4)): compression
collapses from ~50–160× to ~4×. **The τp ≤ τa boundary is a regime change in
description complexity** — above it, persistence is a relative-phase property;
below it, it is pinned to absolute refractory depth, and no loop-local
relative description works.

**Correction to §5's argument** (the finding stands, the reasoning didn't):
"impure class count rose" does not show something is not a refinement —
appending data always refines, and an impure class can split into several
impure pieces. The monotone quantity is *configurations in impure classes*,
used throughout this section (and asserted monotone down the ladder by the
script itself).

Scope: 2×2 core, θ = 1. The ladder at 3×3 (does absolute-state augmentation
also close the (1,2)/(2,3) failures there?) is untested.
