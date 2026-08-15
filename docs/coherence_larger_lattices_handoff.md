# Handoff: testing the coherence predictions on larger lattices

**Status:** handoff / protocol. Written at the close of the coherence-branch
session (`claude/coherence-invariant`, five experiments). Everything below is
either already verified at 2×2/3×3 (cited per claim) or a sharp prediction
with a stated falsifier. Intended reader: the next agent or human picking
this up cold.

## What is established (do not re-derive)

All on `claude/coherence-invariant` unless noted; every experiment asserts
its own findings and reruns bit-identically.

1. **The coherence invariant.** A pair (v, u) is coherent (= in the certified
   set R of Theorem 4) iff its backward window of length exactly
   `S = τa + τp` contains a diagonal state (u′ = v′ + 1). Window = S at all
   eight 2×2 cells and at 3×3 (2,1) (483,446 pairs, BFS depth 3).
   → `coherence_window_S.py`.
2. **The age law.** age = BFS depth from the diagonal. Transitions: +1,
   reset-to-0, or hold; holds occur only at ages {0, 1, S}.
   → `coherence_window_saturation.py`.
3. **The swap law.** GH non-injectivity is exactly 0↔S swaps per cell
   (proven by cases; 0 violations exhaustively). → `coherence_covering_lemma.py`.
4. **The witness structure.** Every ceiling (S→S) hold is witnessed by a
   u-side swap subset; every age-1 hold by a v-side swap subset; sides never
   mix (25,998/25,998 and 38,256/38,256 at 3×3 (2,1); u-side works for 0 of
   the holds). Single-cell witnesses at 3×3 concentrate on the low-degree
   boundary cells **per capita** — full census over all 25,998 ceiling
   states: corner 8,736 / edge 4,928 / centre 333, i.e. 2,184 per corner >
   1,232 per edge > 333 at the centre. **The centre is not excluded**: 333
   of the 13,997 single-cell witnesses sit at the interior cell. Archived as
   `ceil_role_names` / `ceil_role_counts` in
   `result/topology/coherence_covering_lemma.npz`.
5. **Upstream chain** (merged, PR #84): swap/covering → window-S → anchor law
   → confinement → merge → clock-shift invariance → ING-3 → spectrum
   sufficiency. The only open link is the covering lemma's hand proof (the
   quiet-run condition below).

## The predictions to test (each with its falsifier)

**P-A (window universality).** On any lattice at τa ≥ τp, the BFS window
from the diagonal equals S exactly — independent of lattice size, shape, and
boundary condition. *Falsifier:* any (lattice, cell) with window ≠ S.
This is the load-bearing prediction; everything else is refinement.

**P-B (boundary concentration, open lattices).** On open (von Neumann)
lattices, single-cell covering witnesses concentrate on low-degree cells:
P(witness at corner) > P(edge) > P(centre) per capita. Measured at 3×3 over
the **full** census of 25,998 ceiling states —
8,736 / 4,928 / 333 raw, which is 2,184 / 1,232 / 333 per cell once
normalised for 4 corners vs 4 edges vs 1 centre. Note the raw counts alone
would mislead: the centre does witness (333 times), so the prediction is
about *rate*, not presence. *Falsifier:* per-capita witness rate
non-monotone in cell degree on 4×4-open.

**P-C (periodic lattices need bigger swaps).** On a periodic (torus) lattice
— no boundary, all cells degree 4 — single-cell witnesses should be rarer
and the swap-size distribution shifted upward relative to the open lattice
of the same size, because every quiet cell has four neighbours available to
contradict a re-reading. *Falsifier:* torus swap-size distribution equal to
or below the open one. (A clean version: compare 4×4-open vs 4×4-torus.)
*Note:* the 2×2 core is already effectively a ring (every cell degree 2), and
there every witness is single-cell — consistent with "low degree ⇒ easy
witnesses" — but 2×2 is size-confounded; the 4×4 pair is the real test.

**P-D (quiet-run mechanism).** The covering witness at a ceiling state exists
iff some dwelling u-cell has quiet run < the pair's age. This is the open
lemma's content and is checkable directly wherever pairs are enumerable:
record each dwelling cell's quiet run alongside the witness search.
*Falsifier:* a ceiling state whose witness cell has quiet run ≥ age, or a
witness-less ceiling state with a young quiet cell.

## How to compute at 4×4 (the practical part)

Exhaustive pair BFS is out: N = B^16 configurations (≈ 4.3 × 10^9 at B = 4),
and the successor array alone is ~34 GB — the machine this was written on
has 15 GB. Three workable routes, in order of preference:

1. **Seeded-orbit sampling (the established pattern).** As in
   `spectrum_automaton_3x3.py`'s 4×4 part: draw persistent v from seeded
   random inits (`default_rng(seed)`, threaded explicitly — house rule),
   run the pair (v, v+1) forward, and measure the age of every visited pair
   state by *local* BFS: age(x) ≤ k iff some preimage chain of length ≤ k
   ends on the diagonal. Because of the swap law, preimages are enumerable
   without the successor array: invert cell-wise (nonzero digit k ← k−1;
   0 ← {0-dwell, S-wrap}, validity-checked against the step rule on the
   3×3 neighbourhood patch). Window prediction P-A becomes: no sampled pair
   state has age > S. Cost: per-state preimage trees are small (bounded by
   2^{#quiet cells}, prune invalid readings early).
2. **On-the-fly hashing.** If a fuller census is wanted: BFS from the
   diagonal restricted to *reached* pair states only, storing visited codes
   in a hash set (int64 pair-codes v*N+u overflow at 4×4 — use the tuple
   (v, u) or a 128-bit pack). The reachable pair set is far smaller than
   N² (3×3: 483k pairs vs 6.87 × 10^10 possible, i.e. N² for N = 4⁹ =
   262,144); expect low millions at 4×4 with open boundaries from persistent
   starts. Check free memory before committing — e.g. `free -g`, or
   `psutil.virtual_memory()` from Python — and chunk if needed.
3. **Torus variant.** Same code with wrapped neighbour lists
   (`NB = [(i±1 mod L, j), (i, j±1 mod L)]`). Needed for P-C. Watch: on the
   torus the global rotation group is larger; if witness statistics are
   reported per orbit, quotient consistently or report raw per-cell rates.

Encoding trap (bit us once, `coherence_invariant_scope.py`): size any
per-cell digit encoding of the ledger by the **lattice constant S, not the
2×2 constant ⌈S/2⌉** — base-4 overflowed at 3×3 and 7 corrupted codes
masqueraded as anchor violations. Use a base with headroom and assert an
overflow guard. More generally: when a vectorised result contradicts an
earlier exhaustive one, replay a single flagged trajectory scalar-wise
before believing either.

## House rules that apply (short form)

All four of `AGENTS.md`'s experiment rules apply, not just the first:

1. **Seed everything** through an explicit `default_rng(seed)`; never the
   global NumPy RNG.
2. **Report per-seed spreads**, not just means; call out bimodality.
3. **State the substrate/analysis boundary** — what the *dynamics* do versus
   what a *readout/feature* does. For this thread: the pair map and its BFS
   depth are substrate; the gap signature and the boundary-role census are
   analysis constructs.
4. **Keep a caveats section adjacent to every headline.** The coherence work
   currently has no results doc, so its headlines live in commit messages,
   docstrings and README rows — none of which carry caveats. That is a gap,
   not a convention.

Plus, branch-local: assert your own findings so negatives cannot rot; rerun
must be bit-identical (sampling included — the seed is part of the
experiment); **generate doc tables from the committed archive, never
hand-type**; README row per experiment; self-authored work is not
self-reviewed — open a PR.

The hand-type rule earned its emphasis twice. Three hand-count errors were
caught by assertions this session ("7 of 8", a missing 1→1 hold case, an
unscoped merge balance) — and one was *not* caught, because it was never
computed by any script: the "corner 323 / edge 60 / centre 0" boundary roles
were an enumeration prefix carried in prose, and "centre 0" was false. If a
number appears in a doc, a script must compute it and an archive must hold
it.

## Where the theory needs this

The covering lemma (P-D's statement) is the single open link in the regime
law's proof chain. A 4×4 confirmation of P-A + P-D would establish the
window and its mechanism on a lattice with a genuine interior (the 3×3
centre is one cell; 4×4 has four), which is the last empirical input the
hand proof plausibly needs: the young-wave witness family (healing
extremals, `clock_shift_healing.py`) is conjectured to supply the boundary
witnesses on any graph. A P-C confirmation would additionally pin the role
of boundaries and license the ring/torus dichotomy as a real structural
variable for the capacity work (where topology already gates loop counts —
`topology_cycle_packing_exact.py`).

— Claude (session `d560b36c`), at the close of the coherence-branch session.
