# The coherence invariant: a bounded-memory condition on pairs of trajectories

**Status:** theory. This doc exists because the coherence work previously had
none — six commits established a named invariant and a fourth "theorem" while
touching only `README.md`, one experiment and one archive each, so the claims
lived in commit messages and docstrings where no review pass could see them.
Every number below is read from a committed `.npz`; none is hand-typed.

**Read [Honest scope](#honest-scope) before quoting anything here.** The single
most-misread point: *the window law is a **regime** law*. It holds at τa ≥ τp
and is false outside it.

Setting: 2×2 and 3×3 cores, states 0 (receptive), 1..τa (active), τa+1..S
(passive), `S = τa + τp`, θ = 1, open boundaries. A **pair** is (v, u) of live
configurations; a **diagonal** (clock-shift) state is one with u = v + 1
cell-wise mod (S+1). **Age** = BFS depth of a pair from the diagonal set.

---

## The invariant

> A pair (v, u) is **coherent** — i.e. in the certified set R of the anchor-law
> certificate — iff its backward window of length **exactly S** contains a
> diagonal state.

Equivalently `R = ⋃_{k ≤ S} F^k(diag)`, and the closure saturates at depth
exactly S: nothing older than one cycle matters, and no shorter window suffices.

Source: [`coherence_window_S.py`](../experiments/coherence_window_S.py),
archive `coherence_window_S.npz`.

### Where it holds (τa ≥ τp)

| cell | S | \|R\| | window |
| :---: | ---: | ---: | ---: |
| (1,1) | 2 | 48 | 2 |
| (2,1) | 3 | 276 | 3 |
| (2,2) | 4 | 360 | 4 |
| (3,1) | 4 | 730 | 4 |
| (3,2) | 5 | 1,194 | 5 |
| (3,3) | 6 | 1,344 | 6 |
| (4,3) | 7 | 3,400 | 7 |
| (4,4) | 8 | 3,600 | 8 |
| 3×3 (2,1) | 3 | 483,446 | 3 |
| 4×4 open (2,1) | 3 | sampled (300 inits) | 3 |
| 4×4 torus (2,1) | 3 | sampled (300 inits) | 3 |

### Where it fails (τp > τa) — asserted negative control

| cell | S | window |
| :---: | ---: | ---: |
| (2,3) | 5 | 9 |
| (2,4) | 6 | 12 |
| (3,4) | 7 | 15 |
| (3,5) | 8 | 17 |
| (4,5) | 9 | 20 |

These run as a negative control with `assert w != S`, so the regime boundary
cannot silently widen. (1,2) matches at 2×2 by coincidence and dies at 3×3
— window 10 vs S = 3 — so it is excluded from the negative set.

### What is actually new here

R was *defined* as the forward closure of the diagonal seeds, so "coherent iff
within S steps of a diagonal" is close to a restatement of the construction.
The load-bearing content is two facts that are not definitional:

1. the closure saturates at depth **exactly** S (states first appear at every
   depth 1..S, and BFS terminates at S);
2. **R ⊊ lockstep.** At (2,1): 468 lockstep-bound pairs, 276 coherent, **192
   impostors**. Lockstep is *necessary* (|R − P| = 0) but *not sufficient*.
   The 192 counterexamples are archived as `extras_v`/`extras_u` in
   `coherence_formulate_2x2.npz` — a permanent test set for any future
   candidate invariant.

---

## The supporting laws, by epistemic status

The distinction below is the point of this section. This thread contains
**one** proof.

### PROVEN (a mathematical argument, lattice-free)

**The swap law.** GH non-injectivity is exactly 0↔S swaps per cell. If
`stp(w) = stp(w')` with `w ≠ w'`, both differing entries cannot be nonzero
(nonzero states advance +1 mod B, an injection); so one is 0 and the other
k ≠ 0; equality of images forces k = S. The argument uses only the advance rule
and the two routes to an image-0, so it holds on **any graph**.

Note the direction: non-injectivity *implies* 0↔S differences, one way only —
which is why every witness search must re-check image preservation.

*Caveat:* the accompanying exhaustive 0-violation check is **2×2 only**
(`product(range(B), repeat=4)`). There is no 3×3 exhaustive swap-law check,
even though the 3×3 witness census assumes swaps are the only preimage freedom.

### CERTIFIED (exhaustive over a stated domain — not proofs)

**The age law.** Age transitions are +1, reset-to-0, or hold; holds occur only
at ages {0, 1, S}. Certified at `CENSUS_CELLS = [(2,1), (3,3)]` only, both
τa ≥ τp. Every hold is a younger preimage merging in (12/12 and 40/40 at the
ceiling). — `coherence_window_saturation.npz`.

**The witness structure**, at 3×3 (2,1), full census:

| direction | site | witnessed |
| :--- | :--- | ---: |
| u-side | ceiling (S→S) | 25,998 / 25,998 |
| v-side | age-1 holds | 38,256 / 38,256 |
| u-side | age-1 holds | **0** / 38,256 |
| v-side | ceiling | **0** / 25,998 |

At 4×4 open (2,1) (sampled over 300 random trajectories, seed 12345):
- Ceiling holds: 111 / 115 u-side witnessed (**0** on v-side)
- Age-1 holds: 235 / 246 v-side witnessed (**0** on u-side)
- Swap size distribution at ceiling: 43 / 30 / 25 / 9 / 4 for sizes 1..5

At 4×4 torus (2,1):
- Ceiling holds: 22 / 22 u-side witnessed (**0** on v-side)
- Age-1 holds: 57 / 58 v-side witnessed (**0** on u-side)

So *sides never mix*, confirmed at 2×2, 3×3, and 4×4. — `coherence_covering_lemma.npz`, `coherence_window_4x4.npz`.

**"Theorem 4" (the anchor law) is not a theorem — there is no proof.** It is an
exhaustive check over a 276-state forward closure at (2,1): the four predicates
hold on that finite set. Its extension to 3×3 is a **re-enumeration** over
483,446 pair states, converging in 4 steps — the certificate *method* scales;
no argument transfers the result. The name is historical and the README rows
now say CERTIFIED. — `anchor_law_certificate.npz`, `coherence_invariant_scope.npz`.

**The compression barrier.** Seven local ledger invariants all leak, so no
purely local ledger invariant certifies the anchor law:

| I0 | I1 | I2 | I3 | I4 | I5 | I6 |
| ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 3,728 | 2,316 | 2,348 | 1,860 | 1,676 | 1,776 | 1,216 |

This is a deliberate negative result with a guard that fires if any invariant
ever closes. (I6 carries sig-live endpoints plus the I1 edge constraint only —
not the zero/merged/B4 constraints; an earlier label said "all".)

### OBSERVED (a pattern, not yet a law)

**Per-capita boundary concentration.** Single-cell witnesses concentrate
on low-degree cells *per capita*:
- At 3×3 (full census): corner 8,736 / edge 4,928 / centre 333 = 2,184 per corner > 1,232 per edge > 333 at the centre.
- At 4×4 open (300 inits): corner 14 / edge 20 / centre 9 = **3.50 per corner > 2.50 per edge > 2.25 at the centre**.

**The centre is not excluded**: 333 of the 13,997 single-cell witnesses at 3×3 and 9 of the 43 single-cell witnesses at 4×4 sit at the interior cells.

The count is convention-dependent, and **the ordering survives both
conventions** — worth stating explicitly, since convention ambiguity is what
produced the original error:

| convention | corner | edge | centre | per capita |
| :--- | ---: | ---: | ---: | :--- |
| one role per single-cell witness (archived) | 8,736 | 4,928 | 333 | 2,184 > 1,232 > 333 |
| every cell of every witnessing subset | 24,768 | 16,792 | 2,079 | 6,192 > 4,198 > 2,079 |

(Second row independently re-derived on the `agent-comms-log` response to this
review and reproduced here; not archived, as the first is the one the
prediction is stated over.) An earlier revision reported "corner 323 / edge 60 / centre
0" as a sample; those are the first 2,000 ceiling states in BFS enumeration
order — an order-biased prefix, not a sample — and "centre 0" was false. The
full census costs ~1.3 s. Now archived as `ceil_role_names`/`ceil_role_counts`.

---

## What this does *not* explain

The mechanism behind S is **open**. The covering lemma
([`coherence_covering_lemma.py`](../experiments/coherence_covering_lemma.py))
suggests that quiet cells can disguise their age only up to their quiet run, so
ancestry decays through the refractory pipeline at one tick per step — but that
lemma is explicitly labelled "THE REMAINING HAND OBLIGATION" in its own source.
It is a conjecture, not an account.

An earlier rationale — that a diagonal state recurs on every lockstep cycle, so
S is one excursion+refractory span — was **falsified** twelve minutes after it
was written: on-orbit recurrence gaps are 5/4/6/5/7/9/10/12, strictly above S
at all 8 cells. It is marked superseded in place in the source.

---

## Substrate / analysis boundary

- **Substrate:** the pair map and its BFS depth. Computed by exhaustive
  enumeration of the transition rule. No readout, no learned parameter, no
  feature — this is the dynamics and nothing else.
- **Analysis:** the gap signature, the boundary-role census, the ledger
  invariants (I0–I6) and the age function itself. These are constructs imposed
  on the enumeration to describe it; the substrate does not compute them.
- **Not measured here at all:** anything behavioural. No task, no decoder, no
  accuracy.

This matters for the bridge in [`synthesis.md`](synthesis.md): the coherence
window is pure substrate, while E2's retention number is a readout accuracy.
They share a gating quantity (E2's `τ = act + pas = S`) but **not a span** —
E2's retention is *unbounded* once τ < L, whereas coherence says history older
than S is *irrelevant*. The link is an analogy of structure, is not derived,
and nothing in the proof chain consumes it.

---

## Honest scope

- **The window law is a regime law.** Every cell in the positive table has
  τa ≥ τp. It is false at τp > τa (table above). Commit `06b4946`'s title said
  "every cell, both lattices"; that was wrong and is corrected in the README.
- **One 3×3 cell for most claims.** The witness census, the anchor-law
  re-certification and the 3×3 window all rest on (2,1). The window has since
  been spot-checked at 3×3 (1,1), (2,2) and (3,1) and holds, but those runs are
  not yet committed as assertions.
- **No proof except the swap law.** Everything else is exhaustive enumeration
  over 8 of the project's 12 canonical 2×2 cells plus one or four 3×3 cells.
  A certificate over a finite closure is not a theorem about the family.
- **The swap-law exhaustive check is 2×2 only**, while the 3×3 census assumes it.
- **≤1 merged is a 2×2 artifact** — up to 6 merged at 3×3.
- **P-D is not yet well-posed.** "Age" is a property of a pair *state*; "quiet
  run" is a property of a *trajectory*, and ~5% of reachable pairs at 3×3 admit
  both readings for some cell. See the handoff.
- **Self-falsifications this thread reports against itself**, all with archived
  evidence: on-orbit recurrence (8/8 cells), orbit co-membership (both
  directions), lockstep sufficiency, the compression barrier, the clock-shift
  quotient (60/8,100 mixed keys), ≤1-merged at 3×3, the streak skeleton, the
  live-SFT entropy prediction, and spectrum-automaton soundness. Where a claim
  was superseded, the superseding evidence is asserted in code so it cannot rot.

## Reproduce

```bash
python3 experiments/coherence_window_S.py           # window = S, + negative control
python3 experiments/coherence_window_saturation.py  # age law, on-orbit falsification
python3 experiments/coherence_covering_lemma.py     # swap law, witness census, roles
python3 experiments/coherence_invariant_scope.py    # 3x3 re-certification
python3 experiments/coherence_window_4x4.py         # 4x4 open & torus scaling (P-A, P-B, P-C)
python3 experiments/anchor_law_certificate.py       # "Theorem 4" + compression barrier
python3 experiments/coherence_formulate_2x2.py      # lockstep necessary, not sufficient
```

All are deterministic and RNG-free; each asserts its own findings and rewrites
its archive under `result/topology/`. Total runtime is a few seconds each.
