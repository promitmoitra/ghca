# Pillars filing test — does a vertical cut of this repo actually partition its claims?

**Status:** analysis, not a decision. Run 2026-08-11 on `claude/planning-review-branches-sf4y29`.

The proposal under test: reorganise the project's statement docs into **pillars**
(verticals), each a small tree — `vision → track → tactic → claim` — where tiers
revise at their own cadence and coupling between tiers is by typed edges rather
than containment. See [`process.md`](process.md) for the existing review/plan split
this would absorb.

Rather than argue the taxonomy, this test **files the claims the repo has already
made** into a candidate cut and counts the failures. Method and thresholds are
declared before the counts, below.

---

## Method

**Unit.** A *claim* = one falsifiable statement with a status and an evidence
pointer. The repo already writes these as `## Result N — <statement>` /
`## Negative N — <statement>` headings, so the extraction is mechanical:

```bash
grep -h '^#\{1,3\} \(Result\|Negative\|Finding\) ' docs/*.md            # 29 (E/C-series)
grep '^## \(Result\|Negative\)' docs/lattice_results.md                 # 13 (lattice arc)
```

**n = 42 claims.** 29 from the E-series and C-series results docs, 13 from the
lattice arc. The arc landed on `main` in #76 (merge commit `0580fef`) while this
test was being written; the counts below were taken from
`origin/claude/lattice-timescale-demo` before the merge and are unchanged by it —
the merge preserved all 27 commits, so both the claims and their per-commit audit
trail are now on `main`.

**Thresholds, declared in advance.**

| Failure mode | Meaning | Threshold |
|---|---|---|
| Homeless | claim fits no pillar | any ⇒ cut is missing a pillar |
| Ambiguous | two equally good homes | >10% ⇒ boundary is in the wrong place |
| Thin | pillar receives <3 claims | that's a section, not a pillar; merge upward |

**Tie-break rule** (introduced *after* the 6-pillar count, applied only to the
4-pillar count, and flagged wherever used): *file by what the claim constrains,
not by what produced it.*

---

## Run 1 — six pillars: Substrate / Rule / Anatomy / Interface / Instrumentation / Process

**Result: FAILS.** 0 homeless, **13/42 ambiguous (31%)**, three thin pillars
(Anatomy 2, Interface 2, Process 0).

Where the 13 ambiguities cluster — this is the diagnostic, not the count:

| Boundary | Ambiguous claims | Reading |
|---|---|---|
| Rule ↔ Anatomy | 5 | Plasticity and the structure it grows are not separable in this project |
| Substrate ↔ Instrumentation | 3 | C-series claims are causal *statements about the substrate*; the instrument is provenance, not subject |
| Substrate ↔ Interface | 2 | genuine cross-pillar claims |
| Substrate ↔ Rule | 1 | architecture-necessity claims |
| Substrate ↔ Anatomy | 1 | architecture-necessity claims |
| Rule ↔ Process | 1 | the substrate/analysis boundary is a method rule |

The Rule↔Anatomy cluster is the predicted failure and the largest. The
Substrate↔Instrumentation cluster is a different kind of error: **Instrumentation
is not a subject pillar at all.** It answers *how you know*, which is the same
question Process answers. Both are method.

## Run 2 — four pillars: Substrate / Learner / Interface / Method

Merging Rule+Anatomy → **Learner**, and Instrumentation+Process → **Method**,
resolves 8 of the 13 ambiguities structurally.

**Result: PASSES.** 0 homeless, **5/42 needed the tie-break rule (12%)**, no
pillar under 4 claims.

| Pillar | Claims | Share | What it holds |
|---|---:|---:|---|
| **Substrate** | 12 | 29% | what the excitable medium does: τ-control, nucleation/persistence, timescale nesting, the avoidance theorem, tonic drive |
| **Learner** | 20 | 48% | the τ rule, three-factor gating, value chains, emergent conjunction bases, layers-vs-edges, identity plasticity |
| **Interface** | 4 | 10% | environment↔substrate coupling: transmission-as-action, exogenous timing, actuation, closed-loop interception |
| **Method** | 6 | 14% | causal machinery (`do(τ)`, `do(χ)`, well-posedness, mediation) + reproducibility/process rules |

Two things this exposes that the taxonomy argument would not have:

1. **Learner is 48% of all claims.** Not a defect — the project *is* about a
   learning rule — but it means Learner cannot be flat. It needs the tier-2
   track layer immediately, which it already has (Tracks 1, 3c, 3d, 3e).
2. **Method's 6 claims are all C-series instrumentation. Zero are process
   claims.** The method lessons with the highest reuse value — "n=3 is adequate
   for treatments but misleading for controls", "prefer bounded proportions to
   means of unbounded errors", "derived metrics fail silently", the `*_TAG`
   convention — currently exist only as prose in a comms-log entry that nothing
   points at. They have evidence (the retracted attention-precision claim; the
   blind `D−x` τ metric that reversed a headline read) and no home. **This is the
   single clearest gap the test found.**

### Cross-pillar edges are the payoff, not the failure

7 of 42 claims (17%) bear on a pillar other than their home. The load-bearing
example: **Result 8** (avoidance is not sign-symmetric — refractoriness is a
prefix of time after firing, so transmission is monotone non-decreasing, so
selective avoidance is impossible) is a *Substrate* claim whose consequence
**refutes a whole class of Learner tactics.** A single tree buries that inside
whichever section discovered it. 17% is small enough to author by hand and large
enough to be worth the machinery.

---

## Per-claim filing (auditable)

`S` = Substrate, `L` = Learner, `I` = Interface, `M` = Method.
`†` = tie-break rule applied. `⇢` = cross-pillar edge.

| # | Source | Claim | 6-cut | 4-cut |
|---:|---|---|---|---|
| 1 | c1 | certificate == intervention == ground truth, on every graph | Instr | M |
| 2 | c2 | one macro target, many behaviours | Instr/Sub | S |
| 3 | c5 | chirality is fat-handed *only when read at a fixed locus* | Instr | M |
| 4 | c6A | the generative handle is well-posed for every reader | Instr | M |
| 5 | c6B | necessity: the persistent core is needed for switching, not routing | Sub/Instr | S |
| 6 | c7A | the causal role is (handle, outcome)-relative | Instr | M |
| 7 | c7B | mediation: chirality screens its own generator | Instr | M |
| 8 | e1 | the predicted A-vs-B dissociation holds | Rule | L |
| 9 | e2.1 | memory duration is τ-controlled | Sub | S |
| 10 | e2.2 | Line B holds memory, Line A cannot | Rule | L |
| 11 | e3.1 | latency tracks gate τ | Sub | S |
| 12 | e3.2 | double dissociation (single lines), and A+B interference | Rule | L |
| 13 | e5.1 | options gate routing; switch cost consolidates | Rule/Anat | L |
| 14 | e5.2 | the slow loop is necessary for switching, not for routing | Sub/Rule | L † |
| 15 | e6.1 | three distinct questions, all answered well above chance | Rule | L |
| 16 | e6.2 | genuinely distinct questions (near-orthogonal readouts) | Rule | L |
| 17 | e6.3 | the function is in the probe, not a module | Rule/Proc | M † |
| 18 | e7dr.A | recovers chirality *and* survives core meander | Rule | L |
| 19 | e7dr.B | the learned readout drives routing | Rule | L |
| 20 | e7.1 | nucleation and persistence | Sub | S |
| 21 | e7.2 | rotation direction is readable | Sub | S |
| 22 | e8.1 | prediction tracks predictability; anticipation sharpens | Rule | L |
| 23 | e8.2 | the history window is a `do(τ)` property | Sub/Instr | S ⇢M |
| 24 | e8.3 | surprise is a global scalar, not a per-feature error field | Sub | S |
| 25 | e8.4 | nested two-timescale prediction: slow context gates fast readout | Sub | S |
| 26 | e8.5 | an order-preserving reservoir recovers deep positional history | Sub | S |
| 27 | e8.6 | conditional long-range prediction needs positional memory AND conjunction | Sub/Anat | L † |
| 28 | e9.1 | the conjunction basis self-organises (reward-free) | Rule/Anat | L |
| 29 | e9.2 | reward routing on the emergent basis works; frozen control fails | Rule | L |
| 30 | lat.1 | the rule works on a lattice; the old self-referential rule still fails | Rule | L |
| 31 | lat.2 | exogenous timing does *not* penetrate a recurrent medium | Sub/Iface | I † ⇢S |
| 32 | lat.3 | an attention chain carries a timing reference: a clock, not a filter | Rule/Anat | L |
| 33 | lat.4 | reward as a fourth edge: τ encodes a stimulus–reward interval | Rule/Anat | L |
| 34 | lat.5 | a value chain removes the last hand-set constant | Rule | L |
| 35 | lat.6 | layers instead of edges: a synchronous burst, two hard requirements | Anat | L |
| 36 | lat.7 | the action is transmission, and bootstrapping is free | Iface | I |
| 37 | lat.8 | avoidance is not sign-symmetric with approach; the reason is a theorem | Sub | S ⇢L |
| 38 | lat.9 | coherent rhythmic drive enables autonomous Q-switched actuation | Sub/Iface | I † ⇢S |
| 39 | lat.10 | attention gating resolves channel crosstalk under conflicting reward times | Rule/Anat | L |
| 40 | lat.11 | closed-loop timed interception bootstraps via graded baseline transmission | Iface | I |
| 41 | lat.N1 | plastic cell identity fails structurally | Anat | L |
| 42 | lat.N2 | tonic drive has no window | Sub | S |

---

## Migration cost, measured

The 42 claims above are the ones already written in claim form. The rest of the
corpus is not:

- **14 docs** use the `## Result N —` convention → migrate mechanically.
- **~15 docs** (`stats_sweeps_results.md`, `continual_learning_results.md`,
  `timescale_hierarchy_results.md`, `closed_loop_plasticity_results.md`,
  `e0_topologies.md`, the topology/capacity docs …) are structured by *phase*
  (`## P2`, `## P3`, `## P4`) with claims embedded in prose. Roughly **18
  additional claim clusters** need extraction by hand. **This is the bulk of the
  work.**
- **3 docs** (`process.md`, `core_review.md`, `extensions_review.md`) are not
  claim-bearing; they become Method concepts, not claims.
- `next_steps.md` (795 lines) currently does four tiers at once and holds two
  parallel naming schemes for the same tier (`Track 1–6` plus a later
  `Proposal Track A/B/C (Options 1, 2, 4)`). Splitting it is the migration's
  first step, not an afterthought.

---

## Verdict

The **four-pillar cut** (Substrate / Learner / Interface / Method) survives the
test the six-pillar cut fails. Rule and Anatomy are one pillar in this project;
Instrumentation and Process are one pillar.

Adopt only if the migration also *retires* surfaces rather than adding one.
Today state lives in five overlapping places (`next_steps.md`,
`lattice_results.md`, `lattice_timescale_notes.md`, the two review docs,
`.agents/comms_log.md`) plus `process.md` describing how two of them relate. A
sixth structure beside those makes the problem worse.

## Open items

- **The Method pillar has no claims.** Writing its ~5 process claims — with the two
  that carry real evidence (`n=3`-for-controls, the blind derived metric) wired as
  reproduce contracts against the scripts they came from — is the smallest step that
  both tests the schema and fills the gap this test found.
- **Container format.** [OKF v0.1](https://github.com/GoogleCloudPlatform/knowledge-catalog/blob/main/okf/SPEC.md)
  fits the pillar tree without adaptation: a bundle *is* a directory tree of markdown
  with YAML frontmatter, `type` is the only required field, extension keys are
  explicitly permitted, `tags` carries the cross-cutting tensions, `stale_after` makes
  per-tier cadence a checkable date, and `generated`/`verified` reproduce this repo's
  review/plan decoupling. `type: Attested Computation` (§10 — `runtime`, `parameters`,
  `executor.receipt`, `attester`) is the reproduce contract, aimed squarely at the two
  failures that cost this arc time: a clobbered artifact and a silently-wrong derived
  metric.
  **Gap:** OKF links are deliberately *untyped* (§6.1), so the `supports | strains |
  refutes` up-edge — the piece that makes tiers evolve rather than merely nest — must be
  an extension field with a local validator. `status: draft|stable|deprecated` is
  document lifecycle, not epistemic status; `refuted`/`retracted` needs its own field.
- **Git can carry the ledger, not the graph.** Branch-per-tactic is already the
  convention (`claude/3a-p2-sweeps`, `claude/3e2-cfc`); PR-as-evidence-review and
  merge-as-status-transition formalise existing practice. But up-edges are many-to-many
  over *statements* while merges are tree-shaped over *snapshots*, and threshold
  conditions need something reading frontmatter — that is CI, not git. Note also that
  slow tiers want few reviewed commits on `main`, not long-lived branches: the cadence
  mapping inverts.
