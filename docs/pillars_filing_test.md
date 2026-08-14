# Pillars filing test — does a vertical cut of this repo actually partition its claims?

**Status:** analysis, not a decision. Runs 1–2 on 2026-08-11; **Run 3 and two
corrections added 2026-08-14** after the persistence/topology work landed on `main`.

> ## ⚠ Run 2's verdict does not survive — read the corrections first
>
> Runs 1–2 below concluded that a **four-pillar** cut (Substrate / Learner /
> Interface / Method) passed at 12% ambiguity. **It fails at 23.5%** once the
> claims added by PRs #78 and #79 are filed. Two independent causes, plus two
> arithmetic corrections to the corpus itself — see
> [Corrections](#corrections-to-runs-1-2) and [Run 3](#run-3--the-four-pillar-cut-meets-proof-shaped-work).
> The text of Runs 1–2 is left unedited as the historical record.

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

## Corrections to Runs 1–2

*Added 2026-08-14.*

### C1 — the extraction command was wrong, and so was `n`

The Method section's regex requires a space after the keyword. That silently
drops claims whose heading is a bare `## Result` with the statement on the next
line:

```bash
# what was run (trailing space):                     46 on main
# without the trailing space:                        56 on main
```

Adjudicating the 10-heading difference by hand: **7 are section containers**
(`## Results`, `## Findings about the learning rules (honest caveats)`) and
should not count; **3 are genuine claims** that were missed —
`c3_results.md`, `topology_winding_capacity.md`, and `e0_topologies.md`'s
`## Results — both E0 headlines generalise`.

**Corrected n for Runs 1–2: 42 → 45.** All three file cleanly (Method, Substrate,
Substrate), so Run 1 becomes 13/45 = **29%** and Run 2 becomes 5/45 = **11%**.
Neither verdict flips on this correction alone.

This is the repo's own "derived metrics fail silently" lesson landing on the test
that was written to catalogue it. A trailing space moved a headline by 24%. It is
logged below as a Method claim in its own right.

### C2 — `lattice_timescale_notes.md` holds 4 claim headings that must not be counted

They are *earlier, superseded drafts* of claims that appear in final form in
`lattice_results.md` (e.g. "Negative 2 — the rule does not re-tune" and
"Negative 2, resolved: the rule needs a privileged afferent channel"). The
notebook is append-only by design.

Counting them would double-count; ignoring them loses the correction history.
**The corpus therefore requires a `supersedes` / `superseded_by` edge and a
`deprecated` status before any migration** — this is not an optional refinement.
Run 1 excluded them correctly but by accident, not by rule.

---

## Run 3 — the four-pillar cut meets proof-shaped work

PRs #78 and #79 added 2 docs and 5 experiments (`topology_cycle_packing_exact.md`,
`persistent_set_structure.md`). Neither uses the `## Result N —` convention, so
**the extractor returns zero on both** and this run is a hand extraction: 26
claims at bolded-assertion granularity (16 at section granularity; the ambiguity
percentage moves <4 points, so the verdict is not a granularity artifact).

**Result: FAILS.**

| Pillar | Runs 1–2 (45) | New (26) | Combined (71) |
|---|---:|---:|---:|
| Substrate | 14 | 11 | 25 |
| Learner | 20 | 0 | 20 |
| Interface | 4 | 0 | 4 |
| Method | 7 | 15 | 22 |
| **Ambiguous** | 5 (11%) | **11 (42%)** | **16 (23.5%)** |

Two structural breakages, and **neither is a counting artifact**:

**1. Substrate ↔ Method inverts.** Run 1 merged Instrumentation into Method on the
reasoning that "the instrument is provenance, not subject." In the persistence
thread **the instrument *is* the subject** — the headline is literally *which
descriptions of the persistent set exist* (no linear separation; a gap-multiset
invariant; an O(1) outer hull that is necessary but not sufficient). All 11
ambiguities sit on this one seam. That is more than twice the size of the
Rule↔Anatomy cluster that condemned Run 1's cut.

**2. Three claims are homeless, and they name the missing pillar.** The
cycle-packing claims are structure that is **fixed and pre-dynamical** — the doc
says outright *"This is structure, not dynamics. No dynamics are run here."*
Run 2 dissolved Anatomy into Learner on the premise that "plasticity and the
structure it grows are not separable"; these are structure plasticity did **not**
grow, so they fall through. Substrate as defined holds "what the medium *does*."
By this test's own threshold — *any homeless ⇒ the cut is missing a pillar* — that
is a fail before the ambiguity count is read.

The reason Run 2 missed it: **the original 42 contained no topology claim at all.**
The topology/capacity docs sat in the test's own "~15 docs not yet in claim form"
bucket, so the merge was never tested against the case that breaks it.

### What Run 3 got right

**Method is no longer empty of process claims: 0 → 6.** Run 2's headline open item
is closed by the new work, every one carrying evidence:

| Process claim | Evidence |
|---|---|
| A one-word pick rule in a *merged, reviewed* experiment multiplied a headline by up to 7.5× | ring τ=3: 6 published vs 45 certified |
| Solver time limits are part of the result, not an implementation detail | ring τ=6: 24 at 60 s, 25 at 120 s |
| Cross-check a reimplementation against the repo's primitive on *every* transition; abort on mismatch | 7,461 transitions, 9 cells, 0 mismatches |
| Impure *class count* is not a refinement metric; *configurations in impure classes* is | (3,4): class count rose 8→28 under a genuine refinement |
| A "trivial cell" cutoff must be `P == 0` exactly | `P ≤ 0.01` hid (2,5) (P=0.0039) and produced a false iff |
| State a shortcut's validity precondition and check it exhaustively before relying on it | no non-zero dead attractor at θ=1, verified at 2×2 |

Plus C1 above (a trailing space in a grep). Two of these are the same failure
genus as the blind `D−x` τ metric: **a number wrong in a direction nobody checked.**

The irony is worth recording: the new work closed the test's most-cited gap and
broke its pass grade in the same event. Method grew *because* the new claims are
largely about how you know.

### Evidence is not one kind, and `n` is not universal

`persistent_set_structure.md` states it plainly: *"No RNG. Every number is a
complete enumeration, so there is nothing to seed and no per-seed spread to
report."* The existing convention encodes n in the heading (`(n=20)`), with no
slot for "n = the whole space."

But exhaustive work still carries an n — **extent**, not sample size (81 →
40,353,607 configurations; 7,461 validated transitions) — and it still carries
uncertainty, of three kinds the current prose conflates: **coverage** (which cells
were enumerated), **certification** (4 of 16 packing cells are certified; the rest
are brackets), and **representation dependence** (facet counts are
Qhull-triangulation-sensitive; vertex counts are not). And "exhaustive" does not
mean unseeded: `persistent_set_3x3.py:192` seeds a 2000-sample cross-validation.

So `evidence` must be a **discriminated union on `kind`** — `seeded_simulation`
carrying `n`/`spread`/`seed_policy`; `exhaustive_enumeration` carrying
`domain`/`coverage`/`verification` with `spread: not_applicable` and an explicit
reason; `exact_solve` carrying `value_type: exact|bound`, `certified_cells`, and
`solver_sensitivity`. OKF v0.1 permits extension keys and requires only `type`,
so this stays inside the chosen container format.

Two vocabulary items are now demonstrated rather than hypothetical:

- **`supersedes` with a typed effect** (`refutes` / `revises` / `replaces_reasoning`)
  — required by C2, and by the cycle-packing correction, which refutes one
  published finding and revises another *without editing the doc that carries them*.
- **`argument_status` as an axis distinct from `epistemic_status`.** The refinement-metric
  claim is a case where the **finding stands and the argument was retracted**.
  `refuted` and `retracted` cannot express that.

---

## Verdict

**Superseded by Run 3.** The four-pillar cut passed on a corpus that happened to
exclude every static-structure claim in the repo; it fails at 23.5% once they are
included, with three homeless claims naming the pillar it deleted.

What survives from Runs 1–2: merging **Rule + Anatomy → Learner** was right *for
grown structure*, and merging **Instrumentation + Process → Method** was right
*for provenance*. Both merges over-reached by exactly one case each.

**Proposed Run 4 — five pillars, with the rule declared in advance:**

| Pillar | Holds |
|---|---|
| **Substrate** | what the medium *does* — dynamics, τ-control, nucleation, the avoidance theorem |
| **Carrier** | what the medium *is* — fixed topology, cycle space, capacity bounds, the persistent set as an object |
| **Learner** | plasticity and the structure it *grows* |
| **Interface** | environment ↔ substrate coupling |
| **Method** | instrumentation, tractability/description results, process |

**Declared before counting** (this is the repo's own house rule from `AGENTS.md`,
not a new invention): *file by the substrate/analysis boundary the doc itself
declares.* Both new docs state theirs verbatim.

**This is a proposal, not a result. It has not been run.** Predicted ~7% ambiguity,
but that prediction is worth little until the ~18 unextracted clusters in
`stats_sweeps_results.md`, `continual_learning_results.md` and the topology docs
are extracted — and those will land in the same Carrier-shaped hole, so they test
the fix directly. Do not adopt a cut on a predicted number; that is the mistake
Run 2 made.

Adopt only if the migration also *retires* surfaces rather than adding one.

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
