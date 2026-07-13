# Synthesis — one honest narrative

*How the learning experiments (E0–E9) and the causal/constitution experiments
(C0–C7), the outward-facing artifacts (spiral→data predictions, the causal
testbed), and the two review passes form **one argument**. This is the connective
tissue and the honest ledger; each claim links to its results doc. Positioning:
an **illustrative computational study** of an inside-out, reward-shaped excitable
substrate, plus **one reusable external artifact** (a synthetic-SCM causal
testbed). It illustrates mechanisms and yields falsifiable predictions; it does
not claim to prove them of cortex.*

## The two arcs

- **E-series (learning).** On a Greenberg–Hastings network, a strict scalar
  reward (Sutton) selects and stabilises intrinsic dynamical patterns (Buzsáki
  inside-out), through two plasticity lines: **Line A** learns conduction weights
  (spatial routing), **Line B** learns local timescales `τ` (temporal structure).
- **C-series (causality).** The same substrate, where the wave `W = f(S)` is an
  explicit deterministic coarse-graining of spikes, is a synthetic-SCM testbed for
  the spike–wave causal question (Jalaldoust & Zabeh, arXiv:2511.06602).

They look separate. They are the same object viewed twice.

## The one idea that ties them: `θ` is both the learned variable and the causal handle

The C-series' central positive result (C3) is that the well-posed causal handle
is not the wave aggregate `W` but the **generating parameters `θ`** — timescales
and couplings — because `do(θ)` is modular and unique while `do(W)` is
fat-handed (C2: a 33 σ realization band vs 0.014 σ). The E-series' two plasticity
lines learn **exactly those parameters**: Line A = couplings, Line B =
timescales. So:

> The variable the learner adapts (`θ`) is precisely the variable on which
> intervention is causally well-defined. Learning substrate and causal handle
> coincide.

This is not engineered in; it falls out of both analyses independently. It
reframes the paper's programme: instead of "is the wave causal?" (ill-posed under
constitution), ask "can we intervene on what organises the wave?" — and that
organiser is what plasticity already targets.

## A shared dissociation, learned in E and causally grounded in C

| | spatial / "spike-like" | temporal / "wave-like" |
|---|---|---|
| **coding (C0)** | labeled-line: wave uninformative | collective: wave informative |
| **learned by (E1–E3)** | Line A (weights) → identity | Line B (timescales) → timing / memory |
| **causal handle (C3–C4)** | `do(g_route)` → identity | `do(τ)` → timing |
| **outcome-relativity (C4)** | wave epiphenomenal for identity | wave causal for timing |

Read across any row: the spatial/identity channel and the temporal/timing channel
are distinct all the way down — in how behaviour codes them (C0), which plasticity
line learns them (E1 identity, E2 memory, E3 both), and which `do(θ)` handle
controls them (C4). E3's behavioural double dissociation and C4's diagonal causal
matrix are the **same fact** at two levels of description.

## What C taught E

1. **Why learning `θ` (not states) is the right design.** C2/C3 show state/aggregate
   intervention is fat-handed; parameter-level is well-posed. Making weights and
   timescales plastic — rather than clamping activity — is that well-posedness in
   learning form.
2. **A sharper reading of the E3 A+B interference.** C4's outcome-relativity says
   identity and timing are *different outcomes with different causal handles*, so a
   single scalar reward conflates two credit-assignment problems; the fix is to
   treat them as the distinct causal channels C4 shows them to be.
3. **The wave is a reader's variable, not a controller's.** Read order parameters
   as features/critic signals; drive plasticity through `θ`, never by setting the
   wave.

## The combined thesis

A GH network, driven inside-out and shaped by a strict reward, learns cognitive
functions by adapting its **generative parameters** `θ` — couplings for
identity/routing, timescales for timing/memory. Those same parameters are the only
context-free, well-posed causal handles on the system's behaviour; the collective
**waves** they produce are informative, emergent descriptive variables whose causal
status is real but contingent on observation model, graph, constitution, and
outcome. Spikes and waves are not rivals for causal primacy — they are two readouts
of one parameterised dynamics, and the parameters are where both learning and
causation live.

## From two lines to a full repertoire (E4–E6)

The two-line story (routing vs timing) is the engine; E4–E6 show it composing into
recognisable faculties on the *one* substrate, without adding machinery.

- **[E4](e4_results.md) — attention**: the substrate's native refractory
  annihilation is a winner-take-all (colliding waves suppress each other), so
  *selection* is free once waves exist; a top-down bias moves the annihilation locus.
- **[E5](e5_results.md) — executive control**: E2's persistent loop repurposed as an
  **option** — a slow *timescale* structure (`τ < L`) gating *spatial* routing across
  a block. Ablating persistence kills switching but not single-rule routing.
- **[E6](e6_results.md) — emergent readout**: freeze the shaped substrate and the
  "faculties" are *questions* — memory, attention, executive control are three GVF
  cumulants read off one phase stream by one linear probe, near-orthogonally, no
  dedicated sub-network. The reader/controller split made concrete.

## Onto the native medium and forward in time (E7, C5–C7, E8)

- **[E7](e7_results.md) + the spiral causal split ([C5](c5_results.md)–[C6](c6_results.md)–[C7](c7_results.md)).**
  E7 lifts the E5 option from a 1-D ring to a genuine 2-D spiral core whose
  **rotation direction** is the rule. C5–C7 answer the field's spike-vs-wave question
  for a real spiral: rotation direction is a causal **mediator** (`θ_seed → χ → B`;
  ablating the core kills switching) — *not epiphenomenal* — but contingent on the
  **reader** (C5: fat-handed at a fixed locus 6.2σ; well-posed read topologically
  1.0σ) and the **outcome** (C7: causal for the rule, inert for content). The clean,
  reader-robust handle is the generative **nucleation** `do(θ_χ)` (C6: 0σ). The
  scalar C2→C3→C4 arc, re-run on a topological wave variable — and *refining* C2:
  well-posedness is a property of the *(variable, reader)* pair.
- **[E8](e8_results.md) — predictive dynamics.** The substrate makes time-forward
  predictions of tone sequences as **forward dynamics + a `do(τ)` history window**
  (0→7 tones as τ grows), with a **global** surprise. So "the causal variable is the
  generative timescale" reappears E-side as "the predictive memory depth is `τ`."

## Closing the afforded → learned gap (E9, 1c) — tension 1

The programme's most-cited honest limitation (audit tension 1) was that several
headline capabilities were *afforded* by hand-built readouts/features rather than
*learned on* the substrate. Two follow-ups convert the two biggest instances:

- **[E9](e9_results.md) — emergent conjunction cells.** E5 (and E8.5/8.7) *wired*
  the `(stimulus × rule)` conjunction basis. E9 removes the wiring: unbiased fan-in
  + a reward-free **competitive Hebbian** rule (k-WTA lateral inhibition + DeSieno
  conscience, keeping E5's AND-gate) grows the basis from the input statistics —
  conjunction selectivity **0.00 → 1.00**, all four classes tiled, in ~3
  exposures/combo. Reward routing on the **emergent** basis reaches **0.84** (wired
  0.92); the no-self-organisation **frozen** control fails at chance (0.25), so the
  self-organisation is load-bearing. The representation is now *learned*, not wired
  — answering E1's oldest deferred item.
- **[1c](e8_hardening_results.md) — online, reward-free prediction.** E8's offline
  ridge becomes a bank of online **GVF demons** (E6's TD rule): a single incremental
  reward-free pass matches the ridge on every sequence, and γ>0 demons give
  multi-step predictive knowledge (r ≈ 0.97/0.85) the ridge doesn't. Prediction is
  now learned online and intrinsically.

Honest remainder: both are still *linear readouts of a fixed substrate* — E9's
competition is an imposed k-WTA, and 1c makes the *learning* online but not the
*dynamics* plastic. The direction-selective E7 readout (1b, parked) and fully
plastic dynamics are the frontier this points at but does not cross.

## The predictive-coding foil, and a correction it forced (2a) — tension 3

E8 *asserted* "prediction without predictive coding". [2a](e8_hardening_results.md)
built a minimal untied Rao–Ballard PC model on the same tones to make the contrast
empirical — and, in the spirit of the audits, it **corrected an overclaim** rather
than confirming the framing:

- The clean PC signature reproduces (top-down lesion → prediction 0.37→0.12 chance,
  representation 0.95 spared).
- But **dissociability and error-locality — the two observables the E8 doc leaned
  on — do NOT distinguish the accounts**: E8 is *also* dissociable (its two readouts
  tap different parts of the medium), and both models' errors localise on a deviant.
  The `e8_results.md` "not dissociable" line was wrong and is corrected.
- The real discriminators are **architectural**, now stated as empirical/testable:
  intrinsic per-channel *error units* + a dedicated *generative pathway* (PC) vs a
  global scalar surprise + a passive-medium readout + a `do(τ)` window (E8).

This self-correction is part of the argument's credibility, not an aside: the same
scrutiny that caught E3's bimodal-mean composition headline and the `perturb_tau`
RNG bug also caught this framing, and fixed it in place.

## Outward-facing artifacts (2b, 4b) — tension 3

The toy substrates *illustrate* mechanisms; two artifacts carry the work outward,
where the claims can meet data or other methods.

- **[Spiral → data predictions](spiral_predictions.md) (2b).** Six falsifiable
  predictions for cortical spiral-wave recordings (Gong 2023; Ye/Steinmetz 2026),
  each with an observable, a discriminator against an epiphenomenal account, and a
  falsifier — e.g. P2: disrupting persistent cores should selectively impair
  cognitive *flexibility* while sparing fixed mappings; P4 (with a bridge figure from
  C5): a fixed-ROI rule decoder should collapse as the spiral core drifts while a
  topology-aware one holds. Model → testable neural hypothesis.
- **[The causal testbed](causal_testbed.md) (4b).** `ghca_testbed` packages C0–C7 as
  a reusable synthetic-SCM benchmark: a substrate with an explicit constitution, the
  three `do`-operators, and — unlike in vivo — a **known** graph, so every task ships
  a ground-truth answer. The C1 core runs end-to-end (front-door scored **causal
  0.257** by the correct interventional-vs-observational Definition-1 test, which a
  naïve `do(1)`-vs-`do(0)` test would miss); a `score(method_fn)` harness lets others
  test their epiphenomenality / causal-discovery methods against ground truth. The
  externally most-reusable output of the programme.

## State of the evidence — the honest ledger

The three tensions the audits raised, as they stand now:

1. **Afforded vs learned** — *substantially retired.* E9 (conjunctions) and 1c
   (prediction) convert the two biggest afforded capabilities to learned. Open: the
   E7 direction readout (1b, parked — C5's readout-locality makes it genuinely
   risky), and the deeper step of making the *dynamics* (not just the readout)
   plastic.
2. **Narrow evidence** — *still open.* Small n (3–5 seeds), a single substrate
   family, hand-chosen operating points. Results are robust *at those points* (tight,
   unimodal per-seed spreads — the specific E3 failure mode was checked for and is
   absent in the extensions), but not shown *across the regime*. The cheap fix
   (3a: seeds/CIs/sweeps) and generality checks (3b: other topologies) are not done.
3. **Illustrative, not proof** — *addressed by framing + artifacts.* The study is
   positioned as illustrative; 2b turns it into falsifiable data predictions and 4b
   into a reusable ground-truth benchmark. The toy substrates still illustrate rather
   than prove the cognitive-function claims.

**Corrections made (the culture, not just the results).** E3's A+B composition
headline was walked back from "quadrupling / partly resolved" to the honest
"0.11 → 0.48 ≈ chance → 0.56, joint on 1/5 seeds" (direction supported, not
magnitude); a reproducibility bug (`perturb_tau` using the global RNG) was fixed and
the affected numbers re-run; and 2a corrected E8's "not dissociable" claim. Two
audit passes (an independent hallucination/overreach review of E0–E6/C0–C4 and a
self-audit of the extensions, [`extensions_review.md`](extensions_review.md))
underwrite the reported numbers as reproducible.

## One-paragraph takeaway

Across spikes, waves, spirals, and predictions, one substrate — a reward-shaped
excitable medium — produces a repertoire of cognitive functions, and in every case
the context-free place where both *learning* and *causation* live is the generative
parameters `θ` (couplings, timescales), not the collective wave they emit. The
waves are real, informative, and sometimes causal, but as *readouts* whose status is
contingent on reader and outcome. Where the first pass *afforded* capabilities with
hand-built readouts, the later work *learns* them (E9, 1c); where it *asserted* a
contrast with predictive coding, the later work makes it *empirical* — and, honestly,
corrects it (2a). The programme's most durable outputs are the falsifiable spiral
predictions (2b) and the ground-truth causal testbed (4b).

## Map of the work

- **Plan / substrate:** [`learning_experiments.md`](learning_experiments.md) ·
  [`causal_experiments.md`](causal_experiments.md) ·
  [`causal_spiral_experiments.md`](causal_spiral_experiments.md) ·
  [`predictive_dynamics_experiments.md`](predictive_dynamics_experiments.md)
- **E-series:** [`e0`](e0_results.md) · [`e1`](e1_results.md) · [`e2`](e2_results.md) ·
  [`e3`](e3_results.md) · [`e4`](e4_results.md) · [`e5`](e5_results.md) ·
  [`e6`](e6_results.md) · [`e7`](e7_results.md) · [`e8`](e8_results.md) ·
  [`e8 hardening (1c, 2a)`](e8_hardening_results.md) · [`e9`](e9_results.md)
- **C-series:** [`c0`](c0_results.md) · [`c1`](c1_results.md) · [`c2`](c2_results.md) ·
  [`c3`](c3_results.md) · [`c4`](c4_results.md) · [`c5`](c5_results.md) ·
  [`c6`](c6_results.md) · [`c7`](c7_results.md)
- **Artifacts / outward:** [`spiral_predictions.md`](spiral_predictions.md) (2b) ·
  [`causal_testbed.md`](causal_testbed.md) + `ghca_testbed/` (4b)
- **Reviews / roadmap:** [`extensions_review.md`](extensions_review.md) ·
  [`next_steps.md`](next_steps.md) · 4a parked on branch
  `claude/e10-timescale-hierarchy` (`docs/e10_notes.md` — the τ-hierarchy negative
  result and the new-rule it needs)
