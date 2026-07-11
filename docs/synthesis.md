# Synthesis — Tying the E-series and C-series Together

*How the learning experiments (E0–E3) and the causal/constitution experiments
(C0–C4) form one argument. This is the connective tissue; each claim links to
its results doc.*

## The two arcs

- **E-series (learning).** On a Greenberg–Hastings network, a strict scalar
  reward (Sutton) selects and stabilises intrinsic dynamical patterns (Buzsáki
  inside-out), through two plasticity lines: **Line A** learns conduction weights
  (spatial routing), **Line B** learns local timescales `τ` (temporal structure).
- **C-series (causality).** The same substrate, where the wave `W = f(S)` is an
  explicit deterministic coarse-graining of spikes, is used as a synthetic-SCM
  testbed for the spike–wave causal question (arXiv:2511.06602).

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

This is not a coincidence engineered in; it falls out of both analyses
independently. It also reframes the paper's programme: instead of asking "is the
wave causal?" (ill-posed under constitution), ask "can we intervene on what
organises the wave?" — and that organiser is what plasticity already targets.

## A shared dissociation, learned in E and causally grounded in C

| | spatial / "spike-like" | temporal / "wave-like" |
|---|---|---|
| **coding (C0)** | labeled-line: wave uninformative | collective: wave informative |
| **learned by (E1–E3)** | Line A (weights) → identity | Line B (timescales) → timing / memory |
| **causal handle (C3–C4)** | `do(g_route)` → identity | `do(τ)` → timing |
| **outcome-relativity (C4)** | wave epiphenomenal for identity | wave causal for timing |

Read across any row: the spatial/identity channel and the temporal/timing
channel are distinct all the way down — in how behaviour codes them (C0), in
which plasticity line learns them (E1 identity, E2 memory, E3 both), and in which
`do(θ)` handle controls them (C4). E3's behavioural double dissociation (Line A
→ identity, Line B → timing) and C4's causal matrix (diagonal) are the **same
fact** at two levels of description.

## What C taught E (folded back into the E-series)

1. **Why learning `θ` (not states) is the right design.** C2/C3 show state-level
   / aggregate-level intervention is fat-handed; parameter-level is well-posed.
   The E-series' choice to make weights and timescales plastic — rather than, say,
   clamping activity patterns — is the same well-posedness, in learning form.
2. **A sharper reading of the E3 A+B interference.** E3 found that a single shared
   reward makes Line A and Line B interfere. C4's outcome-relativity says why this
   is expected and how to fix it: identity and timing are *different outcomes with
   different causal handles*, so a single scalar error conflates two credit-
   assignment problems. The principled fix (factored credit / separate
   neuromodulators / a curriculum) is to treat them as the distinct causal
   channels C4 shows them to be — not to hope one reward disentangles them.
3. **The wave is a reader's variable, not a controller's.** C0/C4 show the wave
   carries information (for collective codes) and is a legitimate *descriptive*
   variable, but C2 shows it is not a well-posed *control* variable. For the
   learner this means: read order parameters as features/critic signals (as the
   E-series critic does), but drive plasticity through `θ`, never by trying to set
   the wave.

## The combined thesis

A GH network, driven inside-out and shaped by a strict reward, learns cognitive
functions by adapting its **generative parameters** `θ` — spatial couplings for
identity/routing, local timescales for timing/memory. Those same parameters are
the only context-free, well-posed causal handles on the system's behaviour; the
collective **waves** they produce are informative, emergent descriptive
variables whose causal status is real but contingent on observation model, graph,
constitution, and outcome. Spikes and waves are not rivals for causal primacy —
they are two readouts of one parameterised dynamics, and the parameters are where
both learning and causation actually live.

## From two lines to a full repertoire (E4–E6)

The two-line story (routing vs timing) is the engine; E4–E6 show it composing into
recognisable faculties on the *one* substrate, without adding machinery.

- **[E4](e4_results.md) — attention** uses the substrate's native refractory
  annihilation as a winner-take-all: colliding waves suppress each other, so
  *selection* is free once waves exist (zero inhibitory nodes). A top-down bias
  moves the annihilation locus — competition resolved spatially.
- **[E5](e5_results.md) — executive control** repurposes E2's persistent loop as an
  **option**: a slow *timescale* structure (`τ < L`) that gates *spatial* routing
  across a block. This is the two orthogonal `do(θ)` handles (timescales, couplings)
  composed **hierarchically** rather than fought over one scalar — the constructive
  counterpart to the E3 A+B interference. Ablating the loop's persistence kills
  switching but not single-rule routing, localising the loop's causal role.
- **[E6](e6_results.md) — emergent readout** is the punchline: freeze the shaped
  substrate and the "faculties" are *questions*. Memory, attention and executive
  control are three GVF cumulants read off one phase stream by one identical linear
  probe, near-orthogonally and with no dedicated sub-network (a generic probe
  matches an own-region oracle). This is the reader/controller split made concrete —
  the demons **read** the collective state and never drive it, exactly as C2/C4 say
  the wave should be used.

Read together: the spatial/temporal `do(θ)` dissociation (E1–E3, grounded in C3–C4)
is the mechanism; attention, options and the emergent-category readout (E4–E6) are
what that mechanism *looks like* when asked cognitive questions.

## Onto the native medium and forward in time (E7, C5–C7, E8)

Two extensions, prompted by the rotating-wave and auditory-prediction neuroscience,
push the same thesis further — and, satisfyingly, land on the same conclusion.

- **[E7](e7_results.md) + the spiral causal split ([C5](c5_results.md)–[C6](c6_results.md)–[C7](c7_results.md)).**
  E7 lifts the E5 "option" from a 1-D ring to a genuine 2-D spiral core whose
  **rotation direction** is the rule. C5–C7 then point the C-series machinery at it and
  answer the field's spike-vs-wave question for a real spiral: rotation direction is a
  causal **mediator** (`θ_seed → χ → B`; ablating the core kills switching) — *not
  epiphenomenal* — but its status is **contingent on the reader** (C5: fat-handed at a
  fixed locus, 6.2σ; well-posed read topologically, 1.0σ) and **on the outcome** (C7:
  causal for the rule, inert for content). The clean, reader-robust handle is the
  generative **nucleation** `do(θ_χ)` (C6: 0σ for every reader) — "drive the
  parameters," now for a genuine spiral. This is the scalar C2→C3→C4 arc re-run on a
  topological wave variable, and it *refines* C2: well-posedness is a property of the
  *(variable, reader)* pair, and a topological invariant is better-posed than a scalar
  aggregate.
- **[E8](e8_results.md) — predictive dynamics.** The substrate makes time-forward
  predictions of tone sequences as **forward dynamics + a `do(τ)` history window**
  (window 0→7 tones as τ grows), with a **global** surprise and no per-feature error
  hierarchy — a *predictive* system that is not a *predictive-coding* one. Its slow
  extension (E8.5) reuses the E5 option: a slow persistent context, via a fast×slow
  conjunction, gates prediction across a regime block. So the C-series' "the causal
  variable is the generative timescale" reappears on the E-side as "the predictive
  memory depth is the generative timescale `τ`."

The two extensions rhyme: C5–C7 say *read the wave topologically, control it through
`θ`*; E8 *reads the medium's integrated state and controls its memory through `τ`*.
Across spikes, waves, spirals, and predictions, the parameters `θ` (timescales,
couplings) remain the one context-free place where learning and causation live.

## Map of the work

- Plan / substrate: [`learning_experiments.md`](learning_experiments.md),
  [`causal_experiments.md`](causal_experiments.md)
- E-series results: [`e0`](e0_results.md) · [`e1`](e1_results.md) ·
  [`e2`](e2_results.md) · [`e3`](e3_results.md) · [`e4`](e4_results.md) ·
  [`e5`](e5_results.md) · [`e6`](e6_results.md) · [`e7`](e7_results.md) ·
  [`e8`](e8_results.md)
- C-series results: [`c0`](c0_results.md) · [`c1`](c1_results.md) ·
  [`c2`](c2_results.md) · [`c3`](c3_results.md) · [`c4`](c4_results.md) ·
  [`c5`](c5_results.md) · [`c6`](c6_results.md) · [`c7`](c7_results.md)
- Extension specs: [`causal_spiral_experiments.md`](causal_spiral_experiments.md) ·
  [`predictive_dynamics_experiments.md`](predictive_dynamics_experiments.md)
