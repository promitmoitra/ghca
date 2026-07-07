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

## Map of the work

- Plan / substrate: [`learning_experiments.md`](learning_experiments.md),
  [`causal_experiments.md`](causal_experiments.md)
- E-series results: [`e0`](e0_results.md) · [`e1`](e1_results.md) ·
  [`e2`](e2_results.md) · [`e3`](e3_results.md)
- C-series results: [`c0`](c0_results.md) · [`c1`](c1_results.md) ·
  [`c2`](c2_results.md) · [`c3`](c3_results.md) · [`c4`](c4_results.md)
