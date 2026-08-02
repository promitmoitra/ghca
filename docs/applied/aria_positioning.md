# Applied direction — positioning against ARIA's olfaction programme

> **Status: internal strategy draft, not published.** Not in the MkDocs nav, not on
> the public site. Companion to [`aria_olfaction_feedback.md`](aria_olfaction_feedback.md).

## The one-line pitch

**Self-calibrating temporal receptive fields for chemical sensing:** a local,
gradient-free rule that lets a sensor front-end *learn and continuously re-learn* the
timescales of its own input, so that discriminative temporal structure is preserved
as sensor kinetics drift — without labelled recalibration.

Derived from the GHCA timescale arc (3d–3e.3), extracted out of the Greenberg–Hastings
substrate and into a standard spiking/LIF front-end.

## Where it fits the programme (honest mapping)

| Workstream | Fit | Reality check |
|---|---|---|
| **A** — open hardware + dataset | **Indirect but high-leverage.** Our §1 feedback (preserve raw time-resolved transients) is an ORO *standard* contribution, not a funded build. | We are not a hardware or data-collection team. The exception is `[[ Dognosis canine/clinical data — if real, this is a direct A contribution and the strongest card ]]`. |
| **B** — learning representations | **Partial and orthogonal.** B as written means *chemical* representation (molecular embeddings, odour space dimensionality). We offer *temporal* representation. | Do not position as competing with B. Position as the missing second axis. |
| **C** — novel sensing capability | **Best fit.** C explicitly names "neuromorphic olfactory hardware" and "AI-driven drift compensation," and Table 1 marks drift red for every deployable sensor class. | C is hardware-flavoured; we are software over existing arrays. Needs a hardware partner to be a credible C bid. |

**Realistic posture:** a component contributor inside a consortium, not a prime. The
programme funds systems; a mechanism needs a host.

## The preliminary result to build (simulation only)

One experiment, three nested claims, with an explicit kill condition at each step.
This is the evidence a proposal needs, and it is all public-data + CPU work.

**H1 — time matters at all (precondition).** Classifiers using the *raw response
transient* beat classifiers using the standard reduced steady-state feature set, on
the same samples and splits.
*Kill: if reduced features match raw traces on the chosen dataset, the temporal thesis
is dead for that dataset — report it and switch datasets or stop.*

**H2 — timescale diversity helps.** A front-end with a *tiled set* of temporal
receptive fields (heterogeneous τ) beats a single fixed timescale and beats random-τ
at matched parameter count.
*Kill: if random heterogeneous τ matches tiled τ, the "learned tiling" claim adds
nothing over known heterogeneity results.*

**H3 — self-calibration under longitudinal drift (the differentiator).** Under the
recalibration-free protocol (train on early months, test months later, no target
labels), the local τ rule tracks changing sensor kinetics and degrades more slowly
than (a) fixed τ, (b) BPTT-learned τ frozen at deployment.
*This is the claim nobody else can make, because BPTT-learned time constants are
frozen once training ends. It is also the riskiest.*

**Baselines, non-negotiable:** fixed homogeneous τ; random heterogeneous τ;
BPTT/surrogate-gradient-learned τ (PLIF-style); and at least one classical
drift-compensation method from the e-nose literature (most use target-domain data —
note that asymmetry explicitly rather than beating them unfairly).

### Dataset constraint discovered while planning — and why it validates the feedback

The field's canonical drift benchmark (**Vergara/UCSD, 36 months**) is
**feature-reduced**: 8 hand-engineered scalars per sensor, raw traces discarded. So
*the standard drift benchmark cannot be used to test a temporal hypothesis at all.*

That is a real constraint on this plan — and it is direct evidence for §1 of the
feedback. Worth stating in both documents: the field's most-used longitudinal resource
already foreclosed this question by a preprocessing decision, which is exactly the
mistake ORO should avoid.

**Candidate datasets that may retain raw traces (verify before committing):**
- One-year MOX gas sensor array dataset (*Scientific Data*, 2025) — longitudinal, the
  best candidate for H3 if traces are raw.
- UCI *gas sensor array under dynamic gas mixtures* (Fonollosa et al.) — high-rate raw
  time series; good for H1/H2, but not longitudinal, so it cannot test H3.
- **SmellNet** — cited in ARIA's own bibliography (ref 32); using a dataset from their
  reference list is worth something in a proposal.

Likely split: **H1/H2 on a high-rate dynamic dataset; H3 on the longitudinal one.**
If no public dataset has raw longitudinal traces, H3 becomes a *synthetic* drift study
(inject slow kinetic drift into real transients) — weaker, and it must be labelled as
such.

## Sequencing

1. **Feedback response to ARIA** — days. Time-sensitive: thesis is v1.0, programme
   reportedly launching summer 2026. Do this first; it also forces the positioning to
   be sharp and puts us on the PD's radar before the call.
2. **Dataset triage** — confirm which public datasets retain raw traces. Cheap, and it
   gates everything else.
3. **H1/H2** — the precondition and the heterogeneity result. If these fail, stop.
4. **H3** — the differentiator, and the figure the proposal is built around.
5. **Proposal** — assemble, with a hardware and a clinical/data partner identified.

## Honest risks

- **The mechanism may not transfer.** Our τ rule was validated on synthetic timing
  tasks where the teaching signal (input event times) is clean. Chemosensor transients
  are slow, smooth, and noisy — "input arrival time" is not obviously defined. **This
  is the single biggest technical risk** and should be de-risked first, in H1/H2.
- **Sensor drift ≠ input-timescale drift.** Our re-tiling result is about input timing
  statistics changing; gas-sensor drift is baseline/gain change from aging. The
  defensible claim is narrower: *response kinetics are temporal and do change with
  aging, and a self-calibrating layer can track them.* Do not let a proposal blur this.
- **We are not an olfaction group.** Credibility comes from partnering, and from the
  feedback contribution being genuinely useful rather than self-serving.
- **Emergent < wired.** Our own 3d result: the grown basis recovered ~70% of the
  hand-set optimum. The pitch must be "adapts without supervision," not "beats tuned
  baselines."

## Open questions for the author

1. **Does Dognosis hold canine scent-detection data, protocols, or clinical sampling
   access?** ARIA adopts trained canines as the biological parity benchmark (sub-goal
   2) but gives no methodology for it. If yes, that is a stronger and more defensible
   position than the algorithm contribution — and it is a Workstream A card.
2. Consortium: do we have routes to a **sensor hardware** group and a **clinical
   sampling** partner? A C bid needs both.
3. Appetite: feedback-only (low cost, shapes the call), or full bid preparation?
