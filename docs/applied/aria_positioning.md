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

## Constraints that define the position

- **Independent.** No institutional affiliation, no consortium, no lab.
- **Open source only.** This is a personal open-source repo; no proprietary data is
  available or usable. The programme requires open outputs anyway, so this is
  alignment rather than sacrifice.
- **No wet lab, no hardware, no clinical access.** Public datasets and CPU.

These rule out the obvious "bring your own data" play. They do not rule out the two
things below.

## Where it fits the programme (honest mapping)

| Workstream | Fit | Reality check |
|---|---|---|
| **A** — open hardware + dataset | **Standards input only.** The §1 argument (preserve raw time-resolved transients) is worth making regardless of funding. | Not a hardware or data-collection team. No data to contribute. |
| **B** — learning representations | **Partial and orthogonal.** B means *chemical* representation (molecular embeddings, odour-space dimensionality); we offer *temporal* representation. | Do not position as competing with B. Position as the missing second axis. |
| **C** — novel sensing capability | **Topically the best fit.** C explicitly names "neuromorphic olfactory hardware" and "AI-driven drift compensation," and Table 1 marks drift red for every deployable sensor class. | C is hardware-flavoured and we are software over existing arrays. As an independent with no hardware partner, a standalone C bid is not credible — this is a component inside someone else's bid. |
| **Benchmarks / competitions** (p.13, cross-cutting) | **Strongest fit, and the one that *requires* what we are.** The thesis states ARIA will "fund an independent team to design competitions, set benchmarks and assess progress against the programme metrics." | Independence is an eligibility requirement here, not a handicap. Needs credibility in *methodology*, which is demonstrable, rather than in chemistry, which is not. |

**Revised posture.** The default assumption — component contributor inside a
consortium — is right for the *mechanism*. But the benchmark/competition role is a
named, funded, explicitly-independent slot, and it is the one where being an
unaffiliated open-source researcher with a rigorous public track record is the
qualification rather than the obstacle.

What that role actually rewards is methodology, and it is the thing this repo has
most of: kill conditions declared before running, adversarial controls that overturned
the authors' own headline results (the 3c/P2 causal-credit null; the E7 integer-θ
artifact), published negative results kept rather than buried (the e10 ratchet), fully
seeded and reproducible artifacts, and caveats stated adjacent to claims. A programme
whose central failure mode is "overpromising and underdelivering, particularly about
generalisation and in-the-wild usage" — the thesis's own diagnosis of DARPA Real Nose
on p.6 — has an obvious use for a benchmark team with that disposition.

The candidate benchmarks are already drafted in the feedback: recalibration-free
longitudinal generalisation, cross-instrument transfer on shared physical samples, and
the capacity measurement for cross-domain generality.

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
5. **Decide the line** (see Open Questions) — benchmark-design bid, which needs the
   feedback plus the methodological track record and no partner; or mechanism
   contribution, which needs H1–H3 plus a host team. Assemble accordingly.

## Honest risks

- **The mechanism may not transfer.** Our τ rule was validated on synthetic timing
  tasks where the teaching signal (input event times) is clean. Chemosensor transients
  are slow, smooth, and noisy — "input arrival time" is not obviously defined. **This
  is the single biggest technical risk** and should be de-risked first, in H1/H2.
- **Sensor drift ≠ input-timescale drift.** Our re-tiling result is about input timing
  statistics changing; gas-sensor drift is baseline/gain change from aging. The
  defensible claim is narrower: *response kinetics are temporal and do change with
  aging, and a self-calibrating layer can track them.* Do not let a proposal blur this.
- **Not an olfaction group, and no partner.** Credibility for the *mechanism* line
  depends on finding a host team; the benchmark line does not, which is part of why it
  ranks higher.
- **Emergent < wired.** Our own 3d result: the grown basis recovered ~70% of the
  hand-set optimum. The pitch must be "adapts without supervision," not "beats tuned
  baselines."
- **Benchmark roles are usually awarded to institutions.** The thesis says
  "independent team," and an individual is the limiting case of that. Realistic
  mitigations: partner with one other independent/academic for credibility, or scope
  the contribution as *benchmark design input* to whoever wins the role rather than
  bidding for it outright.

## Open questions for the author

1. **Which line do you actually want?** They diverge:
   - *Benchmark/competition design* — plays to the repo's real strength (methodology),
     needs no data, no hardware, no partner, and independence is an eligibility
     requirement. Lower ceiling, much higher probability.
   - *Mechanism contribution (self-calibrating timescales)* — more interesting
     technically, but needs a host team's bid to sit inside, and needs H1–H3 to survive
     contact with real chemosensor data first.
   They are not exclusive: the preliminary experiment strengthens both, and the
   feedback response serves both.
2. Are you willing to approach a **UK academic group** in chemical sensing or
   neuromorphic engineering as a collaborator? That is the single change that most
   improves the mechanism line, and it is a prerequisite for any Workstream C bid.
3. Appetite: feedback-only (days, shapes the call, no commitment), or feedback plus
   the H1–H3 preliminary experiment (weeks, and the thing a real bid needs)?
