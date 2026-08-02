# Feedback to ARIA — *Hypersensory Intelligence: Olfactory Perception* (thesis v1.0)

> **Status: internal draft, not published.** This is a working document for the
> feedback ARIA invited on their programme thesis v1.0 (Claire Donoghue, PD;
> *Extending Our Perception*). It is deliberately **not** wired into the MkDocs nav
> and is not on the public site. Placeholders marked `[[ ]]` need the author's input
> before sending.

---

## 0. Who this is from, and the basis for the feedback

`[[ Name ]]` — independent researcher, no institutional affiliation. The work this
draws on is a personal open-source project: `[[ repo URL ]]`, with all code, seeded
experiment scripts, committed result artifacts, figures, and an adversarial internal
review record public from the start.

The technical basis is a computational study of **learning in excitable media** —
how a homogeneous spiking substrate can *learn its own temporal structure* from
input statistics using local, gradient-free rules. It is an illustrative
computational study on toy substrates, not an olfaction result, and we say so
plainly below. We are offering it because it bears on one specific and, we think,
under-weighted axis of the thesis: **time**.

Being independent and fully open-source is relevant to two things in the thesis: it
means any contribution we make satisfies sub-goal 3 (publicly accessible, usable by
independent groups) by construction, and it makes us eligible for the *independent*
benchmark-and-competition role described on p.13 — see §6.

The feedback below is organised as: one point about dataset design that we think is
the most consequential and least reversible decision in Workstream A (§1); the
general argument about temporal structure (§2); a concrete competition proposal
(§3); a framing for the cross-domain generality goal (§4); and short answers to
several of the specific questions posed on pp. 13–14 (§5).

---

## 1. The most consequential irreversible decision: **record and release the raw
time-resolved transient, not reduced features**

This is our strongest recommendation and it is independent of anything we would
propose to build.

Most existing chemosensor datasets — including the field's canonical drift benchmark
(Vergara et al., UCSD, 16 sensors × 36 months) — reduce each exposure to a handful of
hand-engineered scalar features per sensor (steady-state response, rise/decay fits,
integrals). The raw response curve is discarded.

We think that is the wrong default for ORO, for three reasons:

1. **The transient is not a nuisance around the steady state — it carries
   discriminative information.** Sensor response and recovery *kinetics* differ
   between analytes even where the equilibrium response does not, because adsorption
   / desorption / surface-reaction rates differ. Feature reduction chosen in 2026
   permanently caps what any 2030 model can extract.
2. **The reduction is hardware-specific, which directly undermines the programme's
   own cross-instrument transfer goal.** Hand-picked features encode assumptions about
   one sensor class's response shape. Raw time series are the common denominator that
   makes calibration transfer between instruments *possible* rather than merely hoped
   for. The thesis identifies lack of cross-instrument comparability as a core failure
   of the field (p.6); feature reduction is one of its mechanisms.
3. **It is cheap now and impossible later.** Storage for time-resolved traces is
   trivial at programme scale; re-collecting 36 months of longitudinal breath samples
   is not. Reduced features can always be *derived* from raw traces. The converse is
   false.

**Concrete ask:** make "raw, time-resolved, timestamped sensor response including
baseline, exposure, and recovery phases, at a documented and uniform sampling rate"
a hard requirement of the ORO standard — alongside the sampling protocol metadata
(flow rate, temperature, humidity, sniff/exposure duration) needed to interpret it.
Publish derived feature sets as a *convenience layer* on top, never as the archive.

This also directly serves the thesis's question *"which instruments should we use to
collect the publicly accessible data"* (p.13): our answer is less about which
instrument than about **what is preserved from it**.

---

## 2. Time is an implicit axis in this thesis; we think it should be explicit

The thesis treats olfactory perception largely as a problem in *chemical* space —
molecular embeddings, the Principal Odour Map, the intrinsic dimensionality of odour
space, compound identity. That framing is well-supported and we are not disputing it.
But it leaves **temporal structure almost entirely implicit**, and the biological
systems the programme benchmarks against are temporally organised at every level:

- **Active sampling.** Sniffing is not incidental to mammalian olfaction; it is the
  sampling clock. Odour information is acquired and processed in sniff-locked packets.
- **Plume dynamics.** In the wild, odour arrives as turbulent intermittent bursts. The
  *temporal statistics* of encounters carry information (notably about source
  distance and direction) that no static concentration reading contains. This matters
  specifically for the canine benchmark the thesis adopts — a dog locating a source is
  solving a temporal problem, not only a discrimination problem.
- **Transduction kinetics.** As in §1, response/recovery time constants are
  analyte-dependent.

The sharpest illustration is in the programme's own neighbourhood. The most cited
neuromorphic olfaction result — Imam & Cleland's olfactory-bulb circuit on Intel
Loihi (*Nature Machine Intelligence*, 2020) — operates by iterating over **sequential
gamma-frequency packets**. That gamma period is a **hand-set designer constant**.
The same is true of essentially every temporal parameter in deployed spiking
olfaction models: someone picks the time constants, offline, once.

This is the specific gap our work speaks to. In a toy excitable substrate we showed
that a **local, reward-free rule can learn the relevant timescales from input timing
statistics** rather than having them specified: the population tiles the delay
statistics it is actually exposed to, re-allocates when those statistics shift, and
self-organises into a fast/slow hierarchy under two-rhythm drive. The mechanism is
local and gradient-free (no backpropagation through time), so it is compatible with
on-chip plasticity engines. We are explicit that this is demonstrated on a toy
substrate and on synthetic timing tasks — **it is not an olfaction result** — but the
mechanism is domain-general and the failure mode we characterised (why the obvious
rule *doesn't* work) is, as far as we can find, unreported.

**Suggested change to the thesis:** add temporal structure as an explicit axis
alongside chemical representation — in the ORO capture standard (§1), in Workstream
B's notion of "representation" (a representation of *odour dynamics*, not only odour
identity), and in Workstream C's evaluation (sensors judged on kinetic fidelity and
recovery, not only sensitivity/selectivity).

---

## 3. A competition proposal: **recalibration-free longitudinal generalisation**

The thesis asks for input on competitions (p.13). We suggest one that targets the
hurdle Table 1 marks red for every deployable sensor class — **drift** — and that
tests the programme's deployability claim directly:

> **Task.** Train on early time-blocks of a longitudinal chemosensor dataset; evaluate
> on blocks months later. **No target-domain labels, no recalibration samples, no
> access to later-block statistics at training time.** Score = accuracy at month *k*
> as a function of *k*, plus the slope of degradation.

Why this specific framing:

- It is **honest about deployment**: a sensor that needs periodic re-calibration
  against labelled reference samples is not a general-purpose sensing capability, it
  is a lab instrument with extra steps.
- It cleanly separates **domain adaptation** (most existing drift-compensation work
  uses target-domain data, often labelled) from **self-calibration** (the system
  maintains itself unsupervised). Only the latter survives contact with a phone.
- Baselines already exist, so it is cheap to stand up: the Vergara UCSD benchmark
  supports it today, and the newer one-year MOX array dataset (*Scientific Data*,
  2025) extends it.
- It is **hardware-neutral** and therefore a fair cross-team benchmark.

A second competition worth considering: **cross-instrument transfer** — train on
instrument A, test on instrument B measuring the same physical samples (which the
proposed sample-banking would make possible). This directly metricises the
comparability failure the thesis identifies as the field's core problem.

---

## 4. On cross-domain generality: it has a measurable capacity limit, and that is
worth designing around

Programme sub-goal 1 is "one system, ≥3 distinct VOC application domains, no
application-specific hardware, via learned representations." We suggest treating
*how many domains one shared representation can hold* as a quantity to be **measured
early**, not assumed.

In our (toy, synthetic) continual-learning experiments we found a consistent and, we
suspect, general pattern:

- A **fixed-size shared representation has a finite capacity**. Adding sequential
  tasks works well up to a point, then degrades to the interference floor — we could
  locate the crossover precisely and it scaled with the representation's size, not
  with the quality of the learning signal.
- **Better credit assignment did not help.** We tested this directly, including a
  low-variance causal estimator, and got a clean null: all learning rules traced the
  same stability–plasticity frontier. The limit was *representational*, not
  algorithmic.
- **What did help was conjunctive capacity** — a representation that tiles
  (stimulus × context) rather than stimulus alone. And such a basis can be *learned*
  unsupervised rather than hand-designed.

The transferable claim for this programme: if three domains are to share one
representation, the binding constraint is likely to be **representational capacity and
whether the representation is context-conjunctive**, not the sophistication of the
training objective. That is measurable with a small early experiment — sequentially
train one representation on 2, 3, 4… domains and measure backward transfer — and the
answer should inform how much capacity the ORO representation targets need, before
large-scale collection is committed.

We offer this as a *framing* with toy-scale evidence behind it, not as an established
result about VOC space.

---

## 5. Short answers to specific questions posed (pp. 13–14)

**On candidate applications and scope.** We support the microbial-VOC-adjacency
rationale for expecting transfer. We would add one selection criterion the thesis
does not state: **prefer at least one application whose signal is intrinsically
temporal** (e.g. spoilage *progression*, or longitudinal health *trajectory* rather
than cross-sectional status). If all three chosen domains are static
discrimination tasks, the dataset cannot support learning temporal representations,
and the "cross-domain transfer" test is weaker for it. Longitudinal health signals
already satisfies this; food spoilage naturally does too if sampled as a time course
rather than at a decision point.

**On dataset size.** We are wary of a single headline number. A more useful framing
for planning: size requirements scale very differently for the two things the dataset
must support. (a) *Learning the representation* — plausibly modest, if the intrinsic
dimensionality claims hold; the fruit-fly existence proof cited on p.8 argues that
discriminative structure is low-dimensional. (b) *Establishing that a signal is real
and clinically actionable* — dominated by confounder control and subject count, not
sample count, and this is where the breath-biomarker literature has repeatedly
failed. We would plan capacity around (b), and expect (a) to come cheaply once the
capture standard is right. Concretely: **more subjects × more timepoints, at lower
per-sample richness, will likely beat fewer subjects × exhaustive analysis** for the
longitudinal health goal.

**On instruments.** See §1 — our substantive view is about what is preserved, not
which device. A secondary point: pairing every deployable-sensor recording with
near-simultaneous GC-MS on the *same* physical sample is what makes the dataset
useful as a *bridge* rather than two disconnected corpora.

**On open platform vs IP.** One suggestion: separate the layers explicitly. Make the
**capture standard, raw data, and benchmarks** unambiguously open (these are the
public good and the flywheel), and allow **models and sensor designs** to carry
creator IP with a requirement to publish evaluation results against the open
benchmarks. The field's problem is non-comparability, which is fixed by open data
and open evaluation — not necessarily by open models.

---

## 6. What we would be interested in contributing — and what we do not have

**Could contribute:**

- **Benchmark and competition design (p.13).** The thesis says ARIA will "fund an
  independent team to design competitions, set benchmarks and assess progress against
  the programme metrics." Independence is a *requirement* for that role rather than a
  handicap, and it is where we think we are strongest. The relevant track record is
  methodological rather than chemical: pre-registered kill conditions, adversarial
  controls that overturn the authors' own headline results, published negative
  findings, seeded and reproducible artifacts, and caveats stated next to claims
  rather than buried. Concretely we would bring the recalibration-free longitudinal
  protocol (§3), the cross-instrument transfer protocol (§3), and the capacity
  measurement (§4) as candidate programme benchmarks.
- A **local, gradient-free timescale-adaptation layer** for chemosensor front-ends —
  self-calibrating temporal receptive fields that track changing response kinetics
  without labelled recalibration data. Relevant to Workstream C's "AI-driven drift
  compensation" and "neuromorphic olfactory hardware" lines, and deliverable as
  open-source software over existing sensor arrays. As an independent contributor this
  is most realistically a component inside another team's bid, not a standalone one.
- **The temporal-axis and dataset-standard arguments** in §1–§2, offered as ORO
  standards input whether or not we are funded for anything.

**Do not have, and will not claim:**
- Any validated olfaction result. Our evidence is a computational study on toy
  excitable substrates and synthetic timing tasks.
- Sensor hardware or materials expertise, or a wet-lab / clinical sampling capability.
- Any proprietary dataset. We work only with public data — which is a constraint we
  are happy with, since the programme requires open outputs anyway.
- Evidence that our timescale mechanism improves chemosensor performance. That is
  a hypothesis we would need public (or ORO) data to test — and the honest position
  is that it should be *tested*, not assumed.

**One clarifying question back to ARIA:** how will *biological parity* (sub-goal 2)
be operationalised against the canine benchmark — matched-sample head-to-head
trials, or literature-derived sensitivity targets? This materially affects how teams
should design evaluation, and it is the sub-goal with the least methodological
detail in v1.0.
