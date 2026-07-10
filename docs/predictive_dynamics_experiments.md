# E-Series Extension: Predictive Dynamics (naturalistic sequences)

*A continuation of the [learning E-series](learning_experiments.md). E0–E7 shaped
the substrate by reward and read cognitive functions off it; **E8** asks whether the
same substrate makes genuine **time-forward predictions** of upcoming sensory input
for progressively complex stimuli — and how that prediction is (dis)similar to
predictive coding. Motivated by the auditory-sequence prediction literature
(Chait/Bianco — *Neural mechanisms of time-forward predictions for naturalistic
auditory tone sequences*; and the predecessor showing neural activity integrates the
last ~7 tones to predict the next pitch, with sharper anticipatory tuning for more
predictable sequences) and the wave-based prediction line (Muller / Benigno et al. —
traveling waves ignite short-term predictions of sensory input).*

---

## 0. Status and scope

Specifies **what to build and measure**, not the implementation; status **PLANNED**.
Reuses the substrate (`ghca_net`), the learner (`ghca_learn`: Line A conduction =
transition graph, Line B timescales = integration window), the GVF/Horde readout
([E6](e6_results.md)), and the `do(θ)` causal machinery (the
[C-series](causal_experiments.md) and its
[spiral extension](causal_spiral_experiments.md)). Headline: **prediction as forward
dynamics + tunable history integration, disambiguated from hierarchical predictive
coding.**

---

## 1. The thesis and its disambiguation from predictive coding

- **Predictive coding** (Rao & Ballard 1999; Friston 2005): an explicit
  *hierarchical generative model*; dedicated *prediction-error units* computing
  (input − top-down prediction) at each level; learning = minimize precision-weighted
  error (approximate inference); representation = the explained-away residual.
- **This substrate:** no error-unit hierarchy; prediction *is* the forward evolution
  of a recurrent excitable medium whose `τ` and `W` were shaped by experience; the
  only error is the **global scalar** TD/surprise `δ`; representation = emergent
  dynamics selected by action (inside-out).
- **The honest overlap** (why disambiguation is non-trivial): the TD error `δ` is a
  *value* prediction error (dopamine-as-TD-error); E6's GVF demons already predict
  future signals (Horde = "predictive knowledge"); a traveling wave already
  forward-models moving input. We are a *predictive* system — just not a
  *predictive-coding* one.
- **Two model-distinguishing predictions this program tests:**
  1. **The integration window is a `do(τ)` property** (Line B / C-series), not a
     hierarchy-depth property — sweeping `τ` slides the "~7 tone" temporal receptive
     field. A causal statement the MEG correlation cannot make directly.
  2. **Prediction and representation share one medium** — there is no separable
     top-down pathway to lesion; contrast PC, where lesioning top-down abolishes
     prediction while sparing feedforward.

---

## 2. Substrate / framework deltas

- **Tonotopic sensory map.** `M` tone channels on a 1-D ordered feature axis (a
  "cochleotopic" line); a tone = a population code on its channel. Generalizes E1's
  `K` unordered stimuli to an ordered axis (so "near in pitch" = "near on the map").
- **History integration = recurrent medium + `τ`.** The hidden medium's fading
  memory integrates recent tones; the effective window is set by `τ` (Line B) — the
  object `do(τ)` sweeps.
- **Transition learning = Line A.** Reward- or self-supervised plasticity on the
  `S → H → (next-tone)` edges learns the sequence's transition graph (the design
  doc's own name for Line A's credit).
- **Predictive readout = GVF demon (E6).** A next-tone-prediction demon reads the
  medium and predicts the upcoming channel; **prediction accuracy vs sequence
  predictability** is the primary metric. A self-supervised next-step objective
  (predict-then-compare, no reward) is the intrinsic-prediction variant.
- **Forward-wave option.** On the tonotopic map, a traveling wave whose speed = the
  tone rate pre-activates the expected next channel (the Muller mechanism); wave
  speed is a `do(θ)` (`τ`, coupling) property — so predictive accuracy is
  parameter-controlled.
- **Surprise.** A violation response = deflection of a global order parameter / `δ`
  (**not** a per-channel error field) — the key PC-disambiguation observable.

---

## 3. The stimulus ladder (progressive complexity)

### E8.0 — Instrument the predictive readout *(gate)*
Confirm the tonotopic map, the next-tone GVF readout, and a predictability metric
are well-defined, and that a periodic sequence yields above-chance next-tone
prediction. *Gate:* if the medium can't anticipate even a periodic sequence, revisit
the `τ`/coupling operating point before proceeding.

### E8.1 — Periodic sequences (pure predictability)
Fully predictable repeats. *Hypothesis → prediction:* anticipatory **pre-activation**
of the next channel emerges; prediction accuracy → high. *Discriminator:*
pre-stimulus activity in the *expected* channel exceeds the others *before* onset.

### E8.2 — First-order Markov (transition structure) → Line A
Tones drawn from a transition matrix. *Prediction:* Line A learns the transition
graph; next-tone prediction tracks transition probabilities (better for
higher-probability transitions). *Discriminator:* accuracy scales with transition
entropy; ablating Line A plasticity removes the learned structure while leaving the
periodic (E8.1) case intact.

### E8.3 — Naturalistic pitch random-walk (history integration) → `do(τ)` *(headline)*
A slowly-evolving pitch (random walk on the tone axis, à la Bianco). *Prediction:*
the medium integrates the last ~`N` tones, and **`N` is set by `τ`** — `do(τ)` slides
the integration window (the causal headline that turns the MEG correlation into an
intervention). *Discriminator:* the number of past tones the prediction depends on
moves monotonically with `τ`; at fixed `τ` the window is stable across sequences.

### E8.4 — Violations / oddballs (surprise; PC disambiguation)
Insert low-probability tones. *Prediction:* a **global** surprise signal
(order-parameter / `δ` deflection) scaled by improbability — *not* a localized
per-channel error unit; and prediction vs representation are **not dissociable**
(no separable top-down pathway to lesion). *Discriminator:* the surprise is
aggregate/global, contrasting the per-feature, hierarchically-localized error of PC.

### E8.5 — Nested / hierarchical regularities *(optional stretch)*
Structure at two timescales (local transitions within a slowly-changing context).
*Prediction:* a slow loop (the [E5](e5_results.md) option) supplies the context/regime
while fast dynamics predict within it — multi-timescale prediction, echoing the
theta–gamma nesting literature. *Discriminator:* ablating the slow sub-population
abolishes context-conditioned prediction while sparing within-context prediction.

---

## 4. Analyses that cut across the ladder

- **Prediction vs `do(τ)`.** The integration-window ↔ `τ` law (E8.3) — resting on the
  causal backbone the C-series spiral extension (C5–C7) validates for `do(θ)`.
- **Reward vs self-supervision.** Does prediction acquire *without* reward
  (intrinsic / inside-out) or does it need it? Distinguishes a forward model in the
  dynamics from a purely value-driven one.
- **PC contrast table.** For each rung, PC's prediction vs the substrate's, and the
  distinguishing observable (error locality, dissociability, window control).
- **GVF panel.** Next-tone demon accuracy / calibration vs predictability (reuses E6),
  and representational-similarity of predictive vs stimulus readouts.

---

## 5. Reference anchors

- Chait/Bianco — *Neural mechanisms of time-forward predictions for naturalistic
  auditory tone sequences*; and *Neural Integration of Stimulus History Underlies
  Prediction for Naturalistically Evolving Sequences* (history integration over
  ~7 tones; predictability → anticipatory tuning) — the target phenomenon.
- Muller; Benigno et al., *Nat. Commun.* 2023 — traveling waves ignite short-term
  predictions (the wave forward-model).
- Rao & Ballard 1999; Friston 2005 — predictive coding (the foil to disambiguate from).
- Sutton, Modayil et al. 2011 — Horde / GVFs (predictive knowledge); Sutton 1988 —
  TD learning (value prediction error).
- [`learning_experiments.md`](learning_experiments.md) ·
  [`causal_experiments.md`](causal_experiments.md) ·
  [`causal_spiral_experiments.md`](causal_spiral_experiments.md) ·
  [`synthesis.md`](synthesis.md).

---

## 6. One-line spine

E8.0 grounds the tonotopic predictive readout → E8.1 periodic anticipation →
E8.2 Line A learns transitions → E8.3 `do(τ)` sets the history-integration window
(causal headline) → E8.4 global surprise, no error-unit hierarchy (PC disambiguation)
→ E8.5 multi-timescale nested prediction.
