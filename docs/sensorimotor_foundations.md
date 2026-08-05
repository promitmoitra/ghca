# Foundations for a body: what the motor layer and the environment must be

*Design document, not a result. Every constraint below is pinned to a measurement in
[`lattice_timescale_notes.md`](lattice_timescale_notes.md); where something is a preference
rather than a finding, it says so. Written before building the closed loop, so that the loop
is not designed around whatever happens to be convenient.*

The arc so far learned **content** (a stimulus–reward interval), then **placement** (where to
gate), then produced an **observable** (a synchronous excitability window). It never had an
action, so "reward" never meant what the word normally means. This is the interface that would
change that, and the point of writing it down first is that the next three questions — growing
task-specific microcircuits, generalising across tasks, appetitive value — all depend on it and
none of them can be answered without it.

## 1. What the action is

**The action is transmission, not emission.** A Greenberg–Hastings cell returning to rest emits
nothing; refractoriness gating propagation is the one thing this substrate does natively. So a
learned τ field is a learned schedule of *when input can get through*, and the medium acts by
passing or blocking.

Measured (`lattice_sensorimotor.py`): the transmission edge sits at the learned interval at
every delay tested — 31.5 / 51.5 / 71.5 for D = 30 / 50 / 70, slope 1.0, constant +1.5 step
offset. Before the edge the motor sheet is exactly silent; after it, saturated.

Consequences for the interface:

- **The environment does not need to read a population code.** It presents an event and observes
  whether the medium passed it. That is a scalar per trial, and it is the cheapest possible
  coupling.
- **The motor layer needs no new machinery**: a sheet of the same cells downstream of H, 3×3
  convergent input, θ_M > 1, no plasticity of its own. It fires when H transmits.
- **Timing is the natural action space.** Not a direction, not a force — a *when*. Any first
  environment should therefore be a timing task, because that is what the substrate can express
  without inventing a second mechanism.

## 2. Bootstrapping is not a problem, and this was the open question

The worry (roadmap item 13) was circular: contingent reward requires a well-timed action, an
untrained medium has no timing, so reward never arrives and nothing learns.

Measured: it dissolves. An untrained sheet has τ scattered across its whole range, so **some**
fraction of cells is excitable at any probe time — transmission is *partial*, not zero. The
untrained curve ramps 0.02 → 0.08 and sits at ~60% of ceiling at the true reward time, against a
trained sheet's all-or-nothing 0.00 → 0.08. Graded credit exists from the first trial.

So the environment needs **no shaping, no reward schedule, no uncontingent warm-up phase.**
Population heterogeneity supplies the gradient for free, and what learning adds is *sharpness* —
it converts a smeared population of independent timers into a synchronised gate. This is the
single most useful thing the foundations pass established, because it removes the one identified
blocker to closing the loop.

## 3. Constraints the environment must respect

These are not stylistic. Each one has a measurement behind it, and violating any of them breaks
a result the arc already depends on.

| constraint | why — the measurement |
|---|---|
| **Reward/value must arrive diffusely, near-simultaneously across the sheet** | A *propagating* value signal destroys population synchrony (0.005–0.041 vs 0.912) while leaving per-cell accuracy fine. Each column gets a different reference, so intervals are individually correct and collectively unalignable. Synchrony needs a common temporal reference. |
| **Value must modulate, never excite** | Adding excitatory A→H on top of the modulatory gate collapses synchrony 0.168 → 0.034 and inflates the attention sheet's activity fivefold. Negative 1 and the failed 3e.2b CFC designs say the same thing. |
| **Distinguishable input streams must be distinguishable *structurally* (or by cell type), never by a property of the signal** | Amplitude fails (coincidence gating: worse than ungated at every threshold, τ→ceiling). Phase fails (gated own-input: fails at every window width). Direction fails (reward as a wave inside the medium floors τ to ~1, because a rightward wave excites cell *x* from the left and then *x−1* sees *x* on its right). Three independent attempts, one conclusion. |
| **Trials need a quiet interval long enough for waves to die** | The medium sustains reentrant activity; overlapping trials reintroduce the self-confirming fixed point that sank the lateral-input rule. |
| **Any contingency must be swept, never tested at one value** | At D=50 alone the unpaired control is *indistinguishable* from trained, because its random reward times average ≈52. A single delay would have shown a control that looked broken. Same shape as the flat τ=38.9 in the reward-edge run. |

## 4. The minimal first environment

A **timed interception**, because timing is the action space:

1. A cue is presented (stimulus patch into H).
2. An event arrives Δ steps later.
3. The medium either transmits it (motor sheet fires) or blocks it.
4. Reward is delivered — diffusely, modulatory — graded by whether the motor sheet fired near
   the event.
5. Crucially, **the action must change the world**: a transmitted event is consumed and the
   trial ends; a blocked event passes and is lost. Without step 5 the loop is still open and
   value stays informational.

Generality then comes from varying Δ across a family (already shown to work: 30/50/70), and
from two stimuli with *conflicting* Δ, which is the condition where attention should finally
beat no attention on a population readout rather than only on per-cell precision.

## 5. Two limits to carry forward, stated now rather than discovered later

**The switch is global, not spatial.** Because every trained cell recovers at the same moment,
the whole sheet becomes transmissive at once. Where the probe enters does not matter — its entry
cells recover at Δ like everything else. So there is currently *one* gate, not a spatial map of
gates, and any claim about localised task-specific circuitry has to survive that.

**The anatomy is still designed.** The hand-drawn 1-D delay line is gone, and what remains
(θ_A, θ_M, receptive-field size, gate window) are homogeneous *cell properties* rather than a
topology. That is a real difference in kind, and it is all it is: that there are several sheets,
and which projects to which, is still drawn by hand.

## 6. What this sets up: growing geometry

The natural next question is whether structure itself can be learned — task-specific
microcircuits rather than a designed stack. The useful framing is a division of labour:

> **Identity supplies the grammar; experience supplies the sentence.** Cell type constrains which
> connections are *permitted*; use decides which permitted ones are *realised*. Geometry then
> emerges inside a type-constrained space rather than being drawn.

There is already a thread of this in the results, under-sold at the time: in the reward-edge run
the *write location* was set by wave kinematics rather than design, and two interleaved delays
spontaneously occupied **different columns**. Functional territory is already allocated by the
dynamics. What is missing is **consolidation** — nothing converts a repeatedly-used functional
locus into a structural one, so nothing persists, nothing is protected from overwriting, and
nothing can be recruited again.

Two predictions to hold the next experiment to:

- **Structural plasticity will ratchet without a conserved budget.** Coupling↑ → excitability↑ →
  coincidence↑ → coupling↑ is the same positive-feedback shape as the e10 τ-ratchet and the
  self-confirming fixed point. Expect runaway unless a resource is conserved.
- **Scattering cells will break the labelled lines unless type replaces location.** Anatomical
  separation is what made every result in this arc work. A scattered mix removes it, so the
  plasticity rule must condition on presynaptic *type*. Note this does not eliminate labelled
  lines — it moves their currency from geometry to identity, which is a gain (identity is a cell
  property, not a topology) but not the same as removing the assumption.

The clean test of whether the designed anatomy was load-bearing at all: make layering a
**continuous parameter** — pure sheets, through a density gradient, to a uniform mix — and sweep
it. If type-tagging carries the labelling, the metrics should be roughly flat in that parameter.
If they collapse as soon as it leaves zero, spatial separation was doing hidden work and the
"what remains is homogeneous cell properties" claim was too generous.
