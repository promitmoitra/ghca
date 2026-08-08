# Lattice Results — learning a timescale, then a reward interval, on an excitable sheet

*Eleven experiments porting the input-timing `τ` rule from the non-spatial pool (3d/3e) onto a
2-D Greenberg–Hastings sheet, then building outward: a localised sensory strip, an attention
chain, reward as a fourth edge, a value chain, 2-D layers, an action primitive, and two
structural negatives. Roadmap items [`3e.5`, `3e.7`–`3e.19`](next_steps.md). Running lab
notebook with full derivations and the mistakes made along the way:
[`lattice_timescale_notes.md`](lattice_timescale_notes.md). Interface spec:
[`sensorimotor_foundations.md`](sensorimotor_foundations.md).*

*Headline conditions re-run at **n=20** (artifacts tagged `_n20`); `tonic` is n=5 and
`attention_gate` remains n=3 — flagged in Caveats.*

## The gap this closes

3d/3e established that a per-cell refractory period `τ` can be learned from **input timing** on
a non-spatial pool, where every unit fanned in from rhythmic drive sources. The natural
demonstration is a 2-D lattice — waves are legible on sight where a histogram is not — but the
rule had never been run on a **recurrent spatial medium**. Porting it did not just visualise the
pool result; it produced different behaviour, and the differences are the content here.

## Mechanism

Each cell is a Greenberg–Hastings unit: rested → fires when active neighbours ≥ `θ` → refractory
for `τ` steps → rested. Only **`τ`** is plastic (`act`, the excited duration, stays fixed
throughout — a deliberate choice, not an oversight). The learning rule is a delta rule on the
interval between input events:

```
τ ← τ + η · (Δt − τ)
```

**One fix was required before anything worked: edges, not levels.** Defining an "input event" as
*suprathreshold input* is nearly always true on a lattice — with `θ=1` some neighbour is almost
always active, and a passing wave holds input high for `act` consecutive steps, so the measured
interval is ~1 and carries no rhythm, and τ collapses to the floor of its [3,34] range. An input
event must be the **rising edge** of suprathreshold drive. *(The specific floor value quoted in
the notebook, 3.1, was observed during development and has no committed artifact — the level-
detection variant was never kept as a condition.)* In the 3e.2 pool the drive sources
fired as discrete pulses, so every input already *was* an edge and the distinction was invisible.

## Result 1 — the rule works on a lattice; the old self-referential rule still fails

Single rhythm, P=6, L=96, 6000 steps, **n=5** (`lattice_timescale_demo_single.json`):

| rule | mean τ | fraction near P | Moran's I |
|---|:--:|:--:|:--:|
| **input-timing** | **7.9** | **0.55** | **+0.484** |
| own-firing (e10) | 33.8 (ceiling) | 0.01 | −0.002 |

The input rule pulls τ toward the drive period and produces genuine spatial structure; the
own-firing rule rails to the ceiling with none — the **e10 ratchet reproduced on a lattice**.

*Provenance note.* This table was **unsourced until now**: the committed artifact was a
*two-rhythm* run (P_F=6, P_S=24, n=2) and the single-rhythm numbers came from an earlier version
of the script that no longer existed. Reproducing it required two fixes — making the periods
configurable, and **decoupling the interval clip from `P_S`**, which had been derived as `P_S × 2`
and so silently capped learned τ at 12 when the slow period was lowered (a single-rhythm run read
τ = 11.95, the clip, rather than the ceiling). Earlier write-ups quoted 0.51 and 0.32 for the last
two columns; the reproducible values are 0.55 and 0.484, so the spatial-structure result is
**stronger** than previously reported.

## Result 2 — exogenous timing does *not* penetrate a recurrent medium

A **global** afferent lets the rule track a drive period exactly (τ = 6.0 / 18.0 / 6.0 across a
6→18→6 schedule, |τ−P| = 0.0). But a global afferent puts every cell effectively *in the sensory
layer*, so the recurrence is not doing the learning. Confining the afferent to a **sensory strip**
along one edge — the repo's own `roles` convention — changes the answer: learning from every input
**never locks anywhere, not even in the adjacent column**. |τ−P| is 1.54 in the column abutting
the strip, jumps to ~2.0 immediately beyond it and drifts to 2.92 at the far wall — but the
profile is **non-monotone** (1.54, 2.07, 2.06, 2.01, 1.94, 1.87, 1.80, 1.82, …), so this is a
failure to lock at any depth rather than a clean gradient with distance. Several wavefronts arrive
per drive cycle and τ settles at **8.9 against a drive period of 6**, i.e. the estimator recovers
intervals *longer* than the period, not shorter. The rhythm is present in the input but not
*recoverable*.

## Result 3 — an attention chain carries a timing reference: a clock, not a filter

A 1-D chain of the same cells, orthogonal to the sensory edge, seeded by the beat and propagating
at one cell per step, locks τ to the drive period at **every** depth (|τ−P| = 0.00, depth 92, at
both P=6 and P=12) — matching the global-afferent bound without a global afferent. Gate arrival
phase rises by exactly one step per column and wraps mod P (a broadcast would be flat), and the
chain's own τ is set from its pulse width alone, never from P.

**But it works as a clock, not a filter.** Using the gate to select *when a cell listens to its
own input* fails at every window width; a window as wide as the period reproduces the ungated
numbers exactly. The timing information lives in the gate, not in what the gate selects.

## Result 4 — reward as a fourth edge: τ encodes a stimulus–reward interval (n=20)

Four edges: stimulus patch left, attention chain along the top, **reward chain along the bottom**
travelling right-to-left. The rule becomes three-factor with three physically distinct factors —
*who* from the cell's own firing time (so the sheet's 2-D wave geometry selects who learns),
*when* from the attention chain, *what* from the reward chain — and τ converges to **that cell's
own stimulus-to-reward interval**. Scored on probe trials as `|(t_rew − t_fire) − τ|`:

| condition | \|Δ\| | within ±2 | cells written | σ_y | τ (D=8→32) |
|---|:--:|:--:|:--:|:--:|:--:|
| **paired, ungated** | **0.16–0.19 ± 0.00** | **0.97–0.98** | 0.49–0.67 | 7.1–7.5 | 31.3 → 43.1 |
| paired + slow gate (K=3) | **0.00 ± 0.00** | **1.00** | 0.016 | 8.9 | 26.8 → 38.8 |
| paired + fast gate (K=1) | 0.80 ± 0.00 | 0.87 | 0.004 | **0.00** | 3.0 → 3.0 |
| unpaired *(the control)* | 24.6 / 14.8 / 8.1, **sd 7.5 / 5.9 / 3.3** | 0.00–0.12 | 0.90 | 7.2 | 37.1 **flat** |
| no reward *(untrained)* | 30.1–30.8 ± 0.35–0.53 | 0.04 | — | 24.9 | 46.5 |

The unpaired control delivers **exactly as many reward events on the same line** and fails; its τ
sits at 37.1 regardless of D. Read the flat τ, not the shrinking error — |Δ| falls from 24.6 to
8.1 only because D=32 approaches the random-delay mean.

**The 2-D medium is load-bearing**: after learning, τ still varies by σ_y ≈ 7 along y at fixed x
(the distance-from-patch term) while the field is accurate to 0.2.

**The fast gate is degenerate** — it meets reward exactly where stimulus and reward arrive
together, so the interval is ≈0, τ pins to the floor (σ_y = 7×10⁻⁶, numerically zero) and only
the *column index* carries the delay. Its 0.87 "within ±2" passes for the wrong reason.

## Result 5 — a value chain removes the last hand-set constant (n=20)

The chain's speed K was the remaining designed constant. A **value chain launched backward from
the write site** carrying τ-saturation (`frac at τ_min − frac at τ_max`) retunes each per-column
delay it passes. It must be an edge, not a scalar: the gate's arrival time at column x is the
*cumulative* delay of every cell upstream, so a badly placed write is the fault of the delays
behind it. All delays start at 1 — the degenerate value.

| condition | d̄ | write column | saturation | \|Δ\| | within ±2 |
|---|:--:|:--:|:--:|:--:|:--:|
| **plastic** | **1.50** | **27.0 ± 0.0** | +0.01 → **+0.02** | **0.00** | **1.00** |
| *(all rows: D=20 only — the n=20 re-run covered one delay)* | | | | | |
| fixed K=1 *(degenerate)* | 1.00 | 41.0 ± 0.0 | +0.01 → **+1.00** | 0.80 | 0.87 |
| fixed K=3 *(hand-set target)* | 3.00 | 20.0 ± 0.0 | +0.00 | 0.00 | 1.00 |
| plastic-unpaired *(control)* | 1.36 | **29.1 ± 3.3** | +0.01 | 16.02 | 0.14 |
| plastic-two *(interleaved delays)* | 1.41 | **27.0 and 33.0** ± 0.0 | +0.01 → +0.00 | 0.00 | 1.00 |

Saturation falls from the +1.00 that fixed-K=1 reaches — *every* written cell pinned at the floor
— to ≈0, and alignment matches the hand-set K=3 condition **without being given the value**. It
settles on d̄ ≈ 1.5 rather than 3: a different valid placement, not a recovery of the designer's
answer. The **predicted lock-in did not occur** (d settles far below its ceiling of 12, no
oscillation) — the feedback is a saturation *boundary*, not a rhythm, and a boundary cannot
confirm itself the way a period can.

The credit-assignment geometry is visible in the learned parameter: the value pulse only travels
left, so upstream delays train further (mean 1.83) than downstream ones (mean 1.26). Downstream
delays did still move — by roughly a quarter of the upstream change — and only the far end sits at
the initial 1.00; the *asymmetry* is the artifact of the mechanism, not a hard boundary.

## Result 6 — layers instead of edges: a synchronous burst, and two hard requirements (n=20)

Replacing both chains with 2-D sheets (value / attention / hidden, plus a motor readout) yields
the arc's one genuinely emergent *output*: if every cell's τ equals its own stimulus-to-reward
interval, cells that fired at **different** times all become excitable again at the **same**
moment. A travelling wave becomes a synchronous population burst.

| condition | sync (D=50) | sync (D=70) | peak time |
|---|:--:|:--:|:--:|
| **ungated** | **0.637 ± 0.000** | **0.912 ± 0.000** | D + 1 |
| A-gated | 0.167 ± 0.005 | 0.346 ± 0.007 | D + 1 |
| **propagating value** | **0.021 ± 0.002** | **0.042 ± 0.003** | D + 49.6 / D + 38.9 |
| **excitatory coupling** | **0.034 ± 0.001** | **0.077 ± 0.001** | D + 1 |
| untrained | 0.047 ± 0.003 | 0.070 ± 0.003 | 88 |
| unpaired | 0.006 ± 0.022 | 0.009 ± 0.026 | 89 |

At D=70, **91% of recovery events land within ±3 steps of the reward time** with zero seed
variance, against 7% untrained. Unpaired scores ~0.01 — its τ converges consistent-but-*wrong*,
so the burst is sharp and in the wrong place, a more informative failure than noise.

Two mechanism requirements, both newly visible only in the layered form:

- **Value must be diffuse.** A propagating value signal collapses synchrony to 0.021–0.042 while
  leaving per-cell accuracy fine (|Δ| 0.12–0.22). Each column gets a different reward reference,
  so intervals are individually correct and collectively unalignable. **Population synchrony
  needs a common temporal reference** — a functional argument for volume transmission.
- **Coupling must be modulatory.** Adding excitatory A→H on top of the gate collapses synchrony
  and inflates the attention sheet's own activity **4.86×** at D=50 and 4.88× at D=70
  (0.00769 → 0.03742). Same cells, same signal,
  same timing; only excitatory instead of modulatory.

## Result 7 — the action is transmission, and bootstrapping is free (n=20)

A GH cell returning to rest **emits nothing**, so the synchronous recovery above is not an output.
What it *is* is a synchronous window of **transmissivity** — refractoriness gating propagation.
The action primitive is therefore "be transmissive at the expected time", which makes a motor
sheet definable with no new machinery and makes the environment's job trivial: present an event,
observe whether it passed.

| | D=30 | D=50 | D=70 |
|---|:--:|:--:|:--:|
| **trained — transmission edge** | **31.5** | **51.5** | **71.5** |
| untrained | 56.2 | 61.7 | 64.3 |
| unpaired | 59.3 | 59.7 | 64.0 |

Slope 1.0 against D with a constant +1.5-step offset; motor response at D is *identically* 0.000
across all 20 seeds in every trained condition. Unpaired's edge barely moves with D (59.3 / 59.7 / 64.0) because its reward times are drawn
uniformly from [5, 130) — mean 67 — so it reflects that distribution rather than the tested delay.
The trained/unpaired gap is 27.8 / 8.2 / 7.5 steps at D=30/50/70, i.e. **narrowest at the longest
delay**. At n=3 the two coincided at D=50 specifically, which is what made the sweep necessary; an
earlier write-up attributed that coincidence to a "~52" distribution mean, which is not a quantity
in this experiment.

**Bootstrapping is free.** An untrained sheet has τ scattered across its range, so some fraction
of cells is excitable at any probe time, so transmission is *partial* rather than zero. But the
amount is strongly delay-dependent: untrained transmission at t=D is **4.0% / 21.7% / 70.4%** of
the trained ceiling at D=30/50/70. Graded credit therefore exists from the first trial without
shaping **at long delays**, and is very weak at short ones — the original blanket "~60%" was the
D=70 figure generalised too far. What learning adds is **sharpness** — a smeared population
of independent timers becomes a synchronised gate.

## Result 8 — avoidance is not sign-symmetric with approach; the reason is a theorem (n=20)

For a **dip** (transmission low at D, high at both D±δ) a cell must be rested at D−δ and
refractory at D. But a rested cell *stays* rested — nothing re-fires it — so every cell
transmitting at D−δ also transmits at D. Population transmission is therefore **monotone
non-decreasing in probe time**. Measured max decrease over all four rules × 2 delays: **+0.0000**. The committed
`mono_viol` metric is computed on the seed-mean curve; recomputing it independently **per seed**
from the stored per-seed arrays also gives +0.0000, so the result does not depend on the
averaging — but the stored metric alone would not have established that.

Refractoriness is a **prefix** of the time after firing, never a window inside it. A cell can
express "pass after T"; it cannot express "block near T, pass either side". Avoidance is only
expressible as **postponed approach**.

- The naive sign flip is **bistable and half wrong**: it moves τ *away* from the interval, so
  cells above climb to the ceiling (correct) and cells below fall to the floor (wrong — rested
  and transmissive at D). At D=50, 51% end at the ceiling and 13% at the floor.
- A prediction of mine that was **wrong**: I expected the failure-gated ratchet to run to the
  ceiling and reproduce the e10 ratchet. It does not (floor 0.00, ceiling 0.00) — gating on
  failure self-limits, because the moment τ exceeds the interval the cell stops updating.
- Both working variants need a **margin constant** that approach does not.

## Negative 1 — plastic cell identity fails structurally (n=20)

Making identity plastic (identity = the firing threshold θ, the propagator↔coincidence axis) with
the minimal one-rate/one-context homeostatic law `θ ← clip(θ + r(ā − a*), 1, 4)`:

| | r=0 | r=30 | r=300 |
|---|:--:|:--:|:--:|
| corr(θ,x) — differentiation *(a\*=0.08)* | +0.00 | −0.31 | **−0.40** |
| propagation reach *(a\*=0.08)* | **0.86 ± 0.00** | **0.21 ± 0.16** | **0.21 ± 0.14** |
| \|τ−P\| *(a\*=0.08)* | 2.3 | 14.8 | 15.1 |
| corr(θ,x) *(a\*=0.16)* | +0.00 | −0.09 | −0.26 |
| propagation reach *(a\*=0.16)* | **0.86 ± 0.00** | **0.44 ± 0.09** | **0.48 ± 0.07** |
| \|τ−P\| *(a\*=0.16)* | 2.3 | 1.9 | 1.7 |

A gradient does form (corr → −0.40) with no wiring. But **any θ motion costs a large fraction of
the propagation reach at either set-point**: 0.86 falls to 0.21 at a\*=0.08 (a **76%** loss) and to
0.44 at a\*=0.16 (**49%**) as soon as r > 0. The set-point only decides whether τ-learning dies
*as well*. (An earlier write-up said "about half" — accurate for a\*=0.16, an understatement for
a\*=0.08.) The failure is **init-independent**:
starting from all-readers (θ=4, a dead sheet) converges to the same fixed point as from
all-propagators, so it is structural, not a tuning miss. The rule's interior fixed point is
uniform activity ā = a\* for every cell, while faithful relaying requires a relay to fire once
per wave — a rate set by the drive period, not by an arbitrary set-point.

**Diagnosis:** own-activity is the wrong context signal, because in an excitable medium **activity
*is* load-bearing relay work**, so activity-homeostasis punishes the cells propagation depends on.

## Negative 2 — tonic drive has no window (n=5)

Two dead ends (the transmission step rather than a band-pass peak; the impossibility of selective
avoidance) both rest on single firing. Weak tonic drive breaks that assumption — and unlocks
neither.

Full sweep, all eight q values (`approach` rule unless noted):

| q | 0 | 1e-7 | 3e-7 | 1e-6 | 3e-6 | 1e-5 | 3e-5 | 1e-4 |
|---|:--:|:--:|:--:|:--:|:--:|:--:|:--:|:--:|
| band-pass peakiness (need > 1) | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | 0.00 | nan† |
| avoidance dip (`avoid-margin`) | nan | nan | nan | nan | nan | nan | nan | nan |
| monotonicity drop, `approach` | +0.000 | +0.000 | +0.003 | +0.003 | +0.023 | +0.031 | +0.026 | +0.024 |
| monotonicity drop, `avoid-margin` | +0.000 | +0.000 | +0.003 | +0.003 | +0.032 | +0.039 | +0.042 | **+0.051** |
| **τ within ±2** | **0.98** | 0.98 | 0.98 | 0.93 | **0.63** | **0.47** | **0.31** | 0.29 |
| reentrant activity | 0.00 | 0.04 | 0.11 | **0.31** | 0.64 | 0.87 | 0.96 | 0.98 |

†`nan` because transmission at the far reference point is itself ~0 there, so the ratio is
undefined — not a peak.

The ordering is the finding: flooding and learning-collapse arrive **before** the monotonicity
violation they were meant to buy, and that violation stays small — at most **+0.051**, for
`avoid-margin` at the largest q. Peakiness never
exceeds 0; the avoidance dip is absent at every q. So the escape hatch is closed by measurement
and Result 8 is *stronger* than when proved.

**Why:** stochastic re-firing happens at *random* times, so excitability windows **decohere**
rather than becoming periodic. A band-pass population response needs excitability periodic *and
phase-aligned* — pointing at a **rhythmic** background, converging with Result 6's finding that a
diffuse near-simultaneous signal preserves synchrony where a propagating one destroys it.

## Interpretation

**A design law for this substrate, from four independent failures.** A raw signal *magnitude or
geometry* cannot serve as a label or driver here:

| attempted label | result |
|---|---|
| **amplitude** (coincidence gating, ≥2 neighbours) | worse than ungated at every threshold; τ → ceiling |
| **phase** (gate the cell's own input) | fails at every window width; window = period reproduces ungated exactly |
| **direction** ("input from my right = reward") | floors τ to ~1: a rightward wave excites x from the left, then x−1 sees x on its right |
| **activity level** (homeostatic identity) | punishes the relays propagation depends on |

Every time, the fix was a **structural or predictive** signal — a separate pathway, or a reference
carried by dedicated cells. This is the single most transferable finding in the arc.

**Deterministic learning, stochastic controls.** Across five experiments at n=20 the pattern is
uniform: conditions whose learned quantity *converges* reproduce to the digit (motor@D sd 0.0000;
write column ± 0.0; sync ± 0.000), while all seed variance lives in conditions where nothing
converges (unpaired sd 3.3–11.3). Practical consequence: **n=3 is adequate for treatments on this
substrate and misleading for controls.**

**Prefer bounded proportions to means of unbounded errors.** Six instrument failures in this arc
were all *derived* metrics. The mechanical reason, visible in Result 6's numbers: `|Δ|` is a mean
over a heavy-tailed per-cell error distribution (gated |Δ| 2.17 ± **2.14** — spread exceeding the
effect), while fraction-within-tolerance is bounded and better-behaved (± 0.07–0.09 at the narrow gate
widths; it degrades to ± 0.15–0.20 at the widest, so "bounded" is the guarantee, not "always
tight").

## Caveats / open items

- **This is not reinforcement learning.** There is no action that changes the world and no policy.
  The value signal in Result 5 is **homeostatic** ("do not saturate your own representation"), not
  appetitive. Result 7 supplies an action primitive but the loop stays open.
- **The anatomy is designed.** Results 4–6 made the *content* and the *placement* emergent, not
  the architecture: how many sheets and which projects to which is still drawn by hand. Negative 1
  is the attempt to fix that, and it failed structurally.
- **The transmission switch is global, not spatial.** Every trained cell recovers at the same
  moment, so there is one gate, not a map of gates.
- **The attention gate is not justified in this task family.** Its cost is certain (population
  readout lost by 8.2 sd at D=70, 12.4 sd at D=50); its benefit is weak and delay-dependent
  (fraction within ±2 separates at D=70, 2.2 sd, but not at D=50, 1.2 sd) and **absent on the
  primary error metric** (0.9 sd at D=70, 0.1 sd at D=50). An earlier "threefold error reduction"
  claim was retracted — see the notes.
- **Coverage gaps.** `tonic` is **n=5**, not 20. `attention_gate` remains **n=3** — deliberately
  deprioritised as least load-bearing, and its result is a step function that n=3 shows cleanly,
  but it is not at the repo's standard. Result 1 is **n=5** (regenerated); Result 2's depth sweep
  is **n=3**; Result 3 is **n=3**. Result 5's n=20 re-run covered **D=20 only**.
- **Negative 1's init-independence claim rests on a weaker run than the rest.** The reverse-init
  check (θ=4 start) is **n=1 at L=48 and 2500 steps**, not the L=64 / 6000-step / n=20 setting of
  the main sweep. Its θ statistics land at the same fixed point, but it should be re-run at the
  main setting before the "structural, not a tuning miss" wording is leaned on further.
- **Selective avoidance is impossible**, not merely unachieved — and now also unreachable via
  re-firing.
- **The lattice is not a free upgrade on the pool.** The pool results (3d/3e, n=20) remain the
  validated science for the *capacity* claims; the lattice has its own, different behaviour.

## Operating point

L=64 sheet (L=96 for Results 1–3), 4-neighbour, walls on all sides except where noted, `act=2`,
`θ=1` (θ_A = θ_M = 3 for coincidence layers), τ ∈ [3, 90] for interval-learning experiments and
[3, 34]–[3, 40] for period-tracking ones, η = 0.15–0.4, 30–50 training trials + 3 probe trials,
ITI 220–260 steps. Every RNG is a `default_rng(seed)` threaded explicitly — no global NumPy RNG.

## Reproduce

```bash
# headline conditions at n=20 (outputs tagged _n20 so variants cannot clobber each other)
SM_N=20  SM_TAG=_n20  python3 experiments/lattice_sensorimotor.py     # Result 7
AV2_N=20 AV2_TAG=_n20 python3 experiments/lattice_avoidance.py        # Result 8
REW_N=20 REW_TAG=_n20 python3 experiments/lattice_reward_edges.py     # Result 4
AV_N=20  AV_TAG=_n20  AV_D=20 python3 experiments/lattice_attention_value.py   # Result 5
LY_N=20  LY_TAG=_n20  LY_D=50 python3 experiments/lattice_layers.py   # Result 6
LY_N=20  LY_TAG=_n20_d70 LY_D=70 python3 experiments/lattice_layers.py
ID_N=20 ID_TAG=_n20_move ID_ASTAR=0.08 ID_R=0,30,300 python3 experiments/lattice_identity.py
ID_N=20 ID_TAG=_n20_keep ID_ASTAR=0.16 ID_R=0,30,300 python3 experiments/lattice_identity.py
TO_N=5  TO_TAG=_n5 python3 experiments/lattice_tonic.py               # Negative 2

# figures (prefers _n20 artifacts, falls back to n=3)
python3 experiments/lattice_arc_figures.py
```

![Transmission edge](figures/lattice_transmission_edge.png)
![Tonic trade](figures/lattice_tonic_trade.png)
![Identity bind](figures/lattice_identity_bind.png)
