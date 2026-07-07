# Learning on a Greenberg–Hastings Network: Proposed Experiments

*A staged experimental program for a reward-driven learning mechanism built on
Greenberg–Hastings (GH) excitable dynamics running on a graph. The framing is
inside-out (Buzsáki): the substrate generates its own repertoire of dynamical
patterns first, and a strict scalar reward (Sutton) selects and stabilises a
subset of them through action. Memory, attention and executive function are not
built-in modules — they are read out, after the fact, as different questions
asked of one homogeneous substrate.*

---

## 0. Status and scope

This document specifies **what to build and what to measure**, not the
implementation. It fixes the substrate, the learning framework, the input /
cue / feedback formats, the hyperparameters and their defaults, and a staged
series of experiments (E0–E6) with well-defined target behaviours and success
criteria. Each experiment names the hypothesis it tests and the prediction that
would falsify it.

The existing code (`ghca_main.py`) already provides the excitable substrate on
a 2D lattice; `ghca_plot.py::plot_param_space` already computes a
**persistence-probability map** over `(active, passive)` space — that map is,
in effect, a first draft of the memory-capacity phase diagram this program
depends on (E0). The work below is a superset of that substrate: generalise the
neighbourhood to a graph, make selected parameters plastic, and close an
action–perception loop with a reward signal.

---

## 1. Substrate specification

### 1.1 Node state

Each node `i` holds an integer phase `φ_i` on a cyclic clock, matching the
current `Population` encoding:

| Phase range        | Label            | Meaning                                    |
|--------------------|------------------|--------------------------------------------|
| `0`                | rest (S)         | quiescent, excitable                       |
| `1 … a`            | active (I)       | excited / infectious; can excite neighbours|
| `a+1 … a+p_i`      | refractory (R)   | cannot be excited; counts down             |

`a` = `active` (excited duration, **fixed** across nodes), `p_i` = `passive`
(refractory duration, **per-node**), and the local timescale is
`τ_i = a + p_i`. Advancing past `τ_i` wraps to `0` (this is the
`if p >= tau0: p = 0` / `if 1 <= p < tau0: p += 1` logic already in
`infect()`).

Refractoriness is the substrate's **native inhibition**: it enforces wave
directionality, causes colliding waves to annihilate (a built-in winner-take-all
— see E4), and prevents immediate re-excitation. It does **not** regulate
recruitment of *rested* cells; that job is given to the homeostatic threshold
(§1.4).

### 1.2 Topology (the graph)

Generalise `Population.nbr(i,j)` from lattice offsets to an explicit adjacency
structure. Neighbours of `i` are `N(i)`, with directed weights `w_{ij} ≥ 0`
for `j ∈ N(i)`.

Supported topologies (chosen per experiment):

- **`ring`** — 1D ring, degree `2k` (nearest `k` each side). Cleanest for
  sequence / timing experiments.
- **`lattice2d`** — the current torus (periodic 2D), range `ρ`. Keeps
  continuity with existing results.
- **`smallworld`** — Watts–Strogatz `(N, k, β_rewire)`. Introduces shortcuts →
  richer reentrant loop structure.
- **`rgg`** — random geometric graph (radius `r_geo`). Biologically flavoured
  distance-dependent connectivity.

### 1.3 Excitation rule (with plastic conduction)

A rested node `i` (`φ_i = 0`) becomes active (`φ_i ← 1`) at the next step iff

```
Σ_{j ∈ N(i), φ_j ∈ [1, a]}  w_{ij}   ≥   θ_i
```

i.e. the weighted count of *active* neighbours meets the node's threshold. This
generalises the current rule, which is the special case `w_{ij} = 1`,
`θ_i = 1` ("any active neighbour excites"). Active/refractory nodes advance
deterministically regardless of neighbours (the GH clock).

### 1.4 Spontaneous firing and homeostasis

- **Spontaneous firing.** Each rested node fires spontaneously with probability
  `p_s` per step (independent). This provides self-generated activity
  (inside-out), doubles as the RL exploration source (no separate ε-greedy —
  `p_s` is effectively the exploration temperature), and keeps the medium in the
  cyclic-CA "soup" regime that nucleates spiral cores without seeding.
- **Homeostatic threshold.** `θ_i` adapts slowly toward a target activity level
  to keep the network in the self-sustaining band and off both the runaway
  ("seizure") and death edges:

  ```
  θ_i ← θ_i + η_θ · ( ρ̄_i − ρ* )
  ```

  where `ρ̄_i` is a running estimate of local activity fraction and `ρ*` the
  target. This is the gain regulator refractoriness cannot provide.

### 1.5 Operating regime (why it must be tuned)

Following Fisch–Gravner–Griffeath threshold-range scaling, an excitable CA of
this form sits in one of three regimes as the threshold-to-neighbourhood ratio
increases: **turbulent/debris** (too excitable) → **self-sustaining spiral
waves** (the useful band) → **death/fixation** (too inhibited). The learning
experiments assume the substrate is homeostatically parked in the spiral band;
**E0 exists to find that band** for each topology before any learning is run.

---

## 2. Learning framework

### 2.1 Strict reward (Sutton)

A single scalar reward `r_t` is emitted **only** by the environment. The agent
maximises the expected discounted return `E[ Σ_k γ^k r_{t+k} ]`. No shaping is
permitted **except** potential-based shaping `F(s,s') = γΦ(s') − Φ(s)`
(Ng–Harada–Russell — the unique form that preserves the optimal policy). Since a
natural choice is `Φ = V`, admissible shaping coincides with the TD signal the
agent already computes; nothing else is fed in.

### 2.2 Critic = order parameter of the medium (with timescale separation)

The value estimate `V_t` is **read from the medium's collective state**, not a
separate weight vector. Candidate order parameters (choose/compare in E0):

- `A_t` — global active fraction (sustained-activity level).
- `C_t` — phase coherence (e.g. Kuramoto-style order parameter over node phases
  `φ_i / τ_i`).
- `S_t` — dominant-loop stability (autocorrelation peak height of `A_t`).

To retain actor–critic stability without a structurally separate critic, run
value estimation on a **faster timescale** than plasticity (two-timescale AC):
the readout `V_t = v · m_t` (linear on the order-parameter vector `m_t`) adapts
with rate `α_v ≫ η_w, η_τ`. Equivalently, implement the critic as a slow GVF
demon (§2.5).

### 2.3 Neuromodulator = TD error

The global learning signal broadcast to every plastic element is the TD error

```
δ_t = r_t + γ V_{t+1} − V_t
```

computed internally. This is the dopamine-as-TD-error identity; it is **not** the
reward.

### 2.4 The two parallel plasticity lines

Both lines share the scaffold: strict reward → internal `δ_t` → eligibility
trace → local update. They differ only in the learned variable and its
eligibility, which lets a single task discriminate them (E3).

- **Line A — Conduction plasticity.** Learn edge weights `w_{ij}` (and let
  `θ_i` float homeostatically). Edge eligibility charges on causal conduction
  (`j` active immediately before `i` fires) and decays:

  ```
  e_{ij} ← λ e_{ij} + 1{ φ_j ∈ [1,a] at t−1  ∧  φ_i : 0→1 at t }
  Δw_{ij} = η_w · δ_t · e_{ij}          (w_{ij} clipped to [0, w_max])
  ```

  Assigns **spatial** credit — *which paths conduct* (routing / transition
  graph).

- **Line B — Timescale plasticity.** Learn per-node refractory length `p_i`
  (hence `τ_i`); weights held fixed/uniform. Eligibility is a period trace: when
  node `i` re-excites after interval `Δt_i`, drive `τ_i` toward the interval of
  *rewarded* loops:

  ```
  Δτ_i = η_τ · δ_t · ( Δt_i − τ_i )     (τ_i clipped to [τ_min, τ_max])
  ```

  Assigns **temporal** credit — *when / at what tempo*; self-organises a
  hierarchy of loop periods.

- **Line A+B — Combined.** Both updates active. Expected to be necessary for
  spatiotemporal sequence tasks. E3 tests whether A and B **compose** and
  whether B-first is a useful curriculum.

Every experiment is run as an **ablation matrix**:
`{A, B, A+B} × {homeostasis on/off} × {p_s on/off}`.

### 2.5 Emergent categories = GVF demons (Horde)

For readout experiments (E6), attach many **General Value Functions**, each a
tuple `⟨ target policy π, cumulant c, continuation γ_c ⟩`, learned **off-policy
via gradient-TD** from the single behaviour stream, in constant time/step. The
claim "memory / attention / executive are emergent categories" is operationalised
as: distinct GVF demons, riding one unchanged substrate, answer these distinct
predictive questions well — no dedicated module per function.

Executive control specifically is framed with **options** `⟨I, π, β⟩`: a slow,
long-`τ` loop (a spiral core) whose initiation set `I` is the context in which it
ignites, whose policy `π` is the fast dynamics it gates, and whose termination
`β` is its decay — an option learned over the medium's own loops.

---

## 3. Common experimental scaffold

### 3.1 Node roles

Partition nodes into disjoint sets:

- **Sensory** `V_S` — cue is injected here.
- **Motor** `V_M` — grouped into `A` action channels `V_M^{(1)} … V_M^{(A)}`.
- **Interneuron / hidden** `V_H` — the bulk; where loops live.
- **Context** `V_C` (E5 only) — receives the task/rule signal.

### 3.2 Input format (cue)

A cue is a **spatial pattern of forced excitation** on sensory nodes for a fixed
onset window:

- **Identity coding.** Stimulus `x ∈ {1…K}` maps to a population code over
  `V_S`: a designated subset `V_S^{(x)}` is clamped to `φ = 1` for
  `Δ_cue` steps (population/one-hot; overlapping codes optional for
  generalisation tests).
- **Temporal coding.** A sequence `(x_1, x_2, …)` is delivered as successive cue
  clamps separated by inter-stimulus intervals `ISI`.
- **Go signal.** A distinct sensory node `g ∈ V_S` (or a brief global drop in
  `θ`) marks when a response is required.

Injection is additive over ongoing spontaneous dynamics — cues **perturb**
self-generated activity rather than overwrite it (action-driven perception).

### 3.3 Output format (action readout)

Over a **response window** `W` steps after the go signal, integrate activity per
motor channel:

```
score(m) = Σ_{t ∈ W} Σ_{i ∈ V_M^{(m)}} 1{ φ_i ∈ [1,a] }
action  = argmax_m score(m)          (ties / sub-threshold ⇒ "no response")
```

A `first-to-fire` variant (channel whose first spike in `W` is earliest) is used
where response *latency*/timing is the dependent variable (E3).

### 3.4 Feedback (reward)

Strict scalar, sparse, delivered at the end of `W`:

```
r = +1   if action == target action
r =  0   otherwise            (optionally r = −c_err for a wrong commit)
```

`r` drives `δ_t` (§2.3), which gates plasticity. No per-node teaching signal is
ever supplied.

### 3.5 Episode structure

```
[ settle: S steps of free dynamics ]
[ cue onset: Δ_cue ]
[ delay: D (0 in E1) ]
[ go signal ]
[ response window: W → evaluate → reward → δ broadcast ]
[ inter-trial interval: ITI ]
```

Trials are drawn i.i.d. from the task distribution; learning curves are reported
over trials.

### 3.6 Metrics

- **Accuracy** — fraction correct, vs trial (learning curve) and vs task
  parameter (e.g. delay `D`, ISI).
- **Latency** — steps from go to correct motor spike.
- **Timing error** (E3) — |produced interval − target interval|.
- **Loop persistence** — does a reentrant loop, detected by recurrence /
  autocorrelation of the activity vector, survive cue removal through the delay?
- **Order parameter trace** — `A_t, C_t, S_t` as the value proxy; correlation of
  `V_t` with realised return (validates the critic).
- **Regime location** — active fraction, spiral-core count, dominant period
  (FFT of `A_t`), to confirm the substrate stays in the spiral band.
- **Sample efficiency** — trials to reach a criterion (e.g. 90% over 100 trials).

---

## 4. Hyperparameters and defaults

Defaults are starting points to be calibrated in E0; ranges indicate sweep
bounds.

### 4.1 Substrate

| Symbol | Meaning | Default | Sweep range |
|--------|---------|---------|-------------|
| `N` | number of nodes | 400 | 100 – 2000 |
| topology | graph type | `smallworld` | ring / lattice2d / smallworld / rgg |
| `k` | base degree (per side / mean) | 6 | 2 – 20 |
| `β_rewire` | small-world rewiring prob | 0.1 | 0 – 0.3 |
| `a` | active (excited) duration | 2 | 1 – 6 |
| `p` (init `p_i`) | refractory duration | 8 | 1 – 17 |
| `τ` (`= a+p`) | local timescale | 10 | 2 – 23 |
| `w_max` | max edge weight | 4.0 | — |
| `θ` (init) | excitation threshold | 2.0 | 1 – k |
| `p_s` | spontaneous firing prob | 5e-3 | 0 – 5e-2 |
| `ρ*` | homeostatic target activity | 0.10 | 0.02 – 0.25 |
| `η_θ` | threshold adaptation rate | 1e-3 | 1e-4 – 1e-2 |

### 4.2 Learning

| Symbol | Meaning | Default | Sweep range |
|--------|---------|---------|-------------|
| `γ` | discount | 0.95 | 0.8 – 0.99 |
| `λ` | eligibility decay | 0.9 | 0.5 – 0.98 |
| `η_w` | conduction learning rate (Line A) | 1e-2 | 1e-3 – 1e-1 |
| `η_τ` | timescale learning rate (Line B) | 5e-2 (in Δτ units) | 1e-2 – 2e-1 |
| `α_v` | critic (order-param readout) rate | 1e-1 | ≫ η_w, η_τ |
| `τ_min, τ_max` | timescale bounds (Line B) | 2, 23 | — |
| `c_err` | wrong-commit penalty | 0 | 0 – 1 |

### 4.3 Task / episode

| Symbol | Meaning | Default | Sweep range |
|--------|---------|---------|-------------|
| `K` | number of distinct stimuli | 2 | 2 – 8 |
| `A` | number of action channels | 2 | 2 – 8 |
| `|V_S^{(x)}|` | sensory nodes per stimulus | 8 | 4 – 32 |
| `|V_M^{(m)}|` | motor nodes per channel | 8 | 4 – 32 |
| `Δ_cue` | cue duration | 3 | 1 – 10 |
| `D` | delay (E2+) | 0 | 0 – 100 |
| `ISI` | inter-stimulus interval (E3) | 6 | 2 – 20 |
| `W` | response window | 5 | 2 – 15 |
| `S` | settle steps | 50 | — |
| `ITI` | inter-trial interval | 20 | — |

---

## 5. Experiment series

Each experiment lists: **purpose**, **setup delta** (relative to the scaffold),
**target action**, **cue/feedback**, **metrics**, **hypothesis → prediction**,
and **discriminator**.

### E0 — Substrate characterisation (no learning)

- **Purpose.** Locate the self-sustaining spiral band and choose the critic's
  order parameter, per topology. Establishes the operating point every later
  experiment assumes.
- **Setup delta.** No plasticity, no reward, no cues. Sweep `(θ, p_s, ρ*)` and
  the timescale `τ`; for lattice2d also sweep range `ρ`.
- **Target action.** None (characterisation only).
- **Metrics.** Sustained active fraction; spiral-core count; dominant period
  (FFT of `A_t`) vs `τ`; lifetime of self-sustained activity after transient
  seeding; the three candidate order parameters `A_t, C_t, S_t` and their
  stationarity.
- **Hypothesis → prediction.** There is an intermediate `θ/⟨deg⟩` band with
  persistent, coherent, non-saturating activity; below it → turbulence/seizure,
  above it → death. Dominant loop period tracks `τ` (validating Line B's control
  variable). Reuse/extend `plot_param_space` to render the persistence map on
  the graph substrate.
- **Discriminator.** If no band supports both persistence *and*
  non-saturation, revisit topology/degree before proceeding.
- **Status: DONE.** See [`e0_results.md`](e0_results.md). Key outcomes: range-1
  (von Neumann) fixates and cannot self-sustain — the discriminator fired and
  moved the operating point to r≥2; the live threshold band widens with range
  (threshold-range scaling); an organised intermediate band sits at
  **r=2, a=6, τ=14, θ≈4**; and the dominant loop period tracks τ almost exactly
  (`period = 1.00·τ + 0.95`, r = 0.9992), validating Line B's control variable.
  Implemented in `ghca_net.py` + `experiments/e0_characterization.py`.

### E1 — Stimulus→response conditioning (identity only, no delay)

- **Purpose.** Simplest closed loop: can strict reward carve a stimulus→action
  mapping into the substrate at all?
- **Setup delta.** `K = A = 2`, `D = 0`. Go signal immediately after cue.
- **Target action.** Stimulus `x` → fire motor channel `m = x` within `W`.
- **Cue/feedback.** Population code on `V_S^{(x)}`; `r = +1` iff `argmax`
  channel `= x`.
- **Metrics.** Accuracy learning curve; trials-to-criterion; correlation of
  `V_t` with return (critic sanity check).
- **Hypothesis → prediction.** Line A alone suffices (a conduction path
  `V_S^{(x)} → V_M^{(x)}` strengthens). Line B alone should be **near chance**,
  because identity mapping needs routing, not timing.
- **Discriminator.** A ≫ B here. If B alone succeeds, the readout is leaking
  spatial information through timing — investigate.
- **Status: DONE.** See [`e1_results.md`](e1_results.md). Over 6 seeds:
  **A = 0.91, B = 0.35 (≤ chance), A+B = 0.86** final accuracy. Line A learns
  the mapping (chance → ~0.91); Line B cannot route identity and stays at/below
  chance; the discriminator passes with no timing leak. Implemented in
  `ghca_learn.py` + `experiments/e1_conditioning.py`. Note: required
  stimulus-selective (channel-biased) hidden wiring — letting the net discover
  selective representations is deferred.

### E2 — Delayed response (working memory)

- **Purpose.** Test memory as a persistent reentrant loop that bridges a gap.
- **Setup delta.** Cue removed after `Δ_cue`; go signal after delay `D > 0`.
  Sweep `D`.
- **Target action.** Reproduce the E1 mapping, but the correct channel must fire
  after the (stimulus-free) delay.
- **Cue/feedback.** As E1; reward at end of `W` after the go.
- **Metrics.** Accuracy vs `D`; **loop-persistence** metric during the delay;
  fraction of trials where a stimulus-specific loop survives cue removal;
  robustness to partial cue (pattern completion) at test.
- **Hypothesis → prediction.** Accuracy holds up to a delay `D_max` set by loop
  stability. **Line B (timescale) is critical**: long-`τ` loops are needed to
  hold activity across `D`; A-only should show accuracy decaying fast with `D`.
- **Discriminator.** B (or A+B) sustains long `D`; A-only degrades. Persistent,
  cue-specific loop = the memory engram.
- **Status: DONE.** See [`e2_results.md`](e2_results.md). The dissociation
  **inverts E1**: with `τ` fixed above the loop transit time, **Line A retains
  only at D=0** (loop dies, retention `[1,0,0,0,0]`), while **Line B learns
  `τ` below `L` (26→~13) and holds memory to D=200** (`[1,1,1,1,1]`); A+B
  matches B. Mechanism sweep confirms memory duration is `τ`-controlled
  (`τ<L` sustains, `τ≥L` dies). Two rule findings: the resonance rule targets
  the death boundary (`τ=L`) so B uses reward-gated perturbation; and per-node
  `τ` hits a weakest-link problem (loop dies at the slowest node), so a **shared
  regional timescale** is what learns robustly. Implemented in `ghca_learn.py`
  (shared-`τ` mode) + `experiments/e2_delayed_response.py`.

### E3 — Temporal sequence reproduction (the A-vs-B discriminator)

- **Purpose.** The central discriminator: separate spatial from temporal credit
  assignment.
- **Setup delta.** Cue = ordered sequence `A → B → C` on distinct sensory
  subsets with fixed `ISI`; after a go, the agent must reproduce the sequence on
  motor channels with correct **order and timing** (or predict the next element).
- **Target action.** Emit motor spikes in order `A,B,C` with inter-emission
  intervals matching the presented `ISI` (within tolerance).
- **Cue/feedback.** Temporal coding (§3.2). Reward decomposed for analysis but
  delivered as one scalar: correct only if both order and timing tolerances met.
- **Metrics.** Order accuracy; timing error; joint success rate; effect of a
  **B-first curriculum** (train timescales, then weights).
- **Hypothesis → prediction.**
  - Line A alone → **order correct, timing distorted** (intervals dictated by
    frozen `τ`).
  - Line B alone → **timing/tempo correct, order/identity confused**.
  - A+B → both; and B-first is expected to speed convergence (clock scaffolds
    routing).
- **Discriminator.** The double dissociation above is the key predicted result of
  the whole program. Its presence/absence validates or refutes the
  two-line decomposition.

### E4 — Selective attention (cue competition)

- **Purpose.** Show attention as biased winner-take-all among competing waves —
  using the refractory annihilation the substrate already has.
- **Setup delta.** Two simultaneous cue streams on disjoint sensory sets, each
  implying a different action. A **top-down bias** (small `θ` reduction, or a
  slow biasing wave from `V_C`) favours one stream.
- **Target action.** Fire the channel implied by the **attended** (biased)
  stream.
- **Cue/feedback.** Two concurrent population codes + a bias vector; `r = +1`
  iff attended-stream action produced.
- **Metrics.** Routing accuracy vs bias magnitude; psychometric curve; margin of
  the winning wave (annihilation locus).
- **Hypothesis → prediction.** Even small bias reliably routes the correct
  stream, because colliding waves annihilate and the biased front wins the race
  to `V_M`. Attention needs **no new inhibitory machinery**.
- **Discriminator.** If accurate routing requires an added inhibitory population,
  the "refractoriness = sufficient competition" claim fails.

### E5 — Executive control / task switching (options)

- **Purpose.** Show a slow loop acting as an option that gates fast
  stimulus→response mappings; test rule reversal.
- **Setup delta.** A context signal on `V_C` selects which of two mappings
  (`x→x` vs `x→¬x`) is rewarded. Context switches in blocks. Encourage a slow
  timescale sub-population (`τ` at the high end) via homeostasis.
- **Target action.** Apply the currently-cued rule to map stimulus → action.
- **Cue/feedback.** Stimulus cue + context cue; reward per current rule.
- **Metrics.** Reversal learning speed; **switch cost** (accuracy dip after a
  switch and recovery time); evidence that a slow-`τ` loop's state predicts the
  active rule.
- **Hypothesis → prediction.** A slow loop self-organises into a context
  variable (an option's initiation/termination structure) that reconfigures fast
  routing; switch cost decreases across blocks as options consolidate.
- **Discriminator.** Ablating the slow (`high-τ`) sub-population abolishes
  flexible switching while leaving single-rule performance intact.

### E6 — Emergent categories (Horde readout)

- **Purpose.** Operationalise "memory / attention / executive are emergent
  categories, not modules."
- **Setup delta.** Freeze a substrate trained through E2–E5. Attach GVF demons
  (off-policy gradient-TD) reading the same node-phase stream, each with a
  different cumulant/policy/continuation:
  - *memory demon*: predict persistence of the current stimulus-specific loop;
  - *attention demon*: predict which stream will capture `V_M`;
  - *executive demon*: predict an imminent context switch.
- **Target action.** None new — demons predict, they do not control (readout
  only).
- **Metrics.** Each demon's prediction accuracy / calibration; representational
  similarity between demons (are they genuinely distinct questions?); do demons
  succeed **without** any function-specific wiring?
- **Hypothesis → prediction.** All three demons achieve well-above-chance
  prediction from one unchanged substrate; the "functions" live in the probes,
  not the machine.
- **Discriminator.** If a demon can only succeed given a dedicated sub-network,
  the strong emergence claim is weakened to a modular one.

---

## 6. Analyses that cut across experiments

- **Ablation matrix** (§2.4) reported for E1–E5: isolates the contribution of
  each plasticity line, homeostasis, and spontaneous firing.
- **Critic validation.** Across all learning experiments, regress realised
  return on the order-parameter value `V_t`; a good critic shows rising
  correlation over training. Compare `A_t` vs `C_t` vs `S_t`.
- **Regime tracking.** Log active fraction / dominant period every trial to
  confirm the substrate never drifts out of the spiral band (homeostasis
  working).
- **Timescale hierarchy.** In A+B runs, inspect the learned `τ_i` distribution
  for emergence of a fast/slow hierarchy and cross-frequency coupling between
  loops of different period.

---

## 7. Implementation notes (mapping to existing code)

Minimal changes to reach E0–E1; the rest are additive.

1. **Graph neighbourhood.** Replace the offset-based `default_nbh` / `nbr()` in
   `ghca_main.py` with an adjacency list + weight lookup. Keep a `lattice2d`
   adapter so existing results and `plot_param_space` still run.
2. **Per-node timescale.** Promote `tau0` to an array `τ_i` (and `p_i`); the
   `infect()` wrap condition becomes per-node (`φ_i >= τ_i → 0`).
3. **Weighted threshold.** `check()` computes the weighted active-neighbour sum
   against `θ_i` instead of "any active neighbour".
4. **Spontaneous firing + homeostasis.** In `run()`, after `infect()`, fire
   rested nodes with prob `p_s` and update `θ_i` toward `ρ*`.
5. **Plasticity module.** Eligibility traces on edges (Line A) and per-node
   period traces (Line B); update on the broadcast `δ_t`.
6. **Environment loop.** Sensory/motor/context node sets, cue injection, action
   readout, reward, and the order-parameter critic — wrap the existing `run()`
   step inside an episode driver.

Keep the current animation/plot utilities for qualitative inspection of loops
and waves during debugging.

---

## 8. Reference anchors

- **Greenberg & Hastings** — excitable-media CA (the substrate).
- **Fisch, Gravner & Griffeath** — *Threshold-range scaling of excitable
  cellular automata* and *Cyclic Cellular Automata in Two Dimensions* (regimes:
  debris → spirals → fixation; the operating-point map for E0).
- **Sutton (1988); Barto, Sutton & Anderson (1983)** — TD learning and
  actor–critic (the learning scaffold, the two-timescale critic).
- **Sutton, Precup & Singh (1999)** — options / temporal abstraction (E5,
  executive control).
- **Sutton, Modayil et al. (2011)** — Horde / GVFs (E6, emergent categories).
- **Ng, Harada & Russell (1999)** — potential-based shaping (the strict-reward
  boundary, §2.1).
- **Buzsáki, *The Brain from the Inside Out*** — action-first, self-generated
  activity, representations as emergent rather than imposed (the framing).
