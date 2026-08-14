# ghca

A study of local timescales in Greenberg–Hastings cellular automata, and an
exploratory program using GH excitable dynamics as the substrate for a
reward-driven **learning mechanism** on a graph.

## Two strands

**1. Lattice GH cellular automata (original).** Excitable-media dynamics on a 2D
lattice with an `(active, passive)` refractory cycle, plus tools to enumerate
which initial configurations self-sustain ("persistence probability" over
`(active, passive)` parameter space).

**2. Learning on a GH network (in progress).** GH dynamics generalised to an
arbitrary weighted graph, with per-node timescales, as the substrate for a
learning mechanism. The framing is inside-out (Buzsáki): the medium generates
its own repertoire of dynamical patterns, and a strict scalar reward (Sutton)
selects and stabilises a subset through action. Memory, attention and executive
function are treated as categories *read out* from one homogeneous substrate,
not as built-in modules.

## Files

| File | Purpose |
|------|---------|
| `ghca_main.py` | Lattice GH substrate: `Population` class, `run`, `plot`, `animate` |
| `ghca_core.py` | Encode/embed integer configurations, run + animate |
| `ghca_plot.py` | Persistence-probability maps over `(active, passive)` space |
| `ghca_net.py` | **GH dynamics on a graph**: per-node timescales, weighted-threshold excitation, spontaneous firing, homeostatic threshold; topology builders and order-parameter observables |
| `ghca_learn.py` | **Reward-modulated learner**: eligibility-trace conduction (Line A) and timescale (Line B) plasticity, order-parameter critic, layered-graph builder |
| `ghca_plasticity.py` | **Multi-axis closed-loop engine**: tri-axis ($\tau$-adaptation, $\theta$-homeostasis, $W$-routing) substrate adaptation |
| `ghca_causal.py` | **Causal instrumentation** (C-series): partial-observation `S_obs`, wave variables `W=f(S)`, and `do(S)` / `do(W)` / `do(θ)` intervention operators |
| `experiments/e0_characterization.py` | E0 — substrate characterisation (find the self-sustaining band) |
| `experiments/e1_conditioning.py` | E1 — stimulus→response conditioning (A-vs-B dissociation) |
| `experiments/e2_delayed_response.py` | E2 — delayed response / working memory (τ-controlled memory) |
| `experiments/e2_information.py` | E2 addendum — memory as a τ-tuned information-destruction rate |
| `experiments/topology_cycle_capacity.py` | E2 addendum — cycle-space (circuit-rank) bound on how many reentrant loops a topology admits |
| `experiments/topology_cycle_packing_exact.py` | **Correction** to the above — its greedy maximises cycle *length*, not *count*; exact ILP packing, certified on the ring |
| `experiments/persistent_set_structure.py` | Does the persistent set have a tractable description? Not linear; combinatorial when `tau_p <= tau_a` (companion to SSRN 4047679) |
| `experiments/phase_space_polytope.py` | Addendum: exact hulls in direct phase-space coords — small saturating outer polytope, necessary-not-sufficient; (1,1) conservation law |
| `experiments/persistent_set_3x3.py` | The gap signature generalises to 3x3 as a plaquette multiset — 40M configs, zero impure at `tau_p <= tau_a`, 2096x compression |
| `experiments/persistent_set_dec.py` | Discrete-calculus ladder: curl (winding) insufficient, ALL relative-phase invariants insufficient at `tau_p > tau_a`; threshold = + absolute state multiset |
| `experiments/spectrum_mechanism_hunt.py` | Theory hunt: 3 mechanism candidates falsified (bridging, lumpability x2), 3 exact laws found (period law `T = k*S + d`, spectrum constancy on attractors, dwell = f(spectrum)) |
| `experiments/topology_winding_capacity.py` | E2 addendum — the GGH (1980) winding number as the exact sustain criterion; calibrates the length gate above |
| `experiments/scaling_capacities.py` | Scaling (Track 3b, size half) — does substrate size buy memory (E2) / attention (E4) / executive control (E5)? |
| `experiments/lattice_capacities.py` | Representation + `lattice2d` port of all three capacity mechanisms; shows E5's hidden layer has **0** recurrent edges |
| `experiments/lattice_animation.py` | GIFs for the three lattice mechanisms (reentry, wave annihilation, held option) |
| `experiments/lattice_timescale_demo.py` | Ports the input-timing `τ` rule to a lattice; two results, two negatives (see [`docs/lattice_timescale_notes.md`](docs/lattice_timescale_notes.md)) |
| `experiments/lattice_afferent_timing.py` | Resolves the lock-in: `τ` learning needs a privileged **afferent** channel |
| `experiments/lattice_afferent_depth.py` | How far exogenous timing penetrates a recurrent medium (answer: it does not) |
| `experiments/lattice_attention_gate.py` | A 1-D attention strip of the same cells gating plasticity — a clock, not a filter |
| `experiments/lattice_reward_edges.py` | Reward as a fourth edge; `τ` encodes a stimulus–reward interval |
| `experiments/lattice_attention_value.py` | A value chain teaches the attention strip *where* to gate |
| `experiments/lattice_layers.py` | 2-D layers instead of 1-D edges; a synchronous burst timed to reward |
| `experiments/lattice_sensorimotor.py` | The action is **transmission**, not emission; bootstrapping is free |
| `experiments/lattice_avoidance.py` | Avoidance is not sign-symmetric with approach (transmission is provably monotone) |
| `experiments/lattice_identity.py` | Plastic cell *identity* via homeostatic θ — a structural negative |
| `experiments/lattice_tonic.py` | Tonic drive: no window; neither stalled thread unlocks |
| `experiments/e3_timed_response.py` | E3 — timed response (identity × latency double dissociation) |
| `experiments/e3_factored_credit.py` | E3 composition study — factored credit + curriculum vs shared reward |
| `experiments/e4_attention.py` | E4 — selective attention as biased WTA by wave annihilation |
| `experiments/e5_executive.py` | E5 — executive control / task switching (a slow-loop option gates fast routing) |
| `experiments/e6_horde.py` | E6 — emergent categories (three GVF demons read one frozen substrate) |
| `experiments/c0_instrumentation.py` | C0 — instrument the causal variables (`W=f(S)`, partial spikes) |
| `experiments/c1_graph_certificates.py` | C1 — validate Theorem-1 epiphenomenality certificate on known SCMs |
| `experiments/c2_fat_handed.py` | C2 — `do(W)` is fat-handed when `W=f(S)` (achievable-band of behaviour) |
| `experiments/c3_do_theta.py` | C3 — `do(θ)` (timescales/couplings) is the well-posed causal handle |
| `experiments/c4_outcome_relativity.py` | C4 — outcome-relativity & degeneracy (causal-emergence cap) |
| `experiments/closed_loop_plasticity.py` | **Closed-loop plasticity**: single-task benchmark evaluating Readout Independence Ratio (RIR) |
| `experiments/sequential_closed_loop.py` | **Sequential learning**: Task A $\to$ Task B $\to$ Task A reversal learning evaluating anti-forgetting |
| `result/` | Saved simulation outputs (`.npy`/`.npz`) and experiment data |

## Documentation

- [`docs/learning_experiments.md`](docs/learning_experiments.md) — the full
  design: substrate spec, strict-reward learning framework, the two parallel
  plasticity lines (conduction weights vs local timescales), input/cue/feedback
  formats, hyperparameters, and the staged experiment series **E0–E6**.
- [`docs/e0_results.md`](docs/e0_results.md) — **findings from E0** (substrate
  characterisation): range-1 fixates, the live threshold band widens with range
  (threshold-range scaling), an organised spiral band at r=2/a=6/θ≈4, and the
  dominant loop period tracking τ (`period = 1.00·τ + 0.95`, r = 0.9992).
- [`docs/e1_results.md`](docs/e1_results.md) — **findings from E1**
  (conditioning): a strict scalar reward carves the stimulus→action mapping;
  the predicted dissociation holds (Line A = 0.91, Line B = 0.35 ≤ chance,
  A+B = 0.86 final accuracy over 6 seeds).
- [`docs/e2_results.md`](docs/e2_results.md) — **findings from E2** (working
  memory): memory is a τ-controlled reentrant loop; the dissociation inverts —
  Line A retains only at zero delay, Line B learns τ below the loop transit time
  and holds memory to D=200. Needs a *shared* regional timescale (per-node τ
  hits a weakest-link problem).
- [`docs/topology_cycle_capacity.md`](docs/topology_cycle_capacity.md) —
  **E2 addendum** (topology analysis, no learning): the circuit rank
  `β₁ = m − N + c` bounds how many independent reentrant loops a substrate
  admits, and the E2 sustain gate (`τ < L`) collapses that ceiling by 20–70× to
  the usable count — which *falls* with τ on lattice/small-world/RGG (the ring is
  a flat exception: its *packed* cycles all exceed the τ range). A
  capacity/duration tradeoff on any Line-B τ policy.
  Capacity is a bound, not a measured count.
  **⚠ Superseded in part — see `topology_cycle_packing_exact.md`: the greedy
  packing maximises cycle length rather than count, so every `K_dyn` there is low
  by 1.4–7.5×, the 20–70× collapse is really 3.9–20×, and the ring is not flat.**
- [`docs/topology_cycle_packing_exact.md`](docs/topology_cycle_packing_exact.md) —
  **correction to the above.** `pack_long_cycles` picks `max(longer, key=len)`;
  maximising the *count* wants `min`. One word costs 1.0–4.7×, an exact set-packing
  ILP up to 7.5×. On `ring(60,k=3)` the optimum is **certified** at all four τ
  tested (45/36/30/25, meeting `min(β₁, ⌊m/(τ+1)⌋)`) versus 6 in the merged doc —
  so ring capacity *does* fall with τ. Best-known K is still monotone decreasing in
  τ on all four topologies, so the capacity/duration tradeoff survives; only the
  numbers and the "flat ring" claim change.
- [`docs/persistent_set_structure.md`](docs/persistent_set_structure.md) —
  **new, external-facing** (structural companion to Moitra & Sen, SSRN 4047679).
  Does the persistent set admit a tractable description? **Geometrically no** — it
  is not linearly separable from the dying set in one-hot coordinates in any of 22
  cells with `P > 0`, so no polytope in those coordinates describes it.
  **Combinatorially yes** whenever `τ_p ≤ τ_a`: persistence is a function of the
  sorted cyclic phase-gap multiset, compressing the space 16–161× into homogeneous
  classes (14641 configs → 91 classes at (5,5)), with compression growing in `S`.
  One-way implication only — `(2,5)` is invariant with `τ_p > τ_a`. Dynamics
  validated against `ghca_net.Network`: 7461 transitions, 0 mismatches.
- [`docs/topology_winding_capacity.md`](docs/topology_winding_capacity.md) —
  **E2 addendum** (calibration of the above): the winding number of a continuous
  cycle, the invariant from Greenberg, Greene & Hastings (1980), is an *exact*
  sustain criterion (45/45) where the length gate scores 40/45 — and every
  length-gate miss is `τ = L`, the marginal death boundary. Winding is a
  *dynamical* quantity, so it calibrates the topological bound rather than
  replacing it; read too early it underperforms the gate.
- [`docs/scaling_capacities.md`](docs/scaling_capacities.md) — **scaling**
  (Track 3b, size half): size helps only where size is the binding constraint.
  E2 memory is size-limited and matches the `τ < L` law 6/6; E4 attention is
  **scale-invariant** (0.08 sensitivity spread over an 8× arena range —
  bias-to-noise sets the collision locus, not room); E5 executive control shows
  **no established size effect** (16× `N_H`, Spearman p=0.51, overlapping CIs,
  per-seed bimodality dominating). A null on the size axis, next to 3c/P4's
  positive result for interference — different capacity, different constraint.
- [`docs/lattice_capacities.md`](docs/lattice_capacities.md) —
  **representation + 2-D port** (with animations): each capacity uses a
  *different* substrate primitive (reentrant loop / colliding waves / held option
  gating a feedforward conjunction), and E5's context ring is E2's loop reused.
  E5's hidden layer has **0 H→H edges** — purely feedforward AND-gates, which is
  *why* the `N_H` null holds. E2's `τ < L` law ports to the 2-D torus **28/28**;
  E4's psychometric curve ports unchanged and stays size-invariant. E5 shows a
  **threshold** at the option's sustain boundary and no graded size dependence
  above it, on either knob.
- [`docs/e3_results.md`](docs/e3_results.md) — **findings from E3** (timed
  response): double dissociation confirmed — Line A learns identity (wrong
  timing), Line B learns timing (not identity). New open problem: naive A+B
  *interferes* (both worse than either alone) under a single shared reward —
  later decomposed: factored credit removes the below-chance collapse (to ≈chance),
  a slow-first curriculum adds a marginal bimodal lift (joint composition on 1/5
  seeds) — direction supported, magnitude not established at n=5.
- [`docs/e4_results.md`](docs/e4_results.md) — **findings from E4** (attention):
  selective attention as biased winner-take-all by wave annihilation — a
  textbook psychometric (accuracy 0.96 at modest bias), the annihilation locus
  linear in the bias, achieved with zero inhibitory nodes.
- [`docs/e5_results.md`](docs/e5_results.md) — **findings from E5** (executive
  control): a persistent reentrant loop (the E2 mechanism) acts as an *option*
  that gates fast routing — switching 0.89 vs 0.20 when the loop is ablated,
  post-switch accuracy consolidating 0.57→0.92, single-rule routing spared by the
  ablation (0.87 vs 0.86). The discriminator localises the loop's role to *holding*
  the rule across a block.
- [`docs/e6_results.md`](docs/e6_results.md) — **findings from E6** (emergent
  categories): three GVF demons on one frozen substrate, reading the same feature
  vector, predict distinct questions well above baseline (memory R²=0.62, attention
  forecast 0.84, executive R²=0.98); their readouts are near-orthogonal and a
  generic probe matches an own-region oracle — memory/attention/executive are
  *questions asked of one machine*, not modules.
- [`docs/lattice_results.md`](docs/lattice_results.md) — **findings from the lattice arc**
  (learned timescales on an excitable sheet, 11 experiments, headlines at n=20): the
  input-timing `τ` rule works on a 2-D medium where the old self-referential rule ratchets to
  its ceiling, but exogenous timing does **not** penetrate a recurrent medium from a localised
  sensory strip (|τ−P| 1.54 at the strip, ~2.0 immediately beyond it, 2.92 at the far wall and
  non-monotone in between — locking nowhere). A 1-D attention chain of the same cells
  carries a timing reference to any depth (0.00 at depth 92) — a *clock, not a filter*. Reward
  as a fourth edge makes `τ` encode each cell's own **stimulus–reward interval** (|Δ| 0.16–0.19,
  97–98% within ±2, against an unpaired control receiving identical reward events that fails
  completely), and a backward value chain removes the last hand-set constant. Layered 2-D sheets
  produce the arc's one emergent *output* — a travelling wave becomes a synchronous burst timed
  to reward (91% within ±3 steps at D=70) — and show that value must arrive **diffusely** and
  couple **modulatorily** or synchrony collapses. The action primitive turns out to be
  **transmission, not emission**: the transmission edge sits at the learned interval (31.5 /
  51.5 / 71.5 for D=30/50/70) and graded credit exists from trial 1, so contingent reward needs
  no shaping. Two structural negatives: selective avoidance is **impossible** (transmission is
  provably monotone in probe time; max per-seed violation +0.0000), and plastic cell identity
  costs 49–76% of the propagation reach at any rate. Across four independent failures — amplitude,
  phase, direction, activity level — the arc's transferable law is that a raw signal *magnitude
  or geometry* cannot serve as a label on this substrate; only a structural or predictive signal
  can. **Not** reinforcement learning: no action changes the world, and the anatomy is still
  designed.
- [`docs/causal_experiments.md`](docs/causal_experiments.md) — **C-series plan**:
  using the substrate (where `W = f(S)` is explicit) as a synthetic-SCM testbed
  for the spike-wave causal question (arXiv:2511.06602) — validate the paper's
  certificates on ground truth, then show `do(W)` is fat-handed under real
  constitution and `do(θ)` is the well-posed handle.
- [`docs/c0_results.md`](docs/c0_results.md) — **findings from C0**: `W=f(S)`
  verified; the wave carries info beyond *partial* spikes for a collective code
  (growing as observation gets sparser) but not for a labeled-line code —
  informativeness is structure-dependent.
- [`docs/c1_results.md`](docs/c1_results.md) — **findings from C1**: on six
  canonical graphs the Theorem-1 certificate matches ground-truth `do(W)` —
  including the confounded case (association without causation) and front-door
  (causal despite an observed mediator).
- [`docs/c2_results.md`](docs/c2_results.md) — **findings from C2** (headline):
  when `W=f(S)` is constituted, one `do(W=w)` admits a huge behavioural band
  (33 σ) for a micro-reading behaviour vs ~0 for a collective one — `do(W)` is
  fat-handed and its causal verdict depends on the realization.
- [`docs/c3_results.md`](docs/c3_results.md) — **findings from C3**: `do(θ)`
  (timescales/couplings) is the well-posed handle — single-valued reproducible
  response, intervention ambiguity 0.014 σ vs `do(W)`'s 33 σ; `θ` is exactly
  what plasticity acts on.
- [`docs/c4_results.md`](docs/c4_results.md) — **findings from C4**: the causal
  role is (handle, outcome)-relative (`do(θ)` matrix is diagonal); the wave is
  the natural causal variable only where behaviour is collective (macro-
  sufficiency 1.03 vs 0.11 — causal emergence).
- [`docs/synthesis.md`](docs/synthesis.md) — **tying note**: the E-series and
  C-series are one argument — `θ` (timescales, couplings) is both the variable
  the learner adapts and the only well-posed causal handle; spikes and waves are
  two readouts of one parameterised dynamics.

### Process & reviews

- [`docs/process.md`](docs/process.md) — **how the project runs its planning and
  review passes**: decoupled in process, linked by a one-directional review→plan
  hand-off, and why. Read before doing either pass. ([`AGENTS.md`](AGENTS.md)
  points agents/humans here.)
- [`docs/core_review.md`](docs/core_review.md) — **independent integrity/overreach
  audit** of the core series (E0–E6, C0–C4): no fabrication; surfaced the E3
  framing overreach and the `perturb_tau` reproducibility bug.
- [`docs/extensions_review.md`](docs/extensions_review.md) — **self-audit** of the
  extensions (E7, C5–C7, E8.x): reproducible and honest; the residual risk is the
  substrate-vs-analysis boundary (afforded vs learned).
- [`docs/next_steps.md`](docs/next_steps.md) — **roadmap / planning**: candidate
  directions scored by which review-surfaced tension they retire.

## Progress

- [x] **E0** — substrate characterisation and operating point (see results)
- [x] **E1** — stimulus→response conditioning (A-vs-B dissociation confirmed)
- [x] **E2** — delayed response / working memory (dissociation inverts: B critical)
- [x] **E3** — timed response (double dissociation confirmed; A+B interference **decomposed**: factored credit removes the below-chance reward-conflation collapse 0.11→0.48, a slow-first curriculum adds a marginal, bimodal lift to 0.56 with genuine joint composition on 1/5 seeds — direction supported, magnitude not established at n=5; substrate-resonance capped)
- [x] **E4** — selective attention as biased WTA by wave annihilation (psychometric accuracy 0.96 at modest bias; **zero inhibitory nodes**)
- [x] **E5** — executive control / task switching: a persistent loop (E2 mechanism) as an *option* gating routing (switching 0.89 vs ablated 0.20; switch cost consolidates 0.57→0.92; single-rule spared 0.87 vs 0.86)
- [x] **E6** — emergent categories (Horde/GVF readout): three demons on one frozen substrate predict distinct questions (memory R²=0.62, attention 0.84, executive R²=0.98), near-orthogonal, no dedicated wiring — **E-series complete**

**C-series** (constitution & causality of spike–wave duality — see [`docs/causal_experiments.md`](docs/causal_experiments.md)):

- [x] **C0** — instrument the causal variables (`W=f(S)`; wave informative beyond partial spikes for a collective code only)
- [x] **C1** — certificate validated on ground truth (all 6 canonical graphs agree; confounded & front-door as key cases)
- [x] **C2** — `do(W)` is fat-handed for a constituted `W=f(S)` (achievable band 33σ vs ~0)
- [x] **C3** — `do(θ)` is the well-posed causal handle (ambiguity 0.014σ vs 33σ; `θ→W→B`)
- [x] **C4** — outcome-relativity (diagonal `do(θ)` matrix) & degeneracy (macro-sufficiency 1.03 vs 0.11) — **C-series complete**

**Closed-Loop Substrate Plasticity** (see [`docs/closed_loop_plasticity_results.md`](docs/closed_loop_plasticity_results.md)):

- [x] **Multi-Axis Plasticity** — tri-axis closed-loop engine ($\tau$-adaptation, $\theta$-homeostasis, $W$-routing) achieves RIR = $0.851 \pm 0.111$ on E1
- [x] **Substrate Credit Assignment & Anti-Forgetting** — Task A $\to$ Task B $\to$ Task A sequential reversal learning yields **70.0% ± 13.8%** retention (vs 29.6% ± 31.2% weight-only) via topological loop protection

Both E/C series and closed-loop substrate plasticity tracks are complete. See [`docs/synthesis.md`](docs/synthesis.md) for how E-series and C-series tie together.

## Reproduce

```bash
uv sync
python3 reproduce_all.py                      # runs 8-step automated test & verification harness
```

Each `experiments/*.py` is self-contained and writes its figures to `docs/figures/`
and data to `result/`.
