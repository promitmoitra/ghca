# C-Series Extension: Causality of the Spiral Option (`do(chirality)`)

*A continuation of the [C-series](causal_experiments.md), pointing the same
synthetic-SCM machinery at the 2-D spiral core built in [E7](e7_results.md). Where
C0–C4 asked whether a scalar wave aggregate `W = f(S)` (coherence, active fraction)
is a causal handle, C5–C7 ask the sharper, field-relevant version: is the **rotation
direction (chirality) of a genuine spiral** — the variable the rotating-wave
neuroscience (Xu/Gong et al. 2023; Ye/Steinmetz et al. 2026) calls task-relevant —
a causal handle on behaviour, or an epiphenomenal reflection of the spiking?*

## Motivation

E7 established that a persistent 2-D spiral core's chirality (a) **carries** the
executive-control rule (decode 1.00), (b) **drives** reward-learned routing
(switching 0.86), and (c) is **necessary** for switching (ablation → 0.49). Those
are strong association + necessity results, but they do not settle the causal
question the spiral-wave literature is fighting over: does the rotating wave *do
computational work*, or is it a readout of the spikes that do the work?

The C-series already built the tools and the answer-shape for the scalar case:
`do(W)` on a constituted aggregate is fat-handed (C2), the well-posed handle is the
generating parameters `do(θ)` (C3), and the verdict is outcome-relative (C4). C5–C7
instantiate that for the spiral, with one genuinely new twist:

- **Chirality is a topological invariant**, not a continuous aggregate. C2's
  fat-handedness came from the many micro-states realizing a given active fraction.
  A winding number `χ` is realization-*invariant* by definition — so C5 asks whether
  the spiral's rotation direction **escapes** C2's fat-handedness (a "better" wave
  variable), or whether the behavioural readout re-introduces it (because the E7
  router reads a *local* winding at a location, and realizations that move the core
  defeat that readout). This is the crux, and C5 is written to *measure* it, not
  assume the C2 result carries over.

## Shared setup (deltas from the C-series)

- **Substrate.** The E7 spiral lattice (`lattice2d`, r=2, a=6, τ=14, θ=4, no-flux)
  + the E7 router, reused. Behaviour `B` = the router's decoded rule / action
  correctness.
- **Wave variable under test.** `χ = winding(f(S))` — the topological charge /
  rotation direction of the collective phase field. Readouts (all already in
  `e7_spiral_option.py`): local winding at centre, local winding at the tracked
  core, global net charge (`signed_charge`, `local_winding`).
- **New intervention `do(χ = c)`** (to add to `ghca_causal.py`): force the phase
  field to a target winding `c` under a **realization policy π** — the multiplicity
  of phase fields with the same `χ` is the whole point. Policies: {core-centred,
  core-displaced by `d`, tight vs loose arm pitch, high vs low activity,
  random-phase constrained to winding `c`}.
- **Generative handle `do(θ_χ)`.** Set chirality via the nucleation seed sign (or a
  pinned chirality-biasing heterogeneity) — modular, no realization freedom. This is
  E7's actual mechanism.
- **Measures (as in C2/C3).** The achievable **band** of `B` across realizations
  (σ units) for `do(χ)`, vs the across-seed band of the generative handle `do(θ_χ)`;
  mediation and an outcome matrix for C7.

## The sequence (C5 → C6 → C7)

### C5 — Is `do(χ)` fat-handed? Topological invariant vs realization *(headline, gate)*
Force `χ = ±1` under ≥4 realization policies; inject each into the router; measure
the band of `B`. Repeat across the three readouts (centre / tracked-core / global).
- *Hypothesis A (fat-handed, like C2).* The local-at-centre readout makes `do(χ)`
  fat-handed: realizations that displace the core carry the same `χ` but yield a
  different/absent routing verdict → wide band.
- *Hypothesis B (topology escapes C2).* Because `χ` is realization-invariant and a
  tracked-core / global readout ignores core position, `do(χ)` pins `B` with a
  *small* band — rotation direction is a better-posed wave variable than a generic
  aggregate.
- *Discriminator.* The band's dependence on readout locality decides A vs B; either
  is a real result. *Gate:* the outcome sets C6's contrast (a C3-style fix if A, or
  a near-zero-band contrast if B).
- *Status: DONE.* See [`c5_results.md`](c5_results.md). **Both, readout-decided.**
  Topological chirality *escapes* C2's 33σ fat-handedness when read topologically
  (tracked-core band **1.0σ**), but a fixed-locus readout re-introduces it (center
  band **6.2σ**; global 2.6σ). Behavioural confirmation on a frozen E7 router at a
  displaced core: routing 0.55 (center) vs 0.78 (tracked). The sharpened lesson:
  well-posedness is a property of the *(variable, reader)* pair — a winding number is
  a better-posed wave variable than a scalar aggregate, but only for a topology-aware
  reader. Implemented in `experiments/c5_do_chirality.py`.

### C6 — `do(θ_χ)` (nucleation) is the well-posed handle + necessity
Set `χ` generatively (seed sign / pinned rotor); measure `B` across seeds → expect
single-valued, tiny band (reuse C5's band for the contrast panel, as C3 reuses C2).
Fold in the E7 ablation as an explicit `do`-style **necessity** test (remove the
persistent core → switching collapses, single-rule spared).
- *Expect.* `θ_χ → χ → B` is modular and reproducible; the generative handle is
  well-posed regardless of C5's A/B outcome — "drive the parameters."
- *Contribution.* The handle E7 actually uses (nucleation) is the causally clean one.
- *Status: DONE.* See [`c6_results.md`](c6_results.md). `do(θ_χ)` (centred
  nucleation) gives a **0.0σ** band for *every* reader (center/tracked/global) —
  well-posed even for C5's fat-handed fixed-centre reader — vs `do(χ)`'s 6.2/1.0/2.6σ.
  Necessity (do-ablation of persistence): switching **0.85 → 0.52** while single-rule
  is spared (**0.90 vs 0.89**) — the persistent core is causally required for
  switching, not routing. Implemented in `experiments/c6_do_theta_chi.py`.

### C7 — Outcome-relativity, mediation, synthesis *(cap)*
- **Mediation `θ→χ→B`.** Conditioning on the decoded `χ` screens off the seed — `χ`
  is a genuine mediator, not a side-effect.
- **Outcome matrix (C4-style).** `do(χ)` moves the *rule* axis (→1) not
  identity-within-rule; `do(g_route)` moves identity (→1) not rule → diagonal. So
  chirality is causal *for the rule outcome specifically*.
- **Synthesis.** The spiral's rotation direction is a genuine causal **mediator**
  (behaviour reads it; `θ` acts through it; ablating it kills switching) — *not
  epiphenomenal* — while its direct `do(wave)` form is (or is not; per C5)
  fat-handed. "Read the wave, drive the parameters," now for a real spiral, and a
  direct entry in the cortical spiral-wave causality debate.
- *Status: PLANNED.*

## Breadth constraints (deliberately excluded)

- One spiral substrate (the E7 lattice); one router (E7's).
- Chirality is the wave variable of interest; not a rescan of coherence / active
  fraction (that was C0–C4).
- No new learning rules; C5–C7 are causal *analysis* of the E7-trained system.
- No real neural data; the contribution is *which query about the spiral is
  well-posed*, not experimental accessibility (per C3's caveat).

## One-line spine

C5 asks whether topological chirality escapes `do(W)`'s fat-handedness → C6 shows
the nucleation handle `do(θ_χ)` is well-posed and the persistent core is necessary →
C7 shows rotation direction is a causal mediator, outcome-relative, not
epiphenomenal.

## Reference anchors

- Jalaldoust & Zabeh, *A Causal Formulation of Spike–Wave Duality* (arXiv:2511.06602)
  — the framework the whole C-series tests.
- Xu/Gong et al., *Nat. Hum. Behav.* 2023; Ye/Steinmetz et al., *Science* 2026 —
  cortical spiral rotation direction as task-relevant (the phenomenon C5–C7 probe).
- Hoel — causal emergence / macro-sufficiency (the C4 lens, reused in C7).
- [`causal_experiments.md`](causal_experiments.md) · [`e7_results.md`](e7_results.md)
  · [`synthesis.md`](synthesis.md).
