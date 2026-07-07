# C-Series: Constitution & Causality of Spike–Wave Duality

*A methodological arc, distinct from the learning E-series, using the GH network
as a substrate where the **constitution `W = f(S)` is explicit and controllable**.
It is a synthetic-SCM testbed for the causal framework of Jalaldoust & Zabeh,
*A Causal Formulation of Spike-Wave Duality* (arXiv:2511.06602, 2025).*

## Motivation

That paper defines a wave `W` as **epiphenomenal** to behaviour `B` given
observed spikes `S` iff the interventional distribution is invariant:
`P(B | S; do(W)) = P(B | S)` (their Definition 1), and gives graphical
certificates (Prop 1, front-door Prop 2, Theorem 1) for when `W` is removable.
Two things are hard in vivo and easy here:

1. **Ground truth.** We *build* the SCM, so we can check whether their
   certificates correctly predict measured (in)variance — impossible with real
   brains, where the graph is unknown (Causal Hierarchy Theorem).
2. **Constitution.** Physically a wave is a deterministic aggregate of spikes,
   `W = f(S)` (the low-pass of the same signal). Treating `W` as an autonomous
   manipulable node makes `do(W)` **fat-handed**: you cannot set the aggregate
   without also setting its constituents. The paper parks this under "recursive
   coupling as a declared assumption." Our substrate makes `W = f(S)` literal,
   so we can test whether `do(W)` is even well-defined, and propose the
   well-posed alternative: intervene on the **wave-generating parameters** `θ`
   (timescales, couplings), not on the wave.

## Shared setup (fixed across the series)

- **`S` — spikes (micro), partially observed.** Node states `φ`; `S_obs` is a
  *subset* of nodes. Partial observation is essential: given the *full* micro
  state, `W = f(S)` adds no information, so the question is only non-trivial
  under partial `S_obs` (their `C/M/U` decomposition).
- **`W` — waves (macro), explicitly `W = f(S)`.** Primary: Kuramoto phase
  coherence `R`. Secondary: active fraction `A`. Deterministic coarse-grainings
  of the *full* node set, recomputed each step.
- **`B` — behaviour.** Motor readout / action / latency (reuses E-series I/O).
- **`θ` — generating parameters.** Per-node `τ`, coupling weights, thresholds.
- **Interventions.**
  - `do(S)`: clamp chosen nodes (existing cue machinery).
  - `do(W=w)`: force `f(S)=w` by projecting `S` onto the level set
    `{S : f(S)=w}` under a **realization policy π** — the multiplicity of valid
    `S` is the whole point.
  - `do(θ)`: set `τ`/couplings and let dynamics run — modular, no state-clamp.
- **Epiphenomenality test (Def 1).** `Z` epiphenomenal to `B` given `X` iff
  `P(B|X; do(Z)) = P(B|X)`; estimated by Monte Carlo + a distributional distance
  and a shuffle test.

## The sequence (decision-gated; stop-early-complete at C2/C3/C4)

### C0 — Instrument variables and operators *(gate)*
Confirm `S_obs, W, B` and the three `do`-operators are well-defined, that
`W = f(S)` is deterministic, and that `W` carries predictive info about `B`
**beyond partial `S_obs`** (the analog of the paper's Eq. 1,
`P(B|S_obs,W) ≠ P(B|S_obs)`), while adding nothing given full `S`.
*Gate:* if `W` adds no info even under partial observation, the substrate can't
host the question.
*Status: DONE.* See [`c0_results.md`](c0_results.md). `W = f(S)` verified
(adds `+0.001` beyond full spikes). Under partial observation the wave adds
predictive info beyond `S_obs` **for a collective code** (up to `+0.11`, growing
as observation gets sparser) but **~0 for a labeled-line code** — so whether the
wave is informative is structure-dependent, decided by how behaviour reads the
population. Gate passed. Implemented in `ghca_causal.py` +
`experiments/c0_instrumentation.py`.

### C1 — Build the paper's graphs; validate the certificates on ground truth
Wire the net into the canonical graphs by construction — fork, mediated
(observed vs unobserved mediator), front-door (mediators observed, confounder
hidden) — controlling which nodes are `S` vs hidden `U`. Treat `W` as a
*designated node* here to isolate graph structure. Compare `P(B|S; do(W))` vs
`P(B|S)` and check against Theorem 1's d-separation prediction.
*Expect:* fork → invariant (epiphenomenal); unobserved-mediator & front-door →
causal. *Contribution:* empirical validation of their theorems on a known SCM.
*Status: DONE.* See [`c1_results.md`](c1_results.md). All six canonical graphs
agree: Theorem-1 certificate == ground-truth `do(W)` == expected verdict. The
`confounded` graph shows strong association (0.64) with zero causal effect
(correlation≠causation); `front-door` is causal (0.26) despite conditioning on
the observed mediator (Prop 2). Implemented in
`experiments/c1_graph_certificates.py`. (Subtlety: Definition 1 must be scored
interventional-vs-*observational*, not `do(W=1)` vs `do(W=0)`, or front-door is
mis-scored.)

### C2 — Constitution / fat-handed `do(W)` *(headline)*
Now `W = f(S)` is the genuine constituted aggregate. Implement `do(W=w)` under
≥3 realization policies `π` (random subset, spatially-clustered, minimal-edit).
Measure `P(B; do_π(W=w))` across `π`.
*Expect:* the causal verdict **varies with `π`** → `do(W)` is ill-posed
(fat-handed) for a constituted `W`. *Gate:* motivates `do(θ)`.

### C3 — `do(θ)`: the generating parameter as the well-posed handle
`do(τ)` and `do(coupling)`; measure effect on `B` and mediation `θ→W→B`.
*Expect:* well-defined, realization-independent, changes `B` via `W`, and
recovers E3 (`do(τ)` moves timing not identity). *Contribution:* in mechanistic
systems the causal question should target generating parameters, not the
aggregate.

### C4 — Outcome-relativity & degeneracy *(synthesis cap)*
Test Def-1 invariance of the *same* `W` under `do(θ)` for two behaviours
(`B_identity`, `B_timing` from E3); quantify `S→W` degeneracy (an
effective-information / invariance measure).
*Expect:* `W` epiphenomenal for one outcome, causal for the other; `W` is the
"natural" causal variable exactly where `S→W` degeneracy is high (causal
emergence, Hoel). Closes the arc.

## Breadth constraints (deliberately excluded)

- One substrate; one primary `W` (coherence), one secondary (active fraction) —
  no topology sweep.
- We always **know** the graph — no identification-from-data (the paper proves
  that's impossible; not our contribution).
- No real neural data; no exhaustive `W`-zoo; no new learning rules.
- C4 is optional; the core result (constitution breaks `do(W)`, `do(θ)` fixes
  it, verdict is outcome-relative) stands after C2–C3.

## One-line spine

C0 grounds the variables → C1 shows the testbed reproduces the theorems → C2
shows `do(W)` is ill-posed under real constitution → C3 shows `do(θ)` is the
well-posed causal handle → C4 shows the verdict is outcome-relative.
