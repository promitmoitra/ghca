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
| `ghca_causal.py` | **Causal instrumentation** (C-series): partial-observation `S_obs`, wave variables `W=f(S)`, and `do(S)` / `do(W)` / `do(θ)` intervention operators |
| `experiments/e0_characterization.py` | E0 — substrate characterisation (find the self-sustaining band) |
| `experiments/e1_conditioning.py` | E1 — stimulus→response conditioning (A-vs-B dissociation) |
| `experiments/e2_delayed_response.py` | E2 — delayed response / working memory (τ-controlled memory) |
| `experiments/e2_information.py` | E2 addendum — memory as a τ-tuned information-destruction rate |
| `experiments/e3_timed_response.py` | E3 — timed response (identity × latency double dissociation) |
| `experiments/e3_factored_credit.py` | E3 composition study — factored credit + curriculum vs shared reward |
| `experiments/c0_instrumentation.py` | C0 — instrument the causal variables (`W=f(S)`, partial spikes) |
| `experiments/c1_graph_certificates.py` | C1 — validate Theorem-1 epiphenomenality certificate on known SCMs |
| `experiments/c2_fat_handed.py` | C2 — `do(W)` is fat-handed when `W=f(S)` (achievable-band of behaviour) |
| `experiments/c3_do_theta.py` | C3 — `do(θ)` (timescales/couplings) is the well-posed causal handle |
| `experiments/c4_outcome_relativity.py` | C4 — outcome-relativity & degeneracy (causal-emergence cap) |
| `result/` | Saved simulation outputs (`.npy`) and experiment data |

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
- [`docs/e3_results.md`](docs/e3_results.md) — **findings from E3** (timed
  response): double dissociation confirmed — Line A learns identity (wrong
  timing), Line B learns timing (not identity). New open problem: naive A+B
  *interferes* (both worse than either alone) under a single shared reward.
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

## Progress

- [x] **E0** — substrate characterisation and operating point (see results)
- [x] **E1** — stimulus→response conditioning (A-vs-B dissociation confirmed)
- [x] **E2** — delayed response / working memory (dissociation inverts: B critical)
- [x] **E3** — timed response (double dissociation confirmed; A+B interference **decomposed & partly resolved**: factored credit + slow-first curriculum lift joint identity 0.20→0.77, residual is a substrate resonance artifact)
- [ ] E4 — selective attention (cue competition)

**C-series** (constitution & causality of spike–wave duality — see [`docs/causal_experiments.md`](docs/causal_experiments.md)):

- [x] **C0** — instrument the causal variables (`W=f(S)`; wave informative beyond partial spikes for a collective code only)
- [x] **C1** — certificate validated on ground truth (all 6 canonical graphs agree; confounded & front-door as key cases)
- [x] **C2** — `do(W)` is fat-handed for a constituted `W=f(S)` (achievable band 33σ vs ~0)
- [x] **C3** — `do(θ)` is the well-posed causal handle (ambiguity 0.014σ vs 33σ; `θ→W→B`)
- [x] **C4** — outcome-relativity (diagonal `do(θ)` matrix) & degeneracy (macro-sufficiency 1.03 vs 0.11) — **C-series complete**

See [`docs/synthesis.md`](docs/synthesis.md) for how the E-series and C-series tie together.
- [ ] E5 — executive control / task switching (options)
- [ ] E6 — emergent categories (Horde/GVF readout)

## Reproduce

```
python3 -m pip install numpy matplotlib scipy
python3 experiments/e0_characterization.py    # writes docs/figures/e0_*.png, result/e0/
```
