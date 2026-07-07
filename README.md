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
| `experiments/e0_characterization.py` | E0 — substrate characterisation (find the self-sustaining band) |
| `experiments/e1_conditioning.py` | E1 — stimulus→response conditioning (A-vs-B dissociation) |
| `experiments/e2_delayed_response.py` | E2 — delayed response / working memory (τ-controlled memory) |
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

## Progress

- [x] **E0** — substrate characterisation and operating point (see results)
- [x] **E1** — stimulus→response conditioning (A-vs-B dissociation confirmed)
- [x] **E2** — delayed response / working memory (dissociation inverts: B critical)
- [ ] E3 — temporal sequence reproduction (the A-vs-B discriminator)
- [ ] E4 — selective attention (cue competition)
- [ ] E5 — executive control / task switching (options)
- [ ] E6 — emergent categories (Horde/GVF readout)

## Reproduce

```
python3 -m pip install numpy matplotlib scipy
python3 experiments/e0_characterization.py    # writes docs/figures/e0_*.png, result/e0/
```
