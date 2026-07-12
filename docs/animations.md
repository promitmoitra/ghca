# Watch the substrate — a narrated animation gallery

The E- and C-series results are reported as static plots (learning curves, bar
charts, phase diagrams). But the whole premise is about *dynamics*: one
homogeneous excitable medium that "generates its own repertoire of patterns,"
from which reward selects useful ones and an observer reads out faculties. This
page lets you **watch** that medium — the same Greenberg–Hastings substrate,
animated for each mechanism.

All animations are rendered by [`ghca_net_viz.py`](../ghca_net_viz.py) from a
recorded phase rollout (`Network.run(...)["phi"]`). Colour = GH state: light grey
rested, bright active crest, darkening refractory tail.

---

## 0. The raw substrate — self-organised spiral waves

At the E0 operating point (`r=2, a=6, θ≈4`) the medium neither dies nor
saturates: it settles into organised, self-sustaining spiral waves. Nothing is
driving it — this is the "repertoire" the later experiments shape and read.

![substrate spirals](figures/demo_lattice.gif)

*(`python3 ghca_net_viz.py`)*

---

## 1. Memory — a τ-tuned reentrant loop (E2)

Two identical directed rings are ignited the same way; they differ **only** in
the local timescale `τ`. Left (`τ=20 < L=24`) sustains its rotating pulse
indefinitely — the stimulus is *held*. Right (`τ=28 ≥ L`) dies within ~`L` steps,
its refractory tail wrapping round to block reentry — the memory is *lost*.
Memory is not a module; it is a loop the substrate has tuned to forget slowly.

![E2 ring memory](figures/e2_ring_memory.gif)

*(`python3 experiments/e2_animation.py`; see [e2_results.md](e2_results.md))*

---

## 2. Attention — biased winner-take-all by wave annihilation (E4)

Two waves ignite opposite ends of an excitable chain (blue = attended left,
orange = right), travel inward, collide, and **annihilate in each other's
refractory wake**. A top-down bias gives the attended stream a head start, so the
collision lands off-centre and the attended wave captures the centre node. This
is winner-take-all with **no inhibitory machinery** — refractoriness *is* the
competition.

![E4 wave annihilation](figures/e4_annihilation.gif)

*(`python3 experiments/e4_animation.py`; see [e4_results.md](e4_results.md))*

---

## 3. Executive control — a slow loop as an option (E5)

The rule in force is literally **which slow loop is rotating**. At each block
start the active rule's context ring is ignited and persists (the E2 mechanism,
`τ < L`); the other ring stays dark. At the next block they flip. This standing
rotation is the held context that gates the fast stimulus→response routing —
executive control as an *option*, not a dedicated controller.

![E5 options](figures/e5_options.gif)

*(`python3 experiments/e5_animation.py`; see [e5_results.md](e5_results.md))*

---

## The point of the gallery

Memory (a sustaining loop), attention (annihilating waves), executive control (a
gating loop) are **three readouts of one excitable medium** — the same
`Network`, the same `(active, refractory, rested)` state, the same colour key.
No mechanism above adds a bespoke module; each is a different question asked of
one substrate. That is the Buzsáki inside-out thesis made watchable, and it is
the E6 "emergent categories" result seen as motion rather than as a table.

## A note on E3 (timed response)

E3 is deliberately **not** given a headline "mechanism" animation. Its solid
finding — response latency tracks the gate timescale — is a one-line relationship
better shown as the static `latency = τ − a` plot. Its *interesting* result is
the **fragility** of composing timing with identity, which is a per-seed scatter,
not a dynamic (see the [hallucination review](hallucination_review.md), E3
deep-dive). Animating a clean E3 "composition" would misrepresent a result that
is genuinely bimodal and seed-dependent — so it isn't done here.

## Reproduce all

```
python3 ghca_net_viz.py                     # docs/figures/demo_lattice.gif
python3 experiments/e2_animation.py         # docs/figures/e2_ring_memory.gif
python3 experiments/e4_animation.py         # docs/figures/e4_annihilation.gif
python3 experiments/e5_animation.py         # docs/figures/e5_options.gif
```

## Future scope / backlog

- **Interactive HTML explorer.** A self-contained web page with live sliders for
  `τ`, `θ`, bias, ring length `L`, etc., re-rendering the substrate in the browser
  — turning these fixed GIFs into something you can drive. Needs either a small JS
  reimplementation of the GH update rule or a pre-rendered parameter sweep. Deferred.
- **Scrollytelling narration.** The same walkthrough as a single scrolling page
  where each animation plays alongside the claim it supports.
- **C-series causal animation.** `do(θ) → W → B` vs the fat-handed `do(W)` as a
  side-by-side dynamic, to give the causal argument the same watchable treatment.
