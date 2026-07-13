# GHCA — a Greenberg–Hastings learning substrate

A study of local timescales in Greenberg–Hastings cellular automata, and an
exploratory program using GH **excitable dynamics** as the substrate for a
reward-driven **learning mechanism** on a graph.

The framing is *inside-out* (Buzsáki): the medium generates its own repertoire of
dynamical patterns, and a strict scalar reward (Sutton) selects and stabilises a
subset through action. Memory, attention and executive function are treated as
categories **read out** from one homogeneous substrate — not as built-in modules.

<div class="grid cards" markdown>

- :material-play-circle: **[Watch the substrate](animations.md)**
  A narrated animation gallery — spiral waves, ring memory, wave-annihilation
  attention, and the executive "option" loop, all on one excitable medium.

- :material-flask: **[E-series](e0_results.md)**
  Reward carves cognition out of the medium: conditioning, working memory, timed
  response, attention, executive control, emergent categories.

- :material-sitemap: **[C-series](causal_experiments.md)**
  A synthetic causal testbed for the spike–wave question: `do(W)` is fat-handed
  under constitution; `do(θ)` is the well-posed handle.

- :material-magnify-scan: **[Core series review](core_review.md)**
  An independent audit of every claim against code and data — what reproduces,
  what is overstated, and one real reproducibility bug.

</div>

## The two strands

**1. Lattice GH cellular automata (original).** Excitable-media dynamics on a 2D
lattice with an `(active, passive)` refractory cycle, plus tools to enumerate
which initial configurations self-sustain.

**2. Learning on a GH network.** GH dynamics generalised to an arbitrary weighted
graph with per-node timescales, as the substrate for a reward-driven learner. Two
parallel plasticity lines — conduction weights (**Line A**) and local timescales
(**Line B**) — and a series of staged experiments (E0–E6) probing what each can
learn, tied to a causal analysis (C0–C4) of which variables are the well-posed
handles.

Start with **[Watch the substrate](animations.md)** to see the dynamics in
motion, then dive into the [design & framework](learning_experiments.md) or the
[synthesis](synthesis.md) tying the two series together.
