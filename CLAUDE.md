# Working in this repo

This is an exploratory computational-neuroscience study: Greenberg–Hastings
excitable dynamics on a graph as the substrate for a reward-driven learner
(E-series), plus a causal-inference instrumentation of the same substrate
(C-series). See [`README.md`](README.md) for the map of files and results.

## Process — read before a planning or review pass

The project runs two recurring meta-passes, **kept decoupled in process but
linked by a one-directional hand-off (review → plan)**. Before doing either,
read [`docs/process.md`](docs/process.md). In short:

- **Review** (adversarial, backward-looking) re-runs experiments, checks every
  number reproduces, and audits for overreach. Artifacts:
  [`docs/core_review.md`](docs/core_review.md) (independent, E0–E6/C0–C4) and
  [`docs/extensions_review.md`](docs/extensions_review.md) (self-audit).
- **Planning** (generative, forward-looking) consumes the review findings and
  lays out the roadmap: [`docs/next_steps.md`](docs/next_steps.md).

Do not fold a review into a plan or vice versa, and never soften a review
finding to justify a roadmap track. Both audits live on `main` with the work.

## House rules (apply to every experiment)

- Seed *everything*: thread `default_rng(seed)` explicitly; never use the global
  NumPy RNG (the `perturb_tau` bug — see `core_review.md`).
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the *dynamics* do vs
  what a *readout/feature* does.
- Keep a caveats section adjacent to every headline.
