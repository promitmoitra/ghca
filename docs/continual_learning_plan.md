# 3c Plan — Continual Learning as Causal Credit Assignment

*Re-scope of Track 3c ([`next_steps.md`](next_steps.md)). The original question —
"can one substrate learn E1→E5 sequentially without catastrophic interference?" —
is sharpened by the continual-learning literature into a question worth running,
and fused with the "learning as causal inference" thread. Nothing built yet.*

## Why the naive version isn't a result

The reservoir-computing literature is blunt about the easy case: **freeze the
recurrent substrate and give each task an independent output head, and
catastrophic interference is a non-issue** — output weights of different heads are
independent ([continual learning with echo state networks, arXiv:2105.07674](https://arxiv.org/pdf/2105.07674)).
But that frozen-substrate/per-head setup **is essentially E6** (Horde demons on a
frozen substrate). So a naive 3c either reproduces E6 (proves nothing) or hits the
known-hard plastic case and interferes (also not a finding on its own). The study
has to target the part that is *both hard and novel*: the plastic shared substrate,
with credit assignment as the lever.

The frontier CL methods handle the hard case with mechanisms **we already have**:
competitive/federated reservoirs gate learning by lateral inhibition so only the
most-predictive head updates ([Bereska & Boix, arXiv:2206.13336](https://arxiv.org/abs/2206.13336));
conceptors carve orthogonal per-task state subspaces ([Jaeger — Overcoming
Catastrophic Interference by Conceptors](https://www.researchgate.net/publication/318489433_Overcoming_Catastrophic_Interference_by_Conceptors)).
E9's k-WTA competition and E4/E5's WTA are the substrate's native version of the
first.

## The question (re-scoped)

**On the hard (plastic-substrate) regime where interference is expected, does
*causal* credit assignment reduce forgetting relative to *correlational* credit —
and how close does it get to the frozen-reservoir no-interference bound?**

This fuses 3c with learning-as-causal-inference. The hypothesis rests on the
C-series: C2/C3 showed `do(θ)` is the well-posed causal handle while `do(W)` is
fat-handed, and Mesnard et al.'s counterfactual credit assignment shows
future-conditional (causal) credit is **provably low-variance** in model-free RL
([ICML 2021](https://proceedings.mlr.press/v139/mesnard21a.html);
[arXiv:2011.09464](https://arxiv.org/abs/2011.09464)). If credit is assigned to the
parameters that *caused* a task's reward rather than those merely correlated with
it, cross-task interference should drop.

## Design (v1)

**Task sequence:** E1 → E2 → E5 on one shared substrate (conditioning → working
memory → executive switching — spans both plasticity lines and reuses clean,
CI-quantified dissociations). E3/E7 deferred (E3 is the fragile one; E7 spiral is
geometry-bound).

**Conditions (2 regimes × 2 credit rules for v1):**

| | correlational credit (current eligibility trace) | causal credit (θ-perturbation `do(θ)`) |
|---|---|---|
| **frozen** substrate + per-task heads | reservoir **upper bound** (≈ E6; expect no interference) | — (credit rule irrelevant when frozen) |
| **plastic** shared substrate | expected **interference floor** | **the hypothesis** |

So three runs: frozen (control/upper-bound), plastic+correlational (floor),
plastic+causal-θ (test). The Mesnard **hindsight/future-conditional** estimator is
a **v2 follow-on**, not v1 — keep the first cut to the contrast the C-series
directly predicts.

**Metrics (with 3a's CI discipline, n≥30, bootstrap CIs, per-seed spreads):** the
standard continual-learning trio measured after the full sequence —
- **average accuracy** across the three tasks,
- **backward transfer / forgetting** (drop on task *k* after training *k+1…*),
- **forward transfer** (head start on later tasks),
plus a per-task learning curve.

**Optional 4th condition (if v1 is inconclusive):** add the substrate's native
competitive gating (E9 k-WTA / E4 WTA) to the plastic regime — tests whether the
*biological* interference solution the substrate already contains does the job,
independent of the credit rule.

## Predictions & honest priors

- **Frozen** wins trivially — framed as a **control / upper bound**, not a result.
- **Plastic + correlational** likely interferes (the expected floor).
- **Plastic + causal-θ** is the falsifiable claim: it should sit *above* the
  correlational floor on backward transfer. **Null path:** it doesn't — which would
  say the interference is *representational* (shared substrate capacity), not a
  credit-assignment artifact. A clean null here is an acceptable, publishable
  outcome and is stated as such up front.

## Scope guards (avoid the E3 failure mode)

- Small task set (3), one substrate size, one operating point — a generality sweep
  is explicitly out of v1.
- The frozen arm is a control, never reported as the headline.
- Report per-seed forgetting distributions and call out bimodality (the E3 lesson);
  do not headline a mean over a split.
- State the substrate/analysis boundary: what adapts (recurrent W/θ) vs what is a
  fixed readout, per condition.

## Phasing

- **P1** — harness: sequential-task driver over E1/E2/E5 sharing one substrate;
  CL metrics (avg acc, BWT, FWT); the frozen control + plastic+correlational floor.
- **P2** — the causal-θ credit rule (node-perturbation `do(θ)` credit) as the third
  condition; the headline contrast.
- **P3** — (optional) native competitive gating condition; and/or the Mesnard
  hindsight estimator (v2).
- **P4** — results doc, figure, and roadmap close-out; flag any result-doc edits.

## Effort & risk

**Effort:** high — this is a new learning harness, not a seed-scaleup. **Risk:**
high and *intended* — the causal-vs-correlational contrast is genuinely uncertain,
and a null is informative. This is the boldest test of the "one homogeneous
machine that genuinely learns" claim, and the concrete vehicle for unifying the
learning and causality arcs.
