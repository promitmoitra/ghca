# 3c Results — Continual Learning (P1 baselines · P2 causal credit · P3 low-variance)

*Track 3c, Phase 1 (see [`continual_learning_plan.md`](continual_learning_plan.md)).
Establishes the continual-learning harness, the metrics, and the two "what adapts"
baselines on one shared substrate, under the current **correlational** credit rule.
The **causal `do(θ)` credit** contrast — the actual hypothesis — is P2.
Experiment: `experiments/continual_learning.py`; figure: `continual_p1_figure.py`;
data: `result/stats/continual_p1.*`.*

## Setup

**One shared substrate** (E1-style `layered_graph`, K=A=3) is trained sequentially
on **three maximally-conflicting remappings** — `T0: a=stim`, `T1: a=stim+1`,
`T2: a=stim+2 (mod 3)` — that share the same readout, so later tasks overwrite
earlier ones (the canonical catastrophic-interference setup). This v1 deliberately
uses a remapping triple rather than the heterogeneous E1→E2→E5, to isolate
*interference* from *architecture-swapping* and reuse E1's validated conditioning
machinery; E1→E2→E5 is v2.

Two "what adapts" regimes, both with the **correlational eligibility-trace** credit:
- **frozen** — the S→H representation is fixed; only the shared H→M head adapts.
- **plastic** — S→H and H→M both adapt.

Metrics after the full T0→T1→T2 sequence, n=30, bootstrap 95% CIs: average accuracy,
**backward transfer** (forgetting of earlier tasks), forward transfer.

![continual P1](figures/continual_p1.png)

## Results

| regime | avg accuracy | backward transfer | forward transfer |
|---|:--:|:--:|:--:|
| frozen | 0.316 [0.305, 0.327] | **−0.258 [−0.309, −0.207]** | −0.236 |
| plastic | 0.291 [0.283, 0.299] | **−0.174 [−0.225, −0.128]** ⚠bimodal | −0.232 |

*(chance = 0.333)*

**1. Catastrophic interference is real and severe.** Average accuracy collapses to
≈chance after the sequence, and backward transfer is strongly negative in both
regimes: a task trained to ~0.7 (e.g. frozen T0, 0.72 just-trained) falls to ~0.32
once the later tasks are learned. This is the necessary precondition for P2 — there
is a large interference gap for a better credit rule to try to close.

**2. Freezing the representation does *not* rescue it — and slightly hurts.** The
reservoir-CL intuition is "freeze the recurrent medium and interference goes away."
Here frozen backward transfer (−0.258) is if anything *worse* than plastic
(−0.174), CIs barely overlapping. The reason is structural: the three tasks
conflict on the **shared readout**, which is plastic in *both* regimes, so freezing
the representation cannot help — and letting the representation adapt gives the
later tasks a little room to re-separate (visible as the plastic matrix's higher
off-diagonal recovery). So the reservoir "freezing helps" claim has a boundary
condition: it needs *per-task* heads, not a shared one.

**3. Forward transfer is negative** (~−0.23 both) — training task *k−1* leaves task
*k* below chance before it is trained, as expected for adversarial relabellings.

## Honest caveats

- **"frozen" here is not the reservoir no-interference upper bound.** That bound
  requires *per-task independent heads* (interference-free by construction, ≈1.0)
  and is not run — it would be a trivial control. Our "frozen" is frozen-
  representation + **shared** head, a genuine interference condition chosen so there
  is something to study.
- **v1 task set is a remapping triple, not E1→E2→E5** (flagged in the plan) — a
  deliberate simplification to isolate interference; the heterogeneous port is v2.
- **plastic backward transfer is bimodal** — some seeds re-separate the tasks, some
  collapse; the mean hides two modes (the E3 lesson — reported, not buried).
- **Absolute accuracy is modest** (single-task ~0.6–0.7 at K=A=3 with a shared
  head); the study is about *relative* interference, not ceiling performance.

## P2 — causal vs correlational credit (and the control that decides it)

**Result: an honest null.** Causal `do(θ)` credit does *not* reduce interference
beyond an effective-learning-rate effect. This is the outcome the plan flagged as
acceptable, and it is well-controlled.

**Setup.** The causal rule (weight-perturbation on the plastic couplings — an
interventional `do(w+ε)` estimator, à la Lansdell–Kording) only learns when
restricted to the low-dimensional **H→M readout** and at **K=2** (perturbing the
full plastic set, or K=3, is too high-variance to learn even one task — the known
weakness of perturbation methods, and the reason the low-variance hindsight
estimator exists). So P2 is run in the cleanest fair setting: **K=2 two-task
reversal** (identity → reversal), **frozen representation** so only the shared H→M
head adapts — both rules act on the *same* parameters, isolating the credit rule.
n=30.

**The apparent effect.** At the fixed default operating point, causal credit looks
like it halves forgetting:

| credit rule | backward transfer | new-task acq / task-0 retention |
|---|:--:|:--:|
| correlational (eligibility) | −0.778 [−0.856, −0.683] | 0.82 / 0.18 |
| causal `do(θ)` (perturbation) | −0.376 [−0.493, −0.261] | 0.71 / 0.29 |

Cohen d = 1.36, CIs cleanly separated. Taken alone this would "support" the
hypothesis.

**The control that kills it.** Sweeping each rule's plasticity knob (correlational
η, causal lr) traces a **stability–plasticity frontier** (retention vs new-task
acquisition). The two rules lie on the **same frontier**: causal @ lr=2 sits at
(acq 0.71, retention 0.29) exactly where a *slow* correlational learner (η≈0.03)
sits; correlational @ default η is just further along the same curve (acq 0.82,
retention 0.18). At **matched acquisition, retention is identical.**

![P2: apparent effect vs frontier control](figures/continual_p2.png)

**Verdict.** The apparent "causal forgets less" is entirely an *effective-learning-
rate* difference — causal weight-perturbation is simply slower, and any slow learner
forgets less. It does **not** assign credit in a way that interferes less at matched
plasticity. So in this setting the interference is **representational**: a shared
linear readout cannot simultaneously hold two anti-correlated mappings (a capacity
limit), and better credit assignment does not fix a capacity limit. (Causal is if
anything strictly worse — it caps at ~0.72 acquisition; higher lr just adds noise.)

## What this means for the two-arcs unification

The arcs connect *conceptually* — perturbation learning **is** interventional
causal-effect estimation, and the substrate's Line-B rule already performs a
`do(θ)` estimate. But **causal credit ≠ better continual learning here.** The
C-series prediction (well-posed `do(θ)` → better learning signal) does not transfer
to *interference resistance* on a capacity-limited shared readout. An honest, useful
negative: it rules out "just assign credit causally" as a continual-learning fix and
points at capacity/representation as the real lever.

## P3 — does *low-variance* causal credit move the frontier? (No.)

P2 left one escape hatch: maybe causal credit only failed because raw weight-
perturbation is *high-variance*. P3 tests a genuinely low-variance causal estimator
— **antithetic central-difference weight-perturbation** averaged over M=4 pairs
(`(r⁺ − r⁻)·ε`, baseline-free, the tractable stand-in for a Mesnard-style low-
variance estimator). It learns better than P2's single-sided rule (reaches
acquisition 0.79 vs the single-sided cap of ~0.72, no noise collapse). Re-run the
stability–plasticity frontier with all three rules (n=15):

![P3 frontier](figures/continual_p3.png)

**All three credit rules trace one frontier.** At matched new-task acquisition the
retention is identical (±0.02, n=15 noise): low-var causal @ acq 0.74 → retention
0.28 vs correlational-interpolated 0.26; @ 0.79 → 0.20 vs 0.21; @ 0.69 → 0.30 vs
0.31. The low-variance estimator reaches *higher* acquisition and is more robust,
but it **slides along the frontier rather than moving it.**

**Verdict — the P2 null was not a variance artifact.** The stability–plasticity
frontier is a genuine **capacity boundary**: no credit rule (correlational, causal,
or low-variance causal) beats it. Credit quality/variance sets *where on the
frontier* you operate, not *the frontier itself*. Continual-learning interference on
a shared readout is **representational**, and the lever is **capacity**, not credit.

## What closes the 3c arc

The learning-as-causal-inference unification is real *conceptually* (perturbation
learning is interventional estimation; the substrate's Line-B rule already performs
a `do(θ)` estimate) but, across P2 and P3, **causal credit — even low-variance
causal credit — does not buy interference resistance.** That is a clean, useful
negative: it redirects the "one homogeneous machine that genuinely learns
continually" question away from *credit assignment* and toward *representational
capacity* (WTA gating, per-task subspaces, conceptors — the CL literature's
representational fixes).

## Honest caveats

- **The task is a hard capacity case by design.** K=2 anti-correlated reversal
  fundamentally cannot be held on a shared linear head, so the frontier is close to
  a wall. Whether credit quality *ever* matters is best tested on **partially-
  overlapping** tasks that admit a coexisting solution — the sharper open question.
- **Antithetic ≠ Mesnard hindsight literally.** The hindsight (future-conditional)
  baseline is built for *multi-step* credit; in this single-step conditioning task
  there is little "future" to condition on, so the antithetic central-difference
  estimator is the appropriate low-variance test here. A temporally-extended task
  would be needed to exercise the hindsight mechanism specifically.
- **n=15 frontier**, 5 operating points per rule — enough to see the curves overlap,
  not to resolve ±0.02 differences.

## Deferred / next

- **Capacity, not credit.** Give the readout room (per-task subspace / higher-dim /
  conceptor-style / WTA gating) — the direction P2–P3 point at as the real lever.
- **Partially-overlapping tasks** where a coexisting solution exists — the fair test
  of whether credit quality ever matters.
- **Temporally-extended credit** — the setting where a true hindsight estimator
  would have something to condition on.
