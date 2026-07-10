# E8 Results — Predictive Dynamics on Naturalistic Tone Sequences

*Run of `experiments/e8_predictive.py`. First cut of the E-series predictive
extension (see [`predictive_dynamics_experiments.md`](predictive_dynamics_experiments.md)).
Asks whether the substrate makes **time-forward predictions** of upcoming input, and
how that prediction differs from predictive coding. Motivated by the auditory-sequence
literature (Chait/Bianco: neural activity integrates ~the last 7 tones to predict the
next pitch; anticipatory tuning sharpens for more predictable sequences).*

## Setup

A tonotopic map of `M = 8` tone channels. A per-channel **trace bank** on the GH
clock: presenting tone `c` sets trace `c` to phase 1, and the phase then advances
`1…τ` and wraps to 0 — so a tone is represented for ~`τ` steps after it occurs, the
phase encoding time-since-onset. The non-rested traces are therefore a **`τ`-deep
history window** whose depth is a `do(τ)` property (Line B). A next-tone **predictive
readout** (a GVF-style linear/ridge demon, self-supervised, no reward) reads the
*whole* trace state — heeding [C5](c5_results.md): read the medium's integrated state,
not a fixed locus — and predicts the upcoming tone. There is no error-unit hierarchy;
prediction is the forward readout of the medium's own dynamics.

## Result 1 — prediction tracks predictability; anticipation sharpens (E8.1/8.2)

Next-tone accuracy and anticipatory confidence (chance = 0.125), 6000-tone streams,
train/test split:

| sequence | next-tone accuracy | anticipatory confidence |
|----------|:------------------:|:-----------------------:|
| periodic | **1.00** | 0.98 |
| Markov α=0.9 | 0.90 | 0.97 |
| Markov α=0.5 | 0.55 | 0.74 |
| Markov α=0.1 | 0.20 | 0.21 |
| random-walk | 0.41 | 0.45 |

Prediction accuracy and the readout's **anticipatory confidence** both grade
monotonically with the sequence's predictability — the substrate anticipates the next
tone, and its pre-onset confidence *sharpens for more regular sequences*, the
signature reported for auditory cortex.

## Result 2 — the history window is a `do(τ)` property (E8.3, the causal headline)

Integration window = the deepest past-tone lag decodable from the trace state (from an
i.i.d. stream, so memory is isolated), swept over the trace timescale `τ`:

| τ | 4 | 8 | 14 | 20 | 26 |
|---|:-:|:-:|:--:|:--:|:--:|
| **window (tones)** | 0 | 1 | 3 | 5 | 7 |
| long-range (K=4) accuracy | 0.33 | 0.58 | 0.61 | 0.67 | 0.59 |

![E8 predictive](figures/e8_predictive.png)

The memory depth grows monotonically with `τ` — **`do(τ)` sets how many past tones
inform the prediction.** This turns the MEG *correlation* ("activity integrates ~7
tones") into an *intervention*: the "~7 tones" is not a fixed architectural depth but
a timescale parameter one can dial. On a long-range recurrence
(`x[t] = x[t−1] + x[t−K] mod M`, `K=4`) that needs history, accuracy rises from 0.33
(no window) to ~0.67 as the window opens — history availability helps prediction, up to
a representational ceiling (see caveats).

## Result 3 — surprise is a global scalar, not a per-feature error field (E8.4, PC disambiguation)

On an oddball stream (periodic standard, 12% random deviants), the prediction-error
magnitude `1 − p(actual)` (a single global scalar):

| | standard | deviant |
|---|:--------:|:-------:|
| surprise `1 − p(actual)` | **0.16** | **0.64** |

Surprise spikes ~4× on deviants — a genuine violation response — but it is a **single
aggregate mismatch signal**, not a hierarchy of per-channel error units. Combined with
Result 1's anticipatory sharpening, this is the predictive-coding disambiguation made
concrete: the substrate *predicts and is surprised*, yet has **no dedicated
prediction-vs-error pathway** — prediction is the forward evolution of one medium, and
the only error is a global TD/surprise scalar (the same `δ` the E-series critic uses).

## Interpretation — prediction without predictive coding

E8 shows the substrate is a genuinely *predictive* system that is not a
*predictive-coding* one:

- **Prediction = forward dynamics + history integration**, read off one medium, rather
  than the top-down output of an inverted generative model.
- **The integration window is a `do(τ)` (Line B) property** — the causal handle the
  [C-series](causal_experiments.md) isolates — not a hierarchy depth. This is the
  clean intervention the biological correlations cannot make, and it is the E-side
  echo of the C-series thesis: *the causally meaningful variable is the generative
  timescale.*
- **Surprise is global, and prediction is not dissociable from representation** (there
  is no separate top-down pathway to lesion) — two observables that distinguish this
  account from hierarchical predictive coding (Rao–Ballard/Friston).

It also carries [C5](c5_results.md)'s lesson across the series: the predictive readout
succeeds because it integrates the medium's *whole* trace state; a fixed-locus readout
(the fat-handed reader of C5) would be the wrong probe here too.

## Caveats / open items

- **Per-channel traces store recency, not positional order.** Each channel's trace
  records *how recently that tone last occurred*, so when tones repeat within the
  window the exact tone *at lag K* is not recoverable — which caps long-range recall
  (Result 2's ~0.67 plateau rather than a sharp threshold at `window ≥ K`). A richer
  recurrent reservoir that preserves order (or an explicit positional code) is the
  natural next step (E8.5), and would sharpen the long-range/nested-regularity results.
- **The learner is the predictive readout (a next-tone GVF), not Line A plasticity.**
  Transition structure is learned by the self-supervised readout; folding it into Line A
  conduction plasticity (so the *substrate* learns the transition graph) is the deferred
  variant named in the spec.
- **Prediction here is self-supervised (no reward)** — an intrinsic/inside-out result;
  the reward-driven variant and the E8.5 nested (two-timescale) task are not yet run.
- The trace bank is a deliberately minimal reservoir (independent clocked nodes); the
  qualitative results (predictability grading, `do(τ)` window, global surprise) are the
  robust claims, not the absolute accuracies.

## Operating point

```
map        : M=8 tone channels; ISI=3 steps between tones; act=2
traces     : one GH-clocked node per channel, tau swept in {4,8,14,20,26} (do(tau));
             theta=99 (only clamped), p_s=0, no edges
readout    : ridge next-tone demon on [trace recency (M) + current-tone one-hot (M) + bias];
             softmax temperature beta=6 for confidence/surprise; train/test split 0.5
sequences  : periodic (period 4); Markov (neighbour-bias alpha); random-walk (+-1);
             long-range recurrence x[t]=x[t-1]+x[t-K] mod M (K=4); oddball (12% deviants)
streams    : 6000 tones; window measured on i.i.d. streams; 3-seed average for E8.3
```

## Reproduce

```
python3 experiments/e8_predictive.py
```

Writes `docs/figures/e8_predictive.png` and `result/e8/e8_data.npz`.
