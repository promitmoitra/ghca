# E7 Results — The Option as a 2D Spiral Core (Phase A: mechanism)

*Run of `experiments/e7_spiral_option.py`. E7 extends the program beyond the
original E0–E6 plan, prompted by the rotating-wave neuroscience (Xu/Gong et al.,
*Nat. Hum. Behav.* 2023; Ye/Steinmetz et al., *Science* 2026): cortical activity
rides rotating waves whose **rotation direction is task-relevant**. E5 realised
the executive-control "option" as a persistent 1-D ring; E7 moves it onto the
substrate's **native 2-D medium** — a genuine spiral wave with a phase-singularity
core — and asks whether rotation direction is a controllable, persistent, readable
variable here. **This doc covers Phase A (the mechanism).** Phase B (rotation
direction as the learned rule, reusing E5's routing) and the C-series causal test
are the sequel.*

## Setup

The lattice from [E0](e0_results.md)'s organised spiral band — `lattice2d`,
`r = 2`, `a = 6`, `τ = 14`, `θ = 4`, `p_s = 0` — with **no-flux (non-periodic)
boundaries**. The boundary choice is not cosmetic: on a periodic torus the total
topological charge must be zero, so a lone spiral necessarily breeds a
compensating anti-spiral; no-flux boundaries let a single signed core persist (a
concrete echo of Ye et al.'s point that boundaries / anatomy shape the wave).

A spiral of chosen handedness is nucleated by imposing a polar-angle phase ramp
around a core (`φ ∝ chirality · atan2(y−c_y, x−c_x)`); the sign of the ramp sets
CW vs CCW. Rotation direction is read out two ways: **globally** as the sign of the
net topological charge (a signed phase-singularity count over 2×2 plaquettes), and
**locally** from two phase probes placed 90° apart around the (tracked) core, where
which probe leads in phase flips with handedness.

## Result 1 — nucleation and persistence

| condition | net charge (mean, t ≥ 60) | net charge at t = 399 |
|-----------|:-------------------------:|:---------------------:|
| CCW seed (+1) | **+0.89** | **+1** |
| CW seed (−1) | **−1.11** | **−1** |
| planar (no core, control) | +0.00 | 0 |

![E7 mechanism](figures/e7_mechanism.png)

A single spiral of the seeded handedness nucleates and **persists for the full 400
steps**, holding net charge +1 (CCW) or −1 (CW); the planar-wave control carries no
core and ~zero net charge, confirming the signal is spiral-specific. (The mean net
charge is not exactly ±1 because of brief transient secondary defects the detector
occasionally registers, but the *sign* — the handedness — is stable throughout.)
This persistence is the causal lever the learning phase will use, the 2-D analogue
of E2's "memory duration is τ-controlled" mechanism result.

## Result 2 — rotation direction is readable

Chirality recovered over 20 trials (varying core position ±6 and adding init
jitter):

| readout | accuracy |
|---------|:--------:|
| net topological charge (global) | **1.00** |
| phase-probe lead (local, core-tracked) | **0.90** |

The global charge readout is perfect. The **local** phase-probe readout — the
on-thesis one, where the rotating wave itself delivers a binary context signal with
no explicit charge computation — reaches 0.90 once the probes are anchored to the
*tracked* core rather than the seeded position (0.70 → 0.90), because the spiral tip
meanders. This is exactly the quantity Phase B will feed into routing, and the
mechanistic counterpart of the fMRI finding that rotation direction classifies the
task.

## Interpretation

The object the rotating-wave papers describe — a persistent rotating wave whose
*direction* is a state variable — exists in our substrate and is controllable
(seed the handedness), persistent (holds across a block), and readable (globally
exact; locally 0.90). E7 thus lifts E5's "option" from a hand-built 1-D ring to a
genuine 2-D spiral core on the medium E0 characterised, closing the gap between the
learning program and the spiral-wave neuroscience. It also sets up the sharpest
version of the field's open controversy: once rotation direction *drives* routing
(Phase B), the C-series machinery ([`synthesis.md`](synthesis.md)) can ask whether
the chirality is **causal** (`do(chirality)`) or merely epiphenomenal.

## Caveats / open items

- **Phase A only.** This establishes the mechanism; it does *not* yet show learning.
  Phase B — CW/CCW as the identity-vs-reversal rule, read via the phase probe and
  routed by E5's conjunction-gate + reward-driven Line A, with the switching /
  discriminator / rotation-decodes-the-task results — is the sequel.
- **Tip meander.** The spiral core drifts; the robust readout is the global charge
  (1.00), while the local probe (0.90) needs core-tracking. Phase B's fixed
  downstream sites will see a chirality-dependent phase relationship but with
  meander noise, so mild core-pinning (a small central heterogeneity) may be
  warranted.
- **Operating point is narrow.** Persistence needs the E0 band (θ ≈ 4); θ ≥ 5 lets
  the medium die and lower θ turns turbulent. This is the same narrow organised band
  E0 flagged.
- **Deterministic run** (`p_s = 0`); trial variation comes from core position and a
  small init jitter, not spontaneous dynamics.

## Operating point

```
substrate : lattice2d L=48, range r=2, act=6, tau=14, theta=4.0, no-flux (periodic=False), p_s=0
nucleation: polar-angle phase ramp about the core; sign of ramp = CW / CCW
readout   : (global) sign of net topological charge; (local) phase-probe lead 90deg
            apart around the tracked core, radius 10
run       : 400 steps; readout accuracy over 20 trials (core jitter +/-6, init jitter)
```

## Reproduce

```
python3 experiments/e7_spiral_option.py
```

Writes `docs/figures/e7_mechanism.png` and `result/e7/e7_mechanism.npz`.
