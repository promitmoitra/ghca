# Track Specification: Closed-Loop Substrate Plasticity

## Overview
This track implements unified multi-axis local plasticity ($\tau$-adaptation for refractory timescales, $\theta$-adaptation for homeostatic excitability, $W_{\text{rec}}$ for Hebbian structural routing) on Greenberg-Hastings excitable graphs. The goal is to determine whether closed-loop local plasticity allows the physical substrate itself to learn and retain task representations without relying on trained linear readout weights or external task-context heads.

## Objectives & Research Questions
1. **Substrate vs. Readout Credit Assignment (Track 1a):** Can local $\tau,\theta,W$ update rules self-organize the excitable graph to solve reward routing tasks (E1/E5) with a fixed/identity readout, achieving a Readout Independence Ratio (RIR) $> 0.85$?
2. **Catastrophic Forgetting Elimination (Track 1b):** Does local $\tau$-adaptation topologically protect learned reentrant loops during sequential learning (Task A $\to$ Task B $\to$ Task A) without task-head masks or frozen representations?

## Requirements
- **RNG Seeding:** All plasticity engines and network simulations must explicitly thread `default_rng(seed)`.
- **Local Event-Driven Update Engine:** Updates to $\tau_i, \theta_i, W_{ij}$ must depend only on local node states $s_i(t)$, local event timing, and a global scalar reward signal $R(t)$.
- **Testing & Benchmarking:**
  - `tests/test_closed_loop_plasticity.py` unit test suite verifying update invariants, bounds, and seed consistency.
  - `experiments/closed_loop_plasticity.py` experiment suite generating quantitative results across $n=30$ seeds.
  - `docs/closed_loop_plasticity_results.md` reporting per-seed spreads, RIR metrics, retention scores, and caveats.

## Acceptance Criteria
- [ ] Local update rules implemented cleanly in `ghca_main.py` or `experiments/closed_loop_plasticity.py`.
- [ ] Unit tests pass with 100% seed reproducibility.
- [ ] Single-task performance achieves RIR $> 0.85$ across $n=30$ seeds.
- [ ] Sequential learning maintains $> 80\%$ Task A retention after Task B acquisition.
- [ ] RNG audit cleanly passes via `review_helper.py audit-rng`.
