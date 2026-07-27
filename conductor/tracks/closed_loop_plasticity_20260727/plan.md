# Track Implementation Plan: Closed-Loop Substrate Plasticity

## Phase 1: Engine Architecture & Unit Testing (Days 1–3)
- [x] Task: Design and implement local multi-axis update engine ($\tau$-adaptation, $\theta$-homeostasis, $W$-routing) in `ghca_plasticity.py`.
- [x] Task: Create unit test suite `tests/test_closed_loop_plasticity.py` checking value bounds, event-timing triggers, and RNG bit-identity across seeds.
- [x] Task: Phase Verification & Checkpoint — Run `python3 -m unittest discover -s tests` and `python3 .claude/skills/experiment-review/review_helper.py audit-rng`.

## Phase 2: Single-Task Self-Organization & Readout Independence (Days 4–7)
- [x] Task: Implement single-task benchmark script in `experiments/closed_loop_plasticity.py` targeting E1 sensorimotor routing and E5 contextual switching across $n=30$ seeds.
- [x] Task: Compute Readout Independence Ratio (RIR) comparing fixed/identity readout vs. trained linear readout.
- [x] Task: Phase Verification & Checkpoint — Verify RIR $> 0.85$ and confirm per-seed spreads/bimodality.

## Phase 3: Sequential Learning & Anti-Forgetting Protocol (Days 8–11)
- [x] Task: Implement sequential learning benchmark (Task A $\to$ Task B $\to$ Task A) without task context heads or weight freezing.
- [x] Task: Measure Task A retention after Task B learning and track $\tau_i$ distribution shifts (topological loop protection).
- [x] Task: Phase Verification & Checkpoint — Verify retention $> 80\%$ and compare against standard MLP/Reservoir baselines.

## Phase 4: Documentation, Verification, & Integration (Days 12–14)
- [x] Task: Generate publication figures (`figures/closed_loop_plasticity.png`) and write comprehensive results doc in `docs/closed_loop_plasticity_results.md`.
- [x] Task: Add benchmark runner to `reproduce_all.py` and run full regression suite.
- [x] Task: Phase Verification & Checkpoint — Final RNG audit and Conductor/Beads task sign-off.
