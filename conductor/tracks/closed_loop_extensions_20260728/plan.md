# Implementation Plan: Closed-Loop Substrate Plasticity Extensions

## Phase 1: Long-Horizon Sequential Task Saturation ($K \ge 5$)
- [x] Task: Implement $K$-task sequential learning engine in `experiments/sequential_k_tasks.py` supporting $K \ge 5$ tasks ($A \to B \to C \to D \to E \to A$).
- [x] Task: Execute substrate scaling sweep across $N \in \{30, 50, 100, 200\}$ with $n=30$ seeds per condition (`default_rng(seed)`).
- [x] Task: Compute retention, backward transfer ($BWT$), and capacity bounds relative to circuit rank $\beta_1$.
- [x] Task: Phase Verification & Checkpoint — Verify RNG compliance, confirm per-seed spreads/bimodality, and archive results in `result/closed_loop_plasticity/sequential_k_tasks.npz`.

## Phase 2: Axis $G$ Structural Plasticity (Reward-Modulated Rewiring)
- [x] Task: Extend `ghca_plasticity.py` to support Axis $G$ structural plasticity (edge pruning for sub-threshold $W_{ij}$ and co-activity sprouting).
- [x] Task: Create unit test suite in `tests/test_axis_g_plasticity.py` verifying adjacency matrix updates, weight zeroing, and RNG bit-identity.
- [x] Task: Implement single-task and sequential benchmarks with Axis $G$ in `experiments/closed_loop_structural.py` across $n=30$ seeds.
- [x] Task: Phase Verification & Checkpoint — Verify graph modularity $Q$ increase and check for $RIR$ performance gains.

## Phase 3: Homeostatic $\tau$-Relaxation & Offline Consolidation
- [x] Task: Implement offline $\tau$-relaxation rule in `ghca_plasticity.py` allowing slow refractoriness decay on non-protected pathways.
- [x] Task: Implement consolidation benchmark in `experiments/closed_loop_consolidation.py` testing Task A $\to$ Consolidation Phase $\to$ Task B $\to$ Task A protocol across $n=30$ seeds.
- [x] Task: Phase Verification & Checkpoint — Confirm $\tau$ recycling without Task A retention loss ($\ge 70\%$).

## Phase 4: Verification, Documentation, & Integration
- [x] Task: Write comprehensive extension results document in `docs/closed_loop_plasticity_extensions.md` with publication figures.
- [x] Task: Register new benchmark scripts and verification checks in `reproduce_all.py`.
- [x] Task: Phase Verification & Checkpoint — Final RNG audit (`review_helper.py audit-rng`) and Conductor/Beads task sign-off.

## Phase 5: Remediation & Metric Corrections (TDD Workflow)
- [x] Task: [TDD] Write unit test suite `tests/test_retention_metrics.py` testing chance-corrected retention $R_{\text{chance}} = \frac{\text{acc} - 0.50}{\text{init} - 0.50}$ and motor drive activation.
- [x] Task: Implement chance-corrected retention and task acquisition tracking in `experiments/sequential_k_tasks.py` (Beads: `ghca-6lm.1`).
- [x] Task: Fix sensory/motor drive parameters in `experiments/closed_loop_structural.py` and `experiments/closed_loop_consolidation.py` to ensure Task 1 initial accuracy $> 0.70$ (Beads: `ghca-6lm.2`).
- [x] Task: Create `scripts/print_extensions_table.py` to auto-generate markdown tables from `.npz` archives and update `docs/closed_loop_plasticity_extensions.md` (Beads: `ghca-6lm.4`).
- [x] Task: Re-run full $n=30$ sweeps across Phase 1, Phase 2, and Phase 3, update `reproduce_all.py` assertions, and run RNG audit (Beads: `ghca-6lm.5`).
