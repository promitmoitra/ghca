# Specification: Closed-Loop Substrate Plasticity Extensions

## 1. Overview
This track extends the tri-axis closed-loop substrate plasticity engine ($\tau$-adaptation, $\theta$-homeostasis, $W$-routing) on Greenberg–Hastings excitable graphs. Building on initial findings in [docs/closed_loop_plasticity_results.md](file:///home/dognosis/Documents/ghca/docs/closed_loop_plasticity_results.md), this work investigates three major frontiers:
1. **Long-Horizon Sequential Saturation ($K \ge 5$ Tasks):** Measuring capacity scaling laws and memory interference ceilings.
2. **Axis $G$ Structural Plasticity:** Implementing reward-modulated physical edge pruning and sprouting to achieve topological circuit partitioning.
3. **Homeostatic $\tau$-Relaxation & Consolidation:** Introducing offline $\tau$-decay to recycle refractoriness without unlearning protected memory traces.

## 2. Core Functional Requirements

### Cluster 1: Long-Horizon Sequential Task Saturation
- **Task Sequence Engine:** Support sequential execution of $K \ge 5$ orthogonal or reversal task mappings ($A \to B \to C \to D \to E \to A$) without network resets or task context cues.
- **Substrate Scaling Sweep:** Benchmark across graph sizes $N \in \{30, 50, 100, 200\}$ and graph circuit ranks $\beta_1 = m - N + c$.
- **Metrics:** Compute Readout Independence Ratio ($RIR$), Task A retention %, backward transfer ($BWT$), and seed variance collapse across $n=30$ seeds (`default_rng(seed)`).

### Cluster 2: Axis $G$ Structural Plasticity (Reward-Modulated Rewiring)
- **Plasticity Engine:** Extend `ghca_plasticity.py` to support physical edge pruning for sub-threshold connections ($W_{ij} < w_{\text{thresh}}$) and co-activity sprouting ($a_i \cdot a_j$) modulated by global scalar reward error $\delta$.
- **Graph Updates:** Dynamically update network adjacency matrix and conduction delay structures during runtime.
- **Metrics:** Evaluate graph modularity $Q$, clustering coefficient, and $RIR$ gain vs. weight-only baseline on E1/E5 benchmarks ($n=30$ seeds).

### Cluster 3: Homeostatic $\tau$-Relaxation & Offline Consolidation
- **Consolidation Mechanics:** Implement slow $\tau$-decay during unrewarded sleep/consolidation intervals:
  $$\frac{d\tau_i}{dt} = -\lambda_\tau (\tau_i - \tau_{\text{base}}) \cdot (1 - \tau_{\text{prot}})$$
- **Memory Retention:** Test if offline consolidation recycles high-$\tau$ refractoriness on non-critical nodes while preserving metaplastically protected pathways.
- **Integration:** Incorporate automated regression tests into `reproduce_all.py` and document results in `docs/closed_loop_plasticity_extensions.md`.

## 3. Non-Functional & Reproducibility Requirements
- **100% Seeded RNG Compliance:** All simulations MUST use explicit `default_rng(seed)` instances. Zero calls to global NumPy RNG.
- **Substrate/Analysis Boundary:** Physical learning MUST occur strictly within substrate graph dynamics; decoders are used solely for metric logging ($RIR$).
- **Statistical Rigor:** Report per-seed standard deviations ($\sigma$) and check for bimodality across $n=30$ seeds.

## 4. Acceptance Criteria
- [ ] $K \ge 5$ sequential task benchmark runs cleanly across $n=30$ seeds for $N \in \{30, 50, 100, 200\}$.
- [ ] Axis $G$ structural plasticity engine passes unit tests (`tests/test_axis_g_plasticity.py`) and RNG compliance audit.
- [ ] Offline $\tau$-relaxation mechanism demonstrates non-destructive refractoriness recycling without task retention loss.
- [ ] All new benchmarks are registered in `reproduce_all.py` and verified bit-identical across runs.
