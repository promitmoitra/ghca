# Research Directions & Brainstorming

This document outlines the current state of research in **`ghca`**, summarizes directions pursued across the E-series and C-series, and proposes concrete next steps evaluated through structured ideation frameworks (*Tension Hunting, Cross-Pollination, Composition, and "What Changed?"*).

---

## 1. Summary of Pursued Research Directions

The work in `ghca` sits at the intersection of excitable media dynamics, reward-modulated substrate adaptation, and causal inference:

```
                                  ghca RESEARCH MAP
                                          │
   ┌───────────────────┬──────────────────┼───────────────────┬──────────────────┐
   ▼                   ▼                  ▼                   ▼                  ▼
[Track 1]           [Track 2]          [Track 3]           [Track 4]          [Track 5]
Substrate &         Closed-Loop        Cognitive           Causal Spike-      Timescale
Excitable           Substrate          Readouts            Wave Duality       Hierarchy &
Dynamics            Plasticity         (E4-E8)             (C0-C7)            Plastic Dynamics
(E0, Topologies)    (E1-E3, E9)                            (Testbed)          (3d-3e, 4a)
```

### Track 1: Substrate & Excitable Media Dynamics (E0, Topologies)
* **Scope:** Characterized the operating point of Greenberg–Hastings (GH) dynamics across lattices, small-world networks, and random geometric graphs (RGG). Identified the self-sustaining spiral wave band via threshold-range scaling ($r=2, a=6, \tau=14, \theta\approx 4$) and validated the exact winding number invariant ($G, G, H \text{ 1980}$) for reentrant loop persistence.
* **Key Finding:** The dominant loop period tracks local timescales almost perfectly ($\text{period} = 1.00\tau + 0.95, r = 0.9992$). Excitable dynamics and reward routing generalize off-lattice to small-world and RGG topologies ([`docs/e0_topologies.md`](e0_topologies.md)).

### Track 2: Closed-Loop Substrate Learning & Plasticity (E1–E3, E9, Multi-Axis Engine)
* **Scope:** Isolated two parallel plasticity lines: **Line A** (conduction weights $W$, spatial routing) and **Line B** (per-node timescales $\tau$, temporal credit). Built a multi-axis closed-loop engine ($\tau$-adaptation, $\theta$-homeostasis, $W$-routing) and developed self-organized conjunction cells via competitive Hebbian learning + DeSieno conscience (E9, [`docs/e9_results.md`](e9_results.md)).
* **Key Finding:** Reversal learning on a plastic shared substrate achieves **70.0% ± 13.8% retention** (vs 29.6% ± 31.2% weight-only) via topological loop protection ([`docs/closed_loop_plasticity_results.md`](closed_loop_plasticity_results.md)).

### Track 3: Cognitive Faculties as Substrate Readouts (E4–E8)
* **Scope:** Demonstrated that core cognitive faculties emerge as questions asked of one homogeneous substrate:
  * **E4 Attention:** Selective attention via wave annihilation (**zero inhibitory nodes required**).
  * **E5/E7 Executive Control:** Persistent reentrant loops / 2D spiral cores acting as *options* $\langle I, \pi, \beta \rangle$ that gate fast routing across task blocks.
  * **E6 Emergent Categories:** General Value Function (GVF) Horde demons reading one phase stream to predict distinct cognitive questions near-orthogonally.
  * **E8 Predictive Processing:** Time-forward predictions via passive-medium $do(\tau)$ history windows and global surprise.

### Track 4: Causal Formulation of Spike–Wave Duality (C0–C7, Testbed)
* **Scope:** Used the substrate as a synthetic-SCM testbed for the Jalaldoust & Zabeh framework ($\text{arXiv:2511.06602}$). Verified Definition 1 epiphenomenality on ground-truth graphs (C1), proved that macro aggregate interventions $do(W)$ under constitution $W=f(S)$ are fat-handed ($33\sigma$ ambiguity), and showed that generating parameters $do(\theta)$ are the unique well-posed causal handles ($0.014\sigma$ ambiguity) (C2–C4). Packaged as a standalone Python benchmark [`ghca_testbed`](../ghca_testbed).
* **Core Synthesis:** The variable the learner adapts ($\theta$) is precisely the variable on which intervention is causally well-defined ([`docs/synthesis.md`](synthesis.md)).

### Track 5: Timescale Hierarchy, Cross-Frequency Coupling & Plastic Dynamics (3d–3e, Track 4a)
* **Scope:** Escaped the Line-B ratchet by introducing an input-timing-driven local $\tau$-plasticity rule. Grown an emergent fast/slow timescale hierarchy under two-rhythm drive (3e.2) and demonstrated theta–gamma-style phase-amplitude cross-frequency coupling ($\text{MI} = 0.594$) on the emergent hierarchy (3e.2b). Shown that a travelling wave turns into a synchronous population recovery burst at the learned stimulus-reward interval (3e.11, [`docs/timescale_hierarchy_results.md`](timescale_hierarchy_results.md)).

---

## 2. Proposed Next Research Steps

### Direction 1: The Composite Multi-Capacity Task (Memory × Attention × Control)
* **Ideation Lens:** *Composition & Decomposition (Framework 9)*
* **Core Motivation:** Every E-series capacity (E2 memory, E4 attention, E5 control) was established in isolation. We do not know what happens when a single task demands all three simultaneously on one substrate. Do they compose smoothly, or do they contend for shared topological cycle space?
* **Two-Sentence Pitch:**
  > *Domain:* Current excitable substrate models evaluate cognitive faculties in isolation, ignoring resource contention on shared media.  
  > *Insight:* We design a composite task requiring simultaneous stream selection (E4), delay maintenance (E2), and rule switching (E5), isolating whether independent capacity constraints compose modularly or collide on topological cycle limits $K_{\text{dyn}}$.
* **Concrete Plan:**
  1. Build a delayed match-to-sample trial where the target stimulus must be selected from competing distractor waves (E4) and held across a variable delay (E2) under a block-switched rule (E5).
  2. Measure performance against the cycle-space bound ($\beta_1 = m - N + c$) from [`docs/topology_cycle_capacity.md`](topology_cycle_capacity.md).
  3. Compare a **monolithic substrate** against a **modular partitioned substrate** with dedicated inter-area pathways at matched total $N$.

---

### Direction 2: Empirical Data Validation Pipeline (Testing Cortical Wave Predictions P1–P6)
* **Ideation Lens:** *The "What Changed?" Principle (Framework 5)* & *Bridge to Data*
* **Core Motivation:** High-density Neuropixels and widefield voltage/calcium imaging datasets in awake, behaving animals (e.g., Liang & Gong 2023, Steinmetz/Ye 2026) are now publicly available. The 6 falsifiable predictions derived in [`docs/spiral_predictions.md`](spiral_predictions.md) can be tested directly on real cortex.
* **Two-Sentence Pitch:**
  > *Domain:* Cortical travelling waves are frequently recorded in vivo, but whether phase singular cores actively gate task routing or are passive epiphenomena remains debated.  
  > *Insight:* We apply our topology-aware phase-singularity tracker to open primate/mouse visual cortex recordings to test Prediction P4: whether fixed-ROI decoders collapse under natural core drift while core-relative topological decoders maintain rule decoding.
* **Concrete Plan:**
  1. Download open widefield voltage/LFP recordings from visual/parietal cortex during visual discrimination/switching tasks.
  2. Implement an automated core-tracking pipeline (`p2b_signature_figure.py` port) to detect phase singularities $\chi(t) = (x_c, y_c, \text{chirality})$.
  3. Compare stimulus/rule decoding accuracy between a fixed-grid spatial ROI and a core-centered frame of reference.

---

### Direction 3: Physical Neuromorphic Substrate Mapping (Event-Driven Wave Hardware)
* **Ideation Lens:** *Cross-Pollination & Solution-First (Framework 4)*
* **Core Motivation:** Standard GPUs simulate GH cellular automata sequentially, but event-driven neuromorphic hardware (e.g., SpiNNaker, Intel Loihi 2, or mixed-signal analog arrays) natively executes local propagation delays, refractory periods, and wave collisions in parallel with near-zero energy consumption.
* **Two-Sentence Pitch:**
  > *Domain:* Simulating large-scale excitable wave dynamics on synchronous digital hardware incurs significant computational overhead for sequential delay updates.  
  > *Insight:* We map the Greenberg–Hastings $\tau$-plasticity and conduction routing rules onto event-driven neuromorphic architectures, leveraging physical axonal transmission delays and asynchronous spike refractoriness as native computing primitives.
* **Concrete Plan:**
  1. Formulate the GH state transition equations ($\phi_i \in \{0, 1..a, a+1..a+p_i\}$) in PyNN / Lava for Loihi 2 or SpiNNaker.
  2. Implement local $\tau$-adaptation via event-driven spike-delay modulation.
  3. Benchmark energy efficiency, execution latency, and scale limits ($N=10^5$ nodes) against the CPU implementation.

---

### Direction 4: Autonomous Motor Closed-Loop (World-Changing Action Transmissivity)
* **Ideation Lens:** *Problem-First & Simplicity Test (Framework 1 & 7)*
* **Core Motivation:** In direction 3e.11–13, the substrate learned to align recovery times into a synchronous population burst at the predicted reward interval. However, the action was passive (an observer counted recovery events). To become a true reinforcement learning agent, the substrate's motor transmission must **change the environment state**.
* **Two-Sentence Pitch:**
  > *Domain:* Current substrate learning models evaluate reward-timed responses using passive readouts that do not affect environment dynamics.  
  > *Insight:* We close the physical sensorimotor loop by requiring the substrate's synchronous transmission window to actively gate environmental feedback, where a blocked transmission permanently forfeits reward and alters future input streams.
* **Concrete Plan:**
  1. Connect the motor transmission layer to a dynamic environment (e.g., an obstacle-avoidance or timing-gate task).
  2. Make reward delivery strictly contingent on event transmission within the excitable window $[t_{\text{rec}} - \epsilon, t_{\text{rec}} + \epsilon]$.
  3. Evaluate whether closed-loop environmental feedback accelerates local $\tau$-retuning compared to open-loop reward signals.

---

### Direction 5: Bioelectric Substrate & Morphogenetic Pattern Repair (Basal Cognition)
* **Ideation Lens:** *Tension Hunting & Abstraction Ladder (Framework 2 & 3)*
* **Core Motivation:** In biology, non-neural bioelectric cellular networks (Levin et al., 2021/2023) use excitable wave dynamics and gap-junction couplings to store spatial anatomical target states and repair morphological damage without a central nervous system.
* **Two-Sentence Pitch:**
  > *Domain:* Most artificial substrate learning focuses strictly on neural I/O, ignoring how non-neural cellular media store, maintain, and repair spatial pattern memories.  
  > *Insight:* We extend $do(\theta)$ parameter plasticity to bioelectric graph substrates, demonstrating that local homeostatic thresholds and gap-junction couplings can store a global target geometry and dynamically regenerate it following arbitrary structural lesions.
* **Concrete Plan:**
  1. Represent a 2D bioelectric tissue sheet with local gap-junction conductances $w_{ij}$ and resting potentials $\theta_i$.
  2. Train the medium via $do(\theta)$ closed-loop plasticity to maintain a stable bioelectric spatial pattern (e.g., a target anatomical voltage map).
  3. Sever or lesion $20\text{--}50\%$ of the graph (simulating injury) and measure pattern restoration time and fidelity.

---

## 3. Decision & Ranking Matrix

| Direction | Primary Lens | Feasibility / Effort | Scientific Impact | Addresses Caveat / Gap |
| :--- | :--- | :--- | :--- | :--- |
| **1. Composite Task (Memory × Attention × Control)** | Composition (F9) | **High** (Days) | **High** | Tests multi-capacity resource contention |
| **2. Empirical Data Validation (P1–P6)** | "What Changed?" (F5) | **Medium** (Weeks) | **Very High** | Bridges toy model to in vivo cortical waves |
| **3. Neuromorphic Hardware Mapping** | Cross-Pollination (F4) | **Medium** (Weeks) | **High** | Demonstrates real hardware energy gains |
| **4. Autonomous Motor Closed-Loop** | Problem-First (F1) | **High** (Days) | **Medium–High** | Closes true action-perception loop |
| **5. Bioelectric Pattern Repair** | Abstraction Ladder (F2) | **Medium** (Weeks) | **High** | Expands framing to non-neural substrate computing |
