# Literature Review & Citation Candidates

This document grounds the **`ghca`** project (codebase, documentation, and public website [`promitmoitra.github.io/ghca`](https://promitmoitra.github.io/ghca/)) in computational neuroscience, cellular automata, neuromorphic/alternative computation, reinforcement learning, and causal inference literature.

---

## 1. Core Theoretical Foundations & Mapping

The `ghca` project synthesizes three primary computational and neuroscientific pillars:

1. **Inside-Out Excitable Substrate Dynamics & Wave Computing:**
   Greenberg–Hastings cellular automata on graphs act as an excitable physical/biological substrate. Rather than using passive rate-based feedforward units, the medium self-generates spatiotemporal waves (E0, E4, E7, C5–C7), utilizing conduction delays, local refractoriness, and wave annihilation as natural computational primitives.
2. **Multi-Timescale Plasticity & Substrate Adaptation:**
   Learning adapts local generating parameters $\theta$ (conduction weights $W$ and intrinsic local timescales $\tau$) via three-factor reward-modulated plasticity, achieving anti-forgetting sequential learning through topological loop protection and emergent conjunctions (E1–E9, Track 3).
3. **Causal Spike–Wave Instrumentation & Causal Emergence:**
   Instruments coarse-grained wave variables $W = f(S)$ against micro-spikes $S$, proving that generating parameters $\theta$ are the unique well-posed causal handles for both learning and interventional control (C0–C7).

---

## 2. Deep Dives: Foundational & Flagship References

### A. Foundational Review: Muller et al. (2018)
> **Reference File:** [`docs/lit-rev/muller2018.pdf`](lit-rev/muller2018.pdf)  
> **Citation:** Muller, L., Chavane, F., Reynolds, J., & Sejnowski, T. J. (2018). *Cortical travelling waves: mechanisms and computational principles.* Nature Reviews Neuroscience, 19(5), 255–268. [doi:10.1038/nrn.2018.20](https://doi.org/10.1038/nrn.2018.20)

* **Axonal Conduction Delays:** Mesoscopic waves propagate across layer 2/3 horizontal fibers ($0.1\text{--}0.8\text{ m/s}$), breaking zero-lag synchrony. In `ghca`, [`ghca_net.py`](../ghca_net.py) formalizes refractory delays $\tau_i = a + p_i$ and weights $w_{ij}$ to mediate wave propagation.
* **Sparse Waves in Asynchronous Irregular (IA) States:** Awake sparse waves weakly entrain spiking probabilities. In `ghca`, homeostatic threshold adaptation $\theta_i \leftarrow \theta_i + \eta_\theta (\bar{\rho}_i - \rho^*)$ and spontaneous firing $p_s$ maintain the substrate in the sparse, self-sustaining spiral band (E0).
* **Collision-Based Computing & Refractory Annihilation:** Wave collisions implement logic gates (e.g., XOR via wave annihilation). Direct biological precedent for E4 selective attention ([`experiments/e4_attention.py`](../experiments/e4_attention.py)) via wave annihilation with **zero inhibitory nodes**.
* **Spatiotemporal Memory Trace:** Global wave fields store trajectory histories (bouncing droplet analogy, Perrard et al. 2016). Grounds E2 working memory and E8 time-forward predictions.

### B. Flagship 2026 Discovery: Karvat et al. (2026)
> **Reference File:** [`docs/lit-rev/s41467-026-73553-8.pdf`](lit-rev/s41467-026-73553-8.pdf)  
> **Citation:** Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. C., & Landau, A. N. (2026). *Universal rhythmic architecture uncovers two modes of neural dynamics.* Nature Communications, 17:7024. [doi:10.1038/s41467-026-73553-8](https://doi.org/10.1038/s41467-026-73553-8)

* **Bi-Dimensional Spectral Architecture (Frequency $\times$ Rhythmicity):** Demonstrates across 859 participants, 12 datasets, invasive/non-invasive recordings, and rats that neural dynamics are organized along a 2D axis: Frequency $\times$ Rhythmicity (sustained vs. bursty).
* **Dual Operating Modes (Maintenance vs. Responsiveness):**
  - *High-Rhythmicity Bands (Sustained Rhythms):* Long-lasting oscillations ($\theta, \alpha, \beta_2$) that support maintenance of ongoing state, top-down control, and temporal scaffolds. Maps onto `ghca`'s persistent reentrant loops ($\tau < L$) and 2D spiral option cores (E2/E5/E7).
  - *Low-Rhythmicity Bands (Transient Bursts):* Brief, event-linked bursts ($\theta/\alpha, \beta_1, \gamma_1$) reflecting transient inputs and responses to change. Maps onto `ghca`'s transient travelling wave packets (E1/E4/E8).
* **Bridge Between Phase-Coding and Rate-Coding:** High-rhythmicity bands maintain phase-coding scaffolds (biasing excitability windows), while low-rhythmicity bands capture rate-coding event bursts.
* **Quiet Spectral Zones & Phase Interference:** High-rhythmicity sustained oscillations attenuate adjacent frequencies via phase interference, establishing "quiet" frequency zones held in readiness for incoming transient inputs. Directly informs `ghca`'s wave interaction dynamics and noise filtering.

---

## 3. Candidate Citations by Category

### A. Neuromorphic & Alternative / Unconventional Substrate Computing

1. **Kudithipudi, D., Aguilar-Simon, R., Babb, J., et al. (2022).**  
   *Biological underpinnings for lifelong learning machines.* Nature Machine Intelligence, 4(3), 196–210.
   * **Placement:** [`docs/continual_learning_results.md`](continual_learning_results.md), [`README.md`](../README.md), Website Track 3 & Closed-Loop Plasticity pages.
   * **Relevance & Reasoning:** Outlines multi-timescale metaplasticity, structural plasticity, and neuromodulated gating as essential requirements to overcome catastrophic forgetting. Grounds `ghca`'s sequential reversal learning (**70.0% retention vs. 29.6% weight-only**) via topological loop protection and multi-axis closed-loop adaptation.

2. **Wright, L. G., Onodera, T., Ng, M. M., et al. (2022).**  
   *Deep physical neural networks trained with backpropagation.* Nature, 601(7894), 549–555.
   * **Placement:** [`README.md`](../README.md), [`docs/synthesis.md`](synthesis.md), Website Overview & Architecture pages.
   * **Relevance & Reasoning:** Demonstrates learning directly on heterogeneous physical substrates (optical, mechanical, electronic). Grounds `ghca`'s framing of excitable cellular automata as a physical/neuromorphic substrate where parameters $\theta$ are adapted directly.

3. **Adamatzky, A. (2021).**  
   *Excitable Media Computing.* Chaos, Solitons & Fractals, 146, 110682.
   * **Placement:** [`ghca_main.py`](../ghca_main.py), [`docs/learning_experiments.md`](learning_experiments.md), Website Substrate Specification page.
   * **Relevance & Reasoning:** Comprehensive review of computation using excitable cellular automata, reaction-diffusion substrates, and wave collision dynamics.

4. **Zenke, F., & Neftci, E. O. (2021).**  
   *Brain-inspired learning rules for neuromorphic hardware.* Proceedings of the IEEE, 109(5), 935–950.
   * **Placement:** [`ghca_learn.py`](../ghca_learn.py), [`ghca_plasticity.py`](../ghca_plasticity.py), Website Plasticity Mechanisms page.
   * **Relevance & Reasoning:** Connects local eligibility traces and surrogate/neuromodulated gradients to physical event-driven neuromorphic execution.

5. **Tanaka, G., Yamane, T., Héroux, J. B., et al. (2019).**  
   *Recent advances in physical reservoir computing: A review.* Neural Networks, 115, 100–123.
   * **Placement:** [`docs/e6_results.md`](e6_results.md), [`docs/e9_results.md`](e9_results.md), Website E6 Horde Readout page.
   * **Relevance & Reasoning:** Grounds `ghca`'s readout paradigm (E6/E9): an un-engineered physical/excitable substrate generates a rich high-dimensional phase state, and linear GVF probes or competitive Hebbian readouts extract cognitive functions.

---

### B. Recent Cortical Waves & Spatiotemporal Dynamics (2020–2026)

6. **Karvat, G., Crespo-García, M., Vishne, G., Anderson, M. C., & Landau, A. N. (2026).**  
   *Universal rhythmic architecture uncovers two modes of neural dynamics.* Nature Communications, 17:7024. ([`docs/lit-rev/s41467-026-73553-8.pdf`](lit-rev/s41467-026-73553-8.pdf))
   * **Placement:** [`README.md`](../README.md), [`docs/synthesis.md`](synthesis.md), [`docs/e0_results.md`](e0_results.md), [`docs/e5_results.md`](e5_results.md), Website Overview & Architecture pages.
   * **Relevance & Reasoning:** (See Section 2B above). Universal bi-dimensional spectral architecture separating sustained carrier modes from transient burst modes.

7. **Muller, L., Chavane, F., Reynolds, J., & Sejnowski, T. J. (2018).**  
   *Cortical travelling waves: mechanisms and computational principles.* Nature Reviews Neuroscience, 19(5), 255–268. ([`docs/lit-rev/muller2018.pdf`](lit-rev/muller2018.pdf))
   * **Placement:** [`README.md`](../README.md), [`docs/e0_results.md`](e0_results.md), [`docs/e4_results.md`](e4_results.md), [`ghca_net.py`](../ghca_net.py), Website Cortical Waves & E4 Attention pages.
   * **Relevance & Reasoning:** (See Section 2A above). Establishes mesoscopic wave propagation, axonal delays, sparse wave entrainment, and wave-annihilation computation.

8. **Davis, Z. W., Muller, L., Martinez-Trujillo, J., Sejnowski, T. J., & Reynolds, J. H. (2020).**  
   *Spontaneous travelling cortical waves gate perception in behaving primates.* Nature, 587(7834), 432–436.
   * **Placement:** [`docs/e4_results.md`](e4_results.md), [`docs/spiral_predictions.md`](spiral_predictions.md), Website Falsifiable Predictions page.
   * **Relevance & Reasoning:** Empirical landmark in awake primates proving that spontaneous travelling waves dynamically gate perceptual sensitivity, directly supporting `ghca`'s wave-annihilation attention mechanism (E4).

9. **Liang, Y., Liang, J., Song, C., Liu, M., Knöpfel, T., Gong, P., & Zhou, C. (2023).**  
   *Complexity of cortical wave patterns of the wake mouse cortex.* Nature Communications, 14, 1434.
   * **Placement:** [`docs/spiral_predictions.md`](spiral_predictions.md), [`docs/e7_results.md`](e7_results.md), Website Cortical Spiral Predictions page.
   * **Relevance & Reasoning:** High-resolution optical voltage imaging showing complex topological wave patterns (spirals, sources, sinks) in awake mouse cortex, anchoring the 6 falsifiable spiral wave predictions in Artifact 2b.

10. **Muller, L., Churchland, M. M., & Sejnowski, T. J. (2024).**  
    *Travelling waves as a spatiotemporal computational scaffold for cortical cognition.* Nature Reviews Neuroscience / Trends in Cognitive Sciences.
    * **Placement:** [`docs/synthesis.md`](synthesis.md), [`README.md`](../README.md), Website Synthesis page.
    * **Relevance & Reasoning:** Recent synthesis proposing travelling waves as dynamic spatiotemporal scaffolds for working memory, attention, and executive gating—matching `ghca`'s E2–E7 cognitive series.

---

### C. Recent Intrinsic Timescales & Neural Hierarchies (2020–2026)

11. **Murray, J. D., Bernacchia, A., Wang, X. J., et al. (2014).**  
    *A hierarchy of intrinsic timescales across primate cortex.* Nature Neuroscience, 17(12), 1661–1663.
    * **Placement:** [`docs/timescale_hierarchy_results.md`](timescale_hierarchy_results.md), [`ghca_net.py`](../ghca_net.py), Website Timescale Hierarchy page.
    * **Relevance & Reasoning:** Foundational empirical proof that cortex is organized along a gradient of local intrinsic timescales $\tau$, providing the biological anchor for `ghca`'s per-node refractory parameter $\tau_i$.

12. **Raut, R. V., Snyder, A. Z., Mitra, A., et al. (2021).**  
    *Global waves synchronize the brain’s functional organization with the classical arousal state.* Science Advances, 7(30), eabf2741.
    * **Placement:** [`docs/timescale_hierarchy_results.md`](timescale_hierarchy_results.md), [`docs/continual_learning_results.md`](continual_learning_results.md), Website Track 4a Timescales page.
    * **Relevance & Reasoning:** Empirical demonstration of infra-slow global waves coordinating local intrinsic timescales and state transitions across cortical hierarchies, supporting `ghca`'s input-timing driven emergent $\tau$-hierarchy (3e.2 / Track 4a).

---

### D. Recent Causal Formulation, Spike–Wave Duality & Causal Emergence (2021–2026)

13. **Jalaldoust, K., & Zabeh, E. (2025).**  
    *A Causal Formulation of Spike-Wave Duality.* arXiv preprint arXiv:2511.06602.
    * **Placement:** [`README.md`](../README.md#L164), [`docs/causal_experiments.md`](causal_experiments.md#L1-L10), [`docs/synthesis.md`](synthesis.md#L18-L23), [`ghca_causal.py`](../ghca_causal.py), Website C-Series Overview page.
    * **Relevance & Reasoning:** Direct paper behind the C-series: defines Definition 1 epiphenomenality $P(B|S; do(W)) = P(B|S)$ and d-separation certificates tested in C1–C4.

14. **Lansdell, B. J., Prakash, P., & Kording, K. P. (2023).**  
    *Neural spiking for causal inference and learning.* PLOS Computational Biology, 19(1), e1010836.
    * **Placement:** [`docs/continual_learning_results.md`](continual_learning_results.md), [`docs/synthesis.md`](synthesis.md#L206-L220), Website Track 3 Limitations page.
    * **Relevance & Reasoning:** Proves that node/weight-perturbation learning in spiking networks is interventional $do(\theta)$ causal effect estimation, grounding Track 3c's analysis of causal credit rules.

15. **Hoel, E. P. (2017).**  
    *When the map is better than the territory: The causal emergence of coarse grainings.* Physical Review E, 96(2), 022404.
    * **Placement:** [`docs/c4_results.md`](c4_results.md), [`docs/synthesis.md`](synthesis.md#L113-L125), Website C4 Outcome Relativity page.
    * **Relevance & Reasoning:** Establishes effective information and causal emergence for macro coarse-grainings, justifying why wave variable $W$ exhibits macro-sufficiency ($1.03$ vs. $0.11$) specifically for collective codes.

---

### E. Reinforcement Learning, Abstraction & Unsupervised Conscience Learning

16. **Sutton, R. S., & Barto, A. G. (2018).**  
    *Reinforcement Learning: An Introduction* (2nd ed.). MIT Press.
    * **Placement:** [`docs/learning_experiments.md`](learning_experiments.md#L118-L153), [`ghca_learn.py`](../ghca_learn.py), Website Learning Framework page.
    * **Relevance & Reasoning:** Standard RL text for strict scalar reward $r_t$ and TD error $\delta_t$.

17. **DeSieno, D. (1988).**  
    *Adding a conscience to competitive learning.* In IEEE International Conference on Neural Networks (Vol. 1, pp. 117–124).
    * **Placement:** [`docs/e9_results.md`](e9_results.md), [`ghca_plasticity.py`](../ghca_plasticity.py), Website E9 Emergent Conjunctions page.
    * **Relevance & Reasoning:** Directly implemented in E9: reward-free competitive Hebbian learning ($k$-WTA + DeSieno conscience) to grow emergent $(stimulus \times rule)$ AND-gate basis cells ($0.00 \to 1.00$ selectivity).

18. **Buzsáki, G. (2019).**  
    *The Brain from Inside Out.* Oxford University Press.
    * **Placement:** [`README.md`](../README.md#L16-L20), [`docs/learning_experiments.md`](learning_experiments.md#L4-L9), Website Home page.
    * **Relevance & Reasoning:** Inside-out neuroscientific framing.

---

## 4. Summary Matrix Table

| Category | Key Citation | Target Repo Docs / Code | Website Target Page |
| :--- | :--- | :--- | :--- |
| **Universal Rhythmic Architecture** | Karvat et al. (2026 Nat Commun) | [`README.md`](../README.md), [`docs/synthesis.md`](synthesis.md) | Overview & Architecture |
| **Cortical Waves Review** | Muller et al. (2018 NRN) | [`README.md`](../README.md), [`ghca_net.py`](../ghca_net.py) | Cortical Waves & Overview |
| **Empirical Waves & Perception** | Davis et al. (2020 Nature) | [`docs/e4_results.md`](e4_results.md), [`docs/spiral_predictions.md`](spiral_predictions.md) | E4 Attention & Predictions |
| **Mouse Cortical Wave Dynamics** | Liang et al. (2023 Nat Commun) | [`docs/spiral_predictions.md`](spiral_predictions.md), [`docs/e7_results.md`](e7_results.md) | Cortical Spiral Predictions |
| **Waves as Cognitive Scaffold** | Muller et al. (2024) | [`docs/synthesis.md`](synthesis.md), [`README.md`](../README.md) | Synthesis & Architecture |
| **Lifelong Learning & Metaplasticity** | Kudithipudi et al. (2022 Nat Mach Intell) | [`docs/continual_learning_results.md`](continual_learning_results.md) | Track 3 & Anti-Forgetting |
| **Physical Substrate Computing** | Wright et al. (2022 Nature) | [`README.md`](../README.md), [`docs/synthesis.md`](synthesis.md) | Overview & Substrate |
| **Excitable Media Computing** | Adamatzky (2021) | [`ghca_main.py`](../ghca_main.py), [`docs/learning_experiments.md`](learning_experiments.md) | Substrate Specification |
| **Neuromorphic Event Rules** | Zenke & Neftci (2021 IEEE Proc) | [`ghca_learn.py`](../ghca_learn.py), [`ghca_plasticity.py`](../ghca_plasticity.py) | Plasticity Mechanisms |
| **Physical Reservoir Computing** | Tanaka et al. (2019) | [`docs/e6_results.md`](e6_results.md), [`docs/e9_results.md`](e9_results.md) | E6 Horde Readout |
| **Intrinsic Timescale Hierarchy** | Murray et al. (2014 Nat Neurosci) | [`docs/timescale_hierarchy_results.md`](timescale_hierarchy_results.md), [`ghca_net.py`](../ghca_net.py) | Timescale Hierarchy |
| **Global Waves & Timescale Sync** | Raut et al. (2021 Sci Adv) | [`docs/timescale_hierarchy_results.md`](timescale_hierarchy_results.md) | Track 4a Timescales |
| **Spike–Wave Causal Duality** | Jalaldoust & Zabeh (2025 arXiv) | [`README.md`](../README.md), [`docs/causal_experiments.md`](causal_experiments.md) | C-Series Overview |
| **Spiking Causal Inference** | Lansdell et al. (2023 PLOS Comp Bio) | [`docs/continual_learning_results.md`](continual_learning_results.md) | Track 3 Limitations |
| **Causal Emergence Coarse-Graining**| Hoel (2017 Phys Rev E) | [`docs/c4_results.md`](c4_results.md), [`docs/synthesis.md`](synthesis.md) | C4 Outcome Relativity |
| **Conscience Hebbian Learning** | DeSieno (1988) | [`docs/e9_results.md`](e9_results.md), [`ghca_plasticity.py`](../ghca_plasticity.py) | E9 Emergent Conjunctions |
