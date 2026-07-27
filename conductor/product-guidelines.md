# Product Guidelines & House Rules — GHCA

## Scientific Rigor & House Rules
1. **Explicit Seed Management:** Thread explicit `default_rng(seed)` instances through all experiments and substrate classes. Never use global NumPy RNG (`np.random.*`).
2. **Per-Seed Spreads & Bimodality:** Always report per-seed ranges/spreads, 95% confidence intervals, and explicitly call out bimodal distributions.
3. **Substrate vs. Readout Boundary:** Clearly distinguish between what the physical *substrate/dynamics* accomplish versus what a fitted *readout/decoder* computes.
4. **Decoupled Planning & Review:** Review passes (backward-looking, adversarial) and Planning passes (forward-looking, generative) are strictly decoupled in process and linked only by a one-directional hand-off (Review → Plan).
5. **Adjacent Caveats:** Every headline finding must have an adjacent Caveats section detailing limitations, boundary conditions, and sampling constraints.
