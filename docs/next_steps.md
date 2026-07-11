# Next Steps — Planning / Roadmap

*Candidate directions after E0–E8 and C0–C7. This is a menu, not a commitment:
each option lists what it is, why it matters, what it would take, effort/risk, and
how it connects to existing work and to the honesty caveats the two review passes
surfaced. Nothing here is built yet.*

## Where the project stands

- **E-series (learning).** E0 substrate/operating point → E1 conditioning → E2
  working memory (τ-controlled loop) → E3 timed response (spatial/temporal double
  dissociation; A+B composition only *directionally* supported) → E4 attention
  (WTA by wave annihilation) → E5 executive control (slow-loop option) → E6
  emergent categories (GVF demons) → E7 the option as a 2-D spiral core (rotation
  direction as rule) → E8.0–8.7 predictive dynamics (prediction, `do(τ)` history
  window, global surprise, nested two-timescale gating, order-preserving reservoir,
  conditional long-range).
- **C-series (causality).** C0–C4 spike–wave causality on the substrate → C5–C7 the
  spiral causal split (`do(χ)` fat-handed at a fixed locus; `do(θ_χ)` well-posed;
  rotation direction is a causal mediator, outcome-relative, not epiphenomenal).
- **Reviews.** An independent hallucination/overreach audit (E0–E6/C0–C4; on branch
  `claude/ai-hallucination-review-ya7aq5`) and a self-audit of the extensions
  ([`extensions_review.md`](extensions_review.md)). No fabrication; one reproducibility
  bug (fixed: `perturb_tau` seeding) and one framing overreach (fixed: E3 composition).

## The three tensions that should steer the roadmap

The audits converge on three honest limitations. Good next steps *retire* one of these.

1. **Afforded vs learned.** Several headline capabilities are *afforded* by hand-built
   readouts/features (the E7 winding readout, the E8.5/8.7 conjunction feature) and a
   reused routing learner — not learned end-to-end *on* the substrate. This is the
   inside-out thesis's own promissory note.
2. **Narrow evidence.** Small n (3–5 seeds), single substrate (`lattice2d`), hand-chosen
   operating points. Results are robust *at those points* but not shown *across the
   regime*.
3. **Illustrative, not proof.** Toy Greenberg–Hastings substrates illustrate (do not
   prove) the cognitive-function claims; the external artifact with the most reach is
   the synthetic-SCM causal testbed.

---

## Track 1 — Close the substrate/analysis gap *(retires tension 1; highest credibility gain)*

### 1a. Emergent conjunction cells *(recommended if the goal is credibility)*
- **What.** Let (stimulus × context) conjunction cells *self-organise* via plasticity on
  the hidden recurrence, instead of the hand-wired channel bias (E1) / hand-provided
  outer-product feature (E5, E8.5, E8.7).
- **Why.** Retires the audit's central caveat and the oldest deferred item ("let the net
  *discover* selective hidden representations", E1). Converts "afforded" → "learned".
- **Takes.** Fold plastic H→H into Line A; a reward or Hebbian rule that grows
  conjunction selectivity; re-run E1/E5 showing the conjunctions emerge rather than being
  wired.
- **Effort.** Medium–high. **Risk.** Medium — discovering conjunctions under a scalar
  reward is genuinely hard (E1 hit this); may need a curriculum or auxiliary objective.
- **Depends on / connects to.** E1, E5, E8.5/8.7.

### 1b. A learned direction-selective readout
- **What.** Replace E7's *computed* local-winding readout with a small population that
  *learns* to read rotation direction from the wave (direction-selective cells via
  reward/Hebbian).
- **Why.** Makes E7's "read the wave" genuinely on-substrate; maps onto real cortical
  direction selectivity.
- **Takes.** A motion/phase-gradient-sensitive readout population + a learning rule;
  show it recovers chirality and drives routing without the god's-eye winding.
- **Effort.** Medium. **Risk.** Medium — direction selectivity from a meandering core is
  noisy (C5 showed the readout-locality problem).
- **Connects to.** E7, C5 (readout-relativity).

### 1c. Reward-/GVF-driven prediction
- **What.** Fold E8's offline-ridge predictor into the substrate's own TD/GVF machinery
  (E6 demons, gradient-TD), so prediction is learned online and intrinsically.
- **Why.** Tests "prediction is inside-out, not taught"; unifies E8 with E6; removes the
  "learning lives in the readout weights" caveat.
- **Takes.** GVF next-tone demons with online (gradient-)TD on the tonotopic substrate;
  compare to the ridge baseline.
- **Effort.** Medium. **Risk.** Low–medium.
- **Connects to.** E6, E8.

---

## Track 2 — Make the neuroscience bridge falsifiable *(retires tension 3; scientific value)*

### 2a. Build the predictive-coding foil *(recommended if the goal is a clean disambiguation)*
- **What.** Implement a small Rao–Ballard/Friston predictive-coding model on the *same*
  tone-sequence task E8 uses, and measure the distinguishing observables (global vs
  per-feature prediction error; prediction/representation dissociability under lesion).
- **Why.** E8 currently *asserts* "prediction without predictive coding"; this makes the
  contrast empirical rather than rhetorical.
- **Takes.** A hierarchical PC network + the same sequences + the same probes.
- **Effort.** Medium–high (a second model). **Risk.** Low — clean comparison either way.
- **Connects to.** E8.

### 2b. Model → data predictions
- **What.** Derive falsifiable claims for the Gong/Steinmetz spiral-wave data: e.g.
  reversing a cortical spiral's chirality should reconfigure routing; ablating persistent
  cores should selectively impair switching-like flexibility.
- **Why.** Bridges toy → testable neural hypothesis; the highest-reach scientific output.
- **Takes.** Mostly analysis/writing atop E7 + C5–C7; identify a measurable signature.
- **Effort.** Low–medium. **Risk.** Low (a hypothesis, not a claim).
- **Connects to.** E7, C5–C7; Gong 2023, Steinmetz/Ye 2026.

---

## Track 3 — Generalise / harden *(retires tension 2)*

### 3a. Statistics & operating-point sweeps
- **What.** Many more seeds, CIs, and sweeps around every hand-chosen operating point.
- **Why.** Turns "shown at n=5 at a point" into "robust across the regime" — directly
  answers the audit's scope caveat.
- **Effort.** Low (compute). **Risk.** Low; may soften some headlines (worth knowing).

### 3b. Other topologies
- **What.** Characterise + re-run the two-line story on `smallworld` / `rgg` (E0 only did
  `lattice2d`; the design doc's `smallworld` default is unused).
- **Why.** Tests generality of routing/timing and of the spiral band.
- **Effort.** Medium. **Risk.** Medium — spirals/loops may behave differently off-lattice.

### 3c. Continual learning on one substrate
- **What.** Can *one* substrate learn E1→E5 sequentially without catastrophic
  interference (rather than E6's post-hoc freeze)?
- **Why.** The real test of the "one homogeneous machine" claim.
- **Effort.** High. **Risk.** High — interference is likely; could be a negative result
  (still informative).

---

## Track 4 — New territory *(novelty)*

### 4a. Emergent timescale hierarchy / cross-frequency coupling *(recommended if the goal is novelty)*
- **What.** With per-node τ plastic, do learned τ distributions self-organise into a
  fast/slow hierarchy with theta–gamma-style cross-frequency coupling?
- **Why.** Deepest genuinely-new phenomenon the substrate could show; connects to the
  nested-waves literature (E8.5 only used two hand-set timescales).
- **Takes.** Line B on a task rewarding multi-scale structure; measure the learned τ
  spectrum and cross-frequency coupling.
- **Effort.** High. **Risk.** Medium–high — may not self-organise; the `perturb_tau`
  mechanism is coarse.
- **Connects to.** E2, E5, E8.5.

### 4b. Package the causal testbed
- **What.** Turn C0–C7 + the substrate into a reusable synthetic-SCM benchmark for the
  Jalaldoust spike–wave framework (clean API, ground-truth graphs, the three `do`
  operators).
- **Why.** The most externally-reusable artifact — a ground-truth SCM others can test
  causal-discovery/epiphenomenality methods on. Low new-science, high utility.
- **Effort.** Medium (engineering/packaging). **Risk.** Low.
- **Connects to.** C0–C7, `ghca_causal.py`.

---

## Track 5 — Consolidate

### 5a. Unified write-up
- **What.** Draw E0–E8 + C0–C7 into one honest narrative (synthesis.md is the seed),
  positioned as an illustrative computational study + a reusable causal testbed.
- **Why.** Makes the arc legible; forces the framing to stay at the right altitude.
- **Effort.** Medium (writing). **Risk.** Low.

### 5b. Reproducibility hygiene
- **What.** Audit every experiment for seeded RNG (the `perturb_tau` bug pattern);
  a `reproduce-all` entry point that regenerates every figure/number; optionally CI.
  Also consider landing `hallucination_review.md` on `main` (currently only on its
  review branch) so both audits sit with the work.
- **Effort.** Low–medium. **Risk.** Low.

---

## Decision matrix (by goal)

| Goal | Do first | Then |
|------|----------|------|
| **Credibility** (retire the honesty gap) | **1a** emergent conjunctions | 1b, 3a |
| **Clean disambiguation** | **2a** predictive-coding foil | 1c |
| **Scientific novelty** | **4a** emergent timescale hierarchy | 2b |
| **External impact / reach** | **4b** causal testbed | 5a write-up, 2b |
| **Lowest-risk strengthening** | **3a** stats/sweeps + **5b** hygiene | 2b |

**Overall recommendation.** If one thing: **1a (emergent conjunction cells)** — it
retires the single most-cited caveat and is the truest test of the inside-out thesis.
If two, pair it with **2b (model→data predictions)** for the neuroscience reach at low
cost. Chase **4a** when the goal shifts from "make the existing story solid" to "find a
genuinely new phenomenon."

## Process notes (apply to whatever is chosen)

- Seed *everything* (the `perturb_tau` lesson); prefer `default_rng(seed)` threaded
  explicitly, never the global RNG.
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the dynamics do vs what a
  readout/feature does.
- Keep a caveats section adjacent to every headline.
