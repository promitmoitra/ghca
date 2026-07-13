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

### 1a. Emergent conjunction cells — ✅ **DONE** (see [`e9_results.md`](e9_results.md))
- **What.** Let (stimulus × context) conjunction cells *self-organise* via plasticity on
  the hidden input edges, instead of the hand-wired channel bias (E1) / hand-provided
  outer-product feature (E5, E8.5, E8.7).
- **Why.** Retires the audit's central caveat and the oldest deferred item ("let the net
  *discover* selective hidden representations", E1). Converts "afforded" → "learned".
- **What was built (E9).** Unbiased hidden fan-in + a reward-free **competitive Hebbian**
  rule (k-WTA lateral inhibition + DeSieno conscience, keeping E5's AND-gate) grows the
  `(stim × rule)` conjunction basis from the input statistics: conjunction selectivity
  **0.00 → 1.00**, all four classes perfectly tiled, within ~3 exposures/combo.
  Reward-driven Line A routing on the **emergent** basis reaches **0.84** (vs wired
  **0.92**); the no-self-organisation **frozen** control fails at chance (0.25) — the
  self-organisation is load-bearing. The plan's worry that a scalar reward alone would
  struggle to discover conjunctions (E1's difficulty) was right; the fix was an
  unsupervised, local, label-free representation rule, with reward left to do only the
  routing.
- **Deferred.** A single *concurrent* process (representation + readout plastic together,
  no phase split) and realising the competition as recurrent inhibitory dynamics rather
  than an imposed k-WTA — see the E9 caveats.
- **Connects to.** E1, E4 (WTA), E5, E8.5/8.7.

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

### 1c. Reward-/GVF-driven prediction — ✅ **DONE** (see [`e8_hardening_results.md`](e8_hardening_results.md))
- **What.** Fold E8's offline-ridge predictor into the substrate's own TD/GVF machinery
  (E6 demons, gradient-TD), so prediction is learned online and intrinsically.
- **Why.** Tests "prediction is inside-out, not taught"; unifies E8 with E6; removes the
  "learning lives in the readout weights" caveat.
- **What was built.** A bank of `M` linear GVF demons (E6's TD rule) learned online: a
  single incremental reward-free pass matches the offline ridge on every sequence
  (periodic 1.00, Markov 0.91/0.56, random-walk 0.41), and γ>0 demons give multi-step
  predictive knowledge (r=0.97/0.85) the ridge doesn't. **Softens** caveat 7 (learning
  is now online/intrinsic) but doesn't fully retire it — still a readout of a fixed
  substrate; plastic dynamics deferred.
- **Connects to.** E6, E8.

---

## Track 2 — Make the neuroscience bridge falsifiable *(retires tension 3; scientific value)*

### 2a. Build the predictive-coding foil — ✅ **DONE** (see [`e8_hardening_results.md`](e8_hardening_results.md))
- **What.** Implement a small Rao–Ballard/Friston predictive-coding model on the *same*
  tone-sequence task E8 uses, and measure the distinguishing observables (global vs
  per-feature prediction error; prediction/representation dissociability under lesion).
- **Why.** E8 currently *asserts* "prediction without predictive coding"; this makes the
  contrast empirical rather than rhetorical.
- **What was built + honest finding.** A minimal untied Rao–Ballard PC model on E8's
  tones. The clean PC signature reproduces (top-down lesion: prediction 0.37→0.12 chance,
  representation 0.95 spared). But the two observables one expects to disambiguate — **do
  not**: E8 is *also* dissociable (its two readouts tap different medium parts), and both
  models' errors localise on the deviant. This **corrects the E8-doc "not dissociable"
  overclaim**. The real discriminators are architectural (intrinsic error units +
  generative pathway [PC] vs global scalar surprise + passive-readout + do(τ) window
  [E8]) — now stated as empirical, testable contrasts.
- **Connects to.** E8; complements 1c (E8 hardened on both tension axes).

### 2b. Model → data predictions — ✅ **DONE** (see [`spiral_predictions.md`](spiral_predictions.md))
- **What.** Derive falsifiable claims for the Gong/Steinmetz spiral-wave data: e.g.
  reversing a cortical spiral's chirality should reconfigure routing; ablating persistent
  cores should selectively impair switching-like flexibility.
- **Why.** Bridges toy → testable neural hypothesis; the highest-reach scientific output.
- **What was built.** Six predictions (P1–P6) from E7 + C5–C7, each with observable,
  discriminator vs an epiphenomenal/non-spiral account, and a falsifier: persistent
  rotation-direction rule code (P1); persistent cores necessary for flexibility not fixed
  mappings (P2, highest-reach); outcome-relativity χ→rule-not-content (P3); a fixed-ROI
  decoder collapsing with core drift while a topology-aware one holds (P4, with a bridge
  figure recast from C5); nucleation as the clean handle + χ mediating θ→behaviour (P5);
  persistent lone cores localising to anatomical boundaries (P6). Bridge figure:
  `experiments/p2b_signature_figure.py`.
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

### 4a. Emergent timescale hierarchy / cross-frequency coupling — ⏸ **ATTEMPTED, PAUSED** (see [`e10_notes.md`](e10_notes.md))
- **What.** With per-node τ plastic, do learned τ distributions self-organise into a
  fast/slow hierarchy with theta–gamma-style cross-frequency coupling?
- **Why.** Deepest genuinely-new phenomenon the substrate could show; connects to the
  nested-waves literature (E8.5 only used two hand-set timescales).
- **Status.** The `Medium–high` risk was realised: the existing Line B resonance rule
  **structurally cannot** build the hierarchy — it only ratchets τ *upward* toward
  multiples of a drive period, never down to a fundamental, so no fast (small-τ)
  population forms (confirmed decisively in the ideal isolated-node case; see
  `experiments/e10_diagnostics.py` on branch `claude/e10-timescale-hierarchy`). 4a is a
  *mechanism-design* task (a new bidirectional "τ tracks input period" rule + E9-style
  competition + balanced channel activity), **not** a reuse-E9 task — correcting the
  earlier de-risking claim. Resume from `e10_notes.md`.
- **Connects to.** E2, E5, E8.5; E9 (competition, for the grouping half only).

### 4b. Package the causal testbed — ✅ **CORE DONE** (see [`causal_testbed.md`](causal_testbed.md); spec/plan alongside)
- **What.** Turn C0–C7 + the substrate into a reusable synthetic-SCM benchmark for the
  Jalaldoust spike–wave framework (clean API, ground-truth graphs, the three `do`
  operators).
- **Why.** The most externally-reusable artifact — a ground-truth SCM others can test
  causal-discovery/epiphenomenality methods on. Low new-science, high utility.
- **What was built.** The `ghca_testbed` package: single-source operator/substrate
  re-exports; the six canonical C1 SCM graphs **executable in-process** with the correct
  interventional-vs-observational `epiphenomenality_test` + Theorem-1 certificate
  (validated end-to-end — front-door scored causal 0.257, all six verdicts match);
  the method-agnostic metric library (`fat_hand_band`, `macro_sufficiency`,
  `outcome_matrix`, `mediation`); a 12-scenario `REGISTRY` (C1 executable, C2–C7 ground
  truth carried as verified data) with JSON export; a `score(method_fn)` harness +
  `python -m ghca_testbed`; `docs/causal_testbed.md`; and `test_ctestbed.py` (10 tests,
  including a guard for the C1 subtlety). Planned via a multi-agent workflow
  ([`causal_testbed_plan.md`](causal_testbed_plan.md)).
- **Deferred (plan 4b.0d).** The verbatim refactor of the eight C-scripts into thin
  callers behind a bit-identical gate, and wiring the C2–C7 substrate scenarios to
  execute in-process through `Testbed` subclasses. Left as reference implementations;
  their ground truth is carried and reproduced by running the scripts.
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

**Progress.** **1a done (E9)** — retired the most-cited caveat. **2b done**
([`spiral_predictions.md`](spiral_predictions.md)) — the toy→data bridge, low cost, high
reach. **4a attempted, paused** ([`e10_notes.md`](e10_notes.md)) — the existing τ rule
structurally can't build the hierarchy; it needs a new bidirectional τ-plasticity rule
(mechanism design), and E9 de-risks only the *grouping* half, not the τ-value rule
(earlier claim corrected). Remaining high-value, lower-risk options: **2a** (the
predictive-coding foil, cleanest disambiguation), **4b** (package the causal testbed,
most reusable artifact), **3a** (stats/sweeps). Return to **4a** only with appetite for
mechanism design.

## Process notes (apply to whatever is chosen)

- Seed *everything* (the `perturb_tau` lesson); prefer `default_rng(seed)` threaded
  explicitly, never the global RNG.
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the dynamics do vs what a
  readout/feature does.
- Keep a caveats section adjacent to every headline.
