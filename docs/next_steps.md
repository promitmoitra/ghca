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
- **Reviews.** An independent integrity/overreach audit (E0–E6/C0–C4;
  [`core_review.md`](core_review.md)) and a self-audit of the extensions
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

### 3a. Statistics & operating-point sweeps — ✅ **CORE DONE (P1–P3b)** (see [`stats_sweeps_results.md`](stats_sweeps_results.md))
- **What.** Many more seeds, CIs, and sweeps around every hand-chosen operating point.
- **Why.** Turns "shown at n=5 at a point" into "robust across the regime" — directly
  answers the audit's scope caveat.
- **What was built.** A seeded stats harness (`ghca_stats.py`: bootstrap/Wilson CIs,
  effect size, an automated bimodality screen) and n=50 (publication-grade) re-runs of
  every headline, plus two operating-point sweeps. **Outcome — tension 2 addressed with
  evidence, not asserted away:** E2/E5/E9/E8 and C5 **strengthen** (tight, cleanly
  separated CIs; E8, the thinnest series at n=3, is the *most* robust); **E7 switching
  softens** (0.86 → 0.75 [0.70, 0.79]) but is **robust across the θ band**; and the one
  genuinely fragile claim, **E3 composition**, is confirmed **operating-point-contingent**
  (joint-success 0%–40% by target latency, 32% [21, 46] at the default, tracking the
  substrate resonance map). E7/E3 result-doc headlines updated with n=50 notes.
  **P3b** (`experiments/stats_p3b.py`, n=30) closed the σ-band/outcome-matrix gap:
  C2's do(W) fat-hand gap holds (0.29σ vs 26.8σ); C4's outcome matrix is **exactly**
  reproducible (bit-identical across all 30 seeds); C5's `tracked`/`global` readers
  are tight, `center` is genuinely spread (not bimodal); and one real methodology
  finding — C7's normalized `do(χ)→content` (0.11 at n=1) is noisy at n=30 (mean
  0.56) because column-max normalisation is unstable when both handles have weak
  raw effects, **not** because chirality's raw content effect grew (it stays modest
  and smaller than `do(route)`'s raw effect throughout) — flagged as a normalization
  caveat, not a reversal.
- **Deferred.** P4b — fold the P3b table + the do(χ)→content normalization caveat
  into the individual C-series result docs; the gate-τ axis (τ is learned, not set)
  and a full θ×τ grid for E3.
- **Effort.** Low (compute). **Risk.** Low; softened E7/E3 as expected (worth knowing).

### 3b. Other topologies — ✅ **DONE (dynamics + E1 learning)** (see [`e0_topologies.md`](e0_topologies.md))
- **What.** Characterise + re-run the two-line story on `smallworld` / `rgg` (E0 only did
  `lattice2d`; the design doc's `smallworld` default is unused).
- **Why.** Tests generality of routing/timing and of the spiral band.
- **What was built.** (i) E0's two topology-agnostic observables re-run on `lattice2d`,
  `smallworld`, and `rgg` at matched N=1600 and mean degree 12 (8 graph seeds, CIs):
  a self-sustaining excitable band exists on all three (`A_ss≈0.40` for θ=1–3, death by
  θ≈4–5), and the `period ~ τ` law holds identically (slope 0.98/0.99/1.01, r≈0.999).
  (ii) The **learned** result too: a non-breaking `hh_topo` knob swaps E1's hidden
  reservoir (smallworld→ring→rgg, matched degree); reward routing (Line A) learns on
  **every** medium (0.86–0.93, all ≫ Line-B control; d=1.06–1.75), smallworld modestly
  best. So both the excitable dynamics and the learned routing dissociation are
  properties of a broad class of recurrent media — the substrate-generality half of the
  narrow-evidence tension is retired.
- **Deferred.** The 2-D **spiral** (E7, C5–C7) is geometry-bound, no off-lattice analogue
  (out of scope). Porting **E2/E5** (τ-loop memory/executive) onto non-lattice media and
  a degree/size sweep remain.
- **Effort.** Medium. **Risk.** Medium — realised: spiral is lattice-only; dynamics *and*
  E1 learning generalise.

### 3c. Continual learning as causal credit assignment — 📋 **SCOPED** (see [`continual_learning_plan.md`](continual_learning_plan.md))
- **What.** Can *one* substrate learn E1→E2→E5 sequentially without catastrophic
  interference (rather than E6's post-hoc freeze)? **Re-scoped**: the CL literature
  shows frozen-substrate + per-head is the *easy* case (≈ E6, no interference), so the
  study targets the **plastic shared substrate** and asks whether **causal credit
  assignment** (`do(θ)` node-perturbation) reduces forgetting vs the current
  **correlational** eligibility trace, relative to the frozen upper bound.
- **Why.** The real test of the "one homogeneous machine" claim, *and* the concrete
  vehicle for unifying the learning and causality arcs (learning = causal inference;
  C2/C3 says `do(θ)` is the well-posed handle → causal-θ credit should interfere less).
- **Design.** 2 regimes (frozen / plastic) × 2 credit rules (correlational / causal-θ);
  CL metrics (avg acc, backward/forward transfer) with 3a-style CIs; Mesnard hindsight
  estimator + native WTA gating as v2 conditions.
- **Effort.** High (new learning harness). **Risk.** High and intended — the
  causal-vs-correlational contrast is genuinely uncertain; a clean null is publishable.
- **Connects to.** E6 (frozen baseline), E9/E4 (WTA gating), C2/C3 (`do(θ)` well-posed),
  Line A/B plasticity.

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

### 5a. Unified write-up — ✅ **DONE** (see [`synthesis.md`](synthesis.md))
- **What.** Draw E0–E8 + C0–C7 into one honest narrative (synthesis.md is the seed),
  positioned as an illustrative computational study + a reusable causal testbed.
- **Why.** Makes the arc legible; forces the framing to stay at the right altitude.
- **What was built.** Expanded `synthesis.md` into the capstone: the two arcs + the
  `θ` = learned-variable = causal-handle spine, the E4–E6 repertoire and E7/C5–C7/E8
  extensions, then the new work woven in — the afforded→learned closure (E9 + 1c), the
  predictive-coding foil and the overclaim it corrected (2a), and the outward artifacts
  (2b predictions + 4b testbed) — capped by an **honest ledger** of the three tensions
  as they now stand (1 substantially retired, 2 still open, 3 addressed by framing +
  artifacts) and the corrections the programme made (E3, `perturb_tau`, 2a).
- **Effort.** Medium (writing). **Risk.** Low.

### 5b. Reproducibility hygiene — ✅ **CORE DONE** (verified 2026-07-13)
- **What.** Audit every experiment for seeded RNG (the `perturb_tau` bug pattern);
  a `reproduce-all` entry point that regenerates every figure/number; optionally CI.
- **What is done.**
  - The `perturb_tau` unseeded-global-RNG bug the [`core_review.md`](core_review.md)
    audit found is **fixed**: `perturb_tau` now draws from the seeded `self.rng`
    (`ghca_learn.py:138/143`), and E2/E3 were regenerated to their seeded numbers
    (E3 curriculum composition mean **0.560**, per-seed `[0.50, 0.51, 0.79, 0.00,
    1.00]`).
  - **Independently verified** by a from-scratch re-run of the previously-flagged
    Line-B conditions: **E2, E3 (main), and E3 (factored) reproduce
    bit-identically** (max\|Δ\| = 0 across every array, including the `ret_B` /
    `tau_B` arrays the audit had found non-regenerable and the composition
    headline). The audit's one "not fully verified" residual, the **E0 period-fit
    `r=0.9992`**, was re-derived from `result/e0/e0_data.npz`: `period =
    0.9991·τ + 0.9480, r = 0.9992` on the 8 oscillating points (the τ=20 `inf`
    non-oscillating point dropped) — matches the doc.
  - The independent core-series audit now sits on `main` as
    [`core_review.md`](core_review.md), alongside the extensions self-audit — the
    "both audits with the work" gap is closed; see [`process.md`](process.md).
- **Deferred.** (i) A single `reproduce-all` entry point that regenerates every
  figure/number, and optional CI. (ii) `ghca_main.py` (the **original lattice-CA
  strand**, README strand 1) still uses the global NumPy RNG (`ghca_main.py:100–103`);
  it is outside the E/C experiment path the audits cover, but bringing it under the
  seed-everything house rule would finish the sweep.
- **Effort.** Low–medium. **Risk.** Low.

---

## Track 6 — Extinguishing reentrant activity without a trackable core *(new capability; not audit-driven — doesn't retire tensions 1–3, prompted by a design discussion rather than a review finding)*

*C6's `do(θ)` ablation is a single global parameter change on a substrate with one
seeded, trackable spiral core — `signed_charge`/`local_winding` assume a locatable
singularity to read out. Real excitable media (cardiac fibrillation is the paradigm
case) often have per-cell heterogeneous timescales and multiple, drifting, or no
clean core at all. This track asks: with nothing dominant to localise, how do you
extinguish self-sustaining activity? Four independent, complementary levers.*

*Shared prerequisite: the C5/C6 success criterion (does the winding-number readout
match the true chirality) has no referent once there is no trackable core. Every
thread below needs a cruder, structure-agnostic order parameter instead — e.g.
global active fraction decaying to (and staying at) the spontaneous-noise floor, or
an entropy/coherence measure separating "organised wavelets" from "noise-driven
flicker."*

### 6a. Percolation-threshold lesioning
- **What.** Force a fraction `p` of cells permanently non-excitable (`theta=inf`, no
  code change needed — `inp >= theta` is simply never true), randomly or in a
  pattern, and sweep `p` (and neighbourhood range `r`) for a critical density above
  which no path anywhere supports a full reentrant loop, vs. below which wavelets
  route around the gaps.
- **Why.** Needs no localisation at all — if the lesion pattern disconnects the
  medium into non-percolating clusters, extinction follows from graph structure, not
  from finding anything. A direct generalisation of E0's own threshold-range (FGG)
  scaling into a percolation account of extinction, and the closest analogue of real
  "substrate modification" ablation strategies used when no single rotor is found.
- **Takes.** A lesion-mask generator on `ghca_net.Network` + a sweep over `p` and
  `r`; measure active-fraction decay to the spontaneous floor; compare against a
  percolation-theory prediction for the critical density at a given `r`.
- **Effort.** Medium. **Risk.** Medium — heterogeneous, multi-wavelet activity may
  not have one sharp threshold; wavelets could adaptively reroute around sparse
  lesions unless density is high.
- **Connects to.** E0 (threshold-range scaling), C2 (constituted aggregates), C6
  (the spatial counterpart to its global `do(θ)`).

### 6b. Closed-loop / reactive ablation
- **What.** A feedback controller that reads a live local-activity signal (e.g. a
  sliding-window active fraction) and raises `theta` / forces refractory only where
  and when activity is currently detected — no a priori patch or global constant,
  purely reactive.
- **Why.** Makes no assumption about a trackable or even locatable core. The
  substrate already has a precedent for readout-driven `theta` changes — the
  homeostatic `rho_star` mechanism (`ghca_net.py`) adapts one *global* threshold from
  the *global* active fraction; 6b is the spatially-local, suppression-oriented
  generalisation of that same pattern. The real-world analogue is anti-tachycardia
  pacing: implantable devices sense local activation and respond without knowing
  where the rotor is in advance.
- **Takes.** A live local order-parameter readout + a control law (threshold/dose as
  a function of detected local activity) layered on `Network.step`; test on both the
  trackable single-spiral substrate (contrast with C6's blind global ablation) and a
  genuinely heterogeneous, untrackable multi-wavelet one.
- **Effort.** Medium–high (a new controller layer). **Risk.** Medium — a naive
  reactive controller could chase activity indefinitely without ever catching up if
  wavelets are fast or numerous.
- **Connects to.** `Network.coherence` / `active_mask` (the existing order
  parameters), the homeostatic `rho_star` rule (the pattern this generalises), C3/C6
  (`theta` as the causal handle).

### 6c. Overdrive pacing / entrainment
- **What.** Drive the substrate with a strong, fast, spatially localised (one or a
  few sites) periodic stimulus; ask whether it entrains the tissue away from
  disorganised reentry, so that withdrawing the drive at the right phase leaves the
  substrate quiescent.
- **Why.** The only lever here that *adds* excitation rather than removing or raising
  a threshold — exploits refractory-period competition to let one dominant rhythm
  outcompete existing wavefronts, then simply stops. Doesn't require finding or
  disabling anything.
- **Takes.** A periodic forced drive (reusing the existing `drive` argument to
  `Network.step`) at one or a few sites, swept over drive frequency, phase, and
  withdrawal timing; success = active fraction reaches and stays at zero after
  withdrawal.
- **Effort.** Medium. **Risk.** Medium–high — entraining a genuinely heterogeneous,
  multi-wavelet medium is the least certain of the four; a narrow working window
  would itself be an interesting (negative-adjacent) finding.
- **Connects to.** E1 (clamped sensory drive, the same `drive` mechanic), C6
  (contrast: an excitatory stopping intervention vs. a threshold-raising one).

### 6d. Heterogeneity as the destabilising lever
- **What.** Instead of raising the *mean* `theta`/`tau` (C6's move), raise the
  *variance* of per-cell timescales or thresholds with the mean held fixed, and ask
  whether heterogeneity alone destabilises otherwise-persistent reentry into
  self-terminating chaos.
- **Why.** The counter-intuitive option — more disorder, not less excitability —
  grounded in a real phenomenon (source–sink mismatch / conduction block at steep
  refractoriness gradients causing wavebreak). If it holds, heterogeneity itself
  becomes a causal handle in the `do(theta)` vocabulary, not just an inconvenience
  the other three threads have to route around.
- **Takes.** Sweep the variance of a per-cell `theta`/`tau` distribution (mean fixed)
  on both a nucleated single spiral and a multi-wavelet substrate; measure
  time-to-extinction or steady-state active fraction vs. variance; contrast against
  C6's mean-shift ablation to see whether variance alone (no mean shift) suffices.
- **Effort.** Low–medium (a parameter sweep on existing machinery). **Risk.**
  Low–medium.
- **Connects to.** C6 (contrast: variance vs. mean as the operative parameter), C4
  (outcome-relativity — is heterogeneity causal for extinction specifically, and
  inert for other outcomes?).

---

## Decision matrix (by goal)

| Goal | Do first | Then |
|------|----------|------|
| **Credibility** (retire the honesty gap) | **1a** emergent conjunctions | 1b, 3a ✅ |
| **Clean disambiguation** | **2a** predictive-coding foil | 1c |
| **Scientific novelty** | **4a** emergent timescale hierarchy | 2b |
| **External impact / reach** | **4b** causal testbed | 5a write-up, 2b |
| **Lowest-risk strengthening** | **3a** stats/sweeps ✅ + **5b** hygiene ✅ done | 1b, 3b |

**Progress.** **1a done (E9)** — retired the most-cited caveat. **1c / 2a done**
([`e8_hardening_results.md`](e8_hardening_results.md)) — online GVF prediction and the
predictive-coding foil. **2b done** ([`spiral_predictions.md`](spiral_predictions.md)) —
the toy→data bridge, low cost, high reach. **4b core done**
([`causal_testbed.md`](causal_testbed.md)) — the reusable synthetic-SCM benchmark.
**5b core done** (verified 2026-07-13) — the `perturb_tau` fix is landed and E0/E2/E3
reproduce bit-identically; only a `reproduce-all` entry point and the original
lattice-CA strand's RNG remain. **3a core done (P1–P3b)**
([`stats_sweeps_results.md`](stats_sweeps_results.md)) — n=50 CIs + operating-point
sweeps: most headlines strengthen, E7/E3 soften as expected, E3 composition confirmed
operating-point-contingent; P3b closed the σ-band/outcome-matrix gap (C2/C4/C5 hold
or strengthen; C7's normalized `do(χ)→content` is noisy at n=30 for methodological,
not substantive, reasons — a normalization caveat, not a reversal). Only the P4b
doc fold-in and the τ-axis sweeps remain.
**3b done** ([`e0_topologies.md`](e0_topologies.md)) — the self-sustaining band and
`period~τ` law generalise cleanly to `smallworld`/`rgg` (r≈0.999), *and* E1 reward
routing learns on every hidden medium (ring/smallworld/rgg, d=1.06–1.75); the 2-D
spiral is geometry-bound (out of scope) and only the E2/E5 port remains. **4a attempted, paused** ([`e10_notes.md`](e10_notes.md))
— the existing τ rule structurally can't build the hierarchy; it needs a new
bidirectional τ-plasticity rule (mechanism design), and E9 de-risks only the *grouping*
half, not the τ-value rule (earlier claim corrected). Remaining high-value, lower-risk
options: **1b** (learned direction-selective readout), the **E2/E5 port** onto
non-lattice substrates (3b follow-on; E1 done), and the deferred **3a P3b** σ-band
headlines. Return to **4a** only with appetite for mechanism design.

## Process notes (apply to whatever is chosen)

- Seed *everything* (the `perturb_tau` lesson); prefer `default_rng(seed)` threaded
  explicitly, never the global RNG.
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the dynamics do vs what a
  readout/feature does.
- Keep a caveats section adjacent to every headline.
