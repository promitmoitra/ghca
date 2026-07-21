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

### 1b. A learned direction-selective readout — ✅ **DONE** (see [`e7_direction_readout_results.md`](e7_direction_readout_results.md))
- **What.** Replace E7's *computed* local-winding readout with a small population that
  *learns* to read rotation direction from the wave (direction-selective cells via
  reward/Hebbian).
- **Why.** Makes E7's "read the wave" genuinely on-substrate; maps onto real cortical
  direction selectivity.
- **What was built.** A population of local elementary motion detectors (EMD,
  Hassenstein–Reichardt: delayed coincidence `active[t−1] at i AND active[t] at i+d`,
  coarse-grained) with a **learned** linear pool. It recovers chirality as well as
  `local_winding` on a centred core (1.00 vs 0.99) and — because it pools local
  detectors instead of reading one locus — stays perfect under core displacement
  (d=8: **1.00 vs 0.17**; d=12: **0.98 vs 0.07**), escaping the [C5](c5_results.md)
  fixed-locus collapse. Dropped into E7's switching it drives routing to **0.72**
  (vs computed 0.78). **Retires** the *computed-integral* and *fixed-locus* parts of
  the readout-honesty gap; the EMD primitive itself is still hand-specified and the
  pool is label- (not reward-) trained — both deferred.
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

### 3c. Continual learning as causal credit assignment — ✅ **P1–P4 DONE (capacity, not credit)** (plan [`continual_learning_plan.md`](continual_learning_plan.md); results [`continual_learning_results.md`](continual_learning_results.md))
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
- **P1 done.** Harness + metrics + baselines (n=30, correlational credit, remapping
  triple on one substrate): catastrophic interference is real (backward transfer
  −0.17…−0.26; avg acc ≈ chance), and **freezing the representation does not rescue it**
  (frozen ≈ or worse than plastic) because the tasks conflict on the shared head — the
  reservoir "freezing helps" claim needs per-task heads.
- **P2 done — honest null.** Causal `do(θ)` (weight-perturbation) credit vs the
  correlational eligibility trace on the K=2 reversal (head-only regime, n=30): causal
  *appears* to halve forgetting (BWT −0.38 vs −0.78, d=1.36), but the stability–
  plasticity **frontier control** shows both rules lie on the *same* curve — the effect
  is an effective-learning-rate difference, not better credit. The interference is
  **representational** (a shared linear head can't hold two anti-correlated mappings),
  not a credit-assignment artifact. Rules out "assign credit causally" as a CL fix.
- **P3 done — the null holds.** A low-variance causal estimator (antithetic central-
  difference perturbation) reaches higher acquisition but **also lies on the same
  frontier** — all three credit rules trace one stability–plasticity curve. So the P2
  null was *not* a variance artifact: the frontier is a genuine **capacity boundary**,
  and interference here is representational, not credit-assignment. Credit quality sets
  *where* on the frontier you sit, not the frontier itself.
- **P4 done — capacity, not credit (positive complement).** Holding the correlational
  credit rule fixed and varying *capacity*, interference collapses monotonically:
  backward transfer −0.78 (shared head) → −0.26 (shared head + task-context conjunction)
  → 0.00 (per-task heads); avg acc 0.50 → 0.66 → 0.96 (n=30, CIs separated). The
  single-head context rung ties directly to E9's conjunction mechanism. So 3c's answer
  is unambiguous: interference here is representational-capacity-limited, and the lever
  is conjunctive/contextual representation, not credit.
- **E9↔3c bridge done.** Grew the (stimulus × context) basis with E9's competitive-
  Hebbian rule and ran the sequential reversal on E9's substrate: the **learned
  (emergent) conjunction basis gives ≈zero forgetting** (backward transfer 0.001, avg
  0.78), approaching the hand-**wired** ceiling (0.89), while a non-conjunctive
  **frozen** basis can't even represent the tasks (0.24). Learned conjunctive tiling
  is the lever — afforded→learned (E9) and continual learning (3c) are the *same*
  mechanism, not credit.
- **P5 done — a fixed conjunction basis saturates.** Swept the number of sequential
  tasks T (balanced-dichotomy tasks, K_stim=4, binary action, frozen substrate,
  correlational credit, n=20). A fixed (stimulus × context) basis (n_h=50) matches
  full capacity at T=2 (avg 0.64 ≈ per-task 0.61) but its average accuracy falls
  **monotonically to the shared interference floor by T=6** (0.452 vs 0.451), while
  per-task heads (capacity ∝ T) stay flat at ~0.61. Crossover at T≈3 (≈12 conjunctions
  on 50 units) locates the basis's capacity. Closes the argument in both directions:
  capacity is the lever (P4/bridge), learnable but **finite** — enough tasks restore
  the interference. (Read avg accuracy, not backward transfer: as the basis saturates
  there is less left to forget, so BWT drifts vacuously to 0 — the E9-frozen floor
  artifact again.)
- **Deferred.** **Partially-overlapping tasks** as the fair test of whether credit
  quality *ever* matters; temporally-extended credit for a true hindsight estimator.
  Both are now best pursued *through* 3d (below), which subsumes the overlapping-task
  question into a sharper form.
- **Connects to.** E6 (frozen baseline), E9/E4 (WTA gating), C2/C3 (`do(θ)` well-posed),
  Line A/B plasticity; **3d** (raising the P5 ceiling along the timescale axis).

### 3d. Timescale as a continual-learning capacity axis — ✅ **WIRED + EMERGENT ARMS DONE** (results [`continual_learning_results.md`](continual_learning_results.md))
- **Emergent arm done (n=20) — the 4a ratchet is escaped.** A local, reward-free,
  **input-timing-driven** `τ`-plasticity rule (each node nudges `τ` toward the *external*
  probe delays it responds to — not its own inter-fire interval — with E9 k-WTA +
  conscience) *grows* a temporal basis from a near-homogeneous start: delay-decode 0.94
  (vs homog 0.53, wired 1.00), and continual-learning per-task 0.53–0.56 — well above the
  homogeneous floor (0.43), recovering ~70% of the wired basis's capacity (0.59). This is
  the **first result where the substrate's own *dynamics* (not a readout) are shaped by
  experience** → the plastic-dynamics caveat (1c/E9/Line B), addressed; and it unblocks
  Track 4a (timescale plasticity works once the teaching signal is external). Gap to
  wired (grown tiling is coarser — a residual cluster at the `τ` init) is a tuning item,
  not a mechanism problem. Follow-ups in **3e**.
- **Wired arm done (n=20).** A *hand-set* timescale-diverse (graded `τ`) basis buys
  genuine, forgetting-free continual-learning capacity on temporal (delay-keyed) tasks
  that a homogeneous basis cannot represent at *any* head capacity: graded+per-task
  0.58–0.62 (bwt ≈ 0, flat over T) vs homog+per-task ~0.43 (chance floor). The gap is
  purely the `τ` distribution → **representational, not capacity**. A shared head over
  the fixed graded basis still forgets (graded+shared ~0.46–0.50, bwt −0.12→−0.06),
  replicating P5's spatial-axis interference on the temporal axis. **The 4a go/no-go
  gate is therefore green** — the emergent hierarchy has a demonstrated capacity payoff,
  so the mechanism-design work is justified. (Mechanism: stimulus → per-node
  refractoriness for `τ` steps → a thermometer code of elapsed time, readable only when
  `τ` is diverse; delay-decode 1.00 graded vs 0.53 homogeneous.) Caveat: wired `τ` does
  not retire the plastic-dynamics caveat — only the emergent arm would; and the tested
  conjunction is essentially *(time)*, not the full *(stimulus × context × time)*.
- **What.** 3c/P5 showed capacity is the continual-learning lever and that a *fixed
  spatial* (stimulus × context) conjunction basis has a **finite** ceiling (it saturates
  to the interference floor by T≈6 on `n_h=50`). The non-"cheating" way to raise that
  ceiling — i.e. *without* per-task heads, which just staple on capacity task-by-task — is
  a **richer, higher-dimensional shared basis**. The substrate has an untapped
  representational axis the spatial basis ignores: **timescale** (`τ`, `act`). Test
  whether a **(stimulus × context × timescale)** conjunction basis raises the P5
  saturation ceiling — i.e. whether *temporal* diversity buys continual-learning capacity.
- **Why it matters (and why now).** It turns two loose threads into one experiment:
  1. It **gives Track 4a a downstream payoff.** 4a's fast/slow hierarchy stops being
     "nested rhythms as a phenomenon" and becomes *the mechanism that grows a temporal
     capacity axis for continual learning* — a concrete functional benefit, which is
     exactly the validation target `e10_notes.md` asks for ("a hierarchical network
     tracks a two-timescale signal better than a `τ`-homogeneous one", here sharpened to
     "learns more sequential tasks before saturating").
  2. It **subsumes 3c's deferred "partially-overlapping tasks" question.** A purely
     *spatial* conjunction basis cannot separate two tasks that share the same stimulus
     *and* context but differ in **timing / temporal structure** (the E3 timed-response
     regime) — same hidden pattern, conflicting required output. A timescale-diverse
     basis *can* (fast vs slow cells integrate the same input over different windows, so
     their conjunction tiles the temporal dimension too). So **temporally-structured
     tasks are the task family** where the timescale axis should pay off and the spatial
     axis provably can't — a cleaner test than generic partial overlap of "does anything
     beyond spatial capacity help?".
- **Design (the afforded→learned ladder, mirroring the E9 bridge — this is what
  de-risks it).** Reuse the P5 saturation harness (`continual_saturation.py`), swap the
  task family for temporally-structured tasks (share stimulus/context, differ in required
  timing), and compare bases by their timescale content:
  - **homogeneous** — all hidden `τ`/`act` equal (P5's basis). The control: predicts
    early saturation, no temporal separation.
  - **wired timescale-diverse** — ✅ **done** (`continual_temporal_saturation.py`,
    n=20). *Hand-set* graded `τ` across hidden units; decoupled from 4a's blocked
    `τ`-rule. The cheap early kill returned **green**: the graded basis learns temporal
    tasks (per-task 0.58–0.62, bwt ≈ 0) where homogeneous `τ` floors at chance (~0.43)
    at any head capacity — timescale diversity is the representational lever.
  - **emergent timescale-diverse** — ✅ **done** (`continual_temporal_emergent.py`,
    n=20). Grown by the input-timing-driven `τ` rule + E9 grouping (the input-tracked
    rule 4a was blocked on, implemented for the delay-tiling setting). Per-task 0.53–0.56
    vs homog 0.43 vs wired 0.59 — recovers most of the capacity and **retires the
    plastic-dynamics caveat** (dynamics, not just readout, shaped by experience).
  Metric: the P5 sweep (avg accuracy / backward transfer vs number of sequential tasks T),
  with the ceiling/crossover-T of each basis as the headline.
- **Dependency on 4a, made explicit.** The *wired* arm needs nothing new and can run
  immediately; it is the honest way to establish the capacity claim before 4a's hard
  mechanism-design work. The *emergent* arm consumes all three 4a ingredients — E9
  grouping (ready), channel-conditioned `act` (PR #42; here it doubles as a capacity
  lever, giving fast/slow cells distinct temporal receptive fields, not only rebalancing
  the hierarchy-formation competition), and the input-tracked `τ` rule (still 4a's one
  genuinely-required, unbuilt change). The E0 minimum-viable-`act` propagation check from
  PR #42 applies here too.
- **Effort.** Medium for the wired arm (task family + basis variant on an existing
  harness); High for the emergent arm (gated on 4a). **Risk.** Medium and honest: the
  timescale axis only adds *separable* capacity if the tasks carry real temporal
  structure — on purely static tasks it is inert, and that is the explicit kill
  condition. State the substrate/analysis boundary: this is a fixed-dynamics basis with a
  plastic *readout* unless/until the emergent arm makes the dynamics themselves plastic.
- **Connects to.** 3c/P5 (the ceiling it tries to raise), 4a (the mechanism it consumes
  and motivates), E3 (temporal task structure), E8.5 (hand-set two-timescale precedent),
  E9 (the frozen/wired/emergent ladder it copies).

### 3e. What the *emergent* timescale mechanism opens up — 📐 **PROPOSED** (directions now unlocked — 3d's emergent arm has landed)
- **The honest boundary first.** For the narrow capacity claim ("does timescale
  diversity buy continual-learning capacity?"), the *emergent* arm adds nothing the
  *wired* arm didn't already prove — a hand-set graded-`τ` basis suffices. The emergent
  arm's value is elsewhere, and it is what makes these directions possible:
  1. **It retires the deepest open tension (plastic *dynamics*).** Every "learned"
     result so far — E1 routing, E9 conjunctions, E8/1c prediction — is a plastic
     *readout* over *fixed* dynamics; even E9 left the timescales hand-set. The emergent
     arm is the first time the substrate's **own dynamics** (`τ`) are shaped by
     experience: *afforded dynamics + learned readout* → **learned dynamics**. This is
     the "plastic dynamics deferred" caveat (1c, E9, Line B) finally addressed.
  2. **The rule is the scientific content.** It shows the 4a ratchet was never intrinsic
     to GH timescale plasticity — only to the *self-referential* signal (a node reading
     its own inter-fire interval). An **external** teaching signal (the probe delay) +
     a bidirectional update dissolves it. Transferable, and what actually unblocks 4a.
  3. **Self-calibration.** The wired basis needs a designer who knows the delay
     statistics in advance (hand-set `linspace`); the emergent basis *discovers* the
     relevant timescale range from the input stream — the inside-out thesis doing real
     work.
- **Directions it unlocks** (ranked by how much they exploit what wired categorically
  cannot do):
  1. **Re-tiling under a shifting delay distribution — ✅ DONE** (`continual_temporal_retile.py`,
     n=20; results [`continual_learning_results.md`](continual_learning_results.md)). Under
     a SHORT→LONG→SHORT delay schedule the emergent basis **re-tiles** (LONG-decode
     0.50→0.92; `τ` migrates bodily; wired-static frozen at 0.50 — the adaptation a
     hand-set basis cannot do). And the plastic **representation has its own
     stability–plasticity frontier**: adapting to LONG costs SHORT (representation
     backward transfer −0.069 [−0.093, −0.045]) — 3c-style interference one level down,
     in the dynamics not the readout. But it is **graceful, not catastrophic** (−0.07 vs
     the readout's −0.4…−0.8; SHORT recovers, LONG partly retained): distributing the
     load across `n_h` units *accumulates* coverage rather than overwriting — capacity
     again softens interference. Deeper form deferred: a fully-online concurrent (not
     phase-split) version.
  2. **Bimodal hierarchy from two-rhythm drive — ✅ DONE, hierarchy half of 4a closed**
     (`timescale_hierarchy.py`, n=20; results [`timescale_hierarchy_results.md`](timescale_hierarchy_results.md)).
     The input-timing rule clusters `τ` at *both* drive periods (near-P_f/near-P_s
     fractions 0.50/0.50; fully **emergent** — nodes fan in from both sources, a
     population conscience splits them ~50/50 and each group locks to its channel),
     where the old self-referential rule forms **no fast cluster** (near-P_f 0.00 — the
     e10 ratchet, reproduced and overturned). The conscience — not channel-conditioned
     `act` — sufficed to fix diagnostic-2's swamping. **Remaining for a *full* 4a close:**
     theta–gamma-style **cross-frequency coupling**, which needs an added inter-population
     pathway (E8.5 nested-waves direction) this pool lacks — the coupling, not the
     hierarchy, is what's left.
  3. **Concurrent co-adaptation (end-to-end inside-out).** Both 3d-emergent and E9 use a
     phase split (grow representation, freeze, learn readout). Run emergent-`τ` and the
     reward readout *together* — do they co-adapt or fight? The strongest form of "one
     homogeneous machine learning end-to-end"; retires the last E9/3d phase-split caveat.
     **Effort** medium–high. **Risk** high (two plastic loops can destabilise).
  4. **Supply E2/E3's temporal machinery from experience.** The emergent basis is a
     self-organised population of *time cells* tiling a delay — exactly what E2 (working
     memory) and E3 (timed response) hand-built. Could retire *their* afforded components,
     extending afforded→learned beyond conjunction cells into temporal cognition.
     **Effort** medium. **Connects** tension 1 to the memory/timing experiments.
  5. **A falsifiable neuroscience prediction (Track 2 / tension 3).** If delays are
     non-uniform, does the rule allocate more `τ`-resolution where intervals are common?
     That is adaptive/efficient temporal coding — a concrete prediction that biological
     time cells tile *experienced* interval statistics, not uniform time. **Effort** low
     (a non-uniform-delay variant of the existing sweep + a tiling-vs-statistics readout).
- **Conceptual thread (ties the two arcs tighter).** C3 established do(`θ`) — intervening
  on timescales — is the well-posed causal handle. The emergent rule is the substrate
  performing **self-directed do(`θ`)**: adjusting its own `τ` from input timing. So `θ`
  is both *the* causal handle (C-series) and the variable that self-organises (E-series)
  — a sharper form of the "learning as causal inference" framing, and a candidate 5a
  synthesis note once 3d-emergent lands.
- **Connects to.** 3d (the arm that unlocks these), 4a (direction 2 closes it), 3c
  (direction 1 extends it), E2/E3 (direction 4), C3/`do(θ)` (the conceptual thread).

---

## Track 4 — New territory *(novelty)*

### 4a. Emergent timescale hierarchy / cross-frequency coupling — ✅ **HIERARCHY DONE (bimodal split); CFC remaining** (results [`timescale_hierarchy_results.md`](timescale_hierarchy_results.md); orig. diagnosis [`e10_notes.md`](e10_notes.md))
- **What.** With per-node τ plastic, do learned τ distributions self-organise into a
  fast/slow hierarchy with theta–gamma-style cross-frequency coupling?
- **Why.** Deepest genuinely-new phenomenon the substrate could show; connects to the
  nested-waves literature (E8.5 only used two hand-set timescales).
- **Status — hierarchy delivered (3e.2).** The blocker was the `τ`-value rule: the old
  Line B rule only ratchets τ *upward* (self-referential — reads a node's own inter-fire
  interval, corrupted once τ overshoots). 3d-emergent's fix (bidirectional,
  *input-timing-driven*) + a population conscience now **grows the fast/slow hierarchy**
  under two-rhythm drive: τ clusters at *both* drive periods (near-P_f/near-P_s 0.50/0.50,
  n=20), fully emergent (nodes fan in from both sources; conscience splits them ~50/50 —
  the fix for diagnostic-2 swamping), where the old rule forms **no fast cluster**
  (near-P_f 0.00). Read the near-period fractions, not Sarle BC (the old rule clears BC
  too, but with its clusters ratcheted to high τ, not at the fundamentals). **Remaining:**
  theta–gamma **cross-frequency coupling** — this pool has no inter-population pathway, so
  none is expected; adding one (E8.5 nested-waves direction) is the last piece for a
  *full* 4a close.
- **Proposed synthesis (2026-07-19, not yet attempted).** A design discussion
  surfaced a third ingredient, prompted by asking whether making the **active**
  duration `act` per-node/tunable (currently a global scalar; only the passive/
  refractory tail is per-node) would help. It does **not** replace the required
  fix above — the ratchet's root cause is self-referential (a node measuring its
  *own* inter-fire interval, which is corrupted once τ overshoots the true period),
  and that pathology reappears in any total-cycle parameter tuned from the same
  signal, `act` included. But `act` is a genuinely useful **third, complementary**
  lever if scoped narrowly:
  - **Channel-conditioned, not freely learned.** Once E9-style k-WTA + conscience
    assigns a node to a fast/slow channel, set `act` to a fixed small/large value
    *for that channel* rather than inventing a new plasticity rule for it — lower
    risk, and sidesteps the self-referential-measurement problem entirely (nothing
    is learning `act` from its own behaviour).
  - **Directly attacks diagnostic 2's actual failure.** The competitive prototype
    died because the fast channel (33% active) swamped the slow channel (8%) in
    the k-WTA competition. Shrinking `act` for the fast channel shrinks its active
    *footprint* too, naturally rebalancing the competition instead of hand-tuning
    drive amplitude/duty-cycle (the fix the original notes proposed).
  - **A second, correlated separation axis.** Fast vs slow channels would then
    differ in both `τ` *and* `act` (short pulse ↔ fast, long pulse ↔ slow), which
    should make the bimodal clustering more robust than separating on `τ` alone.
  - **New risk to check before building.** `act` also sets wavefront width / how
    much drive a firing node hands its neighbours — this is exactly what E0's
    threshold-range characterisation is about. Shrinking `act` too far for the
    fast channel risks under-driving neighbours below θ and **breaking
    propagation entirely**, not just changing rhythm. Any implementation needs an
    E0-style check for a minimum viable `act` at the substrate's operating point
    before assuming the fast channel can be made arbitrarily brief.
  - Net: three ingredients, not two — (i) E9's grouping [already validated,
    reusable as-is], (ii) channel-conditioned `act` [new, low-risk, fixes the
    diagnostic-2 imbalance], (iii) the input-tracked `τ` rule [still the one
    genuinely required change; unblocks nothing by itself if skipped].
- **Downstream payoff — now demonstrated (gate green).** **3d** (Track 3) *consumes*
  the emergent hierarchy as a **temporal capacity axis for continual learning**. Its
  *wired* (hand-set timescale) arm has now run (n=20) and returned **green**: a
  graded-`τ` basis buys forgetting-free continual-learning capacity on temporal tasks
  that a homogeneous basis cannot represent at any head capacity. So 4a's emergent
  hierarchy has a **demonstrated functional payoff** (sharpening `e10_notes.md`'s own
  validation target), and the go/no-go check that could have killed it *before* the hard
  `τ`-rule work instead cleared it. Building the emergent hierarchy is justified.
- **Connects to.** E2, E5, E8.5; E9 (competition, for the grouping half only); E0
  (threshold-range scaling, for the proposed `act` floor constraint); **3d** (the
  continual-learning capacity consumer that motivates the emergent hierarchy).

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
| **Credibility** (retire the honesty gap) | **1a** emergent conjunctions | 1b ✅, 3a ✅ |
| **Clean disambiguation** | **2a** predictive-coding foil | 1c |
| **Scientific novelty** | **4a** emergent timescale hierarchy | 2b |
| **External impact / reach** | **4b** causal testbed | 5a write-up, 2b |
| **Lowest-risk strengthening** | **3a** stats/sweeps ✅ + **5b** hygiene ✅ done | 1b ✅, 3b ✅ |

**Progress.** **1a done (E9)** — retired the most-cited caveat. **1b done**
([`e7_direction_readout_results.md`](e7_direction_readout_results.md)) — a learned
population of local motion detectors replaces E7's computed winding readout, matching
it on centred cores and staying robust under core displacement where the fixed-locus
integral collapses (retires the computed-integral / fixed-locus honesty gap; C5 escaped).
**1c / 2a done**
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
options: the **E2/E5 port** onto non-lattice substrates (3b follow-on; E1 done), and
the deferred **3a P4b** doc fold-in / τ-axis sweeps. Return to **4a** only with
appetite for mechanism design.

## Process notes (apply to whatever is chosen)

- Seed *everything* (the `perturb_tau` lesson); prefer `default_rng(seed)` threaded
  explicitly, never the global RNG.
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the dynamics do vs what a
  readout/feature does.
- Keep a caveats section adjacent to every headline.
