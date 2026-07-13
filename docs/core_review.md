# Core Series Review — Integrity & Overreach Audit (E0–E6, C0–C4)

*Independent audit of the work in this repository (E0–E6 and C0–C4 series),
checking for the failure modes an AI-assisted programme is prone to: fabricated
numbers, claims not backed by code or data, invented citations, and grandiose
theoretical framing disconnected from what the code actually does.*

Reviewer pass date: 2026-07-11. Branch: `claude/ai-hallucination-review-ya7aq5`.
Coverage: data-vs-doc numbers verified for all 13 experiments; citations
verified; per-seed / code-level audit completed for E3, E5, and the C-series
causal operators; full end-to-end reproduction run completed (see below).

## Bottom line

**No fabricated results were found.** Every headline number I checked
reproduces from the committed `.npz` data, every external citation I could
reach is real and accurately characterized, and the prose is consistently
hedged rather than overclaiming. This is an unusually honest body of
exploratory work — most result docs carry an explicit "Caveats" or "Honest
caveats" section that pre-empts the criticisms an auditor would raise.

The residual risk is **not** fabrication; it is the ordinary risk of an
ambitious framing (Buzsáki "inside-out", spike–wave causality) resting on
small hand-built toy simulations. The docs themselves flag this repeatedly.

**One material exception — E3's composition study.** Its numbers reproduce
exactly and its footnoted caveats are honest, but the *top-line framing*
overstates a fragile result: the headline "joint identity 0.20 → 0.77 (~quadrupling)"
is a mean over n=5 seeds that is actually bimodal (2 seeds at ceiling, 2 at
chance), reports only the identity axis, and rests on a hand-picked favorable
operating point. Genuine *joint* composition (identity **and** timing together —
the thing the study claims to partly achieve) occurs on **1 of 5 seeds**. This
is overreach in altitude, not fabrication; see the dedicated section below.

**A reproduction run also surfaced one genuine bug.** Re-running all 14
experiments from scratch, 12 of 14 reproduce bit-identically. The exceptions —
the Line-B / A+B conditions of E2 and E3 — trace to an **unseeded global RNG** in
the Line-B learning mechanism (`perturb_tau`, `ghca_learn.py:138/143`), so those
conditions' "seeds" are not controlled and their committed per-seed numbers are
not regenerable (E3's composition headline moved 0.77 → 0.65 on rerun). E2's
scientific conclusion survives (retention still 1.00); E3's precise composition
numbers do not. This compounds the E3 framing issue with a concrete
reproducibility defect — details in the "Reproduction run" section.

## Method

For each series I did three checks:

1. **Number match** — loaded every `result/*/*.npz` and compared the saved
   arrays against the specific quantitative claims in the corresponding doc.
2. **Citation check** — verified the external references that the arguments
   lean on actually exist and say what the docs claim.
3. **Prose/altitude** — read the result docs for claims the toy experiments
   cannot support, stated as fact rather than analogy.

Depth caveat: this pass verified **data-vs-doc integrity** and **citations**
directly, and read the C-series, E0, E2, E6 code/docs closely. A full
line-by-line audit of every experiment's code for subtle circular setup was
started but not completed (see "Not fully verified" below). Nothing found so
far suggests a problem there; it is simply not exhaustively cleared.

## Per-series verdicts

| Series | Data reproduces? | Verdict |
|--------|:----------------:|---------|
| E0 substrate characterisation | Yes (range-death bands r1→θ1, r2→θ4, r3→θ7; p_ss monotone in θ) | Supported |
| E1 conditioning | Yes (reward_A≈0.67, reward_B≈0.33≈chance, A+B≈0.68) | Supported |
| E2 working memory | Yes (mech mean 0.667=20/30; ret_A=0.20, ret_B/AB=1.0; learned τ_B≈12.6; info half-lives 161/161/161/20/20) | Supported |
| E3 timed response | Yes (all arrays match) | **Numbers supported; framing overstated** — see E3 deep-dive |
| E4 attention | Yes (headline numbers match saved arrays) | Supported |
| E5 executive control | Yes (headline numbers match saved arrays) | Supported |
| E6 emergent categories | Yes (mem R²=0.62, att=0.84, exe=0.98; cosines ≈0.01/0.03; own/other-region controls exact) | Supported |
| C0 instrumentation | Yes (gain_coll, accS_coll/lab arrays match) | Supported |
| C1 certificates | Yes (`all_ok=True`; certificate==measured==expected on all 6 graphs) | Supported |
| C2 fat-handed `do(W)` | Yes (labeled-line band 33.1σ vs collective 0.24σ) | Supported |
| C3 `do(θ)` well-posed | Yes (doW_band_lab=33.09, doW_band_coll=0.24, doTheta_band=0.014) | Supported |
| C4 outcome-relativity | Yes (macro-sufficiency 1.03 vs 0.11; diagonal `do(θ)` matrix) | Supported |

## Citations — all real, all accurately used

The citations are the classic place a model invents references. The three
load-bearing ones check out:

- **arXiv:2511.06602 — "A Causal Formulation of Spike-Wave Duality"**
  (Jalaldoust et al., Nov 2025). Real. It genuinely formalizes epiphenomenality
  via invariance of interventional distributions in SCMs, derives a do-calculus
  certificate of sufficiency, and invokes the Causal Hierarchy Theorem — exactly
  the framing the entire C-series builds on. The docs even note a real subtlety
  in operationalizing the paper's Definition 1 (interventional-vs-observational,
  not interventional-vs-interventional) which they say tripped up a first pass.
- **Casado Noguerales et al., 2026** — "What Does a Discrete Diffusion Model
  Learn?" (IEEE ISIT 2026). Real. It does frame the irreducible cost as the
  "information-destruction rate," which is precisely how the E2 addendum uses it.
- **Fisch–Gravner–Griffeath threshold-range scaling** (E0). A real, standard
  result in the excitable-CA literature; used correctly.
  Buzsáki (inside-out) and Sutton (reward) are used as framing, not as sources
  of specific quantitative claims.

## E3 deep-dive — honest numbers, overstated framing

E3 was audited to the per-seed level (code: `experiments/e3_timed_response.py`,
`experiments/e3_factored_credit.py`; data: `result/e3/e3_data.npz`,
`result/e3/e3_factored.npz`). No fabrication, no circular setup, no rigged
comparison — but this is the repo's weakest link and its top-line story rounds a
fragile result up to a clean one.

**What holds up.** The mechanism is clean (latency tracks gate τ: τ=[6,10,14,18,
22,26,30] → latency [5,10,13,17,22,25,29] ≈ τ−2). Line A learns identity (1.00 on
all seeds); the A-vs-B dissociation *logic* is sound (Line B has no routing
plasticity, so it structurally cannot learn identity).

**Where the framing outruns the data:**

1. **"Stable" means hide high variance (n=5).** Line B identity reported as
   "innate ~0.70" is per-seed `[0.50, 0.76, 1.00, 0.23, 1.00]` — one seed below
   chance, two at ceiling; it is an artifact of the fixed 85%-channel-biased
   random wiring, not a stable baseline. AB-shared "0.20" is `[0, 0, 0.48, 0.55,
   0]` — total collapse on 3/5 seeds.
2. **"Factored credit restores performance" — yet joint identity stays at
   chance.** AB-factored per-seed `[0.53, 0.48, 0.56, 0.48, 0.43]` (mean 0.50):
   every seed at chance. The 0.20→0.50 "lift" is just the removal of the
   catastrophic zeros, landing at chance — not restoration of the faculty.
3. **The headline 0.77 is bimodal and reports only one axis.** AB-factored+
   curriculum identity `[1.00, 0.51, 0.79, 1.00, 0.57]` (mean 0.77) = 2 seeds at
   ceiling, 2 at chance. "Composition" requires identity **and** timing together;
   the matching latency errors on those seeds are `[20, 9.3, 9.3, 0, 14.6]`
   (tolerance 1). Joint success (identity≈1 **and** timing in tolerance) holds on
   **exactly 1 of 5 seeds** (seed 3). The doc reports only the identity axis in
   the headline.
4. **Favorable operating point is hand-picked.** The factored study sets target
   latency = 16 *because* it forces gate τ into [17,19], which the code's own
   resonance map flags as identity-learnable (τ=17–19 → 1.0; τ=13–15 → 0.54,
   chance). Disclosed as a caveat, but it means the "resolution" is contingent on
   choosing a task whose timing lands in a lucky substrate resonance, and the
   curriculum does not steer τ there — it freezes wherever stochastic exploration
   landed.

**Assessment.** Overreach in altitude, not fabrication. All caveats exist in
the doc (it literally says "on the seeds where B's τ lands in the good zone" and
"residual fragility is a substrate resonance artifact"), but they are
subordinated to a confident headline ("~quadrupling," "fully composing"), and the
README compresses further to "decomposed & partly resolved, 0.20→0.77." Suggested
correction: state that joint composition succeeds on 1/5 seeds at a chosen
operating point, and that n=5 with this variance does not support a "partly
resolved" claim without more seeds and a non-resonant substrate.

## What is genuinely good (honesty markers)

- Result docs consistently distinguish measured fact from interpretation, and
  most end with a caveats section that names the weak points (e.g. C2: "the
  rigorous point is partly definitional — that is the intended message"; E6:
  "executive R² is high partly because switch timing is near-deterministic by
  construction"; E2: the resonance rule "lands on the boundary" and per-node τ
  hits a "weakest-link problem").
- Real baselines and controls are present and their numbers are in the data
  (E1 shuffle/baseline arrays; E6 own-region oracle vs other-region control;
  C1 observational-association contrast column).
- The README Progress table does **not** round hedged results up: E3 is marked
  "partly resolved… residual is a substrate resonance artifact," not "solved."

## Reproduction run — completed (turned up a real bug)

All 14 experiments were re-run from scratch and the regenerated `.npz` files
compared numerically against the committed originals. Result: **12 of 14
reproduce bit-identically** (E0, E1, E4, E5, E6, C0–C4, plus the E2 mechanism and
all Line-A conditions). The exceptions are **exactly and only the Line-B /
A+B conditions of E2 and E3** — and finding that isolated the cause.

**Root cause — an unseeded RNG in the Line-B learning mechanism.**
`GHLearner.perturb_tau` (`ghca_learn.py:138` and `:143`) draws its timescale
exploration noise from `np.random.standard_normal()` — the **global, unseeded**
NumPy generator — whereas everything else in the codebase uses a properly seeded
generator (`self.rng`, `np.random.default_rng(seed)`). That τ-perturbation *is*
the entire Line-B learning signal. Consequences:

- The `seed` argument does **not** control Line-B exploration, so the "5 seeds"
  for these conditions are not independent controlled replicates and the committed
  per-seed numbers are **not regenerable** — a rerun draws a different stream.
- This is the reproducibility half of the E3 problem, now with a concrete cause:
  the headline **AB-factored+curriculum composition moved 0.77 → 0.65 on rerun**,
  with the per-seed pattern reshuffling (`[1.0, 0.51, 0.79, 1.0, 0.57]` →
  `[0.74, 0.51, 1.0, 0.0, 1.0]` — seed 3 flipped from full success to total
  failure). The specific "0.20→0.77 ~quadrupling" rests on one uncontrolled draw.

**What survives the bug.** E2's scientific conclusion is robust: retention
`ret_B`/`ret_AB` reproduce **identically at 1.00** because whatever τ the unseeded
perturbation lands on, it stays below the loop transit time `L=24`, so memory
still holds (the τ values wobble, e.g. `[14.1,10.7,14.3,13.7,10.1]` →
`[11.4,11.1,14.5,12.3,13.1]`, but all `<24`). E3's *qualitative* story also
survives (shared-reward collapses; factored-alone lands at chance; curriculum
lifts it fragilely). What does **not** survive is the precise E3 composition
numbers and their per-seed breakdown.

**Fix (recommended, out of scope for this review):** thread `self.rng` through
`perturb_tau` instead of the global RNG, then re-run E2/E3 and update the numbers.
Everything else in the repo is already deterministic under its seeds.

**Undocumented dependencies (minor):** the README's install line lists only
`numpy matplotlib scipy`, but the suite also needs `networkx` (C1) and
`scikit-learn` (C0, C4) to run. Worth adding.

## Not fully verified (residual items)

1. **E0 period-fit `r=0.9992`** — still not independently re-derived (the
   `periods` array has `inf` entries for non-oscillating sweep points). The E0
   `.npz` reproduces bit-identically, so the underlying run is deterministic;
   only the reported regression constant is unconfirmed. Low-stakes.
2. **Code-level circularity audit — completed.** E3 (→ framing issue above), E5
   (→ holds up, note below), and the C-series `do(W)`/`do(θ)` operators in
   `ghca_causal.py` (→ clean, note below) have all been read. No hidden
   circularity found anywhere. The one presentational residual is in C3 (below),
   not a correctness or fabrication problem.

## E5 note — audited, holds up

E5 was checked to the same per-seed depth as E3 (code `experiments/e5_executive.py`,
data `result/e5/e5_data.npz`). Every number reproduces exactly under the code's own
windowing (switching final 0.891; ablated 0.199; single-rule steady-state
0.872/0.861 = last-150-trial mean; post-switch 0.567→0.917). Unlike E3 the per-seed
spread is **tight and unimodal** — switching `[0.96, 0.82, 0.92, 0.91, 0.85]`,
ablated `[0.21, 0.19, 0.25, 0.15, 0.19]` — so this is a robust effect, not a
lucky-seed mean.

The ablation is a **fair dissociation**, which was the main risk. The single-rule
control applies the *same* lesion (ring τ=18≥L) but re-cues the ring every trial;
it still works ablated (0.86), proving the lesion removes *held context*
specifically, not the routing machinery. Same lesion, one task collapses, the
other is spared — a clean 2×2. Two honest points (both disclosed in the doc, not
defects): the option is **cued, not discovered** (the gating structure is
pre-wired; only the H→M routing is learned), and ablated switching falls to 0.20
**below** chance because the hard-AND gate mutes the medium (silence commits
wrong). The latter makes the "0.89 vs 0.20" gap look more dramatic than a "vs 0.50
chance" framing would, but it is explained openly and the single-rule control
anchors the interpretation. Verdict: **Supported.**

## C-series note — operators clean; one presentational caveat in C3

The causal operators in `ghca_causal.py` were read directly. They are honest and
correct: `wave_coherence` / `wave_active_fraction` are genuine deterministic
coarse-grainings (`W=f(S)`); `do_theta` sets τ / θ / coupling directly (modular,
unique); `do_W` explicitly exposes the realization `policy` (random / clustered /
min_edit) as a free parameter — which *is* the fat-handedness being studied, not a
bug. All C2/C3 numbers reproduce exactly (C2 band 0.243 / 33.09; C3 doθ 0.01409;
E3-bridge latency-vs-τ slope 1.00).

Two things a reader should know about the dramatic σ figures:

- **C2's "33.1σ" is the *extremal* achievable band** — `(hi−lo)/std`, where the
  edges activate the k-smallest vs k-largest-weight nodes (an adversarial
  realization), standardized by that behaviour's own observational std. It is a
  real quantity, not a tiny-denominator artifact (raw labeled-line band 76.6 units
  vs collective 6.6; the extra factor in the 33-vs-0.24 ratio comes from each
  behaviour being standardized by its own std, which is the correct way to ask
  "how many of its natural SDs can the intervention move it"). The doc also
  reports the more modest *typical* random-realization spread (2.09σ), so both the
  maximal and typical readings are disclosed. The C2 doc itself concedes the core
  point is "partly definitional." Fair.
- **C3's "0.014σ" is a *different quantity* from C2's band**, and the "~2400× less
  ambiguous" headline places the two on one axis. C2's band is a realization-DOF
  spread ÷ observational-std; C3's `doTheta_band` is `rep_c / std(Ec)` =
  across-seed *estimation noise* ÷ *effect size* (0.146 / 10.36). These are not
  the same kind of σ, so the single-number 2348× ratio and the shared-axis bar
  chart slightly over-dramatize comparability. **However**, the surrounding prose
  states the real, correct point explicitly: `do(W)`'s band is *irreducible* (a
  free policy choice, does not shrink with samples) while `do(θ)`'s realization
  band is *zero by construction* (setting τ=t is unique), leaving only ordinary
  estimation error that shrinks with samples. That qualitative claim is sound and
  honestly explained; only the headline number dramatizes it.

Verdict: C0–C4 **Supported**; C3 carries a minor presentational caveat (the
cross-quantity "2400×"), which the doc's own prose already defuses. Not
fabrication — a framing choice with the honest explanation adjacent.

## Recommendation

Treat the work as **sound and honestly reported** for what it is: an
exploratory computational-neuroscience study using toy Greenberg–Hastings
substrates to illustrate (not prove) claims about memory/attention/executive
function as readouts, and to instantiate an external causal-inference framework
on synthetic ground truth. No cleanup for fabrication is required. The one
concrete follow-up worth doing is a clean reproduction run of all experiments
to close item 1 above.
