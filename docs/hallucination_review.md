# AI Hallucination / Overreach Review

*Independent audit of the work in this repository (E0–E6 and C0–C4 series),
checking specifically for AI-generated hallucination: fabricated numbers,
claims not backed by code or data, invented citations, and grandiose
theoretical framing disconnected from what the code actually does.*

Reviewer pass date: 2026-07-11. Branch: `claude/ai-hallucination-review-ya7aq5`.

## Bottom line

**No hallucinated results were found.** Every headline number I checked
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

**Assessment.** Overreach in altitude, not hallucination. All caveats exist in
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

## Not fully verified (residual items, not defects)

These are honest gaps in *this review's* coverage, not identified problems:

1. **End-to-end re-runs.** I verified that the committed data matches the docs,
   not that re-running each experiment from scratch regenerates that data. The
   code is present and self-contained; a reproduction pass (`python3
   experiments/*.py`) would close this. Cheap to do and recommended.
2. **E0 period-fit `r=0.9992`.** The `periods` array in `e0_data.npz` contains
   `inf` entries (non-oscillating sweep points), so the exact regression
   `period=1.00·τ+0.95, r=0.9992` was not independently re-derived from the
   saved data. Plausible and low-stakes, but not confirmed here.
3. **Code-level circularity audit — mostly completed.** E3 and E5 have now been
   read to the per-seed level. E3 → framing issue above. **E5 holds up cleanly**
   (see E5 note below): the ablation is a fair, well-controlled dissociation and
   the per-seed spread is tight, not bimodal. The one remaining code read is the
   C-series `do(W)`/`do(θ)` operators in `ghca_causal.py` — constructed contrasts
   (33σ vs 0.014σ) where the effect could in principle be partly definitional; the
   C2 doc half-concedes this ("partly definitional — that is the intended
   message"). Nothing seen so far indicates hidden circularity.

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

## Recommendation

Treat the work as **sound and honestly reported** for what it is: an
exploratory computational-neuroscience study using toy Greenberg–Hastings
substrates to illustrate (not prove) claims about memory/attention/executive
function as readouts, and to instantiate an external causal-inference framework
on synthetic ground truth. No AI-hallucination cleanup is required. The one
concrete follow-up worth doing is a clean reproduction run of all experiments
to close item 1 above.
