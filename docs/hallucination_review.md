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
| E3 timed response | Yes (headline numbers match saved arrays) | Supported |
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
3. **Code-level circularity audit.** A careful read for subtle circular setup
   (result baked into construction) across every experiment was begun but not
   completed. The C-series `do(W)`/`do(θ)` operators in `ghca_causal.py` and the
   E5 ablation comparison are the two places most worth a second close read,
   since both are constructed contrasts where the conclusion could in principle
   be built into the world script. Nothing seen so far indicates that, and the
   docs openly acknowledge the near-deterministic-by-construction cases.

## Recommendation

Treat the work as **sound and honestly reported** for what it is: an
exploratory computational-neuroscience study using toy Greenberg–Hastings
substrates to illustrate (not prove) claims about memory/attention/executive
function as readouts, and to instantiate an external causal-inference framework
on synthetic ground truth. No AI-hallucination cleanup is required. The one
concrete follow-up worth doing is a clean reproduction run of all experiments
to close item 1 above.
