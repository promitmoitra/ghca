# Spatial computing (Miller lab) vs. this substrate — a literature note

**Status:** literature note, **citations unverified**. Prompted by a Picower
Institute news piece on Earl K. Miller's proposal that cognition and
consciousness arise from *analog computation* by brain waves.

> ## ⚠ READ THIS BEFORE USING ANY CITATION BELOW
>
> **Not one citation in this note was verified against a primary source.** The
> session that wrote it had **every** scholarly host blocked by its egress
> policy — `nature.com`, `cell.com`, `pnas.org`, `sciencedirect.com`,
> `biorxiv.org`, `arxiv.org`, `pubmed.ncbi.nlm.nih.gov`, `pmc.ncbi.nlm.nih.gov`,
> `europepmc.org`, `api.crossref.org`, `api.openalex.org`, `dspace.mit.edu`,
> `ekmillerlab.mit.edu` — and `picower.mit.edu` itself, so **the news article
> that prompted this note was never read**. Every citation derives from web
> *search-result listings* (titles, author lists, venues, URLs) plus
> search-engine summaries of them.
>
> This violates the repo's normal standard ("if a number appears in a doc, a
> script must compute it"; generate tables, never hand-type). It is committed
> anyway, deliberately and flagged, because the *repo-side* analysis is sound
> and independently checkable — but **every external claim needs one
> confirmation pass against the actual paper before it is quoted anywhere
> else**, and any DOI here should be re-derived, not copied.
>
> Confidence is marked per citation. Nothing marked `uncertain` should be
> cited at all.
>
> **Also unrun:** the independent adversarial review. Three critic passes
> (overreach, technical correctness, completeness) were commissioned and all
> three died on a spend limit. §5 is a **self-audit by the note's own author**,
> which the repo's own convention says is not a substitute
> ("self-authored work is not self-reviewed").

Repo-side numbers in §3–§4 **are** verified: each was read out of the cited
results doc or archive in this repo during writing.

---

## 1. What the theory is, and what is actually peer-reviewed

The distinction matters more than usual here.

**The strong claim — "cognition and consciousness *are* analog computation by
brain waves" — appears only in a lecture and in press coverage of it.** The
talk is Miller's SfN 2025 presidential lecture (reported as 15 Nov 2025). There
is no peer-reviewed paper making that claim. The talk synthesises four strands
that *do* have primary papers, and those papers support materially narrower
statements:

| strand | what the papers actually claim |
| :--- | :--- |
| Working memory as **discrete bursts**, not persistent firing | gamma (45–100 Hz) bursts accompany encoding/reactivation at a *minority* of sites; beta (20–35 Hz) bursts are a default state interrupted by encoding/decoding |
| Beta as **top-down control**, with a laminar address | deep-layer alpha/beta modulates superficial-layer gamma; predictions fed back as alpha/beta pre-inhibit predicted feedforward routes |
| **Spatial computing** — the named framework | alpha/beta is organised into low-dimensional *spatial* patterns that reorganise with task, and are *inversely* correlated with where sensory spiking appears — "inhibitory stencils" |
| **Travelling waves** as the medium | PFC oscillations form travelling waves in theta/alpha/beta; most are **rotating**, not planar; direction preference shifts with task |

**The gap between those four and "consciousness emerges from analog wave
computation" is entirely interpretive.** This note treats the narrow claims as
the object of comparison and flags the strong claim as unsupported by the
peer-reviewed record it cites. That is not a criticism of the lecture — a
lecture is allowed to synthesise — but a note in this repo must not launder a
talk's framing into a citation.

### The consciousness link is observational

The consciousness half rests on anesthesia work, and its causal structure is
worth stating plainly: the manipulated variable is a **systemic drug**;
oscillations, phase alignment, traveling-wave organisation and dynamic
stability are all **dependent** variables measured alongside a responsiveness
readout. That is a common-cause structure —

```
    drug ──> waves
      └────> unresponsiveness
```

— and it does not isolate the waves as mediator. **Nothing in that line
intervenes on the waves themselves.** The one genuine in-brain intervention is
electrical stimulation of central/intralaminar thalamus reversing anesthesia,
which manipulates *upstream generating structure*. Hold that thought for §4.1.

## 2. Citations, with confidence

`search-listing` = title, authors and venue were seen in a search result.
`partial` = some fields unconfirmed. **None** is `verified-against-source`.

| # | citation | conf. |
| :-- | :--- | :--- |
| 1 | Lundqvist, M., Rose, J., Herman, P., Brincat, S. L., Buschman, T. J., & Miller, E. K. (2016). Gamma and Beta Bursts Underlie Working Memory. *Neuron*, 90(1), 152–164. | search-listing (DOI not seen) |
| 2 | Lundqvist, M., Herman, P., Warden, M. R., Brincat, S. L., & Miller, E. K. (2018). Gamma and beta bursts during working memory readout suggest roles in its volitional control. *Nature Communications*, 9, 394. | search-listing |
| 3 | Lundqvist, M., Brincat, S. L., Rose, J., Warden, M. R., Buschman, T. J., Miller, E. K., & Herman, P. (2023). Working memory control dynamics follow principles of spatial computing. *Nature Communications*, 14(1), 1429. | search-listing |
| 4 | Chen, Z., Brincat, S. L., Lundqvist, M., Loonis, R. F., Warden, M. R., & Miller, E. K. (2026). Oscillatory control of cortical space as a computational dimension. *Current Biology*, 36(2), 402–414.e5. | partial |
| 5 | Miller, E. K., Lundqvist, M., & Bastos, A. M. (2018). Working Memory 2.0. *Neuron*, 100(2), 463–475. | search-listing |
| 6 | Bastos, A. M., Loonis, R., Kornblith, S., Lundqvist, M., & Miller, E. K. (2018). Laminar recordings in frontal cortex suggest distinct layers for maintenance and control of working memory. *PNAS*. | partial |
| 7 | Bastos, A. M., Lundqvist, M., Waite, A. S., Kopell, N., & Miller, E. K. (2020). Layer and rhythm specificity for predictive routing. *PNAS*, 117(49), 31459–31469. | search-listing |
| 8 | Bhattacharya, S., Brincat, S. L., Lundqvist, M., & Miller, E. K. (2022). Traveling waves in the prefrontal cortex during working memory. *PLoS Computational Biology*, 18(1), e1009827. | search-listing |
| 9 | Bastos, A. M., Donoghue, J. A., Brincat, S. L., et al. (2021). Neural effects of propofol-induced unconsciousness and its reversal using thalamic stimulation. *eLife*, 10, e60824. | search-listing |
| 10 | Xiong, Y. S., Donoghue, J. A., Lundqvist, M., et al. (2024). Propofol-mediated loss of consciousness disrupts predictive routing and local field phase modulation of neural activity. *PNAS*, 121(42), e2315160121. | search-listing |
| 11 | Eisen, A. J., Kozachkov, L., Bastos, A. M., et al. (2024). Propofol anesthesia destabilizes neural dynamics across cortex. *Neuron*, 112, 2799–2813. | search-listing |
| 12 | Lundqvist, M., Miller, E. K., Nordmark, J., Liljefors, J., & Herman, P. (2024). Beta: bursts of cognition. *Trends in Cognitive Sciences*, 28(7), 662–676. | search-listing |
| 13 | Redinbaugh, M. J., Phillips, J. M., Kambi, N. A., et al. (2020). Thalamus Modulates Consciousness via Layer-Specific Control of Cortex. *Neuron*. | partial |

**Attribution warning caught during research:** Liljefors et al. (2024,
*Nature Communications*) tests the beta/alpha gating account in humans but
**does not include Miller as an author** — do not attribute it to his lab.

Not in the citation set, flagged as hypothesis-only rather than data: the
cytoelectric-coupling and ephaptic-field papers (Pinotsis & Miller and
colleagues), and a 2026 bioRxiv preprint. If the field/ephaptic strand becomes
load-bearing for any repo claim it needs its own verified pass.

## 3. Where the repo already stands

Nothing here is new work; it is what this repo has already measured, gathered
so the comparison in §4 is against evidence rather than impression.

**The wave's causal status (C-series).** All numbers from the linked docs.

| result | measurement |
| :--- | :--- |
| `W = f(S)` is a deterministic coarse-graining | wave adds **+0.001** decode accuracy given full spikes ([C0](c0_results.md)) |
| the wave is informative only for a *collective* code | **+0.110** at `\|S_obs\|=6/120`, falling to **+0.001** at 120; labeled-line ≈ 0 — "predictively epiphenomenal" ([C0](c0_results.md)) |
| `do(W)` is **fat-handed**, and how badly depends on the reader | realization band **0.24 σ** (collective) vs **33.1 σ** (labeled-line) at `t=0`; at n=30, **0.287** [0.280, 0.295] vs **26.82** [25.79, 27.91] ([C2](c2_results.md)) |
| `do(θ)` is single-valued | ambiguity **0.014 σ**; across-seed reproducibility SD 0.15 ([C3](c3_results.md)) |
| the causal role is (handle, outcome)-relative | effect matrix diagonal: `do(τ_gate)`→timing 1.00 / identity 0.00; `do(g_route)`→identity 1.00 / timing 0.06 ([C4](c4_results.md)) |
| the wave is the natural variable only where behaviour is collective | macro-sufficiency **1.03** vs **0.11**; at n=30, **1.027** [1.023, 1.032] vs **0.087** [0.068, 0.106] ([C4](c4_results.md)) |
| a *rotating* wave's chirality escapes fat-handedness only if read **topologically** | `do(χ)` band: center **6.2 σ**, tracked **1.0 σ**, global **2.6 σ** ([C5](c5_results.md)) |
| `do(θ_χ)` is well-posed for every reader, and the core is *necessary* | **0.0 σ** for all three readers; ablation collapses switching **0.85 → 0.52** while sparing single-rule routing (0.90 vs 0.89) ([C6](c6_results.md)) |
| chirality is a genuine **mediator** | on `θ_seed → χ → B`; "**Not epiphenomenal** … A rotating wave here does causal work. But its causal status is contingent" ([C7](c7_results.md)) |

> ⚠ Do **not** quote the "`do(θ)` is ~2400× less ambiguous than `do(W)`"
> comparison. [`core_review.md`](core_review.md) already flags 0.014 σ and
> 33.1 σ as **different quantities** — estimation-noise-over-effect-size vs a
> realization band — so that ratio over-dramatizes comparability.

**Slow-gates-fast, three independent ways.**

| mechanism | measurement |
| :--- | :--- |
| learned fast/slow τ hierarchy + phase–amplitude coupling | Tort MI **0.000 → 0.594**, mean fast amplitude **0.200 → 0.018**, `P_f=5`, `P_s=23` non-commensurate, n=20 ([timescale hierarchy](timescale_hierarchy_results.md)) |
| a slow reentrant loop gating fast routing | switching **0.89** vs **0.20** ablated; switch cost consolidates 0.57 → 0.92; single-rule spared 0.87 vs 0.86 ([E5](e5_results.md)) |
| necessity of a persistent rotating core | switching **0.85 → 0.52** on ablation ([C6](c6_results.md)) |

**Waves doing computational work.** E4 implements selective attention as biased
winner-take-all by refractory wave annihilation, accuracy **0.96** at bias
`|b|=4`, annihilation locus linear in bias (slope 0.50, intercept 19.8) with
**zero inhibitory nodes** — "refractoriness *is* the competition"
([E4](e4_results.md)).

## 4. The comparison

### 4.1 Convergence: the field's own best intervention is a `do(θ)`, not a `do(W)`

The strongest interventional result in the consciousness line does not
manipulate waves. It stimulates thalamus [9, 13] — upstream generating
structure — and the waves and the behaviour follow. In this repo's vocabulary
that is a `do(θ)`, and C3/C6 say exactly why it is the intervention that
*works*: `do(θ)` is single-valued (0.014 σ; 0.0 σ band for every chirality
reader) whereas `do(W)` on a constituted aggregate has an irreducible
realization band (up to 33.1 σ).

This is the note's main observation, and it is a *reframing*, not a refutation:
the repo's formalism predicts that a wave-based theory's decisive experiments
will have to intervene on generators rather than on waves — and the literature,
independently, did that.

### 4.2 Convergence: rotating waves

[8] reports PFC travelling waves are predominantly **rotating**, not planar,
with task-dependent direction. The repo's C5–C7 arc is about exactly that
object — a rotating wave, its topological charge, and whether its chirality
does causal work — and concludes it is a genuine mediator, **not**
epiphenomenal ([C7](c7_results.md)). Two independent routes to "the rotating
wave is the interesting variable, and its direction carries information."

### 4.3 Supported: the coarse claim that a slow rhythm gates a fast one

Three ways, above. The repo shows a substrate of this kind *can* do this.

### 4.4 Contradicted: the **per-cell** reading of the stencil

If the stencil is read as *"an individual cell may fire only during the
permitted phase of the slow wave"*, this repo has already refuted it — four
times, and these are the repo's own negatives, not imported doubts:

- phase gating of a cell's own input was **refuted**: selecting *when* a cell
  listens does no better than selecting *how strongly*;
- the mechanism works "**as a clock, not a filter**" — gate-as-filter fails at
  every window width ([lattice results](lattice_results.md), Result 3);
- the gating *benefit* was **retracted** at n=20 while the gating **cost**
  remains among the most certain numbers in the comparison (readout lost by
  **8.2 sd** at D=70, **12.4 sd** at D=50), concluding "nothing in this task
  family justifies the gate";
- phase is the **third of four** raw signals to fail as a label or driver — the
  arc's transferable law being that a raw signal *magnitude or geometry* cannot
  serve as a label on this substrate.

Any survivable version of the stencil here must therefore be **collective and
geometric**, not per-cell. That is a substantive constraint the repo can
contribute.

### 4.5 Untested — and this is the real opening

Miller's actual claim [3, 4] is that the **spatial pattern** of the slow rhythm
is what constrains where fast activity appears. The repo has never tested that.
Its cross-frequency coupling is **one scalar applied to every fast node**
(`net.theta[fid] = TH_HI - DEPTH * np.clip(s_prev / S_REF, 0, 1)` —
[`timescale_cfc.py:118`](../experiments/timescale_cfc.py)) on a 120-node pool with
**no geometry at all**. So the repo has demonstrated *amplitude* gating and
called it cross-frequency coupling; the stencil claim is about *pattern*
gating, and the gap is large.

Two further asymmetries worth recording:

- **Sign.** The repo's coupling is *permissive* — slow activity **lowers** fast
  thresholds, so fast fires where slow is active. Miller's stencil is
  *inhibitory* — alpha/beta patterns are **inversely** correlated with sensory
  spiking [4]. Same mechanism class, opposite sign. Both are one parameter
  apart in the existing code, and nothing tells us they behave symmetrically.
- **Direction.** The theory privileges slow-constrains-fast. This repo's
  [lattice notes](lattice_timescale_notes.md) record the opposite tendency:
  *the fastest pacemaker entrains the whole medium* (near-`P_f` 0.51, near-`P_s`
  **0.00** at θ=1, Moran's I 0.32), and confining fast pacemakers to one half
  changed nothing. So slow-dominates is not a property this substrate hands you
  for free.

### 4.6 Equivocation risk: "analog"

**This substrate is not analog.** `Network.phi` is `np.int64` with states
`0 / 1..act / act+1..tau` ([`ghca_net.py:130`](../ghca_net.py)), and the entire
theory arc in this repo is *exhaustive enumeration of finite state sets*. Any
claim of the form "our substrate supports analog computation" would be a
category error. What the substrate can speak to is **spatially extended,
continuous-valued-in-the-limit wave organisation** — a different and weaker
thing. This note deliberately declines the word.

A second equivocation to guard: the repo's `W` is a scalar coarse-graining
(Kuramoto coherence, active fraction). A cortical "beta wave" is a
band-limited, phase-resolved, spatially extended field. Whether `do(W)`
fat-handedness transfers to that object is **argued, not shown** — the argument
being that any deterministic many-to-one aggregate inherits a realization band,
which is generic; but the specific 33.1 σ number is a property of *this* readout
and must never be quoted as a fact about cortex.

## 5. Self-audit (the independent pass did not run)

Commissioned critic passes on overreach, technical correctness and completeness
all failed on a spend limit. This section is the author auditing their own
note, and by this repo's convention that is **not sufficient**. Treat §4 as
un-reviewed until someone else attacks it.

Problems I can already see:

1. **The §4.1 convergence is the most seductive claim here and the easiest to
   oversell.** "The field's best intervention is a `do(θ)`" is an *analogy of
   structure*. Thalamic stimulation is not literally setting a per-node
   timescale; it is a gross electrical perturbation of a nucleus with diffuse
   cortical projections. The honest statement is that both are interventions on
   *generators rather than aggregates*, which is a real structural parallel and
   nothing more.
2. **§4.4's "contradicted" is scoped to this substrate's task family**, and the
   retraction it leans on was itself an n=20 correction. It should not be read
   as evidence against the per-cell stencil *in cortex*.
3. **Citation integrity is the largest single risk** and is unresolved (banner).
4. **No number in §3 is new.** This note contributes framing and a gap
   analysis, not results. If it reads as if the repo has tested spatial
   computing, it is misleading and should be rewritten.
5. **Selection effects in the literature pass.** Search-result summaries
   preferentially surface a lab's own framing. No systematic search for
   *critiques* of spatial computing survived (that agent completed, but its
   output was not incorporated into §4 before the critic passes died), so §4 has
   no adversarial-literature counterweight.

## 6. What this substrate could actually test

Prioritised by falsification power, not by how illustrative they are. None of
these is run.

1. **PAC-matched shadow foil (cheap, decisive).** Build a *reverse*-causation
   model in which fast refractory load *generates* the slow envelope, matched to
   the real model on Tort MI, then intervene. If a shadow reproduces
   PAC ≈ 0.594, then **PAC cannot license the stencil claim at all**.
   Structurally this is [C1](c1_results.md)'s `confounded` row — observational
   association 0.642 with a `do(W)` effect of 0.001, the textbook trap the
   certificate already gets right. ~100 lines on
   [`timescale_cfc.py`](../experiments/timescale_cfc.py) plus the existing
   epiphenomenality certificates.
2. **Give the stencil geometry.** Port the coupling from one global scalar onto
   a lattice with real geometry, so the slow *pattern* rather than the slow
   *level* sets fast excitability. Then ask the only question that matters: does
   a *structured* mask beat a *level-matched unstructured* mask? If not, the
   stencil is sparsification, not computation.
3. **Report cost before benefit.** Given §4.4, every stencil experiment must
   report what the mask *costs* (synchrony, reach, accuracy) before what it
   buys. The repo's history is that the cost is the certain number and the
   benefit is the fragile one.
4. **Flip the sign.** Run inhibitory-stencil and permissive-gate variants
   head-to-head (one parameter apart) and check whether they are symmetric.
5. **Test the privileged direction.** Given that the fastest pacemaker entrains
   this medium, find what — if anything — makes slow-constrains-fast the stable
   arrangement. A negative here would be informative about the theory's
   substrate requirements.

## 7. Honest scope

- **A cellular automaton cannot confirm or refute a claim about cortex.** This
  note compares mechanisms and locates gaps. Nothing in it is evidence about
  brains.
- **Citations unverified** (banner) and the **independent review did not run**
  (§5).
- The prompting article was **never read** — `picower.mit.edu` was blocked.
- The strong "cognition and consciousness are analog computation" claim is a
  **lecture/press claim**, not a peer-reviewed one, and is treated as such.
- The repo's substrate is **finite-state**, so "analog" does not apply to it
  (§4.6).
- `θ = 1` and the usual substrate caveats apply to every repo number cited.
