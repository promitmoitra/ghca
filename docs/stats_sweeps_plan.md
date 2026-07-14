# 3a Plan — Statistics & Operating-Point Sweeps

*Scope for Track 3a of [`next_steps.md`](next_steps.md). Turns every headline
from "shown at n=3–5 at one hand-chosen operating point" into "mean ± 95% CI,
robust across a neighborhood." Directly retires the **narrow-evidence** tension
(tension 2) the audits raised — **except substrate generality (other topologies),
which is Track 3b and explicitly out of scope here.** No new mechanisms; the only
expected casualty is the precise magnitude of a few fragile headlines (E3, the E8
tails), which is the point of measuring them properly.*

This is an **experiment/hygiene** track, not a review or planning pass — it runs
under the seed-everything house rules ([`AGENTS.md`](AGENTS.md)) and reports
per-seed spreads and caveats per [`process.md`](process.md).

## Current state (why this track exists)

Seed counts and the single inherited operating point, as they stand on `main`:

| Experiment | Seeds now | Headline that needs a CI |
|---|:--:|---|
| E1 conditioning | 6 | A=0.91 vs B≈chance |
| E2 memory | 5 | retention vs delay (Line-B holds to D=200) |
| E3 timed / factored | 5 (resonance map **2**) | double dissociation; **composition 0.560** |
| E4 attention | psychometric bias sweep (confirm n) | accuracy 0.96 at modest bias |
| E5 executive | 5 | switching 0.89 vs ablated 0.20 |
| E6 horde | 5 | mem R²=0.62 / att 0.84 / exe 0.98 |
| E7 spiral | 5 | switching 0.86; persistence |
| E8 predictive/nested/conditional | **3** (some ridge deterministic) | window 0→7; nested 1.00 vs 0.55; L>K threshold |
| E9 conjunction | 5 | selectivity 0→1; routing 0.84 |
| C5 / C6 / C7 | 3 / 5 / 3 | the σ bands; necessity 0.85→0.52 |
| C0–C4 | largely deterministic / ≤3 | the σ/macro-sufficiency figures |

**The operating point** is essentially one point everything inherits: the E0
spiral band **r=2, a=6, p=8 (τ=14), θ≈4** (θ≥5 → death, θ<4 → no cores). E7 sits
at `L=48, θ=4.0`; E8 fixes `τ=14`; E3's flagged composition sets target latency=16
*because* it forces gate τ∈[17,19], a **lucky resonance** (τ=13–15 → chance in the
code's own map). Partial single-axis sweeps already exist (E0 over θ and τ; E3's
τ resonance map at n=2; E8 over τ / τ_ctx / L at n=3) — they need seeds, CIs, and a
neighborhood, not new scaffolding.

## Seed target

**n=50 for every headline dissociation** (publication-grade CIs). With the
observed tight SDs (≈0.02–0.05) this puts the bootstrap 95% CI half-width far
below the effect gaps; for the bimodal E3 composition, n=50 also gives a usable
Wilson interval on the joint-success *rate*, not just the mixture mean. Sweep
grid-points (θ/τ neighborhood, E8 sweeps) run at n=50 per cell as well, which is
the bulk of the compute below.

## Deliverables

1. **`ghca_stats.py`** — a small reusable harness: bootstrap 95% CI, an
   effect-size gap (intact−ablated, standardised), an automated **bimodality
   flag** (the E3 discipline, so no future headline hides a two-mode mean), and a
   standard strip+CI plot helper. Seeded like everything else.
2. **Per-experiment**: raise n to target, save per-seed arrays, plot mean ± CI +
   the per-seed strip.
3. **`docs/stats_sweeps_results.md`** — one master table (headline · n · mean ±
   95% CI · effect size · holds-across-sweep? · bimodal?), plus sweep panels for
   the two priority cases below.
4. **Targeted doc/README edits** — only where a number moves materially; update
   `next_steps.md` 3a to done with an honest ledger of what softened.

## Priority sweep matrix (highest value first)

1. **E3 composition across target latency / gate τ.** The audit's marquee
   fragility. Sweep target latency across the full resonance map at n=30; report
   joint-success *rate* as a function of τ, not one lucky point. Answers "magnitude
   not established at n=5."
2. **Substrate θ/τ neighborhood.** θ∈{3.5, 4, 4.5}, τ∈{10, 12, 14, 16, 18} at
   r=2, a=6; re-run the E-series headline dissociations across the grid. Retires
   "hand-chosen operating point" for the whole series at once.
3. **Seed scale-up + CI** on E1/E2/E5/E6/E7/E9/C5–C7 (mechanical, low-risk).
4. **E8 series** (thinnest at n=3): n→15 on the τ / τ_ctx / L sweeps, with CIs.

## Phasing (one reviewable PR each, using the new templates)

- **P1** — `ghca_stats.py` + seed/CI scale-up on the E-series dissociations
  (E1/E2/E3/E5/E6/E7/E9). Mechanical, high-confidence; establishes the harness.
- **P2** — the two priority operating-point sweeps (E3 latency; substrate θ/τ).
  The parts most likely to move headlines.
- **P3** — C-series + E8 seed/CI and sweeps.
- **P4** — `stats_sweeps_results.md` + targeted doc/README updates; close 3a.

## Effort & risk

**Effort:** real compute — n=50 is ~10× the current seed counts, and the sweep
grids multiply that per cell; expect several hours of parallelisable runtime
(E8 × 50 seeds × sweep grid is the heaviest single piece), run in the background
in phases. **Risk:** low
to the programme (findings are already hedged); moderate to specific headline
*magnitudes* — E3's composition and the E8 tails are the likeliest to soften. A
softened-but-honest number is the intended outcome, not a regression.

## Open confirmations (resolve in P1)

- E4's exact seed handling (its psychometric already sweeps bias) and whether the
  C0–C4 deterministic panels need any seeded replicates at all.
- Whether the substrate θ/τ grid re-runs *all* E-series headlines or only the
  three most operating-point-sensitive (E5, E7, E3).
