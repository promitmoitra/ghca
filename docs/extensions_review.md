# Extensions Self-Audit — E7, C5–C7, E8.x

*An integrity/overreach audit of the extension work (spiral option E7, the
spiral causal split C5–C7, and predictive dynamics E8.0–E8.7), applying the same
scrutiny the independent review ([`core_review.md`](core_review.md))
applied to E0–E6 / C0–C4 — with particular attention to the failure mode it caught
in E3 (a confident headline over a bimodal, lucky-seed mean).*

Audit date: 2026-07-11. Branch: `claude/e-series-experiments-review-dew74b`.
Method: re-ran all nine extension experiments from scratch and compared stdout /
`.npz` against every quantitative claim in the docs; pulled per-seed arrays to
check for hidden bimodality; verified citations; read the code for
circularity/leakage and for the boundary between what the *substrate* does and
what the *analysis* does.

## Bottom line

**All extension numbers reproduce bit-identically**, and none rest on the E3
failure mode: the headline per-seed spreads are **tight and unimodal**, not
bimodal. There is no fabrication and no lucky-seed rounding.

The honest residual risk is **altitude of a different kind than E3's**: several
extensions credit the *substrate* with a capability whose learning actually
lives in a **hand-built readout or feature** (a computed topological readout, a
reused routing learner, an engineered conjunction feature). The claims are stated
as "what this representation affords," which is accurate, but a reader should keep
the substrate/analysis boundary in view. Details below.

## Reproducibility — all exact

Re-running on this branch (whose `ghca_learn.py` still has the unseeded-RNG bug the
review found) reproduces every extension bit-identically — **because none of the
extensions use `perturb_tau` / Line B**; they are Line-A-only (E5, E7, C5–C7) or
ridge readouts on a seeded substrate (E6, E8). So the reproducibility fix (PR #5)
does not change any extension number; it should still propagate to this branch when
PR #5 merges.

| experiment | reproduces? | headline (re-run == doc) |
|---|:--:|---|
| E7 Phase A (spiral mechanism) | yes | persistence +0.89/−1.11/0; readout charge 1.00, probe 0.90 |
| E7 Phase B (rotation = rule) | yes | switching 0.86 vs 0.49; single 0.90/0.89; rotation→rule 1.00 |
| C5 (`do(χ)` fat-handed) | yes | center 6.2σ, tracked 1.0σ, global 2.6σ; routing 0.55/0.78 |
| C6 (`do(θ_χ)`, necessity) | yes | 0σ all readers; switching 0.85→0.52, single 0.90/0.89 |
| C7 (outcome-relativity) | yes | do(χ)[1.00, 0.11], do(g_route)[0.57, 1.00]; screening follows injected χ |
| E8.0–8.4 (prediction) | yes | 1.00/0.90/0.55/0.20/0.41; window 0→7; surprise 0.16/0.64 |
| E8.5 (nested) | yes | nested 1.00 vs 0.55; within-context 1.00/1.00 |
| E8.6 (order-preserving) | yes | grid recall 1.00 vs trace ≈chance; depth threshold |
| E8.7 (conditional) | yes | grid+conj 0.96, others ≈chance; L>K threshold |

## The E3 check — are the headlines bimodal? (No)

The review's central catch was that E3's "0.20→0.77" was a mean over a bimodal
per-seed set (joint composition on 1/5 seeds). I pulled the per-seed arrays for the
extension headlines:

- **E7 switching** intact `[0.89, 0.84, 0.86, 0.85, 0.87]` (SD ≈ 0.02); ablated
  `[0.51, 0.44, 0.57, 0.49, 0.42]`. Tight, unimodal, well-separated.
- **E7 single-rule** intact/ablated `[0.92,0.85,0.89,0.90,0.93]` / `[0.93,0.73,0.93,0.91,0.92]`.
- **E5 switching** intact `[0.96,0.82,0.92,0.91,0.85]`; ablated `[0.21,0.19,0.25,0.15,0.19]`.
- **C6 necessity** switching 0.85 ± 0.03 vs 0.52 ± 0.01 (SEM over 5 seeds).

So the extension dissociations are robust means, not lucky-seed artifacts — the
specific thing that sank E3's composition headline is absent here.

## Honest caveats — the substrate/analysis boundary

These are the real "altitude" risks. All are disclosed in the individual results
docs; collected here so they are not lost under confident headlines.

1. **Computed readouts, not emergent circuits.** E7's rotation-direction readout
   (local winding / net topological charge) is a *computed* topological operation,
   not a learned or spiking sub-circuit; likewise E8.5/E8.7's fast×slow
   **conjunction** feature and the trace / delay-line **reservoirs** are hand-built
   representations fed to a linear readout. What the substrate supplies is the
   *dynamics* (a persistent spiral core, τ-controlled traces, a shift-register
   delay line); what the analysis supplies is the readout and, in the conditional
   tasks, the interaction feature. The claims are about *what these representations
   afford a linear readout* — true as stated — but "the substrate does X" should be
   read as "the substrate's dynamics make X linearly decodable," not "the substrate
   learned X end-to-end."
2. **E7's rotation→rule decode (1.00) is partly definitional.** The rule *is* the
   chirality and the readout recovers the chirality, so 1.00 is near-tautological —
   it shows rotation direction is cleanly *readable* (qualified by C5's
   fat-handedness result), not that a task-classifier *discovered* rotation as the
   relevant variable (the fMRI framing it alludes to). Stated as an analogy in the
   doc; flagging it as engineered, not emergent.
3. **E7 reuses E5's routing learner.** Phase B's H→M routing and its learning are
   E5's, reused verbatim; E7's novel content is the spiral option + the readout, not
   a new learner. The E7 doc says so; worth keeping in view when reading "switching
   0.86."
4. **Chosen operating points.** Like E3, the extensions evaluate at hand-chosen
   operating points (E0's narrow spiral band θ≈4; E8's τ/L/K/T_ctx). *Unlike* E3,
   the results are robust across seeds at those points (see the E3 check above) — so
   this is a scope limitation ("shown at this operating point"), not a lucky-seed
   inflation.
5. **Small n.** 3–5 seeds throughout; several E8 panels are deterministic single
   runs (ridge) or 3-seed averages. Adequate given the tight spreads, but small.
6. **C7's matrix is diagonal-*dominant*, not clean-diagonal** — `do(g_route)` moves
   `O_rule` by 0.57 (erasing routing also destroys the parity), vs C4's fully
   independent latency/channel outcomes. Disclosed in the C7 doc.
7. **E8 prediction is a self-supervised readout, not reward-driven Line-A
   plasticity**, and "surprise" is the readout's error, not a substrate signal. The
   reservoir-computing framing (fixed substrate + trained linear readout) is honest
   but means the *learning* credited to E8 is in the readout weights.

## Citations

- **Xu/Gong et al., *Nat. Hum. Behav.* 2023** (interacting spiral waves; rotation
  direction task-relevant) — real, accurately used (E7 motivation).
- **Ye/Steinmetz et al., "Brain-wide topographic coordination of traveling spiral
  waves"** (bioRxiv 2023.12.07.570517; *Science* 2026) — real (abstract verified);
  the "axonal architecture shapes the spiral / topographic mirroring" claim is
  accurately characterized (E7, spiral spec).
- **Muller / Benigno et al., *Nat. Commun.* 2023** ("Waves traveling over a map of
  visual space can ignite short-term predictions") — real, accurately used (E8
  wave-forward-model motivation).
- **Jalaldoust & Zabeh, arXiv:2511.06602** — real (confirmed by the prior review);
  underpins C5–C7.
- **Rao & Ballard 1999; Friston 2005; Sutton (TD, Horde)** — standard, used as
  framing/foil, not as sources of numbers.
- **⚠ Chait/Bianco, "Neural mechanisms of time-forward predictions for naturalistic
  auditory tone sequences" (cited for E8)** — **the exact 2026 paper could not be
  independently verified** (the DOI did not resolve in search). Its *predecessor*,
  "Neural Integration of Stimulus History Underlies Prediction for Naturalistically
  Evolving Sequences" (*J. Neurosci.*), is real and its mechanism (history
  integration over ~7 tones → next-pitch prediction; anticipatory sharpening with
  predictability) is what E8 actually leans on and characterizes accurately. **The
  E8 docs should hedge the specific 2026 title as "reported" / cite the verified
  predecessor** until the 2026 reference is confirmed. This is the one
  citation-integrity item in the extensions.

## Verdict & recommendation

The extensions are **reproducible and honestly reported**, and they specifically
avoid the E3 failure mode (no bimodal-mean headlines). The genuine limitation is
that they are **reservoir/readout demonstrations**: the substrate provides the
dynamics, hand-built readouts and features do the decoding, and several headline
capabilities are *afforded* rather than *learned end-to-end*. That is disclosed
per-doc and is legitimate for an illustrative study, but the framing should keep
the substrate/analysis line visible.

Concrete follow-ups: (1) hedge/repoint the **Chait/Bianco 2026 citation** in the E8
docs; (2) when PR #5 merges, propagate the `perturb_tau` seeding fix onto this
branch (no extension numbers change); (3) the deepest open item is to move the
computed readouts/features (winding, conjunction) *into* the substrate as learned
or spiking mechanisms — the honest frontier these experiments point at but do not
yet cross.

*(The independent core-series audit referenced above now lives on `main` as
[`core_review.md`](core_review.md) — renamed from its original branch filename.)*

**Update (2026-07-11): follow-up (3) is now partly crossed for the conjunction
case.** [E9](e9_results.md) replaces the hand-built `(stimulus × rule)` conjunction
feature (E5, and the outer-product features of E8.5/E8.7) with a basis that
**self-organises on the substrate** via a reward-free competitive Hebbian rule —
conjunction selectivity emerges 0.00 → 1.00 and reward routing on the emergent
basis (0.84) approaches the wired upper bound (0.92) while a no-self-organisation
control fails at chance. So the conjunction is now *learned*, not afforded. The
residual computed readouts (E7's topological winding; the E8 delay-line/trace
reservoirs) and the imposed competition (k-WTA) in E9 itself remain on the far side
of the frontier.
