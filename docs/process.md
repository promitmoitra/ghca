# Process — Planning and Review

*How this project runs its two recurring meta-passes. Read this before doing a
review pass or a planning pass; the point is to keep them **decoupled in
process** while **linked by a one-directional hand-off**.*

*This process is **agent-agnostic** (the [`publish-viz`](https://github.com/promitmoitra/ghca/blob/main/.claude/skills/publish-viz/SKILL.md)
convention): any agent — or a human at a terminal — can follow it. It relies on
no proprietary tools.*

## The two passes

**Review** is adversarial and backward-looking. Its job is to re-run the
experiments, check that every headline number reproduces from committed data,
hunt for the failure modes we know about (bimodal-mean headlines dressed as
clean means — the E3 lesson; unseeded RNG — the `perturb_tau` lesson; invented
or mischaracterised citations; capability credited to the *substrate* that
actually lives in a hand-built *readout*), and draw the substrate-vs-analysis
line. A review's virtue is skepticism. It should read as if written by someone
with no stake in the roadmap.

**Planning** is generative and forward-looking. Its job is to lay out candidate
directions as a menu — each with what it is, why it matters, effort/risk, and
how it connects to existing work — and to recommend what to build next. A plan's
virtue is honest ambition anchored to what the reviews actually found.

## Why they are decoupled

1. **Review objectivity is this project's core value.** The whole programme is
   about *not overclaiming*. A review written in the same breath as the roadmap
   has a built-in incentive to go easy — "this caveat is fine, here's the
   follow-up that fixes it." Keeping the passes separate protects the audit from
   the plan's optimism. The independent audit ([`core_review.md`](core_review.md))
   carries more weight precisely because it was not also planning next steps.
2. **Different cadence.** Review fires *after* a batch of experiments lands
   ("audit E7/C5–C7/E8 for overreach"). Planning fires when we need to *choose*
   the next build. Forcing them into one pass means reviewing too early or
   planning too late.
3. **Different mindset.** Skeptical/adversarial vs generative/optimistic.
   Bundling blunts both — you soften the audit to protect the roadmap, or the
   roadmap inherits the audit's pessimism.

## The hand-off (the one link that stays)

Decoupled is not disconnected. The passes are joined by a **one-directional
dependency: review → plan.**

- Review runs first and produces a **standalone artifact** that stands on its
  own as a record (see the review docs below).
- Planning **consumes** the review's findings as its input constraints and
  cites them. In [`next_steps.md`](next_steps.md) this is literal: its spine is
  "the three tensions that should steer the roadmap," and those tensions *are*
  the audit findings; every track is scored by which tension it retires.

Planning never edits a review's findings to make a track look better. If a plan
disagrees with a review, that is a signal to run a fresh review, not to soften
the old one.

## Artifacts and where they live

| Pass | Artifact | Scope |
|------|----------|-------|
| Review | [`core_review.md`](core_review.md) | independent audit, E0–E6 / C0–C4 |
| Review | [`extensions_review.md`](extensions_review.md) | self-audit, E7 / C5–C7 / E8.x |
| Planning | [`next_steps.md`](next_steps.md) | roadmap / menu of directions |

**Both audits land on `main` alongside the work they cover.** A review kept only
on its own branch is decoupling gone too far — the record should sit with the
code it audits, not drift off on a side branch.

## Rules for either pass (agent or human)

- **Doing a review?** Re-run from scratch; compare regenerated output against
  every quantitative claim in the docs. Pull per-seed arrays and check for
  bimodality before trusting a mean. Verify every load-bearing citation
  resolves and says what the doc claims. Read code for circularity/leakage and
  for the substrate/analysis boundary. Keep the tone skeptical; do not propose
  the roadmap in the same document. State what you did **not** verify.
- **Doing a plan?** Start from the current review artifacts. Frame directions as
  a menu, not a commitment. Score each against the open tensions the reviews
  surfaced. Recommend, don't just survey. Do not relitigate or soften a review
  finding to justify a track.
- **Either pass:** seed everything (`default_rng(seed)`, never the global RNG);
  report per-seed spreads, not just means; state the substrate/analysis boundary
  explicitly; keep a caveats section adjacent to every headline.
