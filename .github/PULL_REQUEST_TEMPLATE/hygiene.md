<!--
Reproducibility / infrastructure / fix PR.
The seed-everything house rule is the bar; see AGENTS.md and core_review.md.
-->

## What & why
<!-- The fix / infra change and the problem it addresses. -->

## Verification
- [ ] Affected experiments re-run from scratch
- [ ] Result: bit-identical to committed data — OR numbers updated with the diff explained below
- [ ] No new global/unseeded RNG introduced (`default_rng(seed)` threaded)

## Numbers moved?
<!-- If any committed number changed, show old → new and why it is now correct. -->

## Docs
- [ ] `docs/next_steps.md` §5b (or the affected doc) updated if scope/status changed
