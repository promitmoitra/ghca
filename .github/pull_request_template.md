<!--
Default template — experiment / results PR (the common case).
Other templates (append ?template=NAME.md to the PR-compare URL, or ?expand=1):
  ?template=review.md    — an audit / review pass
  ?template=planning.md  — a roadmap / next_steps update
  ?template=hygiene.md   — reproducibility / infra / fixes
  ?template=publish.md   — publish docs to the GitHub Pages site (deploy-viz-page)
See AGENTS.md (house rules) and docs/process.md (review vs planning).
-->

## What & why
<!-- Which experiment/series (E#/C#) and the question it answers. -->

## Substrate / analysis boundary
<!-- State it explicitly: what the *dynamics* do vs what a *readout/feature* does. -->

## Results
- Seeds (n):
- Per-seed spread (not just the mean — call out any bimodality; the E3 lesson):
- Operating point(s):
- Headline numbers:

## Reproducibility
- [ ] All RNG seeded — `default_rng(seed)` threaded explicitly, no global NumPy RNG (the `perturb_tau` lesson)
- [ ] Headline numbers regenerate from committed `result/**.npz`
- [ ] Figures written to `docs/figures/`

## Caveats
<!-- Adjacent to the headline. What this does NOT show; the scope limits. -->

## Docs
- [ ] `docs/<id>_results.md` added/updated (caveats section present)
- [ ] README progress table updated (do not round hedged results up)
