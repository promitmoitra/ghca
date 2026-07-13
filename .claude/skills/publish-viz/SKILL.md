---
name: publish-viz
description: Publish a doc / result / animation to the project's GitHub Pages site. Cherry-picks the chosen files from a source ref onto the dedicated `deploy-viz-page` branch, wires them into the MkDocs nav, validates the build, and pushes (which auto-redeploys). Use when asked to publish, deploy, ship, or add a result to the public docs site.
---

# publish-viz

Publish curated content to the project's public docs site
(https://promitmoitra.github.io/ghca/).

This skill is **agent-agnostic**: it is plain `git` + `mkdocs`, with an optional
helper script. Any agent — or a human at a terminal — can follow it. It uses no
proprietary tools; pushing to the deploy branch triggers the existing GitHub
Actions workflow that builds and deploys the site.

The skill **aids** a manual write path — it does not replace it. Everything below
can be done by hand (that is the point of keeping it plain git + mkdocs). Two ways
to write to the page, both first-class:
- **Direct push** (fast, for a trusted publish) — the procedure below.
- **Reviewable PR** into `deploy-viz-page` — open a PR with the
  [`publish` template](../../../.github/PULL_REQUEST_TEMPLATE/publish.md) so a
  maintainer sees what goes public before the merge deploys it. The steps are the
  same; you push a publish branch and open the PR instead of pushing the deploy
  branch directly.

## The branch model (read first)

- **`main`** is the research trunk — full history, work in progress.
- **`deploy-viz-page`** is a **dedicated publication branch**. GitHub Pages builds
  only from it, and it owns the publication layer: `mkdocs.yml`, the deploy
  workflow (`.github/workflows/deploy-docs.yml`), `docs/index.md`, and the rendered
  `docs/figures/*.gif`.
- Flow is **one-way into** `deploy-viz-page`: mature results are cherry-picked from
  `main` (or a dev branch) for public display. It may diverge from `main`; that is
  intended. Never open a `deploy-viz-page -> main` PR — that is the wrong direction.

## Prerequisites

```
python3 -m pip install mkdocs-material pillow
```
`pillow` is only needed if you must (re)render GIFs from the `experiments/*_animation.py`
scripts; already-rendered GIFs just need to be copied.

## Procedure

Inputs you need: a **source ref** (e.g. `main` or a commit SHA) and the **paths**
to publish (a result doc, its figures, any new animation GIFs).

1. **Get onto the deploy branch with the new files.** Either run the helper:
   ```
   .claude/skills/publish-viz/publish.sh <source-ref> <path> [<path> ...]
   # e.g. .claude/skills/publish-viz/publish.sh main docs/e8_results.md docs/figures/e8_tone.png
   ```
   …or do it by hand:
   ```
   git fetch origin deploy-viz-page <source-ref>
   git checkout -B deploy-viz-page origin/deploy-viz-page
   git checkout <source-ref> -- <path> [<path> ...]
   ```

2. **Wire it into the site.** Edit `mkdocs.yml` to add a `nav:` entry for any new
   page, and add an image embed (`![alt](figures/xxx.gif)`) in the relevant doc if
   you are showcasing an animation. Keep nav labels consistent with existing ones.

3. **Validate the build (do not skip).**
   ```
   mkdocs build --strict
   ```
   `--strict` fails on broken links or files missing from nav. Fix any warning
   before pushing. (Links to files outside `docs/` — e.g. `../ghca_net_viz.py` —
   will fail; make those plain text or point them at a GitHub URL.)

4. **Commit and push.** Pushing to `deploy-viz-page` auto-triggers the deploy.
   ```
   git add -A
   git commit -m "publish: <what> to the docs site"
   git push origin deploy-viz-page
   ```

5. **Confirm the deploy is green.** Check the **Actions** tab for the
   "Deploy docs (MkDocs -> GitHub Pages)" run on `deploy-viz-page`; both the
   `build` and `deploy` jobs should succeed. The site updates at the URL above
   (allow a minute for CDN propagation).

## First-time setup (only if the site was never deployed)

In the repo's **Settings -> Pages**, set **Source: GitHub Actions**. Then, under
**Settings -> Environments -> github-pages -> Deployment branches**, allow
`deploy-viz-page` (or "All branches"), otherwise the deploy job fails immediately
with an environment/branch-policy error.

## Guardrails

- Only publish content you were asked to publish. The public site is outward-facing.
- Keep binaries lean; prefer GIFs under ~1 MB (tune frame count / `stride` / `dpi`
  in the animation scripts).
- Do not merge `deploy-viz-page` into `main`.
