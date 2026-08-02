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

## Deploy from an isolated worktree (project convention)

**Every publish runs in a dedicated `git worktree`, never by switching the branch
of your main checkout.** The old in-place `git checkout -B deploy-viz-page` flow
moved your working tree onto the deploy branch — stranding in-progress edits and
leaving you on the wrong branch after the push. Instead, the deploy branch is
checked out at a **gitignored** path, `.publish-worktree/`, and removed when you
are done; `main` never moves.

The helper (`publish.sh`) sets this up for you. By hand it is:

```
git fetch origin deploy-viz-page
git worktree add -B deploy-viz-page .publish-worktree origin/deploy-viz-page
cd .publish-worktree          # do ALL publish steps in here
# … stage, wire nav, build --strict, commit, push …
cd - && git worktree remove --force .publish-worktree   # clean up when green
```

`.publish-worktree/` is in `.gitignore`; it is a local staging area only and is
never committed to `main`.

## Prerequisites

Any one MkDocs runner:

```
uvx --with mkdocs-material mkdocs build --strict   # ephemeral, no install (preferred)
# or install it:
python3 -m pip install mkdocs-material pillow
```
`pillow` is only needed if you must (re)render GIFs from the `experiments/*_animation.py`
scripts; already-rendered GIFs just need to be copied. `publish.sh` auto-selects an
installed `mkdocs` if present, otherwise falls back to `uvx --with mkdocs-material`.

## Procedure

Inputs you need: a **source ref** (e.g. `main` or a commit SHA) and the **paths**
to publish (a result doc, its figures, any new animation GIFs).

1. **Stage the new files into the deploy-branch worktree.** Run the helper (it
   creates/refreshes `.publish-worktree/`, brings the paths in, and runs a first
   build check):
   ```
   .claude/skills/publish-viz/publish.sh <source-ref> <path> [<path> ...]
   # e.g. .claude/skills/publish-viz/publish.sh main docs/e8_results.md docs/figures/e8_tone.png
   ```
   …or do it by hand (see "Deploy from an isolated worktree" above):
   ```
   git fetch origin deploy-viz-page <source-ref>
   git worktree add -B deploy-viz-page .publish-worktree origin/deploy-viz-page
   git -C .publish-worktree checkout <source-ref> -- <path> [<path> ...]
   ```
   **All remaining steps run inside `.publish-worktree/`** (`cd .publish-worktree`).

2. **Wire it into the site.** Edit `mkdocs.yml` to add a `nav:` entry for any new
   page, and add an image embed (`![alt](figures/xxx.gif)`) in the relevant doc if
   you are showcasing an animation. Keep nav labels consistent with existing ones.
   **Updating a doc that already exists on the deploy branch?** Prefer a *surgical
   edit* over overwriting it from `<source-ref>` — the deploy branch owns
   publication-layer content (e.g. GIF embeds) that the source ref may not have, and
   a blind overwrite silently drops it.

3. **Validate the build (do not skip).**
   ```
   mkdocs build --strict           # or: uvx --with mkdocs-material mkdocs build --strict
   ```
   `--strict` fails on broken links or files missing from nav. Fix any warning
   before pushing. (Links to files outside `docs/` — e.g. `../ghca_net_viz.py` —
   will fail; make those plain text or point them at a GitHub URL.)

4. **Commit and push from the worktree.** Pushing to `deploy-viz-page` auto-triggers
   the deploy.
   ```
   git add -A
   git commit -m "publish: <what> to the docs site"
   git push origin deploy-viz-page
   ```

5. **Confirm the deploy is green.** Check the **Actions** tab for the
   "Deploy docs (MkDocs -> GitHub Pages)" run on `deploy-viz-page`; both the
   `build` and `deploy` jobs should succeed. The site updates at the URL above
   (allow a minute for CDN propagation).

6. **Clean up the worktree.** Use `--force`: after the push nothing in the
   worktree is at risk, and it sidesteps a spurious "working trees containing
   submodules" refusal that a transient nested `.git` can trigger.
   ```
   cd .. && git worktree remove --force .publish-worktree
   ```

## First-time setup (only if the site was never deployed)

In the repo's **Settings -> Pages**, set **Source: GitHub Actions**. Then, under
**Settings -> Environments -> github-pages -> Deployment branches**, allow
`deploy-viz-page` (or "All branches"), otherwise the deploy job fails immediately
with an environment/branch-policy error.

## Scrollytelling page conventions (`docs/scroll/`)

The `docs/scroll/index.html` walkthrough is **hand-authored HTML copied verbatim**
by MkDocs — it is *not* a Markdown page and is *not* run through MkDocs extensions.
Two consequences bite if forgotten (both have caused real regressions):

- **No Markdown, no MathJax on this page.** `$...$` / `$$...$$` LaTeX renders as
  literal dollar-signs and backslashes here (arithmatex + MathJax only apply to the
  MkDocs *markdown* notebook). Use plain Unicode — `τ θ ρ ± ≥ ≤ → × ≈ ⁻³` — and
  `<code>` for variable-ish terms, matching the page's existing asides. Keep it
  self-contained: inline CSS/JS, no external CDNs.
- **Scrub sprites are clean pixel-art — never bake text into a frame.** Frames are
  `ghca_net_viz.state_colors` output, nearest-upscaled. *All* narration (titles,
  status text, per-node labels, badges) lives in the HTML `.cap`/`.hint` beside the
  `<canvas>`, never burned into the PNG by matplotlib — baked-in labels are tiny,
  and on side-by-side sheets they overlap and garble. A frame needing a legend gets
  it in HTML, not in the image.

Other scroll-page rules:

- **Render scripts live on `main`** (`experiments/scroll_sprites*.py`,
  `experiments/e6_animation.py`); the sprite/GIF **outputs** live on
  `deploy-viz-page`. Always regenerate outputs from the *committed* script — never
  hand-edit a sprite or publish one from an uncommitted variant, or `main` can no
  longer reproduce what is live.
- **When a regenerated sprite's dimensions change, update the `<canvas>`
  attributes** — `data-frames` / `data-cols` / `data-fw` / `data-fh` (and
  `width`/`height`). The scroll scrubber reads its geometry from those; a stale
  value scrubs the wrong region.
- **Lay a substrate out by its structure**, not by convenience: e.g. render a
  layered graph sensory → hidden → motor, not 198 nodes padded into a square that
  hides the architecture.
- **Illustrative ≠ overclaiming.** If an effect is a many-trial statistical average
  (not visible in one rollout), let the animation show the mechanism/response and
  keep the number in the prose — never cherry-pick a seed to fake a visible
  difference. The caption must describe what is actually on screen.

## Guardrails

- Only publish content you were asked to publish. The public site is outward-facing.
- **`mkdocs build --strict` validates links and nav — not that anything *renders*.**
  Before publishing a scroll-page, animation, or math change, open the built page in
  a browser (Chromium + Playwright are preinstalled; drive `_site/` over a local
  `http.server`) and actually look: scrubs advance frame-by-frame, no baked/garbled
  text, no raw `$...$`, math typesets on notebook pages, images load, links resolve.
- Keep binaries lean; prefer GIFs under ~1 MB (tune frame count / `stride` / `dpi`
  in the animation scripts).
- Do not merge `deploy-viz-page` into `main`.
