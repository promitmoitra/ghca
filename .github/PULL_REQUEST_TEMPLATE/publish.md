<!--
Publish to the GitHub Pages site — the MANUAL write path to `deploy-viz-page`.

Base this PR on `deploy-viz-page` (head = your publish branch). This is the
reviewable path: a maintainer sees exactly what goes public before it deploys,
and merging triggers the "Deploy docs" workflow. The `publish-viz` skill AIDS
this — it does the cherry-pick, nav wiring, and `mkdocs build --strict` for you —
but the manual path below stands on its own (plain git + mkdocs, no skill needed).

Flow is one-way INTO `deploy-viz-page`. Never open a `deploy-viz-page → main` PR.
-->

## Publishing
- **Source ref** (where the content is mature — usually `main`):
- **Paths** (result docs, figures, GIFs):

## Manual write path (what this PR did)
```
git fetch origin deploy-viz-page <source-ref>
git checkout -B <publish-branch> origin/deploy-viz-page
git checkout <source-ref> -- <path> [<path> ...]
# wire new pages into mkdocs.yml nav, then:
mkdocs build --strict
```
<!-- Or let the skill do the above: .claude/skills/publish-viz/publish.sh <source-ref> <path>… -->

## Checks
- [ ] Content is mature on the source ref (published downstream of the trunk, not ahead of it)
- [ ] New pages wired into `mkdocs.yml` nav (labels consistent with existing)
- [ ] `mkdocs build --strict` passes (no broken links / files missing from nav)
- [ ] Only content intended to be public — the site is outward-facing
- [ ] Binaries lean (GIFs ≲ 1 MB)
