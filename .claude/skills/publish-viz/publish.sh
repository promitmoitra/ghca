#!/usr/bin/env bash
#
# publish.sh — stage curated files onto the deploy-viz-page branch IN AN ISOLATED
# GIT WORKTREE, then validate the build there. Agent-agnostic: plain git + mkdocs,
# no proprietary tools.
#
# Why a worktree (project convention): a publish must never disturb your main
# working tree. The old flow ran `git checkout -B deploy-viz-page` in place, which
# switched your checked-out branch (usually `main`) and could strand edits or
# leave you on the wrong branch after a push. Deploying from a dedicated worktree
# keeps `main` exactly where it was; the deploy branch is checked out somewhere
# else and cleaned up when you are done. See SKILL.md.
#
# Usage:
#   publish.sh <source-ref> <path> [<path> ...]
#
# Example:
#   publish.sh main docs/e8_results.md docs/figures/e8_tone.png
#
# It creates (or refreshes) a worktree checked out to deploy-viz-page under a
# gitignored path, brings the requested paths in from <source-ref>, and (if a
# mkdocs runner is available) runs a strict build check IN THE WORKTREE. It
# intentionally stops before commit/push so you can wire nav/embeds and review the
# diff first — see SKILL.md steps 2-5. Run git add/commit/push FROM the worktree.

set -euo pipefail

DEPLOY_BRANCH="deploy-viz-page"
REPO_ROOT="$(git rev-parse --show-toplevel)"
WORKTREE_DIR="$REPO_ROOT/.publish-worktree"

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <source-ref> <path> [<path> ...]" >&2
  echo "  e.g. $0 main docs/e8_results.md docs/figures/e8_tone.png" >&2
  exit 2
fi

SRC="$1"; shift

echo ">> fetching $DEPLOY_BRANCH and $SRC"
git fetch origin "$DEPLOY_BRANCH" "$SRC"

# (Re)create the isolated worktree for the deploy branch. Idempotent: if the
# worktree already exists, reset it to the freshly-fetched deploy branch tip.
if git worktree list --porcelain | grep -qx "worktree $WORKTREE_DIR"; then
  echo ">> refreshing existing worktree at $WORKTREE_DIR"
  git -C "$WORKTREE_DIR" checkout -B "$DEPLOY_BRANCH" "origin/$DEPLOY_BRANCH"
else
  echo ">> creating worktree at $WORKTREE_DIR for $DEPLOY_BRANCH"
  git worktree add -B "$DEPLOY_BRANCH" "$WORKTREE_DIR" "origin/$DEPLOY_BRANCH"
fi

echo ">> bringing paths from $SRC into the worktree"
git -C "$WORKTREE_DIR" checkout "$SRC" -- "$@"

echo
echo ">> staged in worktree ($WORKTREE_DIR):"
git -C "$WORKTREE_DIR" status --short

# Choose an mkdocs runner: an installed mkdocs, else an ephemeral uvx invocation.
if command -v mkdocs >/dev/null 2>&1; then
  MKDOCS=(mkdocs)
elif command -v uvx >/dev/null 2>&1; then
  MKDOCS=(uvx --with mkdocs-material mkdocs)
else
  MKDOCS=()
fi

echo
if [ "${#MKDOCS[@]}" -gt 0 ]; then
  echo ">> validating build (${MKDOCS[*]} build --strict) in the worktree"
  # Non-fatal: a brand-new page fails --strict until you add it to nav (step 1
  # below). Report the result but do not abort the staging run.
  if ( cd "$WORKTREE_DIR" && "${MKDOCS[@]}" build --strict --site-dir /tmp/publish-viz-check ); then
    echo ">> build OK"
  else
    echo ">> build --strict failed — expected if a new page is not yet in nav;"
    echo "   wire mkdocs.yml nav (step 1) and rebuild before pushing."
  fi
else
  echo ">> no mkdocs runner found — install uv (for 'uvx') or 'pip install mkdocs-material'."
fi

cat <<EOF

Next (see SKILL.md) — all IN THE WORKTREE:
  cd $WORKTREE_DIR
  1. Edit mkdocs.yml nav + add any figure embeds in the doc. When updating a doc
     that already exists on the deploy branch, prefer a surgical edit over
     overwriting from <source-ref>: the deploy branch owns publication-layer
     content (e.g. GIF embeds) that <source-ref> may not have.
  2. ${MKDOCS[*]:-mkdocs} build --strict     # must pass before you push
  3. git add -A && git commit -m "publish: <what> to the docs site"
     git push origin $DEPLOY_BRANCH
  4. Confirm the "Deploy docs" Actions run is green.
  5. cd "$REPO_ROOT" && git worktree remove --force .publish-worktree   # clean up
     (--force: after the push nothing is at risk, and it sidesteps a spurious
      "working trees containing submodules" refusal from a transient nested .git)
EOF
