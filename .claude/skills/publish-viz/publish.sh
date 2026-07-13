#!/usr/bin/env bash
#
# publish.sh — stage curated files onto the deploy-viz-page branch for the
# GitHub Pages site. Agent-agnostic: plain git + mkdocs, no proprietary tools.
#
# Usage:
#   publish.sh <source-ref> <path> [<path> ...]
#
# Example:
#   publish.sh main docs/e8_results.md docs/figures/e8_tone.png
#
# It fetches the deploy branch, checks it out, brings the requested paths in
# from <source-ref>, and (if mkdocs is available) runs a strict build check.
# It intentionally stops before commit/push so you can wire nav/embeds and
# review the diff first — see SKILL.md steps 2-5.

set -euo pipefail

DEPLOY_BRANCH="deploy-viz-page"

if [ "$#" -lt 2 ]; then
  echo "usage: $0 <source-ref> <path> [<path> ...]" >&2
  echo "  e.g. $0 main docs/e8_results.md docs/figures/e8_tone.png" >&2
  exit 2
fi

SRC="$1"; shift

echo ">> fetching $DEPLOY_BRANCH and $SRC"
git fetch origin "$DEPLOY_BRANCH" "$SRC"

echo ">> checking out $DEPLOY_BRANCH"
git checkout -B "$DEPLOY_BRANCH" "origin/$DEPLOY_BRANCH"

echo ">> bringing paths from $SRC"
git checkout "$SRC" -- "$@"

echo
echo ">> staged from $SRC:"
git status --short

echo
if command -v mkdocs >/dev/null 2>&1; then
  echo ">> validating build (mkdocs build --strict)"
  mkdocs build --strict --site-dir /tmp/publish-viz-check
  echo ">> build OK"
else
  echo ">> mkdocs not found — run 'pip install mkdocs-material' then 'mkdocs build --strict'"
fi

cat <<EOF

Next (see SKILL.md):
  1. Edit mkdocs.yml nav + add any figure embeds in the doc.
  2. mkdocs build --strict          # if you changed nav/docs after this run
  3. git add -A && git commit -m "publish: <what>" && git push origin $DEPLOY_BRANCH
  4. Confirm the "Deploy docs" Actions run is green.
EOF
