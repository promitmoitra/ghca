#!/usr/bin/env bash
# Classify every remote branch as SAFE (content is upstream) or KEEP, and flag
# the two failure modes an ancestry test misses on this repo.
#
# Why this exists: docs/gitops_branch_pruning_plan.md and
# docs/branch_preservation_inventory.md both carry branch LISTS. Lists go stale
# -- the inventory's sole-holder count fell 10 -> 5 in two days. Run this
# instead of reading either table. The rules are durable; the names are not.
#
# Three tests, in increasing order of authority:
#   1. ancestry      -- git merge-base --is-ancestor: fast, but blind to
#                       squash-merges (over-preserves; harmless).
#   2. patch-equality -- git cherry: catches squash-merges.
#   3. per-path novelty -- does any tracked path exist NOWHERE in main's
#                       history? This is the only test that detects a sole
#                       holder of orphaned work, and the only one that is
#                       meaningful when there is no merge-base.
#
# Hard stops, applied regardless of what the tests say:
#   - open PR head       -> never delete
#   - protected branch   -> never delete
#   - NO MERGE-BASE      -> never delete without a human decision; ancestry
#                           results on these branches are meaningless
#   - sole holder        -> never delete unless an archive/* tag covers the tip
#
# Usage:  bash scripts/branch_safety_check.sh [--upstream origin/main]
# Needs:  GITHUB_TOKEN for the open-PR query (degrades to a warning without it).
set -uo pipefail

UPSTREAM="origin/main"
[[ "${1:-}" == "--upstream" ]] && UPSTREAM="${2:-origin/main}"

PROT='origin/HEAD|origin/main|origin/deploy-viz-page|origin/planning-and-review|origin/project-config|origin/agent-comms-log|origin/publish-viz-skill|origin/__dolt_remote_info__'

git fetch origin --prune --quiet 2>/dev/null

OPEN_HEADS=""
if [[ -n "${GITHUB_TOKEN:-}" ]]; then
  OPEN_HEADS=$(curl -s --max-time 25 -H "Authorization: Bearer $GITHUB_TOKEN" \
    "https://api.github.com/repos/promitmoitra/ghca/pulls?state=open&per_page=100" \
    | python3 -c 'import sys,json;print("\n".join(p["head"]["ref"] for p in json.load(sys.stdin)))' 2>/dev/null)
else
  echo "WARNING: no GITHUB_TOKEN -- open-PR heads are NOT protected in this run." >&2
fi

TAGGED=$(git ls-remote --tags origin 2>/dev/null | awk '{print substr($1,1,40)}')

# every path that has ever existed in the upstream history
MAIN_PATHS=$(mktemp)
git log --pretty=format: --name-only "$UPSTREAM" | sort -u > "$MAIN_PATHS"

printf '%-58s %-7s %-7s %-6s %s\n' BRANCH ANCESTOR CHERRY NOVEL VERDICT
for b in $(git for-each-ref --format='%(refname:short)' refs/remotes/origin); do
  [[ "$b" =~ ^($PROT)$ ]] && continue
  short="${b#origin/}"

  mb=$(git merge-base "$b" "$UPSTREAM" 2>/dev/null || true)
  anc=no; [[ -n "$mb" ]] && git merge-base --is-ancestor "$b" "$UPSTREAM" 2>/dev/null && anc=yes
  cherry=$(git cherry "$UPSTREAM" "$b" 2>/dev/null | grep -c '^+' || true)

  novel=0
  while IFS= read -r p; do
    [[ -z "$p" ]] && continue
    grep -qxF "$p" "$MAIN_PATHS" || novel=$((novel + 1))
  done < <(git ls-tree -r --name-only "$b" 2>/dev/null)

  tip=$(git rev-parse "$b")
  verdict="SAFE"
  if grep -qxF "$short" <<< "$OPEN_HEADS"; then verdict="KEEP: open PR"
  elif [[ -z "$mb" ]]; then verdict="KEEP: NO MERGE-BASE (hard stop)"
  elif (( novel > 0 )); then
    if grep -qF "$tip" <<< "$TAGGED"; then verdict="SAFE: sole holder, archive tag covers tip"
    else verdict="KEEP: sole holder of $novel path(s), UNTAGGED"; fi
  elif (( cherry > 0 )) && [[ "$anc" == no ]]; then verdict="KEEP: $cherry patch(es) not upstream"
  fi

  printf '%-58s %-7s %-7s %-6s %s\n' "$short" "$anc" "$cherry" "$novel" "$verdict"
done
rm -f "$MAIN_PATHS"

cat <<'NOTE'

Reading this table:
  ANCESTOR=no with CHERRY=0  -> squash-merged; content is upstream.
  NOVEL>0                    -> tracked paths absent from the whole upstream
                                history. Verify by hand before deleting; an
                                archive/* tag at the tip makes deletion safe.
  A blank merge-base is reported as the NO MERGE-BASE hard stop, not as ANCESTOR=no.
Deleting a ref is irreversible from this side: tag first, delete second.
NOTE
