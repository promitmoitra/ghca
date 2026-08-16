# Branch preservation inventory — read before pruning any ref

**Verified 2026-08-16 against `origin/main` @ `a32802f`. 64 remote heads: `main`
+ 63 branches, every one classified below.** The counts reconcile exactly
(9 + 30 + 24 = 63); if they ever stop reconciling, the table is stale and a
pruning pass must re-run the method in [Reproducing this](#reproducing-this)
rather than trust it.

This exists because **an ancestry-based pruning pass will silently destroy work
in this repo.** The repo has used both squash-merges and merge commits, so
`git merge-base --is-ancestor` reports "unmerged" for branches whose content is
fully landed — **30 of the 63**. Trusting that signal over-preserves, which is
harmless. The real hazard is the other direction, and it is documented under
[Method](#method): the check that decides "safe to delete" has to be about
*content*, not ancestry.

A pre-existing plan, [`gitops_branch_pruning_plan.md`][plan], lives unmerged on
`planning-and-review` and classifies branches by ancestry. **Do not execute it
as written** without reconciling against the tables below.

[plan]: https://github.com/promitmoitra/ghca/blob/planning-and-review/docs/gitops_branch_pruning_plan.md

## Classification

| Class | Count | Safe to delete? |
|---|---:|---|
| (c) ancestor of `main` | 24 | yes |
| (b) content-landed on `main` (ancestry says otherwise) | 30 | yes |
| (a) holds a path absent from all of `main`'s history | 9 | **7 no, 2 yes — see below** |
| `main` itself | 1 | — |
| **total remote heads** | **64** | |

Class (a) is not uniformly "do not delete": 7 must be kept, and 2 are exact
duplicates of branches in that 7.

## Class (a) — holds content `main` has never had

### Sole holders of orphaned work — do not delete

| Branch | Tip | Unique content | Archive tag |
|---|---|---|---|
| `claude/what-you-see-0dgx0h` | `81c20eb` | `RESEARCH_LOG.md`, `chirality_sweep.py`, `ghca_cluster.py`, `random_chirality_sweep.py` | `archive/what-you-see` |
| `claude/ai-hallucination-review-ya7aq5` | `b0ccd48` | `docs/hallucination_review.md`, plus a second copy of the MkDocs site infra (`mkdocs.yml`, `docs/index.md`, `docs/animations.md`, `.github/workflows/deploy-docs.yml`, 4 GIFs) — 9 paths | `archive/hallucination-review` |
| `claude/aria-olfaction-positioning` | `ca3bc68` | `docs/applied/aria_positioning.md`, `aria_olfaction_feedback.md`, `nanopore_collaboration_email.md` — **needs an owner decision, see below** | `archive/aria-olfaction` |
| `planning-and-review` | `9d57f65` | `docs/gitops_branch_pruning_plan.md` (the plan above) | `archive/gitops-pruning-plan` |

All four are covered by an annotated tag, so the content survives even if the
branch is deleted. The branches are still listed as do-not-delete because a tag
is not a branch and nothing has confirmed the work is finished with.

### Unique by design — never delete, never tag

| Branch | Tip | Holds |
|---|---|---|
| `deploy-viz-page` | `cd1ade6` | the live GitHub Pages source: `mkdocs.yml`, `.github/workflows/deploy-docs.yml`, `docs/index.md`, `docs/animations.md`, 8 GIFs, MathJax + scroll assets — 17 paths |
| `agent-comms-log` | `6f2c831` | `.agents/comms_log.md`, the cross-agent channel |
| `__dolt_remote_info__` | `e482d38` | `DOLT_REMOTE.md`, machine-managed by beads |

These are live infrastructure, not archives. Tagging them would imply they are
finished; they are not.

### Redundant duplicates — safe to delete

| Branch | Tip | Duplicate of | Evidence |
|---|---|---|---|
| `claude/comms-claude-dir-tmpfs` | `95735e4` | `agent-comms-log` | ancestor of it; `git diff agent-comms-log ↔ this` is 201 deletions, 0 insertions — a strict earlier subset |
| `publish/k-dyn-correction` | `f071bfe` | `deploy-viz-page` | ancestor of it with an **identical tree** (empty diff) |

Delete the duplicate, never the branch it duplicates.

> **Correction to the previous revision.** It described
> `claude/comms-claude-dir-tmpfs` as holding an "identical `.agents/comms_log.md`
> blob". It does not — the blobs differ (`7b15f31` vs `75e12e2`) and
> `agent-comms-log` carries 201 lines more. The conclusion (safe to delete) is
> unchanged, but the stated reason was wrong; the strict-subset diff above is the
> real one.

## Class (c) — ancestors of `main`, safe to delete (24)

`c5-animation` `c6-animation` `claude/2b-spiral-predictions`
`claude/3a-p3b-sigma-bands` `claude/4b-causal-testbed` `claude/5a-synthesis`
`claude/clock-shift-merge` `claude/coherence-invariant`
`claude/cycle-packing-exact` `claude/e8-hardening`
`claude/e-series-experiments-review-dew74b` `claude/gh-vs-symbolic-viz`
`claude/greenberg-hastings-network-j0yerd` `claude/lattice-timescale-demo`
`claude/persistent-set-structure` `claude/pr-templates`
`claude/reproducibility-framing-fixes` `claude/scroll-animation-scripts`
`claude/why-the-spectrum-decides` `e7-animation`
`feat/experiment-review-and-uv` `project-config` `rephrase-review-terms`
`viz-animator`

## Class (b) — content landed, ancestry says otherwise, safe to delete (30)

`claude/1b-direction-readout` `claude/3a-p2-sweeps` `claude/3a-p3-cseries-e8`
`claude/3a-p4-headline-edits` `claude/3b-learning-port`
`claude/3b-other-topologies` `claude/3c-4a-timescale-capacity`
`claude/3c-capacity` `claude/3c-continual-learning-plan` `claude/3c-e9-bridge`
`claude/3c-p2-causal-credit` `claude/3c-p3-lowvar-credit`
`claude/3c-saturation` `claude/3d-emergent-timescale`
`claude/3d-timescale-capacity` `claude/3e2-cfc` `claude/3e2-hierarchy`
`claude/3e3-concurrent` `claude/3e-retile` `claude/e10-timescale-hierarchy`
`claude/e6-timing` `claude/next-steps-planning`
`claude/planning-review-branches-sf4y29` `claude/publish-viz-scroll-conventions`
`claude/scroll-scrub-rework` `claude/synthesis-refresh-capacity`
`claude/timescale-arc-synthesis` `publish-viz-skill` `track4a-synthesis`
`track/closed_loop_extensions_20260728`

Every one of the 798 distinct paths across these 30 branches was individually
checked against `main`'s history; none is novel.

## Resolved and outstanding

- **`e10_notes.md` — resolved, and the branch has now moved class.** `main`'s
  `docs/next_steps.md` linked to it six times while the file existed only on an
  orphan branch. Both `docs/e10_notes.md` and `experiments/e10_diagnostics.py`
  are now in `main`'s tip, so `claude/e10-timescale-hierarchy` reclassified from
  (a) to (b) — it is now safe to delete. The `archive/e10-timescale-hierarchy`
  tag (`97bc8a8`) is kept as a belt-and-braces record; it is the one archive tag
  that no longer guards a sole holder.
- **`claude/aria-olfaction-positioning` — still needs an owner decision, and it
  is not a technical one.** The three `docs/applied/` files concern olfaction
  positioning, feedback, and a nanopore collaboration approach. This repository
  is maintained as independent open-source work, and material originating from a
  commercial affiliation does not belong in it. **Deliberately not restored to
  `main`.** Decide whether it should be kept on a private remote, kept on the
  branch, or removed — do not land it here by default.

## Method

Class (a) test, per branch: for every tracked path at the branch tip, does that
path appear anywhere in `main`'s history?

```bash
git log --oneline origin/main -- <path> | wc -l   # 0 => never existed in main
```

Any single path returning `0` makes the branch class (a). This is deliberately a
*path*-level test over `main`'s whole history, not a tip-vs-tip content diff: it
errs toward over-preserving, which is harmless, and it is what the earlier
revisions of this file used.

**Ancestry is not the test, and it is not a tiebreaker.** 30 of 63 branches are
fully landed yet report "not an ancestor". The reverse — a branch that *is* an
ancestor but holds novel content — cannot happen, which is why class (c) is
short-circuited.

### Two traps, both of which produced wrong answers in earlier revisions

1. **Counting novel paths with `wc -l` on a variable.** `printf '%s' "$novel" |
   wc -l` counts *newlines*, so a branch with exactly one novel path counts `0`
   and is reclassified as "safe to delete". This silently demoted
   `__dolt_remote_info__`, `agent-comms-log`, and `planning-and-review` — three
   must-keep branches — in a draft of this very table. Count non-empty lines
   (`grep -c .`) instead. The failure is silent and points the dangerous way.
2. **A previous revision claimed "15 branches have no merge-base with `main`,
   they predate a history rewrite", and built a second justification on it. That
   was wrong when written.** Re-checked against that revision's own verification
   point `6f41eca`, exactly **one** branch had no merge-base, and it still does:
   `__dolt_remote_info__`, which is an orphan by design because beads manages it.
   The most likely cause is running `git merge-base` against branches that were
   never fetched, so the command failed and the failure was read as "no common
   ancestor". Fetch everything first:
   `git fetch origin '+refs/heads/*:refs/remotes/origin/*' --prune`.

### Reproducing this

The per-path form above is O(paths x history). The equivalent fast form builds
the set of every path ever in `main` once, then set-differences each branch tree
against it:

```bash
{ git log --pretty=format: --name-only origin/main
  git ls-tree -r --name-only origin/main; } | sed '/^$/d' | sort -u > main_paths
git ls-tree -r --name-only origin/<branch> | sort -u | comm -23 - main_paths
```

Non-empty output => class (a). Both forms were run for this revision and agree
on all 63 branches.

## Archive tags

Annotated `archive/*` tags exist as of 2026-08-16 and are the machine-checkable
backstop; this doc is the reviewable one. Both are kept — a tag only helps
someone who already knows to look, and a doc is what a pruning pass will read.

| Tag | Commit |
|---|---|
| `archive/what-you-see` | `81c20eb` |
| `archive/hallucination-review` | `b0ccd48` |
| `archive/e10-timescale-hierarchy` | `97bc8a8` |
| `archive/aria-olfaction` | `ca3bc68` |
| `archive/gitops-pruning-plan` | `9d57f65` |

All five are annotated tag objects carrying "sole holder" messages, each
resolving to the branch tip recorded above.

### The 403 was never a permissions limit

Worth recording so the next agent does not re-derive it. This `gh` build has no
`gh auth token` subcommand, so a `git push` wired to `$(gh auth token)` sends
*help text* as the password and GitHub answers "Invalid username or token" —
which reads like a rights problem and is not one. The stored `gh` credential has
`repo` scope and is fine for both tag creation and ref deletion; it is reachable
through `gh api` even when the git-over-HTTPS path is not:

```bash
# create an annotated tag: object first, then the ref
SHA=$(gh api -X POST repos/promitmoitra/ghca/git/tags \
        -f tag='archive/NAME' -f message='MSG' -f object='COMMIT' -f type=commit --jq .sha)
gh api -X POST repos/promitmoitra/ghca/git/refs -f ref="refs/tags/archive/NAME" -f sha="$SHA"

# delete a ref
gh api -X DELETE repos/promitmoitra/ghca/git/refs/heads/BRANCH
```

## Housekeeping

`tmp-proxy-probe` (`8a617b0`) was pushed as a write-permission probe. **Deleted
2026-08-16.** Before deleting, it was confirmed to hold nothing unique: `git diff
origin/tmp-proxy-probe origin/main` was pure insertions, i.e. no path existed on
the probe that `main` lacks. Ancestry alone would not have shown this — the probe
was *not* an ancestor of `main`, the same class (b) false negative this whole
file exists to warn about.
