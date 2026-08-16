# Implementation Plan: GitOps Branch Inventory, Protection & Pre-Pruning Plan

> [!CAUTION]
> **Do not execute the step lists below as written.** They were authored
> 2026-07-28 against a repository state that no longer exists, and the safety
> analysis they rely on is incomplete in a way that destroys work. Read this
> header, then run `scripts/branch_safety_check.sh` and act on **its** output.
> The plan's *structure* (protect → land open PRs → harmonize → prune → branch
> out) is still the right sequence; only its lists are stale.
>
> **1. The ancestry test used here is unsafe on this repo — but not for the
> reason the inventory gives.**
> - **Verified and real:** the repo has used both squash-merges and merge
>   commits, so `git merge-base --is-ancestor` reports "unmerged" for branches
>   whose content is fully landed. Reproduced here: 15 branches are `ANCESTOR=no,
>   CHERRY=0`. This over-preserves; harmless. (Recorded in the PR #71 review.)
> - **What actually loses work is per-path novelty, not ancestry.** A branch can
>   be the sole holder of a tracked path that has never existed anywhere in
>   `main`'s history, while looking like an ordinary feature branch. `git cherry`
>   does not catch this; only walking the paths does. The check script below runs
>   that test and refuses to mark such a branch SAFE unless an `archive/*` tag
>   covers its tip.
> - **Correction to [`branch_preservation_inventory.md`](branch_preservation_inventory.md)
>   (on `main`):** that doc states "15 branches have no merge-base with `main` at
>   all… they predate a history rewrite," and names three. **This does not
>   reproduce — including at the inventory's own stated basis, `main` @
>   `6f41eca`.** All three named branches have real merge-bases there and now
>   (`55da552`, `eec7554d`, `857a84c8`; the first is a 2022 commit that is an
>   ancestor of both). Across every remote head (64 when checked; the count
>   moves) exactly **one** ref lacks a merge-base, and it is
>   `__dolt_remote_info__`, a Dolt metadata ref rather than a code branch. The most likely cause is an incomplete local object graph in
>   the session that produced the inventory: `git merge-base` exits non-zero when
>   objects are missing, and that error reads identically to "no common
>   ancestor." **Fetch fully before concluding two refs share no history.**
> - The "no merge-base ⇒ hard stop" rule is still worth keeping — it is one line,
>   and an unfetched-objects error should stop a deletion script rather than be
>   interpreted — but keep it as cheap defence, not as a description of a known
>   15-branch class. The inventory's sole-holder analysis is the part that
>   protects work, and that part reproduces.
>
> **2. The inventory is a snapshot, not a live view.** It was verified against
> `main` @ `6f41eca`. Re-verified 2026-08-16 against `main` @ `cfb60c8`, its
> sole-holder count fell from **10 to 5**, and two of those five dissolve on
> inspection (`claude/comms-claude-dir-tmpfs` is fully landed on
> `agent-comms-log` — `git cherry` shows 0 unlanded patches, it looks novel only
> because `.agents/comms_log.md` does not exist on `main`;
> `publish/k-dyn-correction` is content-identical to `deploy-viz-page`, so its
> "novel" paths are just the site infra). Three genuine sole holders remain, and
> **all three are covered by `archive/*` tags at their exact tips**, so the work
> is protected by immutable refs independent of branch lifetime. Never act on a
> stored table; re-run the check.
>
> **3. Most of this plan is already done.** Step 0's PRs (#70, #71) are merged,
> as are #80 and #84 since. **All nine** of Step 3's remote deletion targets are
> already gone. The `archive/*` tags Step 0 asks for exist (10 on the remote).
> Re-derive the current state before running anything.
>
> — reconciliation added 2026-08-16 (session `d560b36c`), verified against
> `origin/main` @ `cfb60c8`.

## Goal Description
Conduct a comprehensive GitOps pre-pruning audit of all local and remote branches in the `promitmoitra/ghca` repository. Based on current project state (PRs #1–#69), `main` history, and cross-agent coordination entries in `agent-comms-log`, classify every branch into protection, merge/rebase, or pruning categories, and establish a safe step-by-step removal and harmonization sequence.

---

## User Review Required

> [!IMPORTANT]
> **Explicit Protection Constraints:**
> Per user rules and cross-agent research requirements, the following core branches will **NEVER be deleted**:
> 1. `main` / `origin/main` — Production codebase & primary research state.
> 2. `deploy-viz-page` / `origin/deploy-viz-page` — Public scrollytelling site build & deploy branch.
> 3. `planning-and-review` / `origin/planning-and-review` — Process guidelines, core/extensions review templates, and roadmap tracking.
> 4. `project-config` / `origin/project-config` — Project initialization, `.agents` metadata, `GEMINI.md`, and `CLAUDE.md`.
> 5. `agent-comms-log` / `origin/agent-comms-log` — Asynchronous cross-agent communication & coordination log.
> 6. `track/closed_loop_extensions_20260728` / `origin/track/closed_loop_extensions_20260728` — Active track branch containing TDD remediation commit `ffb237c`.
> 7. `feat/experiment-review-and-uv` / `origin/feat/experiment-review-and-uv` — Project skill tooling & `uv` environment configuration.
> 8. `publish-viz-skill` / `origin/publish-viz-skill` — Standardized visualization publishing skill definition.
> 9. `claude/harness-resolve-helper-from-git` / `origin/claude/harness-resolve-helper-from-git` — Active open PR #70 (`12ba4ba`).

---

## Comprehensive Branch Inventory & Classification

Below is the complete classification for all local and remote branches identified across the repository:

### 1. Protected Core, Active PRs & Infrastructure Branches (KEEP / DO NOT DELETE)

| Branch Name | Status | Rationale / Protected Purpose |
| :--- | :--- | :--- |
| `main` / `origin/main` | Active | Primary development branch (`6aa86bb`). |
| `deploy-viz-page` | Protected | Public GitHub Pages site deployment branch (`5d5bb0b`). |
| `planning-and-review` | Protected | Core roadmap & process review branch (`3403858`). |
| `project-config` | Protected | Agent config, `.agents`, `GEMINI.md` initial setup (`2ba9b72`). |
| `agent-comms-log` | Protected | Cross-agent comms log (`0c6e0e7` / `503385c` on origin). |
| `track/closed_loop_extensions_20260728` | Active Track | Closed-loop extensions TDD remediation (`ffb237c`). |
| `feat/experiment-review-and-uv` | Infrastructure | `experiment-review` skill & `uv` lockfile setup (`08a4cb2`). |
| `publish-viz-skill` | Infrastructure | `publish-viz` skill definition (`1d0b8bb`). |
| `claude/harness-resolve-helper-from-git` | Active PR #70 | Open PR #70 (`12ba4ba`); contains `_resolve_review_helper` fix. |
| `__dolt_remote_info__` | Dolt Metadata | Dolt system tracking ref (`e482d38`). |

---

### 2. Local Tracking Branches Ready for Pruning (Fully Merged into Main/Comms)

All local tracking branches below correspond to Pull Requests that have already been reviewed and merged into `origin/main` or `origin/agent-comms-log`:

| Local Branch | HEAD Commit | Target Action | Status / PR Reference |
| :--- | :---: | :--- | :--- |
| `claude/comms-followup-remediation` | `9660e99` | Delete Local (`-d`) | Merged into `agent-comms-log` via PR #69 (`0c6e0e7`). |
| `claude/track-c-composite-task` | `c4266aa` | Delete Local (`-d`) | Merged into `main` via PR #68 (`6aa86bb`). |
| `claude/comms-feedback-closed-loop-ext` | `73bc6a3` | Delete Local (`-D`) | Pre-rebase local orphan (`73bc6a3`); remote ref `4ea9fe2` merged via PR #67. |
| `claude/scroll-narrative-repair` | `feb7ff1` | Delete Local (`-d`) | Merged into `main` & `deploy-viz-page` via PR #66 (`5d5bb0b`). |
| `docs/correct-e5-rir-and-variance-framing` | `7f0092d` | Delete Local (`-d`) | Merged into `main` via PR #60. |
| `claude/phase2-table-from-archive` | `c2d730f` | Delete Local (`-d`) | Merged into `main` via PR #62. |
| `claude/packaging-importable` | `55ff0bc` | Delete Local (`-d`) | Merged into `main` via PR #61. |
| `claude/loop-margin-control` | `e6c6e51` | Delete Local (`-d`) | Merged into `main` via PR #57. |
| `review-54` | `a059144` | Delete Local (`-D`) | Stale review worktree branch; files byte-identical to `main`. |


---

### 3. Remote Feature Branches Ready for Pruning (Merged PRs)

The following remote branches correspond to completed, merged PRs (#1–#65) and are safe for remote cleanup:

- `origin/c5-animation` (PR #15)
- `origin/c6-animation` (PR #20)
- `origin/claude/1b-direction-readout` (PR #40)
- `origin/claude/2b-spiral-predictions` (PR #7/#8)
- `origin/claude/3a-p2-sweeps` (PR #26)
- `origin/claude/3a-p3-cseries-e8` (PR #27)
- `origin/claude/3a-p3b-sigma-bands` (PR #33)
- `origin/claude/3a-p4-headline-edits` (PR #28)
- `origin/claude/3b-learning-port` (PR #31)
- `origin/claude/3b-other-topologies` (PR #29/#32)
- `origin/claude/3c-4a-timescale-capacity` (PR #43)
- `origin/claude/3c-capacity` (PR #38)
- `origin/claude/3c-continual-learning-plan` (PR #34)
- `origin/claude/3c-e9-bridge` (PR #39)
- `origin/claude/3c-p2-causal-credit` (PR #35)
- `origin/claude/3c-p3-lowvar-credit` (PR #36)
- `origin/claude/3c-saturation` (PR #41)
- `origin/claude/3d-emergent-timescale` (PR #45)
- `origin/claude/3d-timescale-capacity` (PR #44)
- `origin/claude/3e-retile` (PR #46)
- `origin/claude/3e2-cfc` (PR #48)
- `origin/claude/3e2-hierarchy` (PR #47)
- `origin/claude/3e3-concurrent` (PR #49)
- `origin/claude/4b-causal-testbed` (PR #9)
- `origin/claude/5a-synthesis` (PR #14)
- `origin/claude/ai-hallucination-review-ya7aq5` (PR #8)
- `origin/claude/e-series-experiments-review-dew74b` (PR #2–#6)
- `origin/claude/e10-timescale-hierarchy` (PR #12/E10 notes)
- `origin/claude/e8-hardening` (PR #12)
- `origin/claude/greenberg-hastings-network-j0yerd` (PR #1)
- `origin/claude/loop-margin-control` (PR #57)
- `origin/claude/next-steps-planning` (PR #18)
- `origin/claude/packaging-importable` (PR #61)
- `origin/claude/phase2-table-from-archive` (PR #62)
- `origin/claude/pr-templates` (PR #19)
- `origin/claude/reproducibility-framing-fixes` (PR #5)
- `origin/claude/scroll-animation-scripts` (PR #59)
- `origin/claude/scroll-narrative-repair` (PR #66)
- `origin/claude/synthesis-refresh-capacity` (PR #37)
- `origin/claude/timescale-arc-synthesis` (PR #50)
- `origin/claude/track-c-composite-task` (PR #68)
- `origin/claude/what-you-see-0dgx0h` (PR #6)
- `origin/e7-animation` (PR #13)
- `origin/rephrase-review-terms` (PR #17)
- `origin/track4a-synthesis` (PR #42)
- `origin/viz-animator` (PR #10)

---

## Step-by-Step Execution Sequence

### Step 0: Open PR Consolidation & Active Track Landing
```bash
# 1. Merge PR #70 (harness: resolve review_helper.py from git) into main
gh pr merge 70 --squash --delete-branch=false

# 2. Merge PR #71 (comms: review of the GitOps pruning plan) into agent-comms-log
gh pr merge 71 --squash --delete-branch=false

# 3. Create PR and merge active track/closed_loop_extensions_20260728 into main
gh pr create --head track/closed_loop_extensions_20260728 --base main \
  --title "track: closed loop extensions motor-drive cue timing remediation & retention metrics" \
  --body "Remediates cue timing, adds chance-corrected retention metrics, unit tests, and automated table generator."
gh pr merge --squash
```

### Step 1: Harmonize & Fast-Forward Core Local Branches
```bash
# 1. Fast-forward local main to origin/main
git checkout main
git merge --ff-only origin/main

# 2. Fast-forward local agent-comms-log to origin/agent-comms-log
git checkout agent-comms-log
git merge --ff-only origin/agent-comms-log
```

### Step 2: Safe Local Branch Cleanup (Merged Feature Branches)
```bash
# Delete fully merged local tracking branches (-d)
git branch -d claude/comms-followup-remediation \
               claude/track-c-composite-task \
               claude/scroll-narrative-repair \
               docs/correct-e5-rir-and-variance-framing \
               claude/phase2-table-from-archive \
               claude/packaging-importable \
               claude/loop-margin-control

# Delete orphan/stale local tracking branches (-D)
git branch -D claude/comms-feedback-closed-loop-ext \
               review-54
```

### Step 3: Remote Stale Branch Pruning
```bash
# Fetch and prune stale remote tracking references
git fetch --prune origin

# Delete merged feature branches on remote (batch clean)
git push origin --delete \
  claude/comms-followup-remediation \
  claude/track-c-composite-task \
  claude/comms-feedback-closed-loop-ext \
  claude/scroll-narrative-repair \
  claude/phase2-table-from-archive \
  claude/packaging-importable \
  claude/loop-margin-control \
  docs/correct-e5-rir-and-variance-framing \
  claude/comms-review-pruning-plan
```

### Step 4: Branch Out for Next Development Track
```bash
# Checkout clean updated main and branch out for next track
git checkout main
git pull --ff-only origin main
git checkout -b track/closed_loop_phase1_3_rerun
```

---

## Verification Plan

### Automated Checks
- **Before deleting anything**, run `bash scripts/branch_safety_check.sh` and act
  on its output rather than on any list in this document or in
  `branch_preservation_inventory.md`. It applies three tests (ancestry, patch
  equality via `git cherry`, per-path novelty) and four hard stops (open-PR head,
  protected branch, missing merge-base, untagged sole holder). Snapshot at the time of writing (2026-08-16, `main` @ cfb60c8, 65 remote
  heads): **41 SAFE / 17 KEEP** — a dated observation, not a standing
  figure. The script is the authority; these numbers move whenever a branch is
  pushed.
- Tag before deleting. A branch the script marks `sole holder … UNTAGGED` becomes
  SAFE once an `archive/*` tag points at its tip; deletion is irreversible from
  the client side, tagging costs nothing.
- Run `git branch -vv` and `git branch -r` to confirm the resulting state matches
  the script's SAFE/KEEP partition.
- Run `python3 reproduce_all.py` to confirm that branch harmonization did not break any reproducibility assertions or unit tests.

### Manual Verification
- Verify that `main`, `deploy-viz-page`, `planning-and-review`, `project-config`, and `agent-comms-log` remain completely intact and protected.

