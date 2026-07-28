# Implementation Plan: GitOps Branch Inventory, Protection & Pre-Pruning Plan

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

---

## Comprehensive Branch Inventory & Classification

Below is the complete classification for all local and remote branches identified across the repository:

### 1. Protected Core & Infrastructure Branches (KEEP / DO NOT DELETE)

| Branch Name | Status | Rationale / Protected Purpose |
| :--- | :--- | :--- |
| `main` / `origin/main` | Active | Primary development branch (`6aa86bb`). |
| `deploy-viz-page` | Protected | Public GitHub Pages site deployment branch (`5d5bb0b`). |
| `planning-and-review` | Protected | Core roadmap & process review branch (`3403858`). |
| `project-config` | Protected | Agent config, `.agents`, `GEMINI.md` initial setup (`2ba9b72`). |
| `agent-comms-log` | Protected | Cross-agent comms log (`0c6e0e7` on origin). |
| `track/closed_loop_extensions_20260728` | Active Track | Closed-loop extensions TDD remediation (`ffb237c`). |
| `feat/experiment-review-and-uv` | Infrastructure | `experiment-review` skill & `uv` lockfile setup (`08a4cb2`). |
| `publish-viz-skill` | Infrastructure | `publish-viz` skill definition (`1d0b8bb`). |

---

### 2. Local Tracking Branches Ready for Pruning (Fully Merged into Main/Comms)

All local tracking branches below correspond to Pull Requests that have already been reviewed and merged into `origin/main` or `origin/agent-comms-log`:

| Local Branch | HEAD Commit | Target Action | Status / PR Reference |
| :--- | :---: | :--- | :--- |
| `claude/comms-followup-remediation` | `9660e99` | Delete Local | Merged into `agent-comms-log` via PR #69 (`0c6e0e7`). |
| `claude/track-c-composite-task` | `c4266aa` | Delete Local | Merged into `main` via PR #68 (`6aa86bb`). |
| `claude/comms-feedback-closed-loop-ext` | `73bc6a3` | Delete Local | Merged into `agent-comms-log` via PR #67 (`4dd01e1`). |
| `claude/scroll-narrative-repair` | `feb7ff1` | Delete Local | Merged into `main` & `deploy-viz-page` via PR #66 (`5d5bb0b`). |
| `claude/harness-resolve-helper-from-git` | `6aa86bb` | Delete Local | Equivalent to `origin/main` (`6aa86bb`). |
| `docs/correct-e5-rir-and-variance-framing` | `7f0092d` | Delete Local | Merged into `main` via PR #60. |
| `claude/phase2-table-from-archive` | `c2d730f` | Delete Local | Merged into `main` via PR #62. |
| `claude/packaging-importable` | `55ff0bc` | Delete Local | Merged into `main` via PR #61. |
| `claude/loop-margin-control` | `e6c6e51` | Delete Local | Merged into `main` via PR #57. |
| `review-54` | `a059144` | Delete Local | Stale review worktree branch (PR #54 handled). |

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

### Step 1: Harmonize & Fast-Forward Core Local Branches
```bash
# 1. Fast-forward local main to origin/main
git checkout main
git merge --ff-only origin/main

# 2. Fast-forward local agent-comms-log to origin/agent-comms-log
git checkout agent-comms-log
git merge --ff-only origin/agent-comms-log

# 3. Return to active track branch
git checkout track/closed_loop_extensions_20260728
```

### Step 2: Safe Local Branch Cleanup (Merged Feature Branches)
```bash
# Delete fully merged local tracking branches
git branch -d claude/comms-followup-remediation \
               claude/track-c-composite-task \
               claude/comms-feedback-closed-loop-ext \
               claude/scroll-narrative-repair \
               claude/harness-resolve-helper-from-git \
               docs/correct-e5-rir-and-variance-framing \
               claude/phase2-table-from-archive \
               claude/packaging-importable \
               claude/loop-margin-control \
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
  docs/correct-e5-rir-and-variance-framing
```

---

## Verification Plan

### Automated Checks
- Run `git branch -vv` to verify only protected core branches and active track branches remain locally.
- Run `git branch -r` to verify remote branch cleanliness.
- Run `python3 reproduce_all.py` to confirm that branch harmonization did not break any reproducibility assertions or unit tests.

### Manual Verification
- Verify that `main`, `deploy-viz-page`, `planning-and-review`, `project-config`, and `agent-comms-log` remain completely intact and protected.
