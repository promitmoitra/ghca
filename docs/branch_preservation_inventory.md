# Branch preservation inventory — read before pruning any ref

**Verified 2026-08-14 against `origin/main` @ `6f41eca`.** 59 remote heads.

This exists because **an ancestry-based pruning pass will silently destroy work in
this repo.** Two independent reasons:

1. The repo has used **both** squash-merges and merge commits, so
   `git merge-base --is-ancestor` reports "unmerged" for branches whose content is
   fully landed — 36 of them. Trusting that signal over-preserves (harmless).
2. **15 branches have no merge-base with `main` at all** — they predate a history
   rewrite. Ancestry checks on them are meaningless, and a pruning script will
   classify them as ordinary unmerged feature branches with no warning. Three of
   them are sole holders of real content. Trusting that signal **destroys work**.

A pre-existing plan, [`gitops_branch_pruning_plan.md`][plan], lives unmerged on
`planning-and-review` and classifies branches by ancestry. **Do not execute it
as written** without reconciling against the table below.

[plan]: https://github.com/promitmoitra/ghca/blob/planning-and-review/docs/gitops_branch_pruning_plan.md

## Classification

| Class | Count | Safe to delete? |
|---|---:|---|
| (c) ancestor of `main` | 12 | yes |
| (b) content-landed on `main` (ancestry says otherwise) | 36 | yes |
| (a) sole holder of content absent from `main` | 10 | **no** |
| `main` itself | 1 | — |

Method for (a): for every tracked path on the branch, `git log --oneline
origin/main -- <path> | wc -l` → `0` means that path has never existed anywhere in
`main`'s history.

## Class (a) — do not delete

### Sole holders of orphaned work

| Branch | Tip | Unique content |
|---|---|---|
| `claude/what-you-see-0dgx0h` ⚠ no merge-base | `81c20eb` | `RESEARCH_LOG.md`, `chirality_sweep.py`, `ghca_cluster.py`, `random_chirality_sweep.py` |
| `claude/ai-hallucination-review-ya7aq5` ⚠ no merge-base | `b0ccd48` | `docs/hallucination_review.md`, plus a second copy of the MkDocs site infra (`mkdocs.yml`, `docs/index.md`, `docs/animations.md`, `.github/workflows/deploy-docs.yml`, 4 GIFs) |
| `claude/e10-timescale-hierarchy` ⚠ no merge-base | `97bc8a8` | `docs/e10_notes.md`, `experiments/e10_diagnostics.py` — **resolved**: both restored to `main`, see below |
| `planning-and-review` | `9d57f65` | `docs/gitops_branch_pruning_plan.md` (the plan above) |
| `claude/aria-olfaction-positioning` | `ca3bc68` | `docs/applied/aria_positioning.md`, `aria_olfaction_feedback.md`, `nanopore_collaboration_email.md` — **needs an owner decision, see below** |

### Unique by design — never delete

| Branch | Holds |
|---|---|
| `deploy-viz-page` | the live GitHub Pages source: `mkdocs.yml`, `.github/workflows/deploy-docs.yml`, `docs/index.md`, 8 GIFs, scroll assets |
| `agent-comms-log` | `.agents/comms_log.md`, the cross-agent channel |
| `__dolt_remote_info__` | `DOLT_REMOTE.md`, machine-managed by beads |

### Redundant, safe once the above is settled

`claude/comms-claude-dir-tmpfs` (`95735e4`) is an ancestor of `agent-comms-log`
holding an identical `.agents/comms_log.md` blob. Delete the redundant one, never
`agent-comms-log`.

## Resolved and outstanding

- **`e10_notes.md` — resolved.** `main`'s `docs/next_steps.md` linked to it six
  times while the file existed only on an orphan branch. Both it and
  `experiments/e10_diagnostics.py` are now restored to `main`; the RNG audit
  passes with them in the tree. The branch is now redundant, but is left tagged
  in this table until someone confirms nothing else on it is wanted.
- **`claude/aria-olfaction-positioning` — needs an owner decision, and it is not
  a technical one.** The three `docs/applied/` files concern olfaction
  positioning, feedback, and a nanopore collaboration approach. This repository
  is maintained as independent open-source work, and material originating from a
  commercial affiliation does not belong in it. **Deliberately not restored to
  `main`.** Decide whether it should be kept on a private remote, kept on the
  branch, or removed — do not land it here by default.

## Why this file exists instead of tags

Annotated `archive/*` tags were the intended mechanism. **Tag pushes and ref
deletions are both rejected (HTTP 403) by this session's credentials**, while
branch pushes succeed — so the preservation had to take a form that could
actually be written. A doc on `main` is arguably the better artifact anyway: it
is reviewable, diffable, and a pruning pass will read it, whereas a tag only
helps someone who already knows to look.

To create the tags from a machine with full push rights:

```bash
git tag -a archive/what-you-see            81c20eb -m 'sole holder: RESEARCH_LOG.md + 3 scripts'
git tag -a archive/hallucination-review    b0ccd48 -m 'sole holder: docs/hallucination_review.md'
git tag -a archive/e10-timescale-hierarchy 97bc8a8 -m 'sole holder: e10 notes + diagnostics'
git tag -a archive/aria-olfaction          ca3bc68 -m 'sole holder: docs/applied/ — see affiliation note'
git tag -a archive/gitops-pruning-plan     9d57f65 -m 'sole holder: gitops_branch_pruning_plan.md'
git push origin --tags
```

## Housekeeping

`tmp-proxy-probe` was pushed from this session as a write-permission probe and
**could not be deleted** (the same 403 above). It points at `8a617b0` and carries
nothing unique. Delete it:

```bash
git push origin --delete tmp-proxy-probe
```
