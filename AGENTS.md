# Working in this repo

This is an exploratory computational-neuroscience study: Greenberg–Hastings
excitable dynamics on a graph as the substrate for a reward-driven learner
(E-series), plus a causal-inference instrumentation of the same substrate
(C-series). See [`README.md`](README.md) for the map of files and results.

This file is **agent-agnostic** (the [`publish-viz`](.claude/skills/publish-viz/SKILL.md)
convention): it is plain guidance any agent — or a human — can follow, with no
proprietary tools. `AGENTS.md` is the vendor-neutral name coding agents look for.

## Process — read before a planning or review pass

The project runs two recurring meta-passes, **kept decoupled in process but
linked by a one-directional hand-off (review → plan)**. Before doing either,
read [`docs/process.md`](docs/process.md). In short:

- **Review** (adversarial, backward-looking) re-runs experiments, checks every
  number reproduces, and audits for overreach. Artifacts:
  [`docs/core_review.md`](docs/core_review.md) (independent, E0–E6/C0–C4) and
  [`docs/extensions_review.md`](docs/extensions_review.md) (self-audit).
- **Planning** (generative, forward-looking) consumes the review findings and
  lays out the roadmap: [`docs/next_steps.md`](docs/next_steps.md).

Do not fold a review into a plan or vice versa, and never soften a review
finding to justify a roadmap track. Both audits live on `main` with the work.

## House rules (apply to every experiment)

- Seed *everything*: thread `default_rng(seed)` explicitly; never use the global
  NumPy RNG (the `perturb_tau` bug — see `docs/core_review.md`).
- Report per-seed spreads, not just means; call out bimodality (the E3 lesson).
- State the substrate/analysis boundary explicitly: what the *dynamics* do vs
  what a *readout/feature* does.
- Keep a caveats section adjacent to every headline.

## Available local skills

The repository defines standardized, agent-agnostic skills under `.claude/skills/` to automate key tasks:

1. **[`publish-viz`](.claude/skills/publish-viz/SKILL.md):**
   - **Purpose:** Publish curated result documents, figures, and animation GIFs to the project's public GitHub Pages site.
   - **How to run:** Cherry-picks files onto the dedicated `deploy-viz-page` branch and validates with `mkdocs build --strict`. Includes a helper (`publish.sh`).
   - **Convention:** every deploy runs in an **isolated git worktree** (gitignored `.publish-worktree/`), never by switching your main checkout's branch — `main` never moves during a publish. The helper sets this up; see the skill's "Deploy from an isolated worktree" section.
2. **[`experiment-review`](.claude/skills/experiment-review/SKILL.md):**
   - **Purpose:** Automate the dual-track review and planning passes (scanning for global RNG usage, scaffolding structured core/extensions reviews, and roadmap planning templates).
   - **How to run:** Run the executable Python CLI helper:
     ```bash
     # Scan codebase for global NumPy/stdlib RNG usage
     python3 .claude/skills/experiment-review/review_helper.py audit-rng

     # Scaffold core or extensions review templates
     python3 .claude/skills/experiment-review/review_helper.py scaffold-review --type [core|extensions] --output [path]

     # Scaffold a fresh roadmap/planning template
     python3 .claude/skills/experiment-review/review_helper.py scaffold-plan --output docs/next_steps.md
     ```


<!-- BEGIN BEADS INTEGRATION v:1 profile:minimal hash:970c3bf2 -->
## Beads Issue Tracker

This project uses **bd (beads)** for issue tracking. Run `bd prime` to see full workflow context and commands.

### Quick Reference

```bash
bd ready              # Find available work
bd show <id>          # View issue details
bd update <id> --claim  # Claim work
bd close <id>         # Complete work
```

### Rules

- Use `bd` for ALL task tracking — do NOT use TodoWrite, TaskCreate, or markdown TODO lists
- Run `bd prime` for detailed command reference and session close protocol
- Use `bd remember` for persistent knowledge — do NOT use MEMORY.md files

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.

## Agent Context Profiles

The managed Beads block is task-tracking guidance, not permission to override repository, user, or orchestrator instructions.

- **Conservative (default)**: Use `bd` for task tracking. Do not run git commits, git pushes, or Dolt remote sync unless explicitly asked. At handoff, report changed files, validation, and suggested next commands.
- **Minimal**: Keep tool instruction files as pointers to `bd prime`; use the same conservative git policy unless active instructions say otherwise.
- **Team-maintainer**: Only when the repository explicitly opts in, agents may close beads, run quality gates, commit, and push as part of session close. A current "do not commit" or "do not push" instruction still wins.

## Session Completion

This protocol applies when ending a Beads implementation workflow. It is subordinate to explicit user, repository, and orchestrator instructions.

1. **File issues for remaining work** - Create beads for anything that needs follow-up
2. **Run quality gates** (if code changed) - Tests, linters, builds
3. **Update issue status** - Close finished work, update in-progress items
4. **Handle git/sync by active profile**:
   ```bash
   # Conservative/minimal/default: report status and proposed commands; wait for approval.
   git status

   # Team-maintainer opt-in only, unless current instructions forbid it:
   git pull --rebase
   bd dolt push
   git push
   git status
   ```
5. **Hand off** - Summarize changes, validation, issue status, and any blocked sync/commit/push step

**Critical rules:**
- Explicit user or orchestrator instructions override this Beads block.
- Do not commit or push without clear authority from the active profile or the current user request.
- If a required sync or push is blocked, stop and report the exact command and error.
<!-- END BEADS INTEGRATION -->

<!-- BEGIN BEADS CODEX SETUP: generated by bd setup codex -->
## Beads Issue Tracker

Use Beads (`bd`) for durable task tracking in repositories that include it. Use the `beads` skill at `.agents/skills/beads/SKILL.md` (project install) or `~/.agents/skills/beads/SKILL.md` (global install) for Beads workflow guidance, then use the `bd` CLI for issue operations.

### Quick Reference

```bash
bd ready                # Find available work
bd show <id>            # View issue details
bd update <id> --claim  # Claim work
bd close <id>           # Complete work
bd prime                # Refresh Beads context
```

### Rules

- Use `bd` for all task tracking; do not create markdown TODO lists.
- Run `bd prime` when Beads context is missing or stale. Codex 0.129.0+ can load Beads context automatically through native hooks; use `/hooks` to inspect or toggle them.
- Keep persistent project memory in Beads via `bd remember`; do not create ad hoc memory files.

**Architecture in one line:** issues live in a local Dolt DB; sync uses `refs/dolt/data` on your git remote; `.beads/issues.jsonl` is a passive export. See https://github.com/gastownhall/beads/blob/main/docs/SYNC_CONCEPTS.md for details and anti-patterns.
<!-- END BEADS CODEX SETUP -->
