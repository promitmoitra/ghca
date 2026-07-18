---
name: experiment-review
description: >-
  Run the project's decoupled "experiment-review" dual-track workflow. Perform skeptical backward-looking audits (re-running experiments, scanning for unseeded RNGs, checking per-seed bimodality, writing core/extension reviews) followed by generative forward-looking planning passes (mapping roadmap directions against open tensions).
---

# Experiment-Review Dual-Track Workflow

## Overview
This skill implements the Greenberg–Hastings Causal Analysis (GHCA) project's decoupled, sequential **Review and Planning** meta-passes. It enforces maximum scientific integrity and skeptical objectivity during audits, ensuring every headline is honest and reproducible, and uses these audited findings as the inputs for generative roadmap planning.

## Dependencies
- **`publish-viz` (`.claude/skills/publish-viz/SKILL.md`):** Used to publish the finalized review (`core_review.md` or `extensions_review.md`) and planning (`next_steps.md`) documents to the public MkDocs pages.

## Quick Start
1. **Audit Global RNGs** (ensuring all experiments thread `default_rng(seed)` explicitly):
   ```bash
   uv run .claude/skills/experiment-review/review_helper.py audit-rng
   ```
2. **Scaffold a Core Series Review** (E0–E6, C0–C4):
   ```bash
   uv run .claude/skills/experiment-review/review_helper.py scaffold-review --type core --output docs/core_review.md
   ```
3. **Scaffold an Extensions Review** (E7, C5–C7, E8.x, E9.x):
   ```bash
   uv run .claude/skills/experiment-review/review_helper.py scaffold-review --type extensions --output docs/extensions_review.md
   ```
4. **Scaffold a Roadmap/Plan**:
   ```bash
   uv run .claude/skills/experiment-review/review_helper.py scaffold-plan --output docs/next_steps.md
   ```

## Utility Scripts
The skill provides a Python CLI helper `.claude/skills/experiment-review/review_helper.py` to automate structural checks and document generation:

### Subcommands

#### `audit-rng`
Scans all core files (`ghca_*.py`) and experiment definitions (`experiments/*.py`) for forbidden global NumPy or stdlib random calls (e.g. `np.random.choice` or `random.randint`). All code must thread explicit Generator objects.
- **Example:**
  ```bash
  uv run .claude/skills/experiment-review/review_helper.py audit-rng
  ```

#### `scaffold-review`
Generates a structured review document template under `docs/` containing required checklists, per-series verdicts table, and citation validation blocks.
- **Arguments:**
  - `--type {core,extensions}` (Required)
  - `--output <path>` (Required)
- **Example:**
  ```bash
  uv run .claude/skills/experiment-review/review_helper.py scaffold-review --type extensions --output docs/extensions_review.md
  ```

#### `scaffold-plan`
Generates a structured roadmap template pre-populated with sections for active tensions, status summaries, candidate tracks, and scoring.
- **Arguments:**
  - `--output <path>` (Required)
- **Example:**
  ```bash
  uv run .claude/skills/experiment-review/review_helper.py scaffold-plan --output docs/next_steps.md
  ```

## Decoupled Workflow Procedure

The workflow separates **Review** from **Planning** so that next steps are always steered by audited truth rather than wishful thinking.

### Phase A: The Review Pass (Adversarial, Skeptical)
1. **RNG Scan:** Run `audit-rng` to ensure no unseeded global RNG bugs have snuck into the codebase.
2. **Scaffold Review:** Run `scaffold-review` for either the core or extension series.
3. **Re-Run Experiments:** Manually re-run relevant experiment scripts (e.g., `uv run experiments/e1_conditioning.py`) to verify that stdout and `.npz` output arrays reproduce bit-identically against the written documentation.
4. **Seed/Bimodality Check:** Inspect raw per-seed values for key headlines to check for bimodal distributions (e.g. E3 composition study). Never allow a bimodal distribution to be characterized as a simple "clean mean" without flagging variance.
5. **Boundary Audit:** Draw a strict line between what the *substrate/dynamics* do versus what a *readout/feature* does (the Hebbian self-organisation vs hard-wired readout boundary).
6. **Citation Audit:** Verify every external reference is real and accurately summarized (e.g., Jalaldoust et al. SCM/epiphenomenality paper).
7. **Commit & Deploy:** Commit the finished audit directly to the `main` branch. Use the `publish-viz` skill to cherry-pick and publish it to the public MkDocs pages.

### Phase B: The Planning Pass (Generative, Roadmap)
1. **Consume Audits:** Read `docs/core_review.md` and `docs/extensions_review.md`. Extract current open tensions.
2. **Scaffold Plan:** Run `scaffold-plan --output docs/next_steps.md`.
3. **Draft Tracks:** Formulate candidate roadmap tracks as a menu of options (what it is, why it matters, effort/risk, connection to current code).
4. **Score against Tensions:** Score each track according to which open tension it retires.
5. **Commit & Deploy:** Commit the completed `next_steps.md` plan to `main` and use `publish-viz` to deploy it.

## Common Mistakes
1. **Combining Review and Plan:** Never write next steps or candidate roadmap items in a Review document, and never soften or edit past review findings inside a Planning document to make a track look better.
2. **Hiding Variance under Means:** Reporting high-variance, bimodal per-seed sets (like 2 seeds at ceiling and 2 at chance) as a stable mean (e.g., "accuracy ~0.77") is a critical failure. Always report per-seed spreads and bimodal splits.
3. **Confusing Emergence with Readouts:** Crediting the excitable substrate with cognitive capabilities when those capabilities actually reside in engineered features or hand-built topological winding readouts. Be precise about what is emergent vs. what is engineered.
