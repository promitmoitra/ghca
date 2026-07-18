# Agent Communication Log

This file is a bridge for cross-agent collaboration (Antigravity/Gemini & Claude Code). Use this log to coordinate development, document workflows, and sync state.

## Update Note for Claude (2026-07-18)

Hello Claude! I have completed two major upgrades to the Greenwich-Hastings Causal Analysis (GHCA) project. All implementation files are currently staged in PR [#30](https://github.com/promitmoitra/ghca/pull/30) on the `feat/experiment-review-and-uv` branch.

Below is a detailed handoff of the new systems:

### 1. 🛠️ New Local Skill: `experiment-review`
We have successfully implemented and packaged a project-specific skill under `.claude/skills/experiment-review/` to automate the project's decoupled "experiment-review" dual-track workflow.

*   **[`SKILL.md`](../.claude/skills/experiment-review/SKILL.md):** Defines standard instructions, checklists, and guardrails for running independent review and planning passes.
*   **[`review_helper.py`](../.claude/skills/experiment-review/review_helper.py):** An executable Python CLI tool containing three subcommands:
    *   `audit-rng`: A static analyzer that scans the repository for unseeded global RNG calls (preventing the `perturb_tau` bug).
        > [!NOTE]
        > The static scanner has been refined with negative lookbehinds `(?<!np\.)` to prevent false positives on valid local generators initialized via `np.random.default_rng(seed)`.
    *   `scaffold-review --type [core|extensions] --output [path]`: Generates structured markdown review files with standard checklists and tables under `docs/`.
    *   `scaffold-plan --output [path]`: Generates a structured roadmap template pre-populated with active tensions and tracks.

### 2. ⚡ Modern Package Management with `uv`
The repository has been successfully transitioned to a `uv` managed system to guarantee reproducible environments across runs.
*   **[`pyproject.toml`](../pyproject.toml):** Configured as a Python application with minimal Python version set to `>=3.13`.
*   **Locked Dependencies:** Added the official dependencies listed in `README.md` and local skills (`numpy`, `matplotlib`, `scipy`, `networkx`, `scikit-learn`, `pillow`, `mkdocs-material`) to the dependencies block and generated the locked dependency graph in [`uv.lock`](../uv.lock).
*   **Execution Pattern:** The new review helper script can be cleanly invoked via:
    ```bash
    uv run .claude/skills/experiment-review/review_helper.py audit-rng
    ```

### 3. 🎯 Discoverability
*   **[`AGENTS.md`](../AGENTS.md):** Updated the main project guide to include a new **`## Available local skills`** section. Both `publish-viz` and `experiment-review` are now linked and documented directly so any arriving agent can easily discover and execute them.

---

### 🚀 Next Steps / Handoff
*   The active Pull Request [#30](https://github.com/promitmoitra/ghca/pull/30) contains the complete implementation. Please review and merge it.
*   When conducting future Review or Planning passes, remember to run `uv run .claude/skills/experiment-review/review_helper.py audit-rng` first to catch any unseeded random calls in development.
