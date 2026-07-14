# GEMINI.md

This is the configuration and reference bridge for the Antigravity (Gemini) coding assistant.

The project's full development conventions, house rules, and meta-passes are vendor-neutral and live in:
* [`AGENTS.md`](AGENTS.md) — Main conventions and rules. Please read this file first.
* [`docs/process.md`](docs/process.md) — Detailed description of the Planning and Review meta-passes.

## Project Guidelines Recap for Antigravity

1. **Seed Everything**: Always thread `default_rng(seed)` explicitly. Never use global NumPy RNG to avoid the `perturb_tau` bug.
2. **Report Spreads**: Always report per-seed spreads, not just means. Watch out for and report bimodality.
3. **Boundary Clarity**: Clearly distinguish between what the *substrate/dynamics* do versus what a *readout/feature* does.
4. **Decoupled Passes**: Keep Review and Planning passes completely separate.
5. **Publishing Visualizations**: Follow the steps in [`.claude/skills/publish-viz/SKILL.md`](.claude/skills/publish-viz/SKILL.md) to publish results and animations to the public docs site via the `deploy-viz-page` branch.

@AGENTS.md
