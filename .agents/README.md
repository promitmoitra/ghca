# Project-scoped Antigravity (.agents) Configuration

This directory contains project-specific configurations, rules, and custom skills for the Google Antigravity (AGY) agentic development platform.

## Standard Layout
- `skills/`: Custom, project-scoped skills.
- `mcp_config.json`: Workspace-specific Model Context Protocol (MCP) servers.

## Project Rules & Guidelines (Greenberg-Hastings study)
When executing or running tasks in this repository, always ensure the following:
1. **Seed Everything:** Pass `default_rng(seed)` explicitly to avoid the global NumPy RNG.
2. **Report Spreads:** Report per-seed spreads and check for bimodality.
3. **Boundary Clarity:** Always distinguish between substrate dynamics and readout features.
4. **Decoupled Passes:** Keep adversarial Review and Planning passes completely separate.
