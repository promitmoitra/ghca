# Development Workflow — GHCA

## Operational Rules & TDD Standard
1. **Spec-First Planning:** Every feature or experiment track must have an approved `spec.md` and `plan.md` in `conductor/tracks/<track_id>/` before code implementation begins.
2. **Test-Driven Development (TDD):**
   - Write or update unit tests (`tests/test_*.py`) before or alongside substrate modifications.
   - Run `pytest` to confirm tests pass cleanly.
3. **RNG Compliance Verification:**
   - Execute `python3 .claude/skills/experiment-review/review_helper.py audit-rng` before finalizing any PR or experiment module.
4. **Data Persistence & Reproducibility:**
   - Save experimental data to `result/<track>/<experiment>.npz`.
   - Ensure `reproduce_all.py` includes new experiments if they form core baseline capabilities.
5. **Documentation & Review:**
   - Document results in `docs/<track>_results.md` including per-seed spreads, substrate vs readout boundary, and caveats.
