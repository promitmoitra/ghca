<!--
Review / audit pass — adversarial, backward-looking.
See docs/process.md: a review is DECOUPLED from planning. Do NOT fold a roadmap
into a review, and never soften a finding to justify future work.
-->

## Scope & date
<!-- Which experiments/series audited; audit date; branch/refs re-run. -->

## Bottom line
<!-- Fabrication? Overreach in altitude? Reproducibility? One honest paragraph. -->

## Reproduction
- [ ] Re-ran the covered experiments from scratch; regenerated output compared to committed `result/**.npz`
- [ ] Pulled per-seed arrays and checked for bimodality behind a clean mean (the E3 lesson)
- [ ] Scanned for global/unseeded RNG use (the `perturb_tau` lesson)

## Citations
- [ ] Every load-bearing citation resolves and says what the docs claim
- [ ] Unverifiable references hedged as "reported" / repointed to a verified predecessor

## Findings
<!-- Each: what / severity / evidence (file:line, array, number). -->

## Not verified
<!-- Be explicit about residual items you did not clear. -->
