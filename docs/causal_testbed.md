# `ghca_testbed` — a synthetic-SCM benchmark for spike–wave causality

*Usage doc for the Track-4b package (spec [`causal_testbed_spec.md`](causal_testbed_spec.md),
plan [`causal_testbed_plan.md`](causal_testbed_plan.md)). Packages the C-series
(C0–C7) as a reusable benchmark: a substrate where the constitution `W = f(S)` is
explicit, the three `do`-operators are available, the causal graph is **known**,
and every task ships a ground-truth answer. Others can score their
epiphenomenality / causal-discovery / causal-emergence methods against it —
impossible in vivo, where the graph is unknown (Causal Hierarchy Theorem).*

## Install

```
pip install numpy networkx scikit-learn   # networkx: Theorem-1 certificate; sklearn: decoders
```

## Quick start

```python
import ghca_testbed as tb

tb.reproduce()                 # print ground truth + the reference results
print(tb.score())              # score the built-in reference impls over the registry
```

`score()` runs the **executable** C1 scenarios in-process and reports match-vs-
ground-truth; the substrate scenarios (C2–C7) are *carried* (ground truth as data)
and reproduced by running their C-script.

## What's executable now vs carried

- **C1 (T1) — executable & validated in-process.** The six canonical spike–wave
  SCM graphs run end-to-end (`build_scm` → `epiphenomenality_test` +
  `theorem1_certificate`); `test_ctestbed.py` asserts all six verdicts and effects.
- **C2–C7 (T2–T5) — carried.** These need the GH substrate and minutes of compute;
  their ground truth lives in `REGISTRY` (verified against the C-series results
  docs) and full reproduction runs `experiments/cN_*.py`. Wiring them to execute
  in-process through `Testbed` subclasses is the spec's deferred 4b.0d step.

## API surface

- **Operators / substrate** (`tb.do_S`, `tb.do_W`, `tb.do_theta`, `tb.wave_coherence`,
  `tb.wave_active_fraction`, `tb.make_observation_mask`, `tb.read_spikes`,
  `tb.Network`, topology builders) — re-exported single-source from `ghca_causal.py`
  / `ghca_net.py`.
- **SCM** (`tb.build_scm(name, n, seed)`, `tb.SCM_GRAPHS`, `tb.EXPECTED_VERDICT`) —
  the executable C1 graphs.
- **Certificates** (`tb.epiphenomenality_test`, `tb.theorem1_certificate`,
  `tb.strata_means`, `tb.observational_association`) — method-agnostic, take
  **samples** not a live net.
- **Metrics** (`tb.fat_hand_band`, `tb.macro_sufficiency`, `tb.outcome_matrix`,
  `tb.mediation`, `tb.binarize_median`).
- **Registry** (`tb.REGISTRY`, `tb.export_ground_truth`) — 12 `Scenario`s (6 C1 +
  6 substrate); export writes `result/testbed/{scenarios,ground_truth}.json`.
- **Harness** (`tb.score(method_fn, scenarios, tasks)`, `tb.reproduce(mode)`).

### The one load-bearing subtlety

`epiphenomenality_test` scores each interventional stratum mean against the
**observational** stratum mean (`samples_do` is a dict `{w: samples}` so the
contrast cannot degenerate to `do(1)`-vs-`do(0)`). A `do(1)`-vs-`do(0)` test reads
the **front-door** graph as epiphenomenal (effect ≈ 0) when it is causal
(effect ≈ 0.26): intervening severs the latent confounder `C→W` and shifts
`P(B|M; do(W))` off the observational baseline without the effect varying with the
set value. `test_frontdoor_needs_interventional_vs_observational` guards this.

## Score your own method (≈20 lines)

```python
import ghca_testbed as tb
from ghca_testbed.scm import build_scm

def my_method(scenario):
    """Return measured quantities keyed like scenario.ground_truth.scalars."""
    if scenario.kind != "scm":
        return {}                      # skip substrate scenarios for now
    name = scenario.name[len("c1_"):]
    obs, do, edges, S = build_scm(name, n=100_000, seed=0)
    # ... your epiphenomenality estimator here; must NOT do do(1)-vs-do(0) ...
    epi = tb.epiphenomenality_test(obs, do, cond=tuple(S))
    return {"verdict": epi.verdict, "effect": epi.effect}

report = tb.score(my_method, tasks=("T1",))
print(report)                          # per-scenario match-vs-ground-truth
assert report.ok
```

## Ground-truth table

Verdicts and C1 effects are exact-portable; substrate numbers (T2–T5) are checked
by **ordering + tolerance**, not exact decimals (small-`n` spiral panels are
±0.02–0.05). Full detail in `result/testbed/ground_truth.json`.

| task | scenario | quantity | value |
|------|----------|----------|-------|
| T1 | c1_fork | effect / assoc / verdict | 0.001 / 0.004 / epiphenomenal |
| T1 | c1_confounded | effect / assoc / verdict | 0.001 / 0.642 / epiphenomenal |
| T1 | c1_direct | effect / assoc / verdict | 0.400 / 0.799 / causal |
| T1 | c1_mediator_unobs | effect / assoc / verdict | 0.322 / 0.640 / causal |
| T1 | c1_mediator_obs | effect / assoc / verdict | 0.008 / 0.001 / epiphenomenal |
| T1 | c1_frontdoor | effect / assoc / verdict | 0.258 / 0.634 / causal (interventional-vs-obs only) |
| T2 | c2_cnet | fat-hand band, labeled vs collective | 33.1σ vs 0.24σ |
| T3 | c3_cnet_e3 | do(W)_lab / do(W)_coll / do(θ) | 33.1σ / 0.24σ / 0.014σ; E3 latency slope 1.00 |
| T4 | c4_e3_cnet | outcome matrix do(τ)/do(route); macro-suff | [1.00,0.00]/[0.06,1.00]; 1.03 vs 0.11 |
| T5 | c5_spiral | do(χ) band center/tracked/global; routing | 6.2/1.0/2.6σ; 0.55 vs 0.78 |
| T5 | c6_spiral | do(θ_χ) band; switching intact→ablated | 0.0σ; 0.85→0.52 (single-rule 0.90/0.89) |
| T5 | c7_spiral | outcome matrix do(χ)/do(route); screening | [1.00,0.11]/[0.57,1.00]; O_rule tracks injected χ |

## Reproduce

```
python -m ghca_testbed              # reproduce report (C1 executed + substrate carried)
python -m ghca_testbed regression   # nonzero exit on any scored failure
pytest test_ctestbed.py             # the regression suite
```
