# 4b Plan — Implementing the Spike–Wave Causal Testbed

*Implementation plan for Track 4b, spec [`causal_testbed_spec.md`](causal_testbed_spec.md)
(approved). This turns the component designs into an ordered, file-level build a
single engineer can execute step by step. It follows the spec's build sequence:
**4b.0** extract the API behind a bit-identical gate → **4b.1** scenario registry
+ export → **4b.2** test/metric library + score harness → **4b.3** docs + minimal
example → **4b.4** regression suite. No new science; the only risk is a refactor
silently moving a C-series number, and the whole shape of the plan is built to
catch that.*

## Preamble — the one load-bearing principle

The C-series numbers come from RNG-consuming code in `ghca_causal.py`,
`ghca_net.py`, `e3_timed_response.py`, `e7_spiral_option.py`, and the eight C
scripts. **The testbed relocates that code; it never rewrites it.** Every shared
unit is lifted *verbatim* (copy-then-delete-original, not re-express), and the
`Testbed` object model is a thin OO adapter that delegates to those functions in
the same call order. Bit-identity — the 4b.0 gate — holds by construction, not by
luck. Any "cleanup" that reorders `np.random` draws breaks it.

**Deliberate deviation from the spec.** The spec says "a single module
`ghca_testbed.py`." At ~1500 lines spanning two substrates (scalar GH population +
2-D no-flux spiral lattice) plus a pure-SCM sibling, that is unwieldy. This plan
ships a **package `ghca_testbed/`** whose `__init__.py` re-exports the entire
surface, so `import ghca_testbed` and every name in the spec resolve unchanged.
Reversible; flagged here so a reviewer isn't surprised. (§4b.2/§4b.4 signatures
that name `ghca_testbed.py` map onto the package facade.)

---

## 4b.0 — Extract the API behind a bit-identical gate *(gate)*

Move shared machinery into the package; rewrite C0–C7 as thin callers. **Do this
step before any registry/harness work** — the later steps are gated on it.

### Files to create

- **`ghca_testbed/__init__.py`** — public facade. Re-exports `Scenario`,
  `GraphSpec`/`GroundTruth`, `Testbed` (+ subclasses), the operator/wave functions,
  the lifted factories, the certificate/metric library, `REGISTRY`, `score`. Keeps
  `import ghca_testbed` behaving as the spec's single module.
- **`ghca_testbed/operators.py`** — single-source intervention/observation layer.
  Re-exports (does **not** reimplement) `do_S`, `do_W`, `do_theta`,
  `wave_coherence`, `wave_active_fraction`, `WAVE_FNS`, `make_observation_mask`,
  `read_spikes` from `ghca_causal.py`. The dead exports (`do_S`,
  `make_observation_mask`, `read_spikes`) stay alive per spec 4b.0.
- **`ghca_testbed/substrate.py`** — re-exports `Network`, `ring`, `lattice2d`,
  `rgg`, `smallworld`, `count_phase_singularities` from `ghca_net.py`; and
  `local_winding`, `signed_charge`, `seed_spiral`, `probe_chirality` from
  `e7_spiral_option.py`. Keeps `e3_timed_response` and `e7_learning` importable
  (C3/C4 need `build`/`ACT`/`TAU0`/roles; C5–C7 need `make_router`/`make_spiral`/
  `spiral_decode`/`trial`/`run_switching`/`run_single_rule`).
- **`ghca_testbed/factories.py`** — the lifted duplication anchors, copied
  **verbatim** to preserve RNG draw order: `make_cnet`, `behaviour_from_readout`,
  `collect_trials`, `decode`, `build_e3_net`, `make_spiral_net`, `realize`,
  `SPIRAL_READERS`, `decode_chirality`.
- **`ghca_testbed/scenario.py`** — `Scenario` + `GraphSpec` dataclasses and the
  four scenario builders (`make_cnet_scenario`, `make_e3_scenario`,
  `make_spiral_scenario`, `make_scm_scenario`).
- **`ghca_testbed/testbed.py`** — `Testbed` ABC + `Testbed.for_scenario` dispatch
  and the four concrete subclasses. Thin adapter; guards the two operator gotchas
  (coupling `W0` snapshot; `do_W` wave≠`active_fraction` → `NotImplementedError`).

### Files to edit (thin callers — machinery deleted, imports added)

`experiments/c0_instrumentation.py`, `c1_graph_certificates.py`,
`c2_fat_handed.py`, `c3_do_theta.py`, `c4_outcome_relativity.py`,
`c5_do_chirality.py`, `c6_do_theta_chi.py`, `c7_outcome_relativity.py`. Each drops
its inlined `build_cnet`/`collect`/`decode`/`make_spiral`/`realize`/`READERS` and
imports the lifted version; each drops its ad-hoc `sys.path.insert`. **What each
caller must keep verbatim** (these are the silent-regression traps):

- **C0** keeps its own partial-obs column selection (`default_rng(0)`,
  `rng.choice(N_H, k)`); routing it through `read_spikes` would move the RNG.
- **C1** keeps module-level `rng = default_rng(0)`, `N = 200_000`, and `GRAPHS`
  iteration order (one global stream consumed sequentially across all six graphs).
- **C2** `build → make_cnet(p_s=0.01)`; prop rng `default_rng(11)`; baseline via
  `net.seed_random(0.1, 0.1)` 800 steps; `N_TRIALS=250` stays unused.
- **C3** per-`(seed,t)` `default_rng(1000*s + int(t))`; net rebuilt each call
  (so `do_theta` must **not** double-reset — see gotcha); loads
  `result/c2/c2_data.npz` unchanged.
- **C4** `E3Testbed.run` resets `net.W[:] = W0` before each `do`; macro-sufficiency
  reuses `collect_trials` with the readout rng `default_rng(seed+1)` **continued**
  (historical quirk — do not hand it a fresh seed); identity = `argmax(tot)` over
  integrated activity, not first-spike.
- **C5** `default_rng(1000*pi + k)`, eval `default_rng(9000+k)`; `seed_spiral`
  collapses into `realize('centered'/'noisy')`.
- **C6** `default_rng(4000+k)`; keeps the `max==min → exactly 0.0` short-circuit;
  loads `result/c5/c5_data.npz` unchanged.
- **C7** `default_rng(s+31)`; retrains its own router (stands alone).

### Key signatures

```
def make_cnet(seed=0, N_H=120, act=2, pas=8, theta=1.0, p_s=0.02, k=6, beta=0.15) -> CNet
def behaviour_from_readout(active, w_vec, ref=None) -> int | np.ndarray   # median split
def collect_trials(cnet, rng, n, settle=40, seed_frac=0.1, wave_keys=('active_fraction',)) -> dict
def decode(X, y, C=1.0, max_iter=2000, cv=5) -> float                     # 5-fold CV logreg acc
def build_e3_net(seed=0, w_hm=0.15, w_gate=4.0, theta_m=5.0) -> tuple[Network, dict, np.ndarray]  # (net, roles, W0)
def make_spiral_net(seed=0) -> Network
def realize(net, chir, policy, rng, d=0, pitch=0.0, jitter=0.0) -> tuple[int, int]
SPIRAL_READERS: dict[str, Callable[[Network], float]]  # center / tracked / global / phase_probe
def decode_chirality(net, reader, steps=4) -> int

class Testbed(abc.ABC):
    scenario: Scenario
    @classmethod
    def for_scenario(cls, scenario, seed=None) -> 'Testbed'   # dispatch on scenario.kind
    def do_S(self, nodes, values) -> None
    def do_theta(self, *, tau=None, coupling_scale=None, theta=None, nodes=None, reset=True) -> None
    def observe(self) -> np.ndarray
    def wave(self, name=None) -> float | int
    def behaviour(self) -> int
    def rollout(self, interventions=None, n=1, seed=None, **kw) -> dict   # accepts Generator OR int seed

class CNetTestbed(Testbed):
    def do_W(self, target, policy='random', wave='active_fraction', nodes=None, rng=None) -> None
class E3Testbed(Testbed):
    def do_theta(self, *, tau_gate=None, g_route=None, reset=True) -> None   # resets W[:]=W0 then scales
    def run(self, n=40) -> tuple[float, float]   # (mean first-fire latency, P(ch0 wins by argmax(tot)))
class SpiralTestbed(Testbed):
    def do_chi(self, chir, policy='centered', d=0, pitch=0.0, jitter=0.0, rng=None) -> tuple[int, int]
    def do_theta_chi(self, chir, rng=None) -> tuple[int, int]   # nucleate: realize('centered')
    def wave(self, reader='center', steps=4) -> int
class SCMTestbed(Testbed):
    def rollout(self, do_w=None, n=None, seed=0) -> dict   # do_w=None is observational
```

**RNG threading rule (mandatory):** `rollout`/`collect_trials` accept a
`Generator` **or** an int seed, never seed-only. C4's macro-sufficiency reuses the
readout `default_rng(seed+1)` *after* it drew `w_coll`/`w_lab` (continuing that
stream); a seed-only API cannot reproduce it.

**Per-caller seed map (freeze this):** c-net `net=seed`, `readout=seed+1`;
C0 trials `seed+7`, cols `default_rng(0)`; C1 module `default_rng(0)`, `N=200000`,
fixed `GRAPHS` order; C2 `build=0`, prop `default_rng(11)`, baseline `net.rng` 800
steps; C3 `default_rng(1000*s+int(t))`; C4 e3 `seed=0`, macro `seed+1`;
C5 `default_rng(1000*pi+k)`, eval `default_rng(9000+k)`; C6 `default_rng(4000+k)`;
C7 `default_rng(s+31)`.

### Gate / acceptance

Freeze **reference artifacts first**: from the pre-refactor commit, on the pinned
CI image, run each `cX.py`, byte-copy `result/cX/cX_data.npz` →
`result/testbed/refs/cX_data.npz` and save `normalize_stdout(stdout)` →
`result/testbed/refs/cX.stdout.txt`; commit them. Then, after the refactor, for
`c0..c7` in dependency order (C2 before C3, C5 before C6), run the thin caller into
an isolated tmp outdir and require **both**: `assert_bit_identical(new_npz, ref) ==
[]` (identical keys, dtype, shape, `np.array_equal(..., equal_nan=True)` — exact,
not `allclose`) **and** `normalize_stdout(new) == ref_stdout`. Any diff fails the
gate. Single-threaded BLAS + pinned versions are prerequisites (see Risks). This
step alone converts eight scripts into one library and is worth shipping on its own.

---

## 4b.1 — Scenario registry + export

Encode the ground-truth scenarios and dump them as language-agnostic data.
`Scenario` **owns** its ground truth (data, not code that re-derives it) so the
harness and the C scripts read the answer, never recompute it.

### Files to create

- **`ghca_testbed/registry.py`** — `REGISTRY: dict[str, Scenario]` for the six C1
  graphs + c-net (C2) + E3 (C3/C4) + spiral (C5/C6/C7); `export_ground_truth()`.
  Preserves `GRAPHS` iteration order and records cross-experiment run-order
  constraints (C2→C3, C5→C6) as manifest-level `depends_on`.
- **`result/testbed/scenarios.json`** — manifest: one record per scenario
  (name/track/task/kind/seed, build_kwargs, roles as int index lists, C1 graph
  edges as source/target pairs, conditioning set `S`, wave specs, scalar answers).
- **`result/testbed/ground_truth.json`** — flat scoring table:
  `{task, scenario, quantity, value, tolerance, assert_kind}` rows plus the
  qualitative ordering assertions. Mirrors the §3 T1–T5 table; single source of
  expected values for both the harness and the regression suite.
- **`result/testbed/<name>.npz`** — per-scenario arrays that don't fit JSON scalars
  (adjacency matrix + node-order, role masks, `w_coll`/`w_lab`, expected effect
  matrices). e.g. `c1_frontdoor.npz`, `c2_constituted.npz`, `c4_e3net.npz`,
  `c5_spiral.npz`, `c7_spiral.npz`.
- **`result/testbed/registry.npz`** — optional single-file bundle of all arrays
  (name-prefixed keys) for consumers preferring one download.

### Key signatures

```
@dataclass(frozen=True)
class GraphSpec:
    edges: tuple[tuple[str, str], ...]; S: tuple[str, ...]; expected_verdict: str  # 'epiphenomenal'|'causal'

@dataclass(frozen=True)
class GroundTruth:
    task: str; verdict: str | None; scalars: dict[str, float]; matrices: dict[str, list]
    tol: dict[str, float]; ordering: list[str]; exact: bool; notes: str

@dataclass(frozen=True)
class Scenario:
    name: str; kind: str; seed: int; params: Mapping[str, Any]; build: Callable[..., Any]
    waves: Mapping[str, str]; roles: Mapping[str, Any]; graph: GraphSpec | None
    ground_truth: GroundTruth   # kind in {'cnet','e3','spiral','scm'}

def export_ground_truth(path='result/testbed') -> None   # scenarios.json + ground_truth.json + *.npz; idempotent
def scenario_to_json(sc: Scenario) -> dict
def adjacency_matrix(edges, node_order) -> tuple[np.ndarray, list[str]]
def load_registry_json(path='result/testbed/scenarios.json') -> list[dict]
```

- `exact=True` only for C1 (categorical / `N=200000` MC / analytic). Everything
  else `exact=False`: `ground_truth.ordering` carries the robust targets
  (`center>global>tracked`, `do(theta)<<do(W)`, diagonal-dominance).
- `GroundTruth.notes` for each C1 graph records the interventional-vs-observational
  contrast and flags that `do1`-vs-`do0` mis-scores front-door — encode the
  **as-built** test (0.03 threshold, no shuffle null) so 4b.2 cannot silently move
  T1. `frontdoor_certificate` is marked PLANNED/optional.

### Acceptance

`export_ground_truth()` is deterministic and idempotent (seeds threaded); a
re-export byte-matches. Each C1 `Scenario` ships edges + `S` + both interventional
and observational sample paths (so the harness can *detect* a `do1`-vs-`do0`
estimator, not just the verdict). Spiral scenarios wire to
`e7_spiral_option.local_winding`/`signed_charge` (signed, no-flux), **not**
`count_phase_singularities` (unsigned, periodic — not interchangeable).

---

## 4b.2 — Test/metric library + score harness

Package Definition-1, the certificates, and the metrics as method-agnostic,
ground-truth-free functions that take **samples, not a live net** (the
substrate/analysis boundary). Then the harness.

### Files to create

- **`ghca_testbed/certificates.py`** — `epiphenomenality_test`,
  `theorem1_certificate`, `frontdoor_certificate` (new per Prop 2, optional),
  `strata_means`, `observational_association`. Extracted from
  `c1_graph_certificates.py:main` (which currently inlines the Def-1 test).
- **`ghca_testbed/metrics.py`** — `fat_hand_band`, `achievable_band_t0` (kept
  **separate** — analytic vs empirical band are different formulas),
  `macro_sufficiency`, `outcome_matrix`, `mediation`, `binarize_median`.
- **`ghca_testbed/harness.py`** — `score` + `ScoreRow`/`ScoreReport` and the T1–T5
  pass-condition rubric encoded as per-task comparators.

### Key signatures

```
def epiphenomenality_test(samples_obs, samples_do, *, target='B', cond=('S',),
                          distance='mean_gap', min_stratum=200, n_shuffle=0,
                          threshold=0.03, rng=None) -> EpiResult
    # effect = max over set-values z in samples_do AND strata c of
    #          dist(P(B|c; do(Z=z)), P_obs(B|c))   -- interventional-vs-OBSERVATIONAL, never do1-vs-do0
    # verdict 'epiphenomenal' iff effect < threshold
def theorem1_certificate(graph, *, W='W', B='B', S=()) -> CertResult   # d-separation in G_{W(S)}
def frontdoor_certificate(graph, *, W='W', B='B', S=(), M=None) -> FrontdoorResult   # NEW, optional

def fat_hand_band(effect_by_policy, *, baseline_std=None) -> float
    # baseline_std given -> C2/C3 observational-sigma units; None -> C5 pooled within-policy std
    # returns exactly 0.0 when all policy means equal (C6 short-circuit)
def achievable_band_t0(w_vec, w, N) -> tuple[float, float, float]   # analytic: sort weights, sum k smallest/largest
def macro_sufficiency(B, W, S, *, decoder=None, chance=0.5, cv=5) -> float   # I(B;W)/I(B;S), may exceed 1, NOT clipped
def outcome_matrix(handles, outcomes, *, normalize='column', diag_tol=0.2) -> OutcomeResult
def mediation(theta, W, B, *, min_stratum=1) -> MediationResult   # screening: theta -> W -> B

def score(method_fn, scenarios=REGISTRY, tasks=('T1','T2','T3','T4','T5')) -> ScoreReport
    # method_fn(task, testbed, scenario) -> task answer; method_fn=None runs the reference impls (also 4b.3)
```

- Backward-compat defaults (`distance='mean_gap'`, `n_shuffle=0`,
  `threshold=0.03`) make `epiphenomenality_test` reduce to C1 exactly. The
  distributional-distance + stratified shuffle null layer is off by default so
  headline numbers do not move.
- `epiphenomenality_test` **requires** `samples_obs` and compares each
  interventional stratum mean to the *observational* stratum mean — it structurally
  cannot do `do1`-vs-`do0`. `samples_do` is a dict `{w: samples}` to force the
  correct contrast.
- Metrics take samples so a user can plug in their own estimator against the same
  registry without touching the substrate. Operator gotchas are guarded in
  `Testbed`, not here.

### Acceptance

`score(method_fn=None)` reproduces the full ground-truth table. T1 must assert the
front-door row (`effect ≈ 0.258`, `causal`) and the confounded row (`assoc ≈ 0.642`,
`effect < 0.03`, `epiphenomenal`) — the two contrasts that catch the wrong Def-1
scoring. Numeric tasks score by ordering + threshold + tolerance, never exact
decimals. `outcome_matrix` tested for diagonal-**dominance**, not perfect diagonal
(C7 off-diagonals 0.57/0.11 are honest).

---

## 4b.3 — Docs + minimal example

### Files to create

- **`docs/causal_testbed.md`** — usage doc: the API surface, the T1–T5 table, a
  ~20-line "score your method" example, and the ground-truth table (below) inline.
- **`reproduce`/`benchmark` entry point** in `ghca_testbed/harness.py` (invokable
  as `python -m ghca_testbed reproduce`) — prints, per experiment, three blocks:
  **GROUND TRUTH** (from `ground_truth.json` + committed refs), **MODEL RESULTS**
  (a fresh run, with deltas), **VERDICT** (Tier-2 PASS/FAIL per headline +
  Tier-1 bit-identical OK/DIFF per experiment). It reuses the same `HEADLINES`
  registry and tolerance rules as the regression suite, so "what reproduce prints"
  and "what the test asserts" cannot drift.

### Acceptance

The 20-line example runs end-to-end against `REGISTRY` and prints a `ScoreReport`.
`reproduce(mode='report')` exits 0 (human inspection); `mode='gate'` exits nonzero
on any bit-identical DIFF; `mode='regression'` exits nonzero on any Tier-2 FAIL.

---

## 4b.4 — Regression suite

### Files to create

- **`test_ctestbed.py`** — pytest suite. **Tier-1** (parametrized `c0..c7`):
  bit-identical vs committed refs (npz + normalized stdout). **Tier-2**: portable
  scientific assertions — exact for categorical/analytic/integer results,
  tolerance-banded for sklearn- and small-N-spiral-derived results, ordering/
  threshold for the robust qualitative claims.
- **`conftest.py`** — pins determinism (`OMP/OPENBLAS/MKL_NUM_THREADS=1`,
  `PYTHONHASHSEED=0`) **before** numpy import; session fixture runs all eight
  experiments once in dependency order into isolated tmp outdirs, caching
  `(npz-dict, normalized-stdout)`; loader for committed refs.
- **`result/testbed/refs/`** — the frozen pre-refactor `cX_data.npz` +
  `cX.stdout.txt` (from 4b.0), plus an environment lockfile
  (numpy/scipy/scikit-learn/networkx versions).

### Key signatures

```
@dataclass
class Headline:
    label: str; source: str; key: str; derive: Callable | None
    expected: Any; check: str; tol: float   # check in {exact, atol, lt, gt, side_of_0.03, order, array_equal, categorical}
HEADLINES: dict[str, list[Headline]]         # single-source registry, shared with reproduce()
def extract_headlines(exp, data) -> dict[str, Any]
def run_experiment(exp, outdir) -> tuple[dict, str]
def assert_bit_identical(new, ref) -> list[str]   # keys where new != ref by exact match; [] == pass
def normalize_stdout(s) -> str                     # drop 'wrote <path>' lines + abs-path prefixes; keep every number
def committed_ref(exp) -> tuple[dict, str]
```

Representative Tier-2 tests: `test_c1_frontdoor_scored_causal` (guards the
interventional-vs-observational contrast against a `do1`-vs-`do0` regression),
`test_c2_fat_hand_ratio` (`band_l/band_c > 100`), `test_c4_diagonal_and_sufficiency`
(`effn[0,0]==1` exact, `suff_coll/suff_lab > 3`), `test_c5_band_order_and_failure_modes`
(`center>global>tracked`, three distinct failure axes),
`test_c6_wellposed_and_necessity` (`theta_chi_band` all zeros; switching collapses
to chance, single-rule spared), `test_c7_diagonal_dominant_and_screening`
(`Mn[1,0] in [0.3,0.75]` — do not "fix" the honest off-diagonal).

### Acceptance

Full suite green on the pinned image. Tier-1 = "did the refactor change anything";
Tier-2 = "is the science still there." Figures are regenerated but **not** diffed
(non-deterministic matplotlib metadata).

---

## Ground-truth table (reviewer's checklist)

Exact (`exact=True`) unless marked *(ordering/tolerance)*. Full spiral detail lives
in `ground_truth.json`.

| task | scenario | quantity | value |
|------|----------|----------|-------|
| T1 | C1 fork | do(W) effect / obs assoc / verdict | 0.001 / 0.004 / epiphenomenal |
| T1 | C1 confounded | effect / assoc / verdict | 0.001 / 0.642 / epiphenomenal (correlation≠causation) |
| T1 | C1 direct | effect / assoc / verdict | 0.400 / 0.799 / causal |
| T1 | C1 mediator-unobs | effect / assoc / verdict | 0.322 / 0.640 / causal |
| T1 | C1 mediator-obs | effect / assoc / verdict | 0.008 / 0.001 / epiphenomenal |
| T1 | C1 front-door | effect / assoc / verdict | 0.258 / 0.634 / causal (only under interventional-vs-obs) |
| T1 | C1 overall | all_ok (cert == intervention == expected, 6/6) | True (N=200000, seed 0, threshold 0.03) |
| T2 | C2 c-net | achievable band t=0, collective vs labeled *(ratio)* | 0.24σ vs 33.1σ (~140×) |
| T2 | C2 c-net | random-realization spread t=0 | collective 0.015σ vs labeled 2.09σ |
| T2 | C2 c-net | propagation washout | labeled large@t0 → pinned by t≈8 |
| T3 | C3 c-net | intervention ambiguity, do(W)_lab / do(W)_coll / do(θ) *(ordering)* | 33.1σ / 0.24σ / 0.014σ |
| T3 | C3 c-net | do(θ) across-seed reproducibility std | 0.15 |
| T3 | C3 E3 | latency-vs-τ_gate fit slope | 1.00 (tol 0.02) |
| T4 | C4 E3 | outcome matrix (normalized) do(τ_gate) / do(g_route) | [1.00, 0.00] / [0.06, 1.00] (diagonal) |
| T4 | C4 E3 | raw do(τ) latency range | 18 steps |
| T4 | C4 c-net | macro-sufficiency I(B;W)/I(B;S), coll vs lab *(ordering)* | 1.03 vs 0.11 |
| T5 | C5 spiral | do(χ) bands center / tracked / global *(ordering)* | 6.2σ / 1.0σ / 2.6σ |
| T5 | C5 spiral | center acc: centered/d4/d8/d12/pitch/noisy | 1.00/1.00/0.10/0.00/1.00/0.90 |
| T5 | C5 spiral | tracked acc (same policies) | 1.00/1.00/1.00/1.00/1.00/0.88 |
| T5 | C5 spiral | global acc (same policies) | 1.00/1.00/1.00/1.00/1.00/0.50 |
| T5 | C5 spiral | routing (frozen router, d=8) center vs tracked | 0.55 vs 0.78 |
| T5 | C6 spiral | do(θ_χ) nucleation bands center/tracked/global | 0.0σ / 0.0σ / 0.0σ (exact, by construction) |
| T5 | C6 spiral | necessity: switching intact→ablated | 0.85 → 0.52 (collapses to chance) |
| T5 | C6 spiral | necessity: single-rule intact / ablated | 0.90 / 0.89 (spared) |
| T5 | C7 spiral | outcome matrix (norm) do(χ) / do(g_route) | [1.00, 0.11] / [0.57, 1.00] (diagonal-dominant) |
| T5 | C7 spiral | screening O_rule, seed+ [inj+,inj-] / seed- [inj+,inj-] | [0.86, 0.19] / [0.82, 0.14] (columns ≈+0.84/−0.16) |

Small-N spiral numbers (T5, `n_trials=40`, `n_seeds=3–5`) are ±0.02–0.05: check the
**orderings and failure modes**, not the decimals. Only T1 is bit-exact-portable.

---

## Risks / gotchas

**The C1 interventional-vs-observational subtlety (highest value, highest risk).**
`epiphenomenality_test` scores each interventional stratum mean against the
*observational* stratum mean, maxed over set-values and strata — never `do(W=1)` vs
`do(W=0)`. `do1`-vs-`do0` detects only effects that vary with the set value and is
blind to the confounding-breaking effect that makes front-door causal: conditioning
on the observed mediator `M` blocks `W→M→B` so `B` does not vary with `w`, but
`do(W)` severs the latent `C→W` arrow and shifts `P(B|M; do(W))` off observational.
`do1`-vs-`do0` reads front-door as ≈0 (epiphenomenal); the correct contrast yields
0.258 (causal). The first C1 pass made this mistake. Mitigations: the API forbids
omitting `samples_obs`; `samples_do` is dict-shaped to force the contrast; T1
ground truth asserts the front-door row; the regression suite has a dedicated
`test_c1_frontdoor_scored_causal`. The dict shape mitigates but does not *prevent*
a future refactor re-introducing the bug — hence the explicit test.

**Refactor-changes-a-number (the dominant hazard) and how the gate catches it.**
Any reordering of `np.random` draws during the lift silently moves a C-series
number. Defenses, in order: (1) copy-verbatim rule for `factories.py`
(copy-then-delete, never re-express); (2) `Generator`-not-seed threading in
`rollout`/`collect_trials`; (3) the **4b.0 bit-identical gate** — `np.array_equal`
(not `allclose`) on every npz key plus byte-identical normalized stdout, against
refs frozen from the pre-refactor commit. Because the reference is regenerated on
the same pinned single-threaded-BLAS image, a moved draw produces a different float
and the exact compare fails immediately, pinpointing the experiment. Tier-2 then
guards the science across version bumps.

**`do_theta` reset-default mismatch.** `do_theta(coupling_scale=)` multiplies
`net.W` in place and is cumulative. A single reset policy cannot serve both C3
(rebuilds the net every call — must **not** double-reset) and C4 (shares the net —
**must** reset `W[:]=W0`). Defaults are chosen per subclass: `E3Testbed.run` resets;
C3's caller uses a fresh testbed / `reset=False`. Wrong defaults regress C3/C4 with
no error — guarded by `test_c4_diagonal_and_sufficiency` (diagonal `[1,0]/[0.06,1]`)
and the C3 latency-slope-1.00 assertion.

**Seeding / order is part of the fixture.** C1 consumes one module-level
`default_rng(0)` sequentially across all six graphs, so `GRAPHS` iteration order is
load-bearing — reorder it and every C1 number moves. `p_s` is a fixed regime knob
(0.02 for C0/C3/C4, 0.01 for C2), not a free parameter. C4's readout-rng reuse
(`seed+1` continued into the trial loop) is a documented caller obligation, not
type-enforced — "cleaning it up" to a fresh seed changes macro-sufficiency.
Determinism prereqs (single-threaded BLAS, `PYTHONHASHSEED=0`, pinned versions)
must be set before numpy import or Tier-1 flaps.

**Attractive cleanups to resist.** Merging `achievable_band_t0` (analytic) with
`fat_hand_band` (empirical) corrupts the C2 33σ or C5 6.2σ headline. Unifying the
three winding kernels erases the signed/no-flux vs unsigned/periodic distinction
(`count_phase_singularities` is not interchangeable with `local_winding`/
`signed_charge`). Making `do_theta` idempotent, changing C4 identity to first-spike,
or swapping C6's threshold-raise ablation for node deletion each change a number
legitimately-looking-wrong. Keep them distinct; document each pinned subtlety inline
in `test_ctestbed.py` pointing at the note it guards.

**Cross-machine bit-identity is not guaranteed.** `Network.step` uses `W@active`
(BLAS matmul, reduction-order-sensitive) and sklearn shares that exposure, so a
borderline `inp>=theta` flip can cascade. The Tier-1 gate is defined as
*same-pinned-environment refactor-invariance*, not cross-machine identity;
portability is delivered by the Tier-2 tolerance layer. If CI cannot pin an image,
relax Tier-1 to tight `allclose` (atol 1e-9) on deterministic keys and skip sklearn
keys.

**Scope creep / dead API.** `Scenario` cannot be defined without its ground-truth
field, so 4b.0 already touches `Scenario`; keep registry/harness signatures sketched
but gated behind the 4b.0 pass. `do_S`, `make_observation_mask`, `read_spikes` stay
exported (spec 4b.0) but the suite asserts only their import-ability. C0 keeps its
own column selection — deliberately not reconciled through `read_spikes`.

---

## Ordered checklist

1. **4b.0a** — On the pre-refactor commit + pinned image: freeze
   `result/testbed/refs/{cX_data.npz, cX.stdout.txt}` for `c0..c7`; commit the
   environment lockfile.
2. **4b.0b** — Create `operators.py`, `substrate.py` (re-export only).
3. **4b.0c** — Create `factories.py` by copy-verbatim from the scripts; create
   `scenario.py`, `testbed.py`, `__init__.py`.
4. **4b.0d** — Rewrite `c0..c7` as thin callers; delete the originals of the lifted
   units; drop ad-hoc `sys.path.insert`.
5. **4b.0 GATE** — Run `c0..c7` in dependency order; `assert_bit_identical == []`
   and stdout matches for all eight. **Do not proceed until green.**
6. **4b.1** — Add `GroundTruth`/`GraphSpec` fields, `registry.py`, the four scenario
   builders; implement `export_ground_truth`; write `scenarios.json`,
   `ground_truth.json`, per-scenario + bundle npz.
7. **4b.2** — Create `certificates.py`, `metrics.py`, `harness.py`; verify
   `score(method_fn=None)` reproduces the ground-truth table.
8. **4b.3** — Write `docs/causal_testbed.md` + the 20-line example; add the
   `reproduce`/`benchmark` entry point sharing `HEADLINES`.
9. **4b.4** — Write `conftest.py` + `test_ctestbed.py` (Tier-1 bit-identical +
   Tier-2 scientific); green on the pinned image.

## One-line spine

Freeze the refs → lift the machinery verbatim so C0–C7 reproduce bit-identically
(4b.0) → register the ground-truth scenarios and export them (4b.1) → package the
certificates/metrics + a scoring harness (4b.2) → document with a minimal example
(4b.3) → lock it with a two-tier regression suite (4b.4): a synthetic spike–wave
SCM others can score causal methods on, where every answer is known and the wrong
Definition-1 contrast fails in a diagnosable way.
