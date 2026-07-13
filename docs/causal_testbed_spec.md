# 4b Spec — Packaging the Spike–Wave Causal Testbed

*Track 4b of [`next_steps.md`](next_steps.md). A **specification**, not an
implementation (status: PLANNED). Turns the C-series (C0–C7,
[`causal_experiments.md`](causal_experiments.md) +
[`causal_spiral_experiments.md`](causal_spiral_experiments.md)) and its shared
machinery ([`ghca_causal.py`](../ghca_causal.py)) into a **reusable synthetic-SCM
benchmark** for the spike–wave causal framework of Jalaldoust & Zabeh
(arXiv:2511.06602). This is the roadmap's "most externally-reusable artifact":
low new-science, high utility — a ground-truth testbed others can score their
causal-discovery / epiphenomenality / causal-emergence methods on.*

---

## 0. Status and scope

**What it is.** Packaging + a scoring harness. The C-series already *did* the
science — validated the paper's certificates on ground truth (C1), showed `do(W)`
is fat-handed under real constitution (C2), `do(θ)` is the well-posed handle (C3),
the verdict is outcome-relative (C4), and carried all of this to a genuine 2-D
spiral (C5–C7). But it lives as eight standalone experiment scripts that recompute
overlapping machinery. 4b refactors that machinery into a **stable API + a registry
of ground-truth scenarios + a harness that scores an arbitrary method against the
ground truth**, plus docs and a regression suite.

**What it is not.** No new causal theory, no new substrate physics, no real neural
data, no identification-from-data (the paper proves that's impossible — not the
contribution). The existing substrate limitations (one medium; `W` = coherence /
active-fraction / winding) carry over. This is not a general SCM simulator; it is
specifically a *spike–wave constitution* testbed where `W = f(S)` is literal.

**Why it's worth doing.** In vivo you never know the causal graph (Causal Hierarchy
Theorem), so a method claiming to detect an epiphenomenal wave, or to pick the
well-posed intervention, cannot be validated. Here the SCM is *built*, so every
verdict has a ground-truth answer. No such synthetic benchmark for spike–wave
causality currently exists.

---

## 1. What already exists (the raw material)

- **`ghca_causal.py`** — the seed of the API:
  - wave variables `W = f(S)`: `wave_coherence`, `wave_active_fraction` (+ the
    spiral `local_winding` / net charge in `ghca_net.py`);
  - observation model: `make_observation_mask`, `read_spikes` (partial `S_obs`);
  - the three operators: `do_S` (clamp), `do_theta` (τ / coupling / threshold),
    `do_W` (force `f(S)=target` under policies `random` / `clustered` / `min_edit`
    — the fat-handedness probe).
- **C0–C4** (`experiments/c0..c4`): the scalar arc — instrumentation, the six
  canonical graphs + Theorem-1 certificate, fat-handed `do(W)`, `do(θ)` handle,
  outcome-relativity matrix + macro-sufficiency `I(B;W)/I(B;S)`.
- **C5–C7** (`experiments/c5..c7`): the spiral arc — `do(χ)` fat-handedness by
  reader, `do(θ_χ)` well-posedness + necessity, outcome-relativity + `θ→χ→B`
  mediation/screening.

4b's job is to lift the *shared* pieces out of these scripts into one API, keep the
scripts as thin callers (so nothing regresses), and add the registry + harness.

---

## 2. The API surface (proposed)

A single module `ghca_testbed.py` exposing:

- **`Scenario`** — a named, seeded ground-truth SCM. Fields: a substrate builder,
  the node roles (`S_obs`, hidden `U`, behaviour readout `B`), the wave variable(s)
  `W = f(S)`, the graph type, and — the crux — the **ground-truth answer** for each
  benchmark task (see §3): the true epiphenomenality verdict, the expected
  intervention effect, the realization-policy band, etc.
- **`Testbed`** — wraps a `Scenario` and exposes the operators as methods:
  `do_S`, `do_W(target, policy)`, `do_theta(...)`, `observe()` → `S_obs`,
  `wave()` → `W`, `behaviour()` → `B`; plus `rollout(interventions, n, seed)` →
  samples of `(S_obs, W, B)`.
- **Certificate / test library** (method-agnostic, ground-truth-free — what a user
  would apply to real data):
  - `epiphenomenality_test(samples_obs, samples_do)` — Definition 1: is
    `P(B|X; do(Z)) = P(B|X)`? (Monte-Carlo + distributional distance + shuffle
    null). *Careful:* score interventional-vs-**observational**, not
    `do(W=1)` vs `do(W=0)` — the C1 subtlety that mis-scores front-door.
  - `theorem1_certificate(graph)` — the d-separation removability certificate;
    `frontdoor_certificate(graph)` (Prop 2).
  - metrics: `fat_hand_band(effect_by_policy)` (σ units), `mediation(θ, W, B)`
    (screening-off), `macro_sufficiency(B, W, S)` = `I(B;W)/I(B;S)`,
    `outcome_matrix(handles, outcomes)`.
- **`score(method_fn, scenarios)`** — the harness: runs a user's method over the
  registry and reports, per scenario, whether its verdict/estimate matches the
  ground truth (and by how much).

Design rules: everything seeded (`default_rng`, threaded explicitly — the
`perturb_tau` lesson); substrate/analysis boundary explicit; the C-scripts must
produce **bit-identical** numbers through the new API (regression gate).

---

## 3. Benchmark tasks (each ships a ground-truth answer)

| task | question a method must answer | ground truth from | pass condition |
|------|-------------------------------|-------------------|----------------|
| **T1 certificate** | is `W` epiphenomenal to `B` given `S_obs`? | C1 six graphs | verdict == Theorem-1 == built graph |
| **T2 fat-handedness** | does `do(W=w)` have a well-defined effect? | C2 | report policy band; flag ill-posed (labeled-line 33σ) vs pinned (collective 0.24σ) |
| **T3 handle selection** | which handle is well-posed, `do(W)` or `do(θ)`? | C3 | `do(θ)` band ≈ 0σ ≪ `do(W)` |
| **T4 outcome-relativity** | for which outcome is `W` causal? | C4 | recover the **diagonal** `do(θ)` matrix + macro-sufficiency 1.03 vs 0.11 |
| **T5 topological handle** | is spiral chirality a causal, reader-robust handle? | C5–C7 | `do(χ)` fat-handed by reader (center 6.2σ / tracked 1.0σ); `do(θ_χ)` 0σ; `θ→χ→B` mediation |

A user plugs in their own epiphenomenality test / effect estimator / discovery
method; the harness scores it against these known answers. The value is that a
method can *fail* here in a diagnosable way (e.g. mis-score front-door by using the
wrong Definition-1 contrast, or miss `do(W)` fat-handedness) before it is ever
trusted on data where no answer exists.

---

## 4. Build sequence (decision-gated)

- **4b.0 — extract the API *(gate).*** Move the shared machinery into
  `ghca_testbed.py`; rewrite C0–C7 as thin callers. *Gate:* every C0–C7 headline
  number reproduces bit-identically through the API (else the refactor changed
  behaviour). This alone is worth doing — it turns eight scripts into one library.
- **4b.1 — scenario registry.** Encode the C1 canonical graphs, the C2 constituted-`W`
  scenario, the C3/C4 E3-net, and the C5–C7 spiral as `Scenario`s with ground-truth
  answers. Export them as data (`JSON`/`npz`: adjacency, roles, expected effects) for
  language-agnostic reuse.
- **4b.2 — test/metric library + harness.** Package Definition-1, the certificates,
  and the metrics as method-agnostic functions; implement `score(method_fn,
  scenarios)`.
- **4b.3 — docs + minimal example.** A `causal_testbed.md` usage doc and a
  ~20-line "score your method" example; a `reproduce`/`benchmark` entry point that
  prints the ground-truth table and the model's own results.
- **4b.4 — regression suite.** Turn the C0–C7 headline numbers into assertions
  (a `test_ctestbed.py`), so the benchmark doubles as the C-series' regression test.

Stop-early-complete: 4b.0 + 4b.1 already deliver a reusable artifact (API + ground-
truth scenarios); the harness (4b.2) and docs (4b.3) make it *usable by others*;
4b.4 is hygiene.

---

## 5. Deliverables

- `ghca_testbed.py` — the API (`Scenario`, `Testbed`, certificate/test library,
  `score`).
- A scenario registry + exported ground-truth artifacts (`result/testbed/*.json`).
- `docs/causal_testbed.md` — usage doc + minimal example + the ground-truth table.
- `experiments/c0..c7` — refactored to thin callers (unchanged outputs).
- `test_ctestbed.py` — regression assertions over the C0–C7 headlines.

---

## 6. Effort / risk

- **Effort.** Medium — engineering/packaging, not research. The physics and the
  results already exist and are validated.
- **Risk.** Low. The main hazard is a refactor silently changing a C-series number;
  the 4b.0 bit-identical gate and the 4b.4 regression suite exist precisely to catch
  that. No scientific risk (nothing new is claimed).
- **Dependencies.** C0–C7 (on `main`), `ghca_causal.py`. Independent of E9 / 2b / 4a.

---

## 7. One-line spine

4b.0 lift the shared machinery into one API (C0–C7 reproduce bit-identically) →
4b.1 register the ground-truth scenarios → 4b.2 package the certificates/tests +
a scoring harness → 4b.3 document it with a minimal example → 4b.4 lock it with a
regression suite: a synthetic SCM others can score spike–wave causal methods on,
where — unlike in vivo — every answer is known.
