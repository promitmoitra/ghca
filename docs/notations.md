# Notations

**Status:** reference. A single place to look up what a symbol means, which
file defines it, and — the part that actually bites — **which symbols are
reused for different things in different threads**.

This doc describes notation only. It makes no claims; every result lives in
its own results doc. Where a symbol's meaning is thread-specific, the thread is
named. If you are adding a symbol, add it here and check
[Collisions](#collisions) first.

---

## 1. The lattice substrate

Defined in [`ghca_main.py`](../ghca_main.py) (`Population`), used by the whole
lattice/theory thread.

| symbol | code name | meaning |
| :---: | :--- | :--- |
| `τa` | `act` | **active** (excited) duration — how long a cell spikes |
| `τp` | `pas` | **passive** (refractory) duration |
| `S` | `tau0` | `S = τa + τp`, the full non-receptive excursion |
| `B` | — | `B = S + 1`, the number of states per cell (the *base*) |
| `L` | `size` | side of the `L×L` box; `n = L²` cells |
| `θ` | `theta` | excitation threshold (number/weight of active neighbours needed) |
| `r` | `r` | neighbourhood radius; `r = 1` is von Neumann |
| `c` | `p` | a **configuration**: the full state array |

A cell's state is an integer in `{0, 1, …, S}`, in three bands:

| band | states | meaning |
| :--- | :---: | :--- |
| **receptive** | `0` | can be excited |
| **active** | `1 … τa` | spiking; excites receptive neighbours |
| **passive** | `τa+1 … S` | refractory; cannot be excited, does not excite |

**The step rule** (`θ = 1`, the setting of every theory result):

```
c_i = 0  and some neighbour is active   ->  1        (fire)
c_i = 0  and no neighbour is active     ->  0        (DWELL)
c_i != 0                                ->  (c_i + 1) mod B   (advance)
```

So a nonzero state always advances by exactly `+1 mod B`. **All** of the
system's freedom lives in what a receptive cell does. That single fact is what
Lemma R, Theorem C and the dwell-debt ledger are all built on.

### Derived terms

| term | meaning |
| :--- | :--- |
| **dwell** | a cell is receptive at `t` *and still receptive* at `t+1` |
| **dwell-free** | a configuration (or whole cycle) in which no cell dwells |
| **clock-shift** | `c → c+1`, adding 1 to every cell mod `B` |
| **live / persistent** | the trajectory's attractor contains an active cell |
| **dead / dying** | the trajectory reaches the all-zero fixed point |
| `P(τa, τp)` | **persistence probability** — the fraction of *initial configurations* that are live |
| `Ω` | the state space, `{0..S}^(L²)`, of size `B^(L²)` |

> ⚠ **`P` is a basin measure, not an order parameter.** It counts initial
> conditions, not a steady-state density. See
> [`coherent_core.md`](coherent_core.md) honest scope.

### The SIR aliasing in the code

`Population.count()` names the three bands with **epidemiological** letters,
which do not match the excitable-media words used everywhere in the docs:

| code | `scount` / `stime` | `icount` / `itime` | `rcount` / `rtime` |
| :--- | :---: | :---: | :---: |
| SIR reading | susceptible | infected | recovered |
| **docs reading** | **receptive** | **active** | **passive** |

The docstring at [`ghca_main.py`](../ghca_main.py) even says "restive, active or
passive". Read `s/i/r` as `receptive/active/passive` and ignore the epidemic
connotation — there is no infection model here.

---

## 2. The graph substrate

Defined in [`ghca_net.py`](../ghca_net.py) (`Network`) — the generalisation of
the above to an arbitrary weighted graph, used by the E-series and C-series.

| symbol | code | meaning |
| :---: | :--- | :--- |
| `W` | `W` | **weight / coupling matrix**, `N×N` — the adjacency of the substrate |
| `N` | `N` | number of nodes |
| `φ` | `phi` | per-node cyclic phase (the graph analogue of a cell's state) |
| `τ` | `tau` | per-node `act + pas`; **`τ` on a graph is `S`, not `B`** |
| `θ` | `theta` | per-node excitation threshold (weighted active-neighbour sum) |
| `p_s` | `p_s` | per-step **spontaneous firing** probability for a rested node |
| `ρ*` | `rho_star` | homeostatic target active fraction (`None` = static threshold) |
| `η_θ` | `eta_theta` | homeostatic threshold adaptation rate |

`p_s` is the only stochastic knob in the project; the lattice `Population` has
none and is fully deterministic.

### Topology quantities

| symbol | meaning | defined in |
| :---: | :--- | :--- |
| `β₁` | circuit rank `= m − N + c` (`m` edges, `N` nodes, `c` components) — the number of independent cycles | [`topology_cycle_capacity.md`](topology_cycle_capacity.md) |
| `K_dyn` | how many of those cycles can actually carry a wave, after the sustain gate | same (⚠ superseded — see [`topology_cycle_packing_exact.md`](topology_cycle_packing_exact.md)) |
| `w` | **GGH winding number** — phase increments summed once around a cycle, `= (S+1)·w`; `w ≠ 0` is the exact sustain criterion | [`topology_winding_capacity.md`](topology_winding_capacity.md) |
| `τ < L` | the **length gate**: a cycle of length `L` sustains only if `τ < L` (an approximation to `w ≠ 0`) | [`e2_results.md`](e2_results.md) |

> ⚠ In the length gate `τ < L`, **`L` is a cycle length**, not the lattice side.
> This is the worst collision in the project — see [Collisions](#collisions).

---

## 3. The theory thread

Symbols introduced by the persistent-set / spectrum / coherence / coherent-core
work. All assume `θ = 1`, von Neumann, open boundary.

| symbol | meaning | defined in |
| :---: | :--- | :--- |
| **gap signature** (gradient spectrum) | sorted multiset of cyclic phase gaps `mod B` around the reentrant loop | [`persistent_set_structure.md`](persistent_set_structure.md) §3 |
| **plaquette signature** | the 3×3 generalisation: sorted multiset of the four 2×2 plaquette codes | same, §9 |
| `z` | the all-equal-phases symbol `(0,0,0,0)` — the dead sink of the spectrum automaton | `spectrum_automaton.py` |
| **ING-1/2/3** | closure / activity / complement-death, the certificate schema ingredients | [`spectrum_sufficiency_proof.md`](spectrum_sufficiency_proof.md) |
| **P4, P5, P6** | the dwell laws and the regime-boundary relocation | same |
| `d_i(t)` | **damage** field: `u_i(t) − v_i(t) mod B` between a trajectory and its clock-shifted copy | `damage_relaxation.py` |
| `Δ_i(t)` | **dwell debt**: `D_i^v(t) − D_i^u(t)`, the cumulative dwell difference; `d = 1 + Δ` exactly | `dwell_debt_confinement.py` |
| **diagonal** | a pair `(v, u)` with `u = v + 1` cell-wise (a clock-shift pair) | [`coherence_invariant.md`](coherence_invariant.md) |
| **age** | BFS depth of a pair from the diagonal set | same |
| **coherence window** | the claim that age saturates at exactly `S` (a **regime** law: `τa ≥ τp` only) | same |
| `R` | the certified coherent pair set | same |
| `C` | the **coherent set** — see below | [`coherent_core.md`](coherent_core.md) |
| `κ` | attractor-entropy density, `ln(#attractors)/L²` | same |
| `c` (cost) | per-cell **coherence cost** `= ln B − κ` | same |

### The coherent set `C`

> `C = { c : for every cell i, some neighbour j has (c_j − c_i) mod B ∈ [1..τa] }`

Read the condition as: **every cell has an upstream neighbour that will be
active at the moment that cell next becomes receptive.** By Theorem C it is
exactly the union of the dwell-free cycles.

The quantity `(c_j − c_i) mod B` is called the **lag** from `i` to `j`.

### Worked example

`(τa, τp) = (2, 1)`, so `S = 3`, `B = 4`; states `0` receptive, `1–2` active,
`3` passive. On the open 3×3 box, cells indexed row-major `0..8`:

```
c =  2 3 1        bands =  active   passive  active
     3 2 0                 passive  active   receptive
     1 3 2                 active   passive  active
```

Every cell has a witness neighbour at lag `∈ {1, 2}`, so `c ∈ C`:

| cell | state | witness | its state | lag |
| ---: | ---: | ---: | ---: | ---: |
| 0 | 2 | 3 | 3 | 1 |
| 1 | 3 | 2 | 1 | 2 |
| 2 | 1 | 1 | 3 | 2 |
| 3 | 3 | 6 | 1 | 2 |
| 4 | 2 | 1 | 3 | 1 |
| 5 | 0 | 2 | 1 | 1 |
| 6 | 1 | 3 | 3 | 2 |
| 7 | 3 | 6 | 1 | 2 |
| 8 | 2 | 5 | 0 | 2 |

Theorem C then predicts the orbit is the pure clock-shift with period exactly
`B = 4`, and it is:

| t | configuration |
| ---: | :--- |
| 0 | `2 3 1 3 2 0 1 3 2` |
| 1 | `3 0 2 0 3 1 2 0 3` |
| 2 | `0 1 3 1 0 2 3 1 0` |
| 3 | `1 2 0 2 1 3 0 2 1` |
| 4 | `2 3 1 3 2 0 1 3 2` ← back to `t = 0` |

Cell 8 shows why the lag condition is stated on *lags* and not on "has an
active neighbour now": at `t = 0` its witness (cell 5) is **receptive**, not
active. But cell 8 only becomes receptive at `t = 2`, and by then cell 5 has
advanced to state 2 — active, exactly on time. Coherence is a statement about
phase *differences*, not about the current instant.

---

## 4. The C-series (causal instrumentation)

Defined in [`ghca_causal.py`](../ghca_causal.py). **This thread reuses `S`, `W`
and `B` for entirely different objects** — see [Collisions](#collisions).

| symbol | meaning |
| :---: | :--- |
| `S` | **spikes** — the node state (not `τa + τp`) |
| `S_obs` | partially observed spikes (a subset of nodes) |
| `W` | a **wave variable**, `W = f(S)`: a deterministic coarse-graining of the full state (e.g. Kuramoto coherence, active fraction) |
| `B` | the **behaviour readout** supplied by an experiment |
| `θ` | here means **generating parameters** broadly — timescales *and* couplings — not just the scalar threshold |
| `do(·)` | the intervention operator: `do(S)`, `do(W)`, `do(θ)` |
| **fat-handed** | one `do(W = w)` is realizable by many underlying states, so its effect on behaviour is a band rather than a value |

---

## 5. The E-series (learning)

| symbol | meaning |
| :--- | :--- |
| **Line A** | conduction-weight plasticity (learning on `W`, the coupling matrix) |
| **Line B** | timescale plasticity (learning on `τ`) |
| **E0 … E10** | the staged experiment series ([`learning_experiments.md`](learning_experiments.md)) |
| **C0 … C7** | the causal series ([`causal_experiments.md`](causal_experiments.md)) |
| `RIR` | **Readout Independence Ratio** — per-seed mean of `fixed_acc / trained_acc` ([`closed_loop_plasticity_results.md`](closed_loop_plasticity_results.md)) |
| `GVF` | general value function (the Horde demons of E6) |

---

## Collisions

These are real and they have caused errors. Always disambiguate by thread.

| symbol | meaning 1 | meaning 2 | meaning 3 |
| :---: | :--- | :--- | :--- |
| **`S`** | `τa + τp`, the excursion length (**lattice/theory**) | **spikes**, the node state (**C-series**, `W = f(S)`) | — |
| **`W`** | the **weight/coupling matrix** (`ghca_net.Network(W, …)`) | a **wave variable** `W = f(S)` (**C-series**) | the transfer-matrix **window** `2L` (`coherent_core.py`, local) |
| **`B`** | the **base** `S + 1` (**theory**) | the **behaviour readout** (**C-series**) | — |
| **`c`** | a **configuration** (**theory**) | the **coherence cost** `ln B − κ` (**coherent-core**) | **components** in `β₁ = m − N + c` (**topology**) |
| **`L`** | the **lattice side** (`L×L`) | a **cycle length**, in the sustain gate `τ < L` | — |
| **`d`** | the **damage** field `u − v` (**damage-relaxation**) | per-cell **dwell steps** `d_i` in the period law `T = k_i·S + d_i` | — |
| **`τ`** | on a **graph**: `act + pas` per node (= `S`) | in **E2**: the memory timescale being learned (same quantity, different role) | — |
| **`r`** | neighbourhood **radius** (`Population(r=…)`) | initial **refractory fraction** (`Population(r_0=…)`) | — |
| **`θ`** | a scalar **threshold** (substrate) | **generating parameters** broadly, in `do(θ)` (**C-series**) | — |
| **`P`** | **persistence probability** | the `Population` **class** | — |
| **DP** | *avoid* — reads as **directed percolation** in this field. Say **TM** (transfer matrix) for the dynamic-programming counter | — | — |

### The two that matter most

- **`S`** is the single worst one: in the theory docs `S = τa + τp` is an
  integer *timescale*, while in the C-series `S` is the *spike state itself*.
  A sentence like "`W = f(S)`" belongs to the C-series and has nothing to do
  with `τa + τp`.
- **`L`** in `τ < L` is a **cycle length**. The lattice-side `L` and the
  cycle-length `L` appear within a few lines of each other in
  [`lattice_capacities.md`](lattice_capacities.md) ("E2's `τ < L` law ports to
  the 2-D torus 28/28"). On a lattice the relevant cycle is the 4-plaquette.

---

## Conventions

- **Indexing.** Lattice cells are row-major, `i·L + j`, `0`-based. The 2×2
  theory thread orders cells around the reentrant loop as `0 → 1 → 3 → 2`, not
  `0,1,2,3` — the loop order, not the array order.
- **Boundary.** Open (non-periodic) unless a doc says torus. The 2×2 "core" is
  therefore a **4-ring**: every cell has degree 2, which is the source of
  several documented smallness artifacts.
- **`θ = 1`** in every theory result. `θ > 1` is untested there.
- **Regime names.** `τa ≥ τp` and `τp > τa` are *the* regime boundary; "strict"
  cells are `τa > τp`, "diagonal" cells are `τa = τp`.
- **Cells.** "Cell" is overloaded by convention: a `(τa, τp)` **parameter cell**
  in a sweep table, and a lattice **site**. Tables headed `(τa, τp)` mean the
  former; anything about neighbours means the latter.
- **Seeding.** Every RNG is an explicitly threaded `default_rng(seed)`; the
  global NumPy RNG is banned and audited
  (`review_helper.py audit-rng`). See [`AGENTS.md`](../AGENTS.md).
