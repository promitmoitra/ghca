# GHCA Research Log

A running record of findings, methods, and open questions from the study of
local timescales in Greenberg-Hastings Cellular Automata.

---

## 2026-06-23 — Multi-core embedding and density transitions

### Question
At what spatial density of embedded 2×2 cores does activity persist in a large
otherwise-resting lattice?

### Setup
- Lattice size L×L (default L=80), absorbing boundary conditions
- Cores placed on a regular grid of 2×2 slots; density = fraction of slots filled
- Two seeding modes: **random** (uniform draw from all states) and **persistent**
  (draw from the pre-computed persistent config set)
- Survival criterion: any active cell remaining at T=2000 steps

### Three regimes identified

**Regime 1 — Pre-selected persistent cores**
Individual cores drawn from the persistent config set are spiral-wave sources.
A single such core survives indefinitely at any non-zero density — unless the
config set is entirely resting-cell-free (see Regime 3).

**Regime 2 — Lower nucleation threshold (random seeding)**
Individually non-persistent random cores can *collectively* nucleate stable
activity once density exceeds a critical threshold d_crit. Scales with act/tau0:
lower active fraction → higher required density.

| pair    | tau0 | d_crit  | notes                        |
|---------|------|---------|------------------------------|
| (1,1)   |  2   | <0.005  | symmetric, persists easily   |
| (2,2)   |  4   | <0.005  |                              |
| (4,4)   |  8   | <0.005  |                              |
| (8,8)   | 16   | <0.005  |                              |
| (4,12)  | 16   | ~0.02   | long refractory              |
| (6,17)  | 23   | ~0.03   |                              |
| (1,7)   |  8   | N/A     | zero persistent configs      |
| (2,10)  | 12   | ~0.05   | phase-wave pair (see below)  |
| (3,14)  | 17   | ~0.09   | phase-wave pair              |

**Regime 3 — Upper extinction threshold (resting-cell starvation)**
At d=1 (every slot filled), if the entire persistent config set contains zero
resting cells, the lattice also contains zero resting cells. The GHCA infection
rule (resting → active) can never fire; all active cells advance to passive and
activity collapses in a single step.

The transition sharpens as d → 1: P(persist) drops from ~1 at d≈0.90 to 0 at
d=1.00, tracking the vanishing supply of resting cells in the few remaining
empty slots.

---

## 2026-06-23 — Full parameter space sweep: resting-cell-free pairs

### Method
Vectorised base-n decoder (`decode_all`) applied to all 289 (act, pas) pairs
in the 1..17 × 1..17 grid. For each pair: compute resting-cell fraction,
number of distinct sorted-state multisets, and whether all configs equal the
single phase-wave multiset `{act, 2·act, 3·act, 4·act}`.

### Result
Five pairs have entirely resting-cell-free persistent config sets:

| pair   | tau0 | n_configs | n_multisets | structure              |
|--------|------|-----------|-------------|------------------------|
| (2,10) |  12  |     8     |      1      | pure phase-wave        |
| (3,14) |  17  |     8     |      1      | pure phase-wave        |
| (3,13) |  16  |    40     |      5      | resting-cell-free, mixed |
| (4,16) |  20  |   120     |     15      | resting-cell-free, mixed |
| (4,17) |  21  |    40     |      5      | resting-cell-free, mixed |

All five sit at or near the extinction boundary — one step further in `pas` and
`n_configs` drops to zero entirely. The multiset counts follow a combinatorial
pattern (1 → 5 → 15 → ...) consistent with choosing 4-element multisets from
an increasingly dense uniform spacing.

**Pure phase-wave structure:** For (2,10) and (3,14), every persistent config
is a permutation of the single multiset `{act, 2·act, 3·act, 4·act}`. The 8
configs are the 8 spatial arrangements of this multiset that sustain activity.

### n_multisets grid (act=row, pas=col, * = all-phase-wave)

```
act/pas    1    2    3    4    5    6    7    8    9   10   11   12   13   14   15   16   17
      1    6   16   14   30   17   36    0    0    0    0    0    0    0    0    0    0    0
      2   26   35   59   66   82   81   86   61   15    1*   0    0    0    0    0    0    0
      3   61   91  115  165  191  228  233  243  205  152   70   28    5    1*   0    0    0
      4  117  175  235  285  377  435  506  527  545  503  430  335  201  114   51   15    5
      5  201  295  400  505  595  749  855  976 1026 1064 1019  937  814  660  455  303  179
      ...
```

(`*` marks pairs where all configs are pure phase-wave)

---

## 2026-06-23 — Spatial structure and chirality (open thread)

### Observation
The multiset representation collapses spatial information. Two configs from the
same multiset but different spatial arrangements can carry opposite wave
chirality — one CW spiral seed, one CCW — which will interact destructively
when adjacent cores collide.

### Why it matters
- **Lower threshold**: collective nucleation depends not just on the *number*
  of persistent configs but on whether neighbouring cores' wave directions
  reinforce or cancel. Chirality-mixed seeding may shift d_crit upward.
- **Upper threshold**: even before d=1, adjacent CW/CCW core pairs may
  annihilate. The sharpness of the upper transition could be partly driven by
  chirality annihilation, not only resting-cell depletion.

### The ring topology
Under Von Neumann neighbourhood on a 2×2 periodic lattice, cells form a 4-cycle:
`(0,0) ↔ (0,1) ↔ (1,1) ↔ (1,0) ↔ (0,0)`. Each config is a discrete phase
field on this ring. The winding number of that phase field is the natural
chirality invariant.

### Proposed approaches

**1. Ring winding number** (fast, vectorisable)
```python
ring_order = [(0,0),(0,1),(1,1),(1,0)]
states = config[ring_order]
diffs = np.diff(states, append=states[0])
wrapped = ((diffs + tau0//2) % (tau0+1)) - tau0//2
winding = wrapped.sum()   # >0 → CCW,  <0 → CW
```
Splits the 8 phase-wave configs into exactly two ±chirality classes of 4.
Generalises to non-phase-wave configs naturally.

**2. D4 orbit classification**
Group configs by their orbit under the dihedral group of the square (4
rotations × 2 reflections). Phase-wave configs should yield 2 orbits of size 4.
Mixed-multiset resting-cell-free configs may reveal additional topological
invariants beyond chirality.

**3. Empirical wave tracking**
Seed each config in a large resting lattice, run forward, track the rotation
direction of the expanding wavefront. Unambiguous, no analytical assumptions,
but costs one simulation per config ID.

### Open questions
- Do the 8 phase-wave configs for (2,10)/(3,14) split evenly into 4 CW + 4 CCW?
- For the mixed resting-cell-free pairs, does chirality still cleanly partition
  configs, or do the additional multisets introduce new topological classes?
- Does chirality-homogeneous seeding (all CW or all CCW) shift the upper
  extinction threshold relative to mixed seeding?
- Is the D4 orbit structure of the config set correlated with d_crit?

---

## Code

| file | description |
|------|-------------|
| `ghca_main.py` | Original GHCA simulator (pure Python loop, matplotlib visualisation) |
| `ghca_cluster.py` | Multi-core embedding experiment; vectorised CA step; density sweep |
| `ghca_cluster.py::decode_all` | Vectorised base-n config decoder (~10,000× faster than ghca_core) |
| `ghca_cluster.py::analyze_config_structure` | Classify a pair's config set: n_configs, resting_frac, n_multisets, is_phase_wave |
| `result/act_config_ids_states-(AA,PP)_core-04.npy` | Pre-computed persistent config IDs for each (act,pas) pair |
