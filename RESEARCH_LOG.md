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

## 2026-06-23 — Chirality, D4 orbits, and the orientation character

### Method
The 2×2 Von Neumann periodic lattice forms a 4-cycle ring:
`(0,0)→(0,1)→(1,1)→(1,0)→(0,0)` (CW).  Each config is a discrete phase field
on this ring.  **Chirality** is computed by voting on ring-step diffs modulo
M = tau0+1: "forward" diffs (< M/2) vs "backward" diffs (> M/2).

```python
diffs = (ring_states[1:] - ring_states[:-1]) % M   # CW ring diffs
chirality = sign(n_forward - n_backward)             # +1 CW, -1 CCW, 0 mixed
```

**D4 canonical form:** the lexicographically smallest config among all 8 images
under the dihedral group of the square (4 rotations × 2 reflections).

### Results

| pair   | n_configs | CW | CCW | mixed | net Q | purity | D4 orbits |
|--------|-----------|----|----|-------|-------|--------|-----------|
| (2,10) |     8     |  4 |  4 |   0   |   0   |  1.00  |     1     |
| (3,14) |     8     |  4 |  4 |   0   |   0   |  1.00  |     1     |
| (3,13) |    40     | 20 | 20 |   0   |   0   |  1.00  |     5     |
| (4,16) |   120     | 60 | 60 |   0   |   0   |  1.00  |    15     |
| (4,17) |    40     | 20 | 20 |   0   |   0   |  1.00  |     5     |
| (1,1)  |    46     | 24 | 12 |  10   |  12   |  0.78  |     9     |
| (4,4)  |  4060     |1500|1140|1420   | 360   |  0.65  |   590     |
| (6,17) |  8940     |3236|2964|2740   | 272   |  0.69  |  1135     |

### Key findings

**1. Theorem (verified for all 289 pairs):**
> resting\_frac == 0  ⟺  chirality\_purity == 1.0

Configs with no resting cells always have definite chirality (±1); configs with
any resting cell can be chirality-ambiguous.  The equivalence is exact — zero
exceptions across the full 17×17 parameter grid.

Mechanistic reason: a resting cell (state 0) creates a phase discontinuity in
the ring field.  The diff across it can land on either side of M/2, making the
chirality vote a tie.  Without resting cells, all diffs are bounded away from
M/2 for arithmetic-progression multisets.

**2. Net chirality Q = n_CW − n_CCW = 0 for all resting-cell-free pairs.**
The config set is perfectly balanced.  No preferred handedness.

**3. D4 orbit structure: orientation character of D4.**
Every D4 orbit for the resting-cell-free pairs has size 8 (= |D4|) and contains
exactly 4 CW + 4 CCW configs.  Within an orbit:
- The 4 **rotations** (det = +1, orientation-preserving) map CW ↔ CW, CCW ↔ CCW.
- The 4 **reflections** (det = −1, orientation-reversing) flip chirality: CW ↔ CCW.

The chirality W ∈ {±1} is therefore the **determinant character** of D4 — the
unique nontrivial 1-dimensional real representation of D4 that assigns +1 to
rotations and −1 to reflections.  This is the natural integer invariant.

**4. Orbit chirality pattern is universal for (3,13).**
All 5 orbits share the identical chirality pattern `[+1,−1,−1,+1,−1,+1,+1,−1]`
regardless of which multiset they belong to.  The CW/CCW assignment is
determined entirely by spatial arrangement within D4, not by state values.

**5. Multisets of (3,13) are near-uniform phase waves.**
All 5 multisets have states with spacing ≈ act=3:
```
{2,5,8,11}   — shifted by -1 from uniform
{3,5,8,11}   — one defect spacing of 2
{3,6,8,11}   — one defect spacing of 2
{3,6,9,11}   — one defect spacing of 2
{3,6,9,12}   — perfectly uniform {act, 2act, 3act, 4act}
```
These are the only multisets with all states in (0, tau0] and uniform or
near-uniform spacing.  The "defect" variants arise because tau0=16 = 4*act+4,
leaving one "slot" short of the next uniform step.

### Integer representation
The single integer that best characterises each config's spatial topology is its
**topological charge W ∈ {−1, 0, +1}** (chirality index).  For the config set:
- **purity** = (n_CW + n_CCW)/n_total — fraction with definite topology
- **net charge** Q = n_CW − n_CCW — handedness bias

Resting-cell-free ↔ (purity=1, Q=0): the "topologically pure and balanced" corner.

### Open questions
- Does chirality-homogeneous seeding (all CW or all CCW) shift d_crit or the
  upper extinction threshold relative to chirality-mixed seeding?
- For generic pairs, does the excess CW (Q > 0) reflect a physical asymmetry or
  is it an artifact of the D4-canonical ring ordering convention?
- Do CW+CW adjacent core pairs nucleate more readily than CW+CCW pairs?
  (Testable by forcing chirality in `make_lattice_persistent`.)
- Is purity correlated with d_crit independently of n_configs?

---

## Code

| file | description |
|------|-------------|
| `ghca_main.py` | Original GHCA simulator (pure Python loop, matplotlib visualisation) |
| `ghca_cluster.py` | Multi-core embedding experiment; vectorised CA step; density sweep |
| `ghca_cluster.py::decode_all` | Vectorised base-n config decoder (~10,000× faster than ghca_core) |
| `ghca_cluster.py::analyze_config_structure` | Classify a pair's config set: n_configs, resting_frac, n_multisets, is_phase_wave |
| `ghca_cluster.py::chirality_batch` | Topological charge W∈{±1,0} for a batch of configs (orientation character of D4) |
| `ghca_cluster.py::d4_canonical` | D4 canonical form (lex-min over 8 symmetry images) |
| `ghca_cluster.py::chirality_summary` | Per-pair chirality distribution: n_cw, n_ccw, n_mix, net_charge, purity, n_d4_orbits |
| `result/act_config_ids_states-(AA,PP)_core-04.npy` | Pre-computed persistent config IDs for each (act,pas) pair |
