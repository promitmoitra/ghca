# The coherent core: the dwell-free attractor set is a static local condition

**Status:** theory. Two hand proofs valid on **any** graph, plus a certified
`L×L` phase-space census. This closes the empirical half of
[`clock_shift_healing.py`](../experiments/clock_shift_healing.py), which
established Lemma R by hand (`step(c) = c+1` ⟺ `c` is dwell-free) and then
recorded, as an *empirical* fact, that "every live attractor with `S+1 ≥ 4` is
dwell-free". Lemma R says what dwell-free attractors *do*; it does not say which
configurations they *are*. That description turns out to be a closed form.

Experiment: [`experiments/coherent_core.py`](../experiments/coherent_core.py) ·
Archive: `result/topology/coherent_core.npz`

Setting: states 0 (receptive), `1..τa` (active), `τa+1..S` (passive),
`S = τa + τp`, `B = S + 1`, `θ = 1`. The census is an `L×L` von Neumann box with
**open** boundaries; the two theorems assume no lattice at all.

**Read [Honest scope](#honest-scope) before quoting anything here.** The single
most-misread point: this doc is about the **attractor set**, which is
regime-*independent*. It says nothing about basins, where the `τa ≥ τp` regime
line lives.

---

## The coherent set

> `C = { c : for every cell i, some neighbour j has (c_j − c_i) mod B ∈ [1..τa] }`

Every cell has an upstream neighbour that will be **active** at the moment that
cell next becomes receptive. The condition is static, local, and evaluable
without running the dynamics.

## Theorem C (proven; any graph, both regimes)

> `C` is exactly the union of the **dwell-free cycles** of the GH step map. On
> `C`, `step` is the global clock-shift `c → c+1`; `C` is `+1`-invariant; every
> dwell-free cycle is **live**; and every such cycle has period exactly `B`.
> Hence **#dwell-free attractors = |C| / B**, and `B` divides `|C|`.

*Proof.* (a) `step(c) = c+1` ⟺ every receptive cell of `c` has an active
neighbour: a nonzero state always advances `+1 mod B`, so the only cells that
can violate `step(c) = c+1` are those at 0, which map to 1 (fire) or 0 (dwell).
This is Lemma R restated cell-wise; call it `C₀(c)`. (b) `C = ⋂ₖ C₀(c+k)`: cell
`i` is receptive in `c+k` exactly when `k = −c_i`, and neighbour `j` is then
active exactly when `c_j − c_i ∈ [1..τa] mod B`, which is the defining condition
of `C`. (c) `c ∈ C` ⟹ the orbit `c → c+1 → … → c+S → c` is a dwell-free cycle,
and `c+k = c` forces `k ≡ 0 mod B`, so the period is exactly `B`. (d) Conversely
a dwell-free cycle has `step = +1` at each point by (a), hence is the `+1`-orbit
of any point, hence lies in `C`. (e) Liveness: pick any cell `i` and the shift
`k = −c_i`; by `C` some neighbour `j` has `c_j + k ∈ [1..τa]`, so `c+k` — a
point of the same cycle — has an active cell. ∎

## Theorem Z (proven; any graph)

> The all-zero fixed point is the **only** dead attractor.

*Proof.* On a cycle every cell's state sequence is periodic. A nonzero state
advances `+1 mod B` each step, so a cell that is ever nonzero reaches `S`, then
0, and to return (periodicity) must fire — passing through state 1, which is
active. So a cycle with no active cell has every cell identically 0. ∎

This is the fact [`persistent_set_3x3.py`](../experiments/persistent_set_3x3.py)'s
boolean label propagation to the all-zero sink relies on. That script checked it
exhaustively at 2×2; it now holds for every graph and every `(τa, τp)`.

---

## The census

Exhaustive, deterministic, no RNG. Every row asserts Theorem C, Theorem Z, and
`|C|/B = ` the independently counted dwell-free attractors.

| lattice | (τa, τp) | regime | configs | P | \|C\| | attractors = \|C\|/B | live attractors | dwelling | max transient |
| :---: | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |
| 2×2 | (1, 1) | τa ≥ τp | 81 | 0.2963 | 0 | 0 | 2 | 2 | 4 |
| 2×2 | (2, 1) | τa ≥ τp | 256 | 0.7031 | 84 | 21 | 21 | 0 | 5 |
| 2×2 | (1, 2) | τa < τp | 256 | 0.1250 | 8 | 2 | 2 | 0 | 5 |
| 2×2 | (2, 2) | τa ≥ τp | 625 | 0.3200 | 40 | 8 | 8 | 0 | 8 |
| 2×2 | (3, 3) | τa ≥ τp | 2,401 | 0.3265 | 224 | 32 | 32 | 0 | 12 |
| 2×2 | (4, 4) | τa ≥ τp | 6,561 | 0.3292 | 720 | 80 | 80 | 0 | 16 |
| 2×2 | (5, 5) | τa ≥ τp | 14,641 | 0.3306 | 1,760 | 160 | 160 | 0 | 20 |
| 2×2 | (2, 5) | τa < τp | 4,096 | 0.0039 | 16 | 2 | 2 | 0 | 12 |
| 3×3 | (1, 1) | τa ≥ τp | 19,683 | 0.7462 | 84 | 28 | 62 | 34 | 7 |
| 3×3 | (2, 1) | τa ≥ τp | 262,144 | 0.9746 | 48,196 | 12,049 | 12,049 | 0 | 9 |
| 3×3 | (1, 2) | τa < τp | 262,144 | 0.4580 | 136 | 34 | 34 | 0 | 8 |
| 3×3 | (2, 2) | τa ≥ τp | 1,953,125 | 0.7752 | 47,520 | 9,504 | 9,504 | 0 | 12 |

The `P` column reproduces
[`persistent_set_structure.md`](persistent_set_structure.md) §3 and §9 cell for
cell — an independent re-derivation, since this pipeline counts cycles directly
rather than propagating labels to the sink.

### Four readings

1. **The attractor architecture is regime-independent.** Theorem C holds
   verbatim at `τp > τa` — (1,2), (2,5), 3×3 (1,2). The `τa ≥ τp` line, which
   governs spectrum-sufficiency, damage relaxation and debt confinement,
   does **not** touch the attractor set; it governs only which initial
   conditions reach it. This corroborates
   [`clock_shift_merge.py`](../experiments/clock_shift_merge.py)'s "the entire
   regime asymmetry lives in the transients" from the attractor side.
2. **(1,1) is the only cell with dwelling attractors**, and the only cell in the
   project's range with `S+1 < 4`. At 2×2 (1,1) `|C| = 0` yet two live
   period-4 attractors exist; at 3×3 (1,1) `C` supplies 28 of the 62 live
   attractors and the other 34 dwell. At every cell with `S+1 ≥ 4`, live
   attractors are exactly `C`, with zero left over and period exactly `B`.
3. **Phase space is wide and shallow.** Transients are `O(S)`, not `O(L²)`, and
   do not grow with `L`: 12 at 3×3 (2,2) exhaustively, 11 at `L=8` (2,2)
   sampled. The latest death observed at any cell or size is `t = 12`. Yet the
   attractor set is a small and rapidly shrinking slice of the space —
   `|C|/B^(L²)` is 2.43% at 3×3 (2,2) and falls exponentially in `L²` (the
   κ table below) — an enormously wide, very shallow funnel.
4. **The dying set is the small-lattice phenomenon.** `P → 1` fast in `L`.

### Sampled: live attractors stay inside C past the exhaustive range

Seeded random inits landed on their attractors, every attractor state checked
against `C` (`default_rng(20260816)`, 150 samples per point).

| L | (2,1) live→in C | (1,2) live→in C | (2,2) live→in C |
| :---: | :---: | :---: | :---: |
| 3 | 142 → 142 | 58 → 58 | 111 → 111 |
| 4 | 150 → 150 | 100 → 100 | 141 → 141 |
| 5 | 150 → 150 | 135 → 135 | 149 → 149 |
| 6 | 150 → 150 | 148 → 148 | 150 → 150 |
| 8 | 150 → 150 | 150 → 150 | 150 → 150 |

Every sampled live attractor has period exactly `B` and lies in `C`. Distinct
attractors per live sample: at (2,1), 132/142 at `L=3` and one-per-sample from
`L=4` on; at (1,2) and (2,2) one-per-sample from `L=6` on (`L=5`: 130/135 and
149/149). So past `L ≈ 5` essentially every initial condition lands on its own
attractor — no basin dominates.

### Sampled P(L)

Horizon-truncated (a run still active at the horizon is counted persistent), so
these are **upper** bounds; the exact values are the `P` column above.

| (τa, τp) | L=2 | L=3 | L=4 | L=5 | L=6 | L=8 |
| :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| (1, 1) | 0.329 | 0.738 | 0.952 | 0.994 | 1.000 | 1.000 |
| (2, 1) | 0.708 | 0.974 | 1.000 | 1.000 | 1.000 | 1.000 |
| (1, 2) | 0.133 | 0.456 | 0.750 | 0.923 | 0.982 | 0.999 |
| (2, 2) | 0.338 | 0.780 | 0.966 | 0.997 | 1.000 | 1.000 |
| (3, 3) | 0.346 | 0.796 | 0.969 | 0.999 | 1.000 | 1.000 |
| (2, 3) | 0.191 | 0.564 | 0.845 | 0.970 | 0.997 | 1.000 |

### What exactly is sampled

Every sampled quantity in this doc draws **raw configurations**, i.i.d. uniform
per cell over `{0..S}` — never the coarse-grained gap-signature space of
[`persistent_set_structure.md`](persistent_set_structure.md). That choice is
load-bearing three times over:

- **It is the right measure for `P`.** Moitra & Sen define `P` as the fraction
  of initial *configurations* that sustain, so config-uniform sampling estimates
  exactly that. Signature classes have wildly unequal sizes (14,641 configs →
  91 classes at 2×2 (5,5)), so a class-uniform draw would estimate a different
  number entirely.
- **Class-uniform sampling would not even be well defined here.** The signature
  is fate-pure only at `τp ≤ τa`, and the sampled sweep deliberately includes
  (1,2) and (2,3), where classes are mixed — a class has no fate to sample.
- **It makes the `κ` estimator a rare-event estimator**, since `C` is an
  exponentially small subset. That is the weak point of this section and is
  handled explicitly below rather than hidden.

The signature compression is therefore *not* used as a variance-reduction
device anywhere here. It remains the obvious lever for pushing the exhaustive
census past 3×3 (§9 of the companion doc), which is a separate job.

### Attractor-entropy density

`#attractors = |C|/B`, so `κ = ln(#attractors)/L²`. Note that **`κ ≤ ln B` by
construction** — it would equal `ln B` if every configuration were coherent — so
`κ` rising with `L` is largely that trivial ceiling. The informative quantity is
the **per-cell coherence cost `c = ln B − κ`**.

`|C|` is taken exactly from the census where available; otherwise `|C|/B^(L²)`
is estimated by Monte Carlo over configurations, with draws escalated until at
least **300** land in `C`. Points that cannot reach 300 hits at 4,000,000 draws
are **dropped, not quoted**.

| (τa, τp) | ln B | L | source | draws | hits | \|C\|/B^(L²) | κ | c = ln B − κ | SE(κ) |
| :---: | ---: | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| (2, 1) | 1.386 | 2 | exact | — | 84 | 3.281e-01 | 0.7611 | 0.6252 | 0.0000 |
| (2, 1) | 1.386 | 3 | exact | — | 48,196 | 1.839e-01 | 1.0441 | 0.3422 | 0.0000 |
| (2, 1) | 1.386 | 4 | MC | 40,000 | 3,661 | 9.152e-02 | 1.1502 | 0.2361 | 0.0010 |
| (2, 1) | 1.386 | 5 | MC | 40,000 | 1,635 | 4.088e-02 | 1.2030 | 0.1833 | 0.0010 |
| (2, 1) | 1.386 | 6 | MC | 40,000 | 671 | 1.678e-02 | 1.2342 | 0.1521 | 0.0011 |
| (2, 1) | 1.386 | 8 | MC | 400,000 | 667 | 1.667e-03 | 1.2647 | 0.1216 | 0.0006 |
| (2, 1) | 1.386 | 10 | MC | 4,000,000 | 399 | 9.975e-05 | 1.2803 | 0.1060 | 0.0005 |
| (2, 2) | 1.609 | 2 | exact | — | 40 | 6.400e-02 | 0.5199 | 1.0896 | 0.0000 |
| (2, 2) | 1.609 | 3 | exact | — | 47,520 | 2.433e-02 | 1.0177 | 0.5917 | 0.0000 |
| (2, 2) | 1.609 | 4 | MC | 400,000 | 1,977 | 4.942e-03 | 1.1770 | 0.4325 | 0.0014 |
| (2, 2) | 1.609 | 5 | MC | 4,000,000 | 2,709 | 6.773e-04 | 1.2532 | 0.3563 | 0.0008 |
| (3, 3) | 1.946 | 2 | exact | — | 224 | 9.329e-02 | 0.8664 | 1.0795 | 0.0000 |
| (3, 3) | 1.946 | 3 | MC | 40,000 | 1,627 | 4.068e-02 | 1.3739 | 0.5720 | 0.0027 |
| (3, 3) | 1.946 | 4 | MC | 40,000 | 403 | 1.008e-02 | 1.5369 | 0.4090 | 0.0031 |
| (3, 3) | 1.946 | 5 | MC | 400,000 | 801 | 2.002e-03 | 1.6195 | 0.3264 | 0.0014 |
| (3, 3) | 1.946 | 6 | MC | 4,000,000 | 1,283 | 3.208e-04 | 1.6684 | 0.2775 | 0.0008 |

Dropped as under-powered: (2, 2) L=6, 8, 10 and (3, 3) L=8, 10.

`κ` is positive at every measured point, so multistability is extensive in
system size. `c` falls monotonically in `L` at all three cells — and **it is
not converged**: `c` is still falling at the largest affordable `L`, so no
limit is claimed for either `c` or `κ`. The shape of the decay (`c` roughly
`c∞ + b/L`, consistent with a boundary correction — corner cells have two
neighbours and so are the hardest to satisfy) is suggestive but is not fitted
or asserted here.

> **⚠ Superseded measurement.** The first revision of this section used a bare
> fraction floor of `5×10⁻⁵` at `n = 40,000`, which admitted `κ` points resting
> on **one to three hits**: (2,2) L=6 rested on 1 hit and (2,1) L=10 on 3, and
> they were tabulated as if they were data. Those are now dropped or
> re-measured. Re-measuring moved (2,2) L=5 from `κ = 1.237` (18 hits) to
> `1.2532` (2,709 hits) — a real shift, not rounding. The estimator now
> escalates on **hit count**, and the script asserts that every quoted point is
> exact or ≥ 300 hits.

---

## Substrate / analysis boundary

- **Substrate:** the step map, its cycles, and dwell events. Computed by
  enumeration of the transition rule; no readout, no learned parameter.
- **Analysis:** `C` itself is a predicate imposed on configurations. Theorem C
  says the two coincide as *sets*; it does **not** say the substrate computes,
  represents, or is sensitive to `C`.
- **Not measured here at all:** anything behavioural. No task, no decoder.

## Honest scope

- **The theorems are about attractors, not basins.** Nothing here decides
  persistence of an arbitrary configuration — that is the open problem of
  [`spectrum_sufficiency_proof.md`](spectrum_sufficiency_proof.md), and it is
  untouched. `C` describes where trajectories *end*, not which ones get there.
- **"Every live attractor is dwell-free when `S+1 ≥ 4`" is NOT proven.** It is
  certified exhaustively at the twelve rows above and sampled to `L = 8`. It is
  the one load-bearing claim here without a proof, and it is what licenses
  reading `C` as *the* attractor set rather than merely the dwell-free part.
- **The girth reading is a conjecture.** 4 is the girth of the square lattice
  and (1,1) is the unique tested cell with `S+1 < 4`, which makes "the
  threshold is the girth" a natural guess — but no second girth is tested, so
  one cell is carrying that entire interpretation. A triangular or hexagonal
  lattice would settle it cheaply.
- **`θ = 1`, von Neumann, open boundary.** As everywhere else in this thread.
  The torus is untested; prediction P-C of
  [`coherence_larger_lattices_handoff.md`](coherence_larger_lattices_handoff.md)
  expects boundary structure to matter.
- **`κ` is not converged** (above), and the sampled `P(L)` columns are upper
  bounds. All sampling is over raw configurations, uniform per cell; the
  gap-signature quotient is not used to reduce variance anywhere.
- **The `κ` estimator is a rare-event estimator.** At the largest `L` it is
  counting a few hundred hits in millions of draws, and its first revision
  quoted points backed by one hit. It is now hit-count-governed, but this
  remains the least robust section of the doc: an exact transfer-matrix count
  of `|C|` (a `B^(2L)`-state DP) would replace it outright and is the right
  next step if these numbers ever have to carry weight.
- **3×3 stops at (2,2).** `B = 6, 7` at 3×3 (10M and 40M configurations) were
  not run in the exhaustive census; the largest exhaustive row is 1,953,125
  configurations.

## Reproduce

```bash
python3 experiments/coherent_core.py
```

Exhaustive sections are RNG-free; sampled sections thread `default_rng(20260816)`
explicitly. The script asserts every finding above and aborts on regression;
reruns are bit-identical. Runtime is a few minutes, dominated by 3×3 (2,2).
