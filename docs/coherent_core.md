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

## Lemma E — when is `C` empty? (proven; any graph)

> If `C` is non-empty then `G` has a cycle of length `ℓ` and an integer `k ≥ 1`
> with `ℓ ≤ kB ≤ ℓ·τa`.

*Proof.* Take `c ∈ C`. Each cell `i` picks a neighbour `σ(i)` whose lag
`(c_σ(i) − c_i) mod B` lies in `[1..τa]`. `σ` is a self-map of a finite set, so
its functional digraph contains a cycle of *distinct* nodes; since `σ(i)` is
always a neighbour, that is a cycle of `G`, of some length `ℓ`. Summing the lags
once around it returns to the starting phase, so the sum is `0 mod B`; the sum
also lies in `[ℓ, ℓ·τa]` and is at least `ℓ ≥ 3 > 0`, so it equals `kB` for some
`k ≥ 1`. ∎

> **⚠ This corrects a conjecture made in the first revision of this doc**, which
> read the (1,1) anomaly as a **girth** effect — "`S+1 ≥ 4`, and 4 is the girth
> of the square lattice". That reading is wrong. The obstruction is arithmetic,
> not metric: the square lattice is **bipartite**, so every cycle length is
> even, while `B = 3` at (1,1) is odd. At 2×2 the only cycle length is 4 and no
> multiple of 3 lies in `[4, 4]`, so `C` is empty; at 3×3 a 6-cycle exists and
> `6 = 2·3` lies in `[6, 6]`, so `C` is not (`|C| = 84`). Girth alone predicts
> neither.

Necessity is asserted at 16 (lattice, cell) pairs. The converse held at all 16
but is **not** claimed.

| lattice | cell | B | cycle lengths | lemma allows C | \|C\| |
| :---: | :---: | ---: | :---: | :---: | ---: |
| 2×2 | (1,1) | 3 | {4} | no | 0 |
| 2×2 | (1,3) | 5 | {4} | no | 0 |
| 2×2 | (1,4) | 6 | {4} | no | 0 |
| 2×2 | (1,5) | 7 | {4} | no | 0 |
| 3×3 | (1,1) | 3 | {4, 6, 8} | yes | 84 |
| 3×3 | (1,3) | 5 | {4, 6, 8} | no | 0 |
| 3×3 | (1,4) | 6 | {4, 6, 8} | yes | 168 |
| 3×3 | (1,5) | 7 | {4, 6, 8} | no | 0 |

(Abridged to the discriminating rows; the other eight — (2,1), (1,2), (2,2),
(2,5) at both lattices — all have `C` allowed and non-empty.)

Two things fall out. Lemma E **accounts for the three 2×2 cells with `P = 0`
exactly**, which [`persistent_set_structure.md`](persistent_set_structure.md) §4
records as "vacuous" without explaining them. And it correctly separates
3×3 (1,4) — where `6 = 1·6` puts `|C| = 168` — from 3×3 (1,3) and (1,5), where
no multiple of `B` lands in any `[ℓ, ℓ]`.

**Lemma E does not settle the open claim.** It says when `C` is *empty*; the
un-proven claim says when *dwelling* attractors are *absent*. These are
different: 3×3 (1,1) has a non-empty `C` (28 attractors) **and** 34 dwelling
ones. See [Honest scope](#honest-scope).

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

`|C|` now comes from three sources, in descending order of authority: the
exhaustive **census**; the exact transfer-matrix **DP** below; and only then
**MC** over configurations, with draws escalated until at least **300** land in
`C`. Points reaching none of these are **dropped, not quoted**.

#### The exact transfer-matrix count

Cell `m`'s coherence test reads `m−L, m−1, m+1, m+L`, so `m` can be finalised
the moment cell `m+L` is placed. Carrying the last `2L` cells as DP state is
therefore sufficient, and transitions branch only `B` ways — a row-by-row DP
would branch `B^L`. Cost is `O(L² · B^(2L+1))`; two guards (state space, and
`int64` overflow of the `B^(L²)` count ceiling) return a *reason* rather than a
wrong number.

It reproduces **all 12 census counts exactly** — an independent check of both,
since the census counts *cycles* and the DP counts a *static predicate*. Reach:
`L ≤ 5` at `B ≤ 5`, `L ≤ 4` at `B = 7`, and `L ≤ 6` at `B = 3`.

**It does not retire sampling entirely** — an earlier draft of this section said
it would, which was too strong. MC still carries (2,1) `L = 6, 8, 10` and
(3,3) `L = 5, 6`. What it does do is convert every point up to `L = 5` to exact,
and **validate the estimator that remains**: over the 15 (cell, L) pairs where
DP and MC overlap, every deviation is under **1.2 SE** (max 1.12). The first
revision's MC numbers were sound estimates — they were merely under-powered in
the tail, exactly as diagnosed.

| (τa, τp) | ln B | L | source | draws | hits | \|C\|/B^(L²) | κ | c = ln B − κ | SE(κ) |
| :---: | ---: | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: |
| (1, 1) | 1.099 | 3 | census | — | 84 | 4.268e-03 | 0.3702 | 0.7284 | 0.0000 |
| (1, 1) | 1.099 | 4 | DP | — | 17,202 | 3.996e-04 | 0.5409 | 0.5577 | 0.0000 |
| (1, 1) | 1.099 | 5 | DP | — | 16,574,562 | 1.956e-05 | 0.6210 | 0.4776 | 0.0000 |
| (1, 1) | 1.099 | 6 | DP | — | 86,148,328,050 | 5.740e-07 | 0.6689 | 0.4297 | 0.0000 |
| (2, 1) | 1.386 | 2 | census | — | 84 | 3.281e-01 | 0.7611 | 0.6252 | 0.0000 |
| (2, 1) | 1.386 | 3 | census | — | 48,196 | 1.839e-01 | 1.0441 | 0.3422 | 0.0000 |
| (2, 1) | 1.386 | 4 | DP | — | 393,419,652 | 9.160e-02 | 1.1503 | 0.2360 | 0.0000 |
| (2, 1) | 1.386 | 5 | DP | — | 45,560,019,060,572 | 4.047e-02 | 1.2026 | 0.1837 | 0.0000 |
| (2, 1) | 1.386 | 6 | MC | 40,000 | 671 hits | 1.678e-02 | 1.2342 | 0.1521 | 0.0011 |
| (2, 1) | 1.386 | 8 | MC | 400,000 | 667 hits | 1.667e-03 | 1.2647 | 0.1216 | 0.0006 |
| (2, 1) | 1.386 | 10 | MC | 4,000,000 | 399 hits | 9.975e-05 | 1.2803 | 0.1060 | 0.0005 |
| (1, 2) | 1.386 | 2 | census | — | 8 | 3.125e-02 | 0.1733 | 1.2130 | 0.0000 |
| (1, 2) | 1.386 | 3 | census | — | 136 | 5.188e-04 | 0.3918 | 0.9945 | 0.0000 |
| (1, 2) | 1.386 | 4 | DP | — | 72,872 | 1.697e-05 | 0.6131 | 0.7732 | 0.0000 |
| (1, 2) | 1.386 | 5 | DP | — | 152,445,152 | 1.354e-07 | 0.6982 | 0.6881 | 0.0000 |
| (2, 2) | 1.609 | 2 | census | — | 40 | 6.400e-02 | 0.5199 | 1.0896 | 0.0000 |
| (2, 2) | 1.609 | 3 | census | — | 47,520 | 2.433e-02 | 1.0177 | 0.5917 | 0.0000 |
| (2, 2) | 1.609 | 4 | DP | — | 740,712,190 | 4.854e-03 | 1.1759 | 0.4336 | 0.0000 |
| (2, 2) | 1.609 | 5 | DP | — | 204,946,279,408,620 | 6.877e-04 | 1.2538 | 0.3557 | 0.0000 |
| (3, 3) | 1.946 | 2 | census | — | 224 | 9.329e-02 | 0.8664 | 1.0795 | 0.0000 |
| (3, 3) | 1.946 | 3 | DP | — | 1,615,376 | 4.003e-02 | 1.3721 | 0.5738 | 0.0000 |
| (3, 3) | 1.946 | 4 | DP | — | 344,337,172,900 | 1.036e-02 | 1.5387 | 0.4072 | 0.0000 |
| (3, 3) | 1.946 | 5 | MC | 400,000 | 801 hits | 2.002e-03 | 1.6195 | 0.3264 | 0.0014 |
| (3, 3) | 1.946 | 6 | MC | 4,000,000 | 1,283 hits | 3.208e-04 | 1.6684 | 0.2775 | 0.0008 |

`|C| = 0` exactly at (1,1) L=2 — no dwell-free attractor exists at all, so κ is
undefined there. That is a **fact explained by Lemma E**, not a measurement
failure, and it is reported separately from the drops.

Dropped as under-powered (< 300 hits at 4,000,000 draws): (1,1) L=8, 10; (1,2) L=6, 8, 10; (2,2) L=6, 8, 10; (3,3) L=8, 10.

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
- **The girth reading was wrong, and its replacement covers only half.** The
  first revision guessed that the `S+1 ≥ 4` threshold was the square lattice's
  girth. Lemma E replaces that with a proven arithmetic condition — and settles
  only the *emptiness* half. Why **dwelling attractors** are absent whenever
  `S+1 ≥ 4` is still unexplained and still rests on (1,1) alone. The cheap
  discriminator is a lattice of different girth **and parity**: triangular
  (girth 3, non-bipartite, odd cycles available) or honeycomb (girth 6,
  bipartite). Not run here.
- **`θ = 1`, von Neumann, open boundary.** As everywhere else in this thread.
  The torus is untested; prediction P-C of
  [`coherence_larger_lattices_handoff.md`](coherence_larger_lattices_handoff.md)
  expects boundary structure to matter.
- **`κ` is not converged** (above), and the sampled `P(L)` columns are upper
  bounds. All sampling is over raw configurations, uniform per cell; the
  gap-signature quotient is not used to reduce variance anywhere.
- **The surviving `κ` points are still Monte-Carlo.** The DP has retired
  sampling up to `L = 5`, but (2,1) `L = 6, 8, 10` and (3,3) `L = 5, 6` remain
  rare-event estimates of a few hundred hits in millions of draws. They are now
  hit-count-governed and cross-validated against the DP where the two overlap
  (all within 1.2 SE), which is as much assurance as sampling can give. Pushing
  the DP further needs either a larger state budget or arbitrary-precision
  counts — `B = 4, L = 6` is blocked by the `int64` ceiling (`4³⁶ ≈ 4.7×10²¹`),
  not by the state space.
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
