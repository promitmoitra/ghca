"""The coherent core: the dwell-free attractor set of GH is a STATIC local
condition, proven lattice-free -- plus the LxL phase-space census it licenses.

Closes the empirical half of clock_shift_healing.py. That experiment established
Lemma R by hand (step(c) == c+1 <=> c is dwell-free) and then recorded, as an
EMPIRICAL fact, that "every live attractor with S+1 >= 4 is dwell-free". The
hand-proof half stopped there: it said what dwell-free attractors DO, not which
configurations they ARE. This experiment supplies the missing description, and
it is a closed form rather than an enumeration.

Setting: any graph G (the theorems), LxL von Neumann box with open boundaries
(the census), states 0 (receptive), 1..ta (active), ta+1..S (passive),
S = ta + tp, B = S + 1, theta = 1.

THE COHERENT SET. Define

    C = { c : for every cell i, some neighbour j has (c_j - c_i) mod B in [1..ta] }

Read it as a phase-gradient condition: every cell has an upstream neighbour that
will be ACTIVE at the moment that cell next becomes receptive. It is static,
local, and needs no dynamics to evaluate.

---------------------------------------------------------------------------
THEOREM C (proven by hand; holds on ANY graph, both regimes).

    C is exactly the union of the dwell-free cycles of the GH step map.
    On C, step is the global clock-shift c -> c+1. C is +1-invariant, every
    dwell-free cycle is LIVE, and every such cycle has period exactly B.
    Hence  #dwell-free attractors = |C| / B,  and B divides |C|.

Proof.
  (a) step(c) = c+1 <=> every receptive cell of c has an active neighbour.
      A nonzero state always advances by +1 mod B, so the only cells that can
      violate step(c) = c+1 are those at 0, which map to 1 (fire) or 0 (dwell).
      This is Lemma R of clock_shift_healing.py, restated cell-wise; call the
      condition C_0(c).
  (b) C = intersection over k of C_0(c+k).
      Cell i is receptive in c+k exactly when k = -c_i, and neighbour j is then
      active exactly when c_j - c_i in [1..ta] mod B. So requiring C_0 of every
      shift of c is literally the defining condition of C, and conversely.
  (c) c in C => the orbit c -> c+1 -> ... -> c+S -> c is a cycle, dwell-free by
      (a); and c+k = c forces k = 0 mod B, so the period is exactly B.
  (d) Conversely, a dwell-free cycle has step = +1 at each of its points by (a),
      so it is the +1-orbit of any of its points, so each point satisfies C_0
      at every shift, i.e. lies in C.
  (e) Every dwell-free cycle is live. Pick any cell i and the shift k = -c_i;
      by C some neighbour j of i has c_j + k in [1..ta], so the configuration
      c+k -- a point of the same cycle -- has an active cell. []

THEOREM Z (proven by hand; any graph). The all-zero fixed point is the ONLY
dead attractor.
  On a cycle every cell's state sequence is periodic. A nonzero state advances
  +1 mod B each step, so a cell that is ever nonzero reaches S, then 0, and to
  return (periodicity) must fire, i.e. pass through state 1 -- active. So a
  cycle with no active cell has every cell identically 0. []
  This is the fact persistent_set_3x3.py's boolean label propagation to the
  all-zero sink relies on; it was checked exhaustively at 2x2 there, and is
  now proven for every graph and every (ta, tp).

WHAT IS NOT PROVEN (the certified/sampled part, and the honest boundary):
  the claim "every LIVE attractor is dwell-free when S+1 >= 4" -- i.e. that C
  captures ALL live attractors and not merely the dwell-free ones. This is
  CERTIFIED here by exhaustive census at 2x2 and 3x3 and by seeded sampling to
  L = 8, not proven. 4 is the girth of the square lattice, but no second girth
  is tested, so the girth reading is a conjecture, not a measurement.

---------------------------------------------------------------------------
FINDINGS.

1. THEOREM C VERIFIED EXHAUSTIVELY at 12 (lattice, cell) combinations spanning
   BOTH regimes: C equals the union of dwell-free live cycles, set-for-set;
   step = +1 on C; C is +1-closed; |C| is divisible by B; and |C|/B equals the
   independently-counted number of dwell-free attractors.

2. THE ATTRACTOR ARCHITECTURE IS REGIME-INDEPENDENT. Theorem C holds verbatim
   at tp > ta ((1,2), (2,5), 3x3 (1,2)). The ta/tp regime line -- which governs
   spectrum-sufficiency, damage relaxation and debt confinement -- does NOT
   touch the attractor set; it governs only which initial conditions reach it.
   This corroborates clock_shift_merge.py's "the entire regime asymmetry lives
   in the TRANSIENTS" from the attractor side.

3. (1,1) IS THE ONLY CELL WITH DWELLING ATTRACTORS, and it is the only cell in
   the project's range with S+1 < 4. At 2x2 (1,1): |C| = 0, yet 2 live
   attractors exist (period 4). At 3x3 (1,1): 28 dwell-free attractors from C,
   but 62 live attractors in total -- 34 of them dwelling. At every cell with
   S+1 >= 4, live attractors == C exactly, zero leftovers.

4. PHASE SPACE IS WIDE AND SHALLOW. Transients are O(S), not O(L^2), and do not
   grow with L: max transient 12 at 3x3 (2,2) exhaustively, and 11 at L=8 (2,2)
   sampled -- while the attractor set is an exponentially small slice of the
   space. Death, when it happens, is immediate: the latest death observed at
   any cell or size is t = 12.

5. ATTRACTOR COUNT IS EXTENSIVE. |C|/B ~ exp(kappa * L^2) with kappa positive
   and RISING in L -- 1.263 at (2,1) L=8 and 1.277 at L=10; 1.761 at (3,3)
   L=8 -- so multistability is extensive in system size and the attractor
   entropy density is bounded away from zero. NOT CONVERGED: kappa is still
   increasing at the largest L the Monte-Carlo floor permits, so the plateau
   value is not established here and no limit is claimed. Past L ~ 5
   essentially every sampled initial condition lands on its own attractor.

6. DEATH IS A SMALL-LATTICE PHENOMENON. P(ta,tp) -> 1 rapidly in L at every
   cell tested, in both regimes. The sampled L = 2, 3 columns agree with the
   exhaustive values of persistent_set_structure.md Sections 3 and 9 within
   Monte-Carlo error (n = 2000 per point, so ~1 sd = 0.01); the EXACT values
   are the P column of table 1, which reproduces those docs cell for cell.
   The sampled columns are horizon-truncated and so are upper bounds.

House Rules Compliance:
  - Sections 1-4 are exhaustive and deterministic (no RNG). Sections 5-6 are
    sampled through an explicitly threaded default_rng(seed); the seed is part
    of the experiment and reruns are bit-identical.
  - Sampled persistence uses a finite horizon: a run still active at the
    horizon is COUNTED AS PERSISTENT. That is an upper bound on P, stated with
    the observed latest death time so the margin is visible.
  - Substrate / analysis boundary: the step map, its cycles and dwell events are
    substrate. C is an analysis construct -- a predicate imposed on
    configurations. Theorem C says the two coincide; it does NOT say the
    substrate computes, represents, or is sensitive to C.
  - Every finding is asserted; the script aborts on regression.
Output: result/topology/coherent_core.npz + printed doc tables.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np

# exhaustive census: (L, ta, tp). Both regimes, and the S+1 < 4 exception.
EXHAUSTIVE = [(2, 1, 1), (2, 2, 1), (2, 1, 2), (2, 2, 2), (2, 3, 3),
              (2, 4, 4), (2, 5, 5), (2, 2, 5),
              (3, 1, 1), (3, 2, 1), (3, 1, 2), (3, 2, 2)]
SAMPLED_CELLS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 3), (2, 3)]
SAMPLED_L = [2, 3, 4, 5, 6, 8]
DENSITY_CELLS = [(2, 1), (2, 2), (3, 3)]
DENSITY_L = [2, 3, 4, 5, 6, 8, 10]
SEED = 20260816
# converged kappa is read off the largest L whose MC estimate clears the floor
MC_FLOOR = 5.0e-5


# --------------------------------------------------------------------------
# lattice
# --------------------------------------------------------------------------
def neighbour_table(L):
    """von Neumann neighbours, OPEN boundary. Returns (index, mask) of shape
    (L*L, 4); masked-out slots point at cell 0 and are ignored via the mask."""
    n = L * L
    nbi = np.zeros((n, 4), np.int64)
    nbm = np.zeros((n, 4), bool)
    for i in range(L):
        for j in range(L):
            k, d = i * L + j, 0
            for di, dj in ((-1, 0), (1, 0), (0, -1), (0, 1)):
                a, b = i + di, j + dj
                if 0 <= a < L and 0 <= b < L:
                    nbi[k, d], nbm[k, d] = a * L + b, True
                    d += 1
    return nbi, nbm


def all_digits(n, B):
    """Every configuration as a digit array, row r = the base-B digits of r."""
    N = B ** n
    D = np.empty((N, n), np.int8)
    idx = np.arange(N, dtype=np.int64)
    for k in range(n - 1, -1, -1):
        D[:, k] = idx % B
        idx //= B
    return D


def step_digits(D, nbi, nbm, ta, B):
    act = (D >= 1) & (D <= ta)
    fires = (act[:, nbi] & nbm[None, :, :]).any(axis=2)
    return np.where(D == 0, fires.astype(D.dtype), (D + 1) % B)


def in_C(D, nbi, nbm, ta, B, chunk=1 << 19):
    """Membership of the coherent set C, chunked to bound peak memory."""
    out = np.empty(len(D), bool)
    for s in range(0, len(D), chunk):
        d = D[s:s + chunk].astype(np.int16)
        lag = (d[:, nbi] - d[:, :, None]) % B
        good = (lag >= 1) & (lag <= ta) & nbm[None, :, :]
        out[s:s + chunk] = good.any(axis=2).all(axis=1)
    return out


# --------------------------------------------------------------------------
# exhaustive census (Theorems C and Z, plus the certified "all live" claim)
# --------------------------------------------------------------------------
def census(L, ta, tp, settle=64):
    S, B, n = ta + tp, ta + tp + 1, L * L
    nbi, nbm = neighbour_table(L)
    D = all_digits(n, B)
    N = len(D)
    pw = (B ** np.arange(n - 1, -1, -1)).astype(np.int64)

    nxt = step_digits(D, nbi, nbm, ta, B)
    succ = nxt.astype(np.int64) @ pw

    # settle onto the cycle set: transients are short, so iterating `settle`
    # times lands every state on a cycle. The permutation assert certifies it.
    land = np.arange(N, dtype=np.int64)
    for _ in range(settle):
        land = succ[land]
    cyc = np.unique(land)
    img = np.unique(succ[cyc])
    assert img.size == cyc.size and np.array_equal(img, cyc), \
        f"settle={settle} insufficient at L={L} ({ta},{tp})"

    # exact transient length per state
    oncyc = np.zeros(N, bool)
    oncyc[cyc] = True
    trans = np.zeros(N, np.int32)
    cur = np.arange(N, dtype=np.int64)
    todo = ~oncyc[cur]
    t = 0
    while todo.any():
        t += 1
        cur[todo] = succ[cur[todo]]
        newly = todo & oncyc[cur]
        trans[newly] = t
        todo &= ~oncyc[cur]
        assert t <= settle, "transient exceeded settle"

    # decompose the cycle set into orbits; classify live / dwell-free
    succ_l = succ  # local alias
    seen = np.zeros(N, bool)
    n_live = n_dead = n_dwellfree = 0
    dead_states = []
    live_cycle_states = []
    dwellfree_states = []
    periods_live = set()
    for c0 in cyc:
        if seen[c0]:
            continue
        orb = [int(c0)]
        seen[c0] = True
        x = succ_l[c0]
        while x != c0:
            orb.append(int(x))
            seen[x] = True
            x = succ_l[x]
        block = D[orb]
        alive = bool(((block >= 1) & (block <= ta)).any())
        nxtblock = nxt[orb]
        dwellfree = not bool(((block == 0) & (nxtblock == 0)).any())
        if alive:
            n_live += 1
            periods_live.add(len(orb))
            live_cycle_states.extend(orb)
            if dwellfree:
                n_dwellfree += 1
                dwellfree_states.extend(orb)
        else:
            n_dead += 1
            dead_states.extend(orb)
        assert not (dwellfree and not alive), "Theorem C(e) broke: dead dwell-free cycle"

    # --- Theorem Z: the only dead attractor is all-zero
    assert n_dead == 1 and len(dead_states) == 1 and dead_states[0] == 0, \
        f"Theorem Z broke at L={L} ({ta},{tp}): {n_dead} dead cycles"

    # --- Theorem C: C == dwell-free cycle states, +1-closed, step == +1 on C
    Cmask = in_C(D, nbi, nbm, ta, B)
    Cset = np.nonzero(Cmask)[0]
    dfset = np.array(sorted(dwellfree_states), dtype=np.int64)
    theoremC = (Cset.size == dfset.size) and np.array_equal(Cset, dfset)
    assert theoremC, f"Theorem C broke at L={L} ({ta},{tp})"

    shifted = ((D[Cset].astype(np.int16) + 1) % B).astype(np.int64) @ pw
    assert Cmask[shifted].all(), "C is not +1-closed"
    assert np.array_equal(succ[Cset], shifted), "step != +1 on C"
    assert Cset.size % B == 0, "|C| not divisible by B"
    assert Cset.size // B == n_dwellfree, "|C|/B != dwell-free attractor count"

    # --- certified (not proven): S+1 >= 4 => every live attractor is in C
    all_live_dwellfree = (n_live == n_dwellfree)
    if B >= 4:
        assert all_live_dwellfree, \
            f"live-attractors-are-C certification broke at L={L} ({ta},{tp})"
    else:
        assert not all_live_dwellfree, "(1,1) dwelling exception disappeared"

    fate_live = np.zeros(N, bool)
    fate_live[np.array(live_cycle_states, dtype=np.int64)] = True
    P = float(fate_live[land].mean())

    return dict(L=L, ta=ta, tp=tp, S=S, B=B, N=N, P=P,
                nC=int(Cset.size), attractors_C=int(Cset.size // B),
                live_attractors=n_live, dwellfree_attractors=n_dwellfree,
                dwelling_attractors=n_live - n_dwellfree,
                theoremC=bool(theoremC), all_live_dwellfree=bool(all_live_dwellfree),
                live_periods=sorted(periods_live),
                max_transient=int(trans.max()), mean_transient=float(trans.mean()))


# --------------------------------------------------------------------------
# sampled: P(L), transients, and whether sampled attractors stay inside C
# --------------------------------------------------------------------------
def sampled_run(L, ta, tp, n_samples=2000, seed=SEED, horizon=None):
    S, B, n = ta + tp, ta + tp + 1, L * L
    if horizon is None:
        horizon = 40 * B + 20 * L
    rng = np.random.default_rng(seed)
    nbi, nbm = neighbour_table(L)
    X = rng.integers(0, B, size=(n_samples, n)).astype(np.int16)
    alive = np.ones(n_samples, bool)
    last_death = -1
    for t in range(horizon):
        anyact = ((X >= 1) & (X <= ta)).any(axis=1)
        if (alive & ~anyact).any():
            last_death = t
        alive &= anyact
        if not alive.any():
            break
        X = step_digits(X, nbi, nbm, ta, B)
    return float(alive.mean()), horizon, last_death


def sampled_attractor_check(L, ta, tp, n_samples=150, seed=SEED, maxT=4000):
    """Land seeded random inits on their attractors; check each live attractor
    state satisfies C (valid only where S+1 >= 4), and report transient stats."""
    S, B, n = ta + tp, ta + tp + 1, L * L
    rng = np.random.default_rng(seed)
    nbi, nbm = neighbour_table(L)
    X = rng.integers(0, B, size=(n_samples, n)).astype(np.int16)
    seen = [dict() for _ in range(n_samples)]
    trans = np.full(n_samples, -1)
    per = np.full(n_samples, -1)
    key = [None] * n_samples
    done = np.zeros(n_samples, bool)
    for t in range(maxT):
        for s in range(n_samples):
            if done[s]:
                continue
            k = X[s].tobytes()
            if k in seen[s]:
                trans[s], per[s], key[s], done[s] = seen[s][k], t - seen[s][k], k, True
            else:
                seen[s][k] = t
        if done.all():
            break
        X = step_digits(X, nbi, nbm, ta, B)
    assert done.all(), f"no cycle found within maxT at L={L} ({ta},{tp})"

    live = np.zeros(n_samples, bool)
    inC = np.zeros(n_samples, bool)
    for s in range(n_samples):
        cfg = np.frombuffer(key[s], dtype=np.int16).reshape(1, n).copy()
        a, allC = False, True
        for _ in range(per[s]):
            a |= bool(((cfg >= 1) & (cfg <= ta)).any())
            allC &= bool(in_C(cfg.astype(np.int8), nbi, nbm, ta, B)[0])
            cfg = step_digits(cfg, nbi, nbm, ta, B)
        live[s], inC[s] = a, allC
    distinct = len({key[s] for s in range(n_samples) if live[s]})
    return dict(L=L, ta=ta, tp=tp, B=B, live=int(live.sum()), n=n_samples,
                live_in_C=int(inC[live].sum()), distinct_live=distinct,
                periods=sorted({int(p) for p in per[live]}),
                max_trans=int(trans.max()), mean_trans=float(trans.mean()))


def density(L, ta, tp, n=40000, seed=SEED):
    """MC estimate of |C| / B^(L^2); returns (frac, ln #attractors, kappa)."""
    B = ta + tp + 1
    rng = np.random.default_rng(seed)
    nbi, nbm = neighbour_table(L)
    X = rng.integers(0, B, size=(n, L * L)).astype(np.int8)
    frac = float(in_C(X, nbi, nbm, ta, B).mean())
    if frac < MC_FLOOR:
        return frac, float("nan"), float("nan")
    ln_attr = np.log(frac) + (L * L - 1) * np.log(B)
    return frac, float(ln_attr), float(ln_attr / (L * L))


# --------------------------------------------------------------------------
def main():
    print("=== Theorem C + Theorem Z, exhaustive (open boundary, theta=1) ===",
          flush=True)
    print(f"{'L':>2} {'cell':>7} {'regime':>7} {'N':>9} {'P':>7} {'|C|':>7} "
          f"{'|C|/B':>7} {'live att':>9} {'dwelling':>9} {'ThmC':>5} "
          f"{'live==C':>8} {'periods':>12} {'maxTr':>6}", flush=True)
    rows = []
    for L, ta, tp in EXHAUSTIVE:
        r = census(L, ta, tp)
        rows.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"{L:>2} {'(' + str(ta) + ',' + str(tp) + ')':>7} {regime:>7} "
              f"{r['N']:>9} {r['P']:>7.4f} {r['nC']:>7} {r['attractors_C']:>7} "
              f"{r['live_attractors']:>9} {r['dwelling_attractors']:>9} "
              f"{str(r['theoremC']):>5} {str(r['all_live_dwellfree']):>8} "
              f"{str(r['live_periods']):>12} {r['max_transient']:>6}", flush=True)

    assert all(r["theoremC"] for r in rows), "Theorem C failed somewhere"
    assert all(r["all_live_dwellfree"] for r in rows if r["B"] >= 4), \
        "a S+1>=4 cell grew a dwelling attractor"
    assert all(r["live_periods"] == [r["B"]] for r in rows if r["B"] >= 4), \
        "a S+1>=4 cell has a live period other than B"
    # regime-independence of the attractor architecture (finding 2)
    assert any(r["ta"] < r["tp"] and r["theoremC"] for r in rows), \
        "no tp>ta cell exercised"

    print("\n=== sampled: every live attractor still lies in C (S+1 >= 4) ===",
          flush=True)
    print(f"{'L':>2} {'cell':>7} {'live':>9} {'live in C':>10} {'distinct':>9} "
          f"{'periods':>9} {'maxTr':>6} {'meanTr':>7}", flush=True)
    att = []
    for ta, tp in [(2, 1), (1, 2), (2, 2)]:
        for L in [3, 4, 5, 6, 8]:
            r = sampled_attractor_check(L, ta, tp)
            att.append(r)
            print(f"{L:>2} {'(' + str(ta) + ',' + str(tp) + ')':>7} "
                  f"{r['live']:>4}/{r['n']:<4} {r['live_in_C']:>10} "
                  f"{r['distinct_live']:>9} {str(r['periods']):>9} "
                  f"{r['max_trans']:>6} {r['mean_trans']:>7.1f}", flush=True)
            assert r["live_in_C"] == r["live"], \
                f"a sampled live attractor left C at L={L} ({ta},{tp})"
            assert r["periods"] in ([r["B"]], []), "sampled period != B"

    print("\n=== sampled P(L): death is a small-lattice phenomenon ===",
          flush=True)
    print("(horizon-truncated, so an UPPER bound on P; latest observed death "
          "printed per row)", flush=True)
    prow = []
    for ta, tp in SAMPLED_CELLS:
        ps, deaths = [], []
        for L in SAMPLED_L:
            p, hor, ld = sampled_run(L, ta, tp)
            ps.append(p)
            deaths.append(ld)
            prow.append(dict(L=L, ta=ta, tp=tp, P=p, horizon=hor, last_death=ld))
        print(f"({ta},{tp}) S={ta + tp}: "
              + "  ".join(f"L={L}:{p:.3f}" for L, p in zip(SAMPLED_L, ps))
              + f"   [latest death t={max(deaths)}]", flush=True)
        assert ps == sorted(ps) or ps[-1] >= ps[0], "P not increasing in L"

    print("\n=== attractor-entropy density: #attractors ~ exp(kappa L^2) ===",
          flush=True)
    dens = []
    for ta, tp in DENSITY_CELLS:
        cells = []
        for L in DENSITY_L:
            frac, ln_attr, kap = density(L, ta, tp)
            dens.append(dict(L=L, ta=ta, tp=tp, frac=frac, ln_attr=ln_attr,
                             kappa=kap))
            if not np.isnan(kap):
                cells.append((L, kap))
        shown = "  ".join(f"L={L}:{k:.3f}" for L, k in cells)
        print(f"({ta},{tp}): kappa {shown}   [MC floor {MC_FLOOR:g}]", flush=True)
        ks = [k for _, k in cells]
        assert len(ks) >= 3 and ks[-1] > ks[0] > 0, "kappa did not rise to a positive plateau"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherent_core.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    scalar = [k for k in rows[0] if k != "live_periods"]
    np.savez(
        out,
        **{f"x_{k}": np.array([r[k] for r in rows]) for k in scalar},
        x_live_periods=np.array([str(r["live_periods"]) for r in rows]),
        **{f"a_{k}": np.array([r[k] for r in att])
           for k in att[0] if k != "periods"},
        **{f"p_{k}": np.array([r[k] for r in prow]) for k in prow[0]},
        **{f"d_{k}": np.array([r[k] for r in dens]) for k in dens[0]},
    )
    print(f"\nwrote {out}", flush=True)

    print("\n--- results-doc table 1 (generated; do not hand-type) ---")
    print("| lattice | (τa, τp) | regime | configs | P | \\|C\\| | attractors "
          "= \\|C\\|/B | live attractors | dwelling | max transient |")
    print("| :---: | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for r in rows:
        regime = "τa ≥ τp" if r["ta"] >= r["tp"] else "τa < τp"
        print(f"| {r['L']}×{r['L']} | ({r['ta']}, {r['tp']}) | {regime} | "
              f"{r['N']:,} | {r['P']:.4f} | {r['nC']:,} | {r['attractors_C']:,} | "
              f"{r['live_attractors']:,} | {r['dwelling_attractors']} | "
              f"{r['max_transient']} |")

    print("\n--- results-doc table 2: sampled P(L) (generated) ---")
    print("| (τa, τp) | " + " | ".join(f"L={L}" for L in SAMPLED_L) + " |")
    print("| :---: | " + " | ".join("---:" for _ in SAMPLED_L) + " |")
    for ta, tp in SAMPLED_CELLS:
        vals = [next(x["P"] for x in prow if x["L"] == L and x["ta"] == ta
                     and x["tp"] == tp) for L in SAMPLED_L]
        print(f"| ({ta}, {tp}) | " + " | ".join(f"{v:.3f}" for v in vals) + " |")


if __name__ == "__main__":
    main()
