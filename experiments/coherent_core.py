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

LEMMA E (proven by hand; any graph). If C is NON-EMPTY then G has a cycle of
length l and an integer k >= 1 with  l <= k*B <= l*ta.  (Proof in the
emptiness_condition docstring below.)

  This CORRECTS a conjecture an earlier revision of this file made. That
  revision read the (1,1) anomaly as a GIRTH effect -- "S+1 >= 4, and 4 is the
  girth of the square lattice". That reading is wrong. The real obstruction is
  arithmetic, not metric: the square lattice is BIPARTITE, so every cycle
  length is even, while B = 3 at (1,1) is odd. At 2x2 the only cycle length is
  4 and no multiple of 3 lies in [4, 4], so C is empty; at 3x3 a 6-cycle exists
  and 6 = 2*3 lies in [6, 6], so C is not (|C| = 84). Girth alone predicts
  neither.

  Lemma E also accounts for the three 2x2 cells with P = 0 exactly --
  (1,3), (1,4), (1,5) -- which persistent_set_structure.md Section 4 records as
  "vacuous" without explaining them, and it correctly separates 3x3 (1,4)
  (|C| = 168, since 6 = 1*6) from 3x3 (1,3) and (1,5) (both |C| = 0).
  Necessity is asserted at 16 (lattice, cell) pairs; the converse held at all
  16 but is NOT claimed.

WHAT IS NOT PROVEN (the certified part, and the honest boundary):
  the claim "every LIVE attractor is dwell-free when S+1 >= 4" -- i.e. that C
  captures ALL live attractors and not merely the dwell-free ones. This is
  CERTIFIED here by exhaustive census at 2x2 and 3x3 and by seeded sampling to
  L = 8, not proven. Note that Lemma E does NOT settle it: Lemma E says when C
  is empty, whereas this claim says when DWELLING attractors are absent, and
  3x3 (1,1) has both a non-empty C (28 attractors) and 34 dwelling ones. A
  lattice of different girth and parity -- triangular (girth 3, non-bipartite)
  or honeycomb (girth 6, bipartite) -- is the cheap discriminator, and is not
  run here.

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
   the project's range with S+1 < 4. At 2x2 (1,1): |C| = 0 -- explained by
   LEMMA E, not by girth -- yet 2 live attractors exist (period 4). At
   3x3 (1,1): 28 dwell-free attractors from C, but 62 live attractors in total,
   34 of them dwelling. At every cell with S+1 >= 4, live attractors == C
   exactly, zero leftovers.

4. PHASE SPACE IS WIDE AND SHALLOW. Transients are O(S), not O(L^2), and do not
   grow with L: max transient 12 at 3x3 (2,2) exhaustively, and 11 at L=8 (2,2)
   sampled -- while the attractor set is an exponentially small slice of the
   space. Death, when it happens, is immediate: the latest death observed at
   any cell or size is t = 12.

5. ATTRACTOR COUNT IS EXTENSIVE. |C|/B ~ exp(kappa * L^2) with kappa positive
   at every measured point, so multistability is extensive in system size.
   READ c = ln B - kappa, NOT kappa. kappa is bounded above by ln B by
   construction (it would equal ln B if every configuration were coherent), so
   "kappa rises with L" is mostly that trivial ceiling. The informative
   quantity is the per-cell coherence COST c, and it falls monotonically:
   (2,1) 0.625 -> 0.106 over L = 2..10; (2,2) 1.090 -> 0.356 over L = 2..5;
   (3,3) 1.080 -> 0.278 over L = 2..6. NOT CONVERGED: c is still falling at
   the largest affordable L, so no limit is claimed for either c or kappa.
   Past L ~ 5 essentially every sampled initial condition lands on its own
   attractor.

   Every quoted point is backed by an exact |C| from the census or the DP, or
   by >= 300 Monte-Carlo hits with draws escalated up to 4,000,000; points that
   reach none of these are DROPPED and listed, not quoted. An earlier revision
   used a bare fraction floor of 5e-5 at n = 40,000 and thereby quoted kappa
   values resting on ONE to THREE hits -- (2,2) L=6 rested on 1 hit, (2,1)
   L=10 on 3.

7. THE TRANSFER-MATRIX DP RETIRES SAMPLING UP TO L = 5 (not entirely -- see
   below). Carrying the last 2L cells as DP state makes |C| exactly computable
   in O(L^2 * B^(2L+1)); it reproduces all 12 census counts EXACTLY (an
   independent check of both, since the census counts cycles and the DP counts
   a static predicate). Reach, set by a state-space cap and an int64 overflow
   guard: L <= 5 at B <= 5, L <= 4 at B = 7. MC still carries (2,1) L = 6, 8,
   10 and (3,3) L = 5, 6, so the earlier "would retire the sampling entirely"
   was too strong.

   Where the two overlap the DP VALIDATES the MC: every deviation is under
   1 SE (max 0.98 over 8 comparisons), so the surviving MC points are sound
   estimates -- the first revision's numbers were right, merely under-powered
   in the tail. Two previously-MC points are now exact and moved:
   (2,2) L=5 from kappa 1.2532 (2,709 hits) to 1.2538 exact, and (3,3) L=3
   from 1.3739 (1,627 hits) to 1.3721 exact.

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
DENSITY_CELLS = [(1, 1), (2, 1), (1, 2), (2, 2), (3, 3)]
DENSITY_L = [2, 3, 4, 5, 6, 8, 10]
# cells for LEMMA E, including the three 2x2 cells with P = 0 exactly
LEMMA_CELLS = [(1, 1), (2, 1), (1, 2), (2, 2), (1, 3), (1, 4), (1, 5), (2, 5)]
SEED = 20260816
# A kappa point is quoted only if backed by >= MIN_HITS configurations landing
# in C (or by an exact count from the census). An earlier revision used a bare
# fraction floor of 5e-5, which at n = 40,000 admitted points resting on ONE to
# THREE hits; those are not measurements and are now dropped.
MIN_HITS = 300
MC_BUDGET = (40_000, 400_000, 4_000_000)


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


def in_C(D, nbi, nbm, ta, B, max_elems=8_000_000):
    """Membership of the coherent set C. Chunked by ELEMENT count (rows x cells
    x 4), not by row count, so peak memory stays bounded as L grows."""
    n_cells = D.shape[1]
    chunk = max(1024, max_elems // (4 * n_cells))
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


def cycle_lengths(L):
    """The set of cycle lengths of the LxL open von Neumann grid graph."""
    nbi, nbm = neighbour_table(L)
    adj = {i: [int(j) for j, m in zip(nbi[i], nbm[i]) if m] for i in range(L * L)}
    lengths = set()

    def walk(start, node, visited, depth):
        for nxt in adj[node]:
            if nxt == start and depth >= 3:
                lengths.add(depth)
            elif nxt > start and nxt not in visited:
                walk(start, nxt, visited | {nxt}, depth + 1)

    for s in range(L * L):
        walk(s, s, {s}, 1)
    return lengths


def emptiness_condition(L, ta, tp, lengths):
    """LEMMA E (proven; any graph): if C is non-empty then G has a cycle of
    length l and an integer k >= 1 with  l <= k*B <= l*tau_a.

    Proof. Take c in C. Each cell i picks a neighbour sigma(i) whose lag
    (c_sigma(i) - c_i) mod B lies in [1..ta]. sigma is a self-map of a finite
    set, so its functional digraph contains a cycle of DISTINCT nodes; since
    sigma(i) is always a neighbour, that is a cycle of G, of some length l.
    Summing the lags once around it returns to the starting phase, so the sum
    is 0 mod B; the sum also lies in [l, l*ta] and is at least l >= 3 > 0, so
    it equals k*B for some k >= 1. Hence l <= k*B <= l*ta. []

    This is what actually explains the 2x2 (1,1) anomaly, and it is NOT the
    girth: the square lattice is bipartite, so every cycle length is even,
    while B = 3 is odd. At 2x2 the only cycle length is 4 and no multiple of 3
    lies in [4, 4], so C is empty. At 3x3 a 6-cycle exists and 6 = 2*3 lies in
    [6, 6], so C is not (|C| = 84).
    """
    B = ta + tp + 1
    for l in lengths:
        kmin = -(-l // B)                    # ceil(l / B)
        if kmin >= 1 and kmin * B <= l * ta:
            return True
    return False


def exact_C_count(L, ta, tp, max_states=12_000_000):
    """EXACT |C| by a cell-by-cell transfer matrix, no sampling.

    Cell m's coherence test reads m-L, m-1, m+1, m+L, so m can be finalised the
    moment cell m+L is placed. Carrying the last 2L cells as the DP state is
    therefore sufficient, and transitions branch only B ways (a row-by-row DP
    would branch B^L). State space is B^(2L); work is O(L^2 * B^(2L+1)).

    Two guards, both returning (None, reason) rather than a wrong number:
      - state space B^(2L) over max_states;
      - the count ceiling B^(L^2) overflowing int64. (The handoff doc's
        encoding trap: size the guard by the lattice constant, and assert it.)
    """
    B = ta + tp + 1
    n, W = L * L, 2 * L
    if B ** W > max_states:
        return None, f"state space {B}^{W} > {max_states:,}"
    if B ** n > np.iinfo(np.int64).max:
        return None, f"count ceiling {B}^{n} overflows int64"

    NS = B ** W
    s = np.arange(NS, dtype=np.int64)
    dig = lambda t: ((s // B ** (W - 1 - t)) % B).astype(np.int16)
    ok = lambda g: (g >= 1) & (g <= ta)

    # prefix: cells 0..W-1 free, then verify row 0 (whose neighbours are all in)
    counts = np.ones(NS, dtype=np.int64)
    for m in range(L):
        sd, acc = dig(m), np.zeros(NS, bool)
        for nb, exists in ((m + L, True), (m - 1, m % L > 0),
                           (m + 1, m % L < L - 1)):
            if exists:
                acc |= ok((dig(nb) - sd) % B)
        counts *= acc

    # sliding phase: verifying cell m = k - L, which always sits at slot L
    self_d = dig(L)
    up_ok, left_ok = ok((dig(0) - self_d) % B), ok((dig(L - 1) - self_d) % B)
    right_ok = ok((dig(L + 1) - self_d) % B)
    del s
    buf = np.empty(NS, dtype=np.int64)
    for k in range(W, n):
        r, c = divmod(k - L, L)
        static = np.zeros(NS, bool)
        if r > 0:
            static |= up_ok
        if c > 0:
            static |= left_ok
        if c < L - 1:
            static |= right_ok
        new = np.zeros(NS, dtype=np.int64)
        new2 = new.reshape(B ** (W - 1), B)
        for a in range(B):
            v = static | ok((a - self_d) % B) if r < L - 1 else static
            np.multiply(counts, v, out=buf)
            new2[:, a] = buf.reshape(B, B ** (W - 1)).sum(axis=0)
        counts = new

    # final: the last row has no down-neighbour, so verify it from the state
    s = np.arange(NS, dtype=np.int64)
    dig = lambda t: ((s // B ** (W - 1 - t)) % B).astype(np.int16)
    for m in range(n - L, n):
        slot, c = m - (n - W), m % L
        sd, acc = dig(slot), np.zeros(NS, bool)
        for nb_slot, exists in ((slot - L, True), (slot - 1, c > 0),
                                (slot + 1, c < L - 1)):
            if exists:
                acc |= ok((dig(nb_slot) - sd) % B)
        counts *= acc
    total = int(counts.sum())
    assert total >= 0, "int64 overflow slipped past the guard"
    return total, None


def _mc_hits(L, ta, tp, B, n, seed, draw_chunk=200_000):
    """Hits of C among n uniform-random CONFIGURATIONS, drawn in chunks so peak
    memory is bounded. Deterministic given (seed, n, draw_chunk)."""
    rng = np.random.default_rng(seed)
    nbi, nbm = neighbour_table(L)
    hits, drawn = 0, 0
    while drawn < n:
        k = min(draw_chunk, n - drawn)
        X = rng.integers(0, B, size=(k, L * L)).astype(np.int8)
        hits += int(in_C(X, nbi, nbm, ta, B).sum())
        drawn += k
    return hits


def _pack(L, ta, tp, nC, source, n=0, hits=None, se=0.0):
    B = ta + tp + 1
    lnB = np.log(B)
    frac = nC / B ** (L * L)
    kappa = (np.log(frac) + (L * L - 1) * lnB) / (L * L)
    return dict(L=L, ta=ta, tp=tp, source=source, n=n,
                hits=nC if hits is None else hits, frac=float(frac),
                kappa=float(kappa), cost=float(lnB - kappa), se=float(se))


def density(L, ta, tp, seed=SEED, exact_nC=None):
    """Attractor-entropy density kappa = ln(#attractors) / L^2, where
    #attractors = |C|/B and |C| = frac * B^(L^2).

    Three sources, in descending order of authority:
      1. the exhaustive census, where it covers this (L, ta, tp);
      2. the exact transfer-matrix DP (exact_C_count), where it is affordable;
      3. Monte-Carlo over uniformly random CONFIGURATIONS, escalating n until
         at least MIN_HITS land in C.
    A point that reaches none of these is REPORTED AS DROPPED rather than
    quoted -- kappa from a handful of hits is not a measurement.

    Returns the string "EMPTY" when |C| = 0 exactly (no dwell-free attractor
    exists, so kappa is undefined -- a FACT, not a measurement failure), and
    None when the point is dropped as under-powered. The two are different and
    the driver reports them separately.
    """
    B = ta + tp + 1
    if exact_nC is not None:
        return "EMPTY" if exact_nC == 0 else _pack(L, ta, tp, exact_nC, "census")
    dp, _why = exact_C_count(L, ta, tp)
    if dp is not None:
        return "EMPTY" if dp == 0 else _pack(L, ta, tp, dp, "DP")
    for n in MC_BUDGET:
        hits = _mc_hits(L, ta, tp, B, n, seed)
        if hits >= MIN_HITS:
            frac = hits / n
            # delta-method SE on ln(frac), divided through by L^2
            se = np.sqrt((1 - frac) / (frac * n)) / (L * L)
            return _pack(L, ta, tp, frac * B ** (L * L), "MC", n=n, hits=hits,
                         se=se)
    return None


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

    print("\n=== transfer-matrix DP: exact |C|, validated against the census ===",
          flush=True)
    print(f"{'lattice':>8} {'cell':>7} {'census |C|':>12} {'DP |C|':>12} {'match':>6}",
          flush=True)
    dpval = []
    for r in rows:
        got, why = exact_C_count(r["L"], r["ta"], r["tp"])
        assert got is not None, f"DP unaffordable at a census row: {why}"
        match = (got == r["nC"])
        dpval.append(dict(L=r["L"], ta=r["ta"], tp=r["tp"], census=r["nC"], dp=got,
                          match=int(match)))
        print(f"{str(r['L']) + 'x' + str(r['L']):>8} "
              f"{'(' + str(r['ta']) + ',' + str(r['tp']) + ')':>7} "
              f"{r['nC']:>12,} {got:>12,} {str(match):>6}", flush=True)
        assert match, f"DP disagrees with the census at L={r['L']} ({r['ta']},{r['tp']})"

    print("\n=== LEMMA E (emptiness): C non-empty => some cycle length l has "
          "l <= kB <= l*ta ===", flush=True)
    print(f"{'lattice':>8} {'cell':>7} {'B':>3} {'cycle lengths':>22} "
          f"{'lemma allows C':>15} {'|C|':>10} {'consistent':>11}", flush=True)
    lem = []
    for L in sorted({r["L"] for r in rows}):
        lengths = cycle_lengths(L)
        for ta, tp in LEMMA_CELLS:
            B = ta + tp + 1
            allowed = emptiness_condition(L, ta, tp, lengths)
            nC, why = exact_C_count(L, ta, tp)
            assert nC is not None, why
            # necessity: |C| > 0 forces the condition. (The converse is not
            # claimed -- it is reported, not asserted.)
            assert not (nC > 0 and not allowed), \
                f"LEMMA E violated at L={L} ({ta},{tp})"
            lem.append(dict(L=L, ta=ta, tp=tp, allowed=int(allowed), nC=nC))
            print(f"{str(L) + 'x' + str(L):>8} {'(' + str(ta) + ',' + str(tp) + ')':>7} "
                  f"{B:>3} {str(sorted(lengths)):>22} {str(allowed):>15} "
                  f"{nC:>10,} {str((nC > 0) == allowed):>11}", flush=True)
    print("Necessity is ASSERTED. Sufficiency is not claimed; the "
          "'consistent' column reports whether the converse also held.",
          flush=True)

    print("\n=== attractor-entropy density: #attractors ~ exp(kappa L^2) ===",
          flush=True)
    print(f"kappa <= ln B by construction, so the informative quantity is the "
          f"per-cell coherence COST c = ln B - kappa.", flush=True)
    print(f"Exact |C| used where the census supplies it; otherwise MC over "
          f"CONFIGURATIONS, escalated to >= {MIN_HITS} hits or dropped.",
          flush=True)
    print(f"{'cell':>7} {'lnB':>6} {'L':>3} {'source':>7} {'draws':>9} {'hits':>8} "
          f"{'frac':>11} {'kappa':>7} {'c=lnB-kappa':>12} {'SE':>7}", flush=True)
    exact_nC = {(r["L"], r["ta"], r["tp"]): r["nC"] for r in rows}
    dens, dropped, empty, xcheck = [], [], [], []
    for ta, tp in DENSITY_CELLS:
        lnB = np.log(ta + tp + 1)
        for L in DENSITY_L:
            d = density(L, ta, tp, exact_nC=exact_nC.get((L, ta, tp)))
            if d == "EMPTY":
                empty.append((ta, tp, L))
                print(f"{'(' + str(ta) + ',' + str(tp) + ')':>7} {lnB:>6.3f} "
                      f"{L:>3} {'EMPTY':>7} {'—':>9} {0:>8} "
                      f"{'0':>11} {'n/a':>7} {'n/a':>12} {'—':>7}", flush=True)
                continue
            if d is None:
                dropped.append((ta, tp, L))
                print(f"{'(' + str(ta) + ',' + str(tp) + ')':>7} {lnB:>6.3f} "
                      f"{L:>3} {'DROPPED':>7} {'—':>9} {'<' + str(MIN_HITS):>8} "
                      f"{'—':>11} {'—':>7} {'—':>12} {'—':>7}", flush=True)
                continue
            dens.append(d)
            print(f"{'(' + str(ta) + ',' + str(tp) + ')':>7} {lnB:>6.3f} {L:>3} "
                  f"{d['source']:>7} {d['n'] or '—':>9} {d['hits']:>8,} "
                  f"{d['frac']:>11.3e} {d['kappa']:>7.4f} {d['cost']:>12.4f} "
                  f"{d['se']:>7.4f}", flush=True)
        got = [d for d in dens if (d["ta"], d["tp"]) == (ta, tp)]
        assert len(got) >= 3, f"too few usable density points at ({ta},{tp})"
        # every quoted point is backed by real evidence
        assert all(d["source"] != "MC" or d["hits"] >= MIN_HITS for d in got), \
            "a kappa point slipped through under-powered"
        # kappa is strictly below its ceiling ln B, and the cost is falling in L
        assert all(0 < d["kappa"] < lnB for d in got), "kappa outside (0, ln B)"
        costs = [d["cost"] for d in got]
        assert costs == sorted(costs, reverse=True), \
            f"coherence cost not monotone decreasing at ({ta},{tp})"
        assert costs[-1] > 0, "coherence cost hit zero"
    if empty:
        print(f"\n|C| = 0 exactly (no dwell-free attractor; kappa undefined) at "
              + ", ".join(f"({ta},{tp}) L={L}" for ta, tp, L in empty)
              + " -- explained by LEMMA E above, not a measurement failure.",
              flush=True)
    if dropped:
        print(f"dropped {len(dropped)} under-powered point(s): "
              + ", ".join(f"({ta},{tp}) L={L}" for ta, tp, L in dropped),
              flush=True)

    print("\n=== validating the surviving MC estimator against the exact DP ===",
          flush=True)
    print("(the MC points that remain sit beyond the DP's reach, so their "
          "estimator is checked where the two overlap)", flush=True)
    print(f"{'cell':>7} {'L':>3} {'MC frac':>12} {'exact frac':>12} "
          f"{'|dev| in SE':>12}", flush=True)
    for ta, tp in DENSITY_CELLS:
        for L in DENSITY_L:
            mc = None
            for n in MC_BUDGET:
                h = _mc_hits(L, ta, tp, ta + tp + 1, n, SEED)
                if h >= MIN_HITS:
                    mc = (h / n, np.sqrt((1 - h / n) / (h / n * n)))
                    break
            ex, _ = exact_C_count(L, ta, tp)
            if mc is None or ex is None or ex == 0:
                continue
            ex_frac = ex / (ta + tp + 1) ** (L * L)
            dev = abs(np.log(mc[0]) - np.log(ex_frac)) / mc[1]
            xcheck.append(dict(L=L, ta=ta, tp=tp, mc=mc[0], exact=ex_frac,
                               dev_se=float(dev)))
            print(f"{'(' + str(ta) + ',' + str(tp) + ')':>7} {L:>3} "
                  f"{mc[0]:>12.5e} {ex_frac:>12.5e} {dev:>12.2f}", flush=True)
            assert dev < 5.0, \
                f"MC estimator off by {dev:.1f} SE at L={L} ({ta},{tp})"

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
        **{f"v_{k}": np.array([r[k] for r in dpval]) for k in dpval[0]},
        **{f"e_{k}": np.array([r[k] for r in lem]) for k in lem[0]},
        **{f"x_c_{k}": np.array([r[k] for r in xcheck]) for k in xcheck[0]},
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

    print("\n--- results-doc table 3: attractor-entropy density (generated) ---")
    print("| (τa, τp) | ln B | L | source | draws | hits | \\|C\\|/B^(L²) | κ | "
          "c = ln B − κ | SE(κ) |")
    print("| :---: | ---: | :---: | :---: | ---: | ---: | ---: | ---: | ---: | ---: |")
    for d in dens:
        lnB = np.log(d["ta"] + d["tp"] + 1)
        draws = f"{d['n']:,}" if d["source"] == "MC" else "—"
        ev = f"{d['hits']:,}" + (" hits" if d["source"] == "MC" else "")
        print(f"| ({d['ta']}, {d['tp']}) | {lnB:.3f} | {d['L']} | {d['source']} | "
              f"{draws} | {ev} | {d['frac']:.3e} | {d['kappa']:.4f} | "
              f"{d['cost']:.4f} | {d['se']:.4f} |")
    if empty:
        print("\n`|C| = 0` exactly (κ undefined, per Lemma E): "
              + ", ".join(f"({ta}, {tp}) L={L}" for ta, tp, L in empty))
    if dropped:
        print("\nDropped as under-powered (< "
              f"{MIN_HITS} hits at {MC_BUDGET[-1]:,} draws): "
              + ", ".join(f"({ta}, {tp}) L={L}" for ta, tp, L in dropped))

    print("\n--- results-doc table 2: sampled P(L) (generated) ---")
    print("| (τa, τp) | " + " | ".join(f"L={L}" for L in SAMPLED_L) + " |")
    print("| :---: | " + " | ".join("---:" for _ in SAMPLED_L) + " |")
    for ta, tp in SAMPLED_CELLS:
        vals = [next(x["P"] for x in prow if x["L"] == L and x["ta"] == ta
                     and x["tp"] == tp) for L in SAMPLED_L]
        print(f"| ({ta}, {tp}) | " + " | ".join(f"{v:.3f}" for v in vals) + " |")


if __name__ == "__main__":
    main()
