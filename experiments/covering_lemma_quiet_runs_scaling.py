"""Scaling the quiet-run bound: the 2x2 closed form is FALSIFIED at 3x3, H1
itself is FALSIFIED (8 violations), and the covering lemma survives anyway --
its witnesses there use only OLD cells. Quiet-run age is not the mechanism.

Scaling test of covering_lemma_quiet_runs.py, which found max quiet run =
min(ceil(S/2)+1, tau_p+2) exactly at 16 2x2 cells and H1 (every ceiling hold
has a dwelling cell with quiet run Q < S) with zero violations.

FINDINGS:

1. THE CLOSED FORM DOES NOT LIFT (the seventh 2x2 artifact in this program).
   Exhaustive at 3x3 over all live configurations:
     cell   S   maxQ(3x3)   2x2 form
     (1,1)  2      6           2
     (2,1)  3      7           3
     (2,2)  4      9           3
     (3,1)  4      7           3
   Quiet runs at 3x3 EXCEED S rather than falling below it. Sampled 4x4 (400
   seeded orbits per condition, rng=default_rng(4444)) agrees: open (2,1) 5,
   (2,2) 8; torus (2,1) 4, (2,2) 7 -- all > S. The 2x2 bound was a
   size artifact: on a 4-cycle the wave returns within half a cycle, so no
   cell can wait long. Larger media give a quiet cell somewhere for the wave
   to be far from.

2. H1 IS FALSIFIED AT 3x3 -- 8 violations out of 58,588 ceiling holds at
   (2,1) (ceiling holds counted only where u has a dwelling cell). Every
   violation has min dwelling quiet run EXACTLY S (=3); the min-Q distribution
   at ceiling is {0: 57780, 1: 716, 2: 84, 3: 8}, so the failures sit precisely
   at the boundary the hypothesis draws. This falsification was only visible
   because maxQ > S here: at 2x2 the bound made H1 near-vacuous, which is why
   it read as universally true.

3. THE COVERING LEMMA SURVIVES, VIA THE OPPOSITE MECHANISM. All 8 violating
   pairs ARE witnessed (8/8), and in every case the witnessing subset uses
   only cells with Q >= S -- old cells. So "young cells supply the witness"
   (the 2x2 finding) is false in general; a witness can be built from cells
   that have waited longer than a full cycle. Example violation at 3x3 (2,1):
   u = [3,1,0,1,0,0,0,0,0] with Q = [0,1,2,1,2,3,2,3,3]; its dwelling cells
   are two edges and a corner, all at Q = 3 = S, and the witness is among them.

4. DEGREE ORDERING IS REAL AND SURVIVES SCALING, unlike the closed form.
   Per-role maxima at 3x3 order strictly by degree at every cell tested:
   corner (deg 2) > edge (deg 3) > centre (deg 4) -- e.g. (2,2): 9 / 6 / 3.
   Low-degree cells can stay quiet longest, which is the same structural fact
   behind the per-capita witness concentration in coherence_covering_lemma.py.
   At 4x4 the ordering is visible but not strict under sampling (open (2,2):
   8 / 8 / 8), and the torus -- where every cell has degree 4 -- shows no
   boundary structure to order, as expected.

CONSEQUENCE FOR THE PROOF. The hand obligation from the 2x2 experiment
("prove max Q <= min(ceil(S/2)+1, tau_p+2)") is WITHDRAWN: the statement is
false beyond 2x2. Quiet-run age is not the covering mechanism, and any proof
of the covering lemma must permit old-cell witnesses. What does transfer is
the swap law (proven, universal) plus the degree ordering; the open question
becomes which dwelling subsets admit a consistent older reading, independent
of how long they have waited.

House Rules Compliance: 3x3 exhaustive/deterministic; 4x4 seeded
(default_rng(4444), threaded explicitly, no global RNG); asserts the
falsification of both the closed form and H1, the 8/8 old-cell witnessing, the
min-Q distribution, and the 3x3 degree ordering; aborts on regression. Caveats:
4x4 numbers are sampled lower bounds on maxQ, not exhaustive maxima; per-role
4x4 ordering is not asserted because sampling does not resolve it.
Output: result/topology/covering_lemma_quiet_runs_scaling.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import combinations

NB3 = [tuple(j for j in range(9) if abs(j // 3 - i // 3) + abs(j % 3 - i % 3) == 1)
       for i in range(9)]
ROLE3 = {i: {2: "corner", 3: "edge", 4: "centre"}[len(NB3[i])] for i in range(9)}
NB4O = [tuple(j for j in range(16) if abs(j // 4 - i // 4) + abs(j % 4 - i % 4) == 1)
        for i in range(16)]
NB4T = [tuple(((i // 4 + dr) % 4) * 4 + ((i % 4 + dc) % 4)
              for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]) for i in range(16)]
ROLE4 = {i: {2: "corner", 3: "edge", 4: "centre"}[len(NB4O[i])] for i in range(16)}
CHUNK = 1_000_000
SEED_4X4 = 4444

CELLS_3 = [(1, 1), (2, 1), (2, 2), (3, 1)]
Q3_EXPECT = {(1, 1): 6, (2, 1): 7, (2, 2): 9, (3, 1): 7}
ROLE3_EXPECT = {(1, 1): (6, 4, 2), (2, 1): (7, 5, 3),
                (2, 2): (9, 6, 3), (3, 1): (7, 5, 3)}
H1_3X3 = dict(cell=(2, 1), n_ceil=58_588, ok=58_580, viol=8,
              minq={0: 57780, 1: 716, 2: 84, 3: 8})
CELLS_4 = [(2, 1), (2, 2)]
Q4_EXPECT = {("open", 2, 1): 5, ("open", 2, 2): 8,
             ("torus", 2, 1): 4, ("torus", 2, 2): 7}


def succ_and_live_3x3(ta, tp):
    B = ta + tp + 1
    N = B ** 9
    pw = (B ** np.arange(9)).astype(np.int64)
    f = np.empty(N, dtype=np.int64)
    for lo in range(0, N, CHUNK):
        hi = min(lo + CHUNK, N)
        dig = ((np.arange(lo, hi, dtype=np.int64)[:, None] // pw[None, :]) % B).astype(np.int8)
        act = (dig >= 1) & (dig <= ta)
        nxt = dig.astype(np.int64)
        adv = dig > 0
        nxt[adv] = (dig[adv].astype(np.int64) + 1) % B
        for i in range(9):
            na = np.zeros(hi - lo, dtype=np.int8)
            for j in NB3[i]:
                na += act[:, j]
            m = dig[:, i] == 0
            nxt[m, i] = (na[m] >= 1)
        f[lo:hi] = nxt @ pw
    dead = f == 0
    dead[0] = True
    while True:
        nd = dead | dead[f]
        if (nd == dead).all():
            break
        dead = nd
    return f, np.nonzero(~dead)[0], B, pw


def maxq_3x3(ta, tp, chunk=400_000):
    S = ta + tp
    f, live, B, pw = succ_and_live_3x3(ta, tp)
    T = 4 * S + 20
    glob = 0
    per_role = {"corner": 0, "edge": 0, "centre": 0}
    for lo in range(0, len(live), chunk):
        V = live[lo:lo + chunk].copy()
        Q = np.zeros((len(V), 9), dtype=np.int8)
        for _ in range(T):
            for i in range(9):
                Q[:, i] = np.where(((V // pw[i]) % B) == 0, Q[:, i] + 1, 0)
            for i in range(9):
                m = int(Q[:, i].max())
                if m > per_role[ROLE3[i]]:
                    per_role[ROLE3[i]] = m
            glob = max(glob, int(Q.max()))
            V = f[V]
        del Q
    return glob, per_role, len(live)


def h1_and_witnesses_3x3(ta=2, tp=1):
    """H1 audit at 3x3 + whether violations are witnessed, and by which cells."""
    S = ta + tp
    f, live, B, pw = succ_and_live_3x3(ta, tp)
    N = B ** 9
    idxall = np.arange(N, dtype=np.int64)
    sh = np.zeros(N, dtype=np.int64)
    for k in range(9):
        sh += (((idxall // pw[k]) % B + 1) % B) * pw[k]

    depth = {}
    V, U = live.copy(), sh[live]
    for c in (V * N + U).tolist():
        depth[c] = 0
    while True:
        nV, nU = f[V], f[U]
        code = nV * N + nU
        mask = np.fromiter((c not in depth for c in code.tolist()),
                           dtype=bool, count=len(code))
        if not mask.any():
            break
        d = max(depth.values()) + 1 if depth else 1
        for c in code[mask].tolist():
            depth[c] = d
        V, U = nV[mask], nU[mask]
    codes = np.fromiter(depth.keys(), dtype=np.int64, count=len(depth))
    ags = np.fromiter(depth.values(), dtype=np.int16, count=len(depth))
    order = np.argsort(codes)
    codes_s, ages_s = codes[order], ags[order]

    def age_of(Vv, Uu):
        c = Vv * N + Uu
        pos = np.clip(np.searchsorted(codes_s, c), 0, len(codes_s) - 1)
        hit = codes_s[pos] == c
        out = np.full(len(c), -1, dtype=np.int16)
        out[hit] = ages_s[pos[hit]]
        return out

    n_ceil = ok = viol = 0
    minq = {}
    viol_pairs = []
    Vv, Uu = live.copy(), sh[live]
    Q = np.zeros((len(Vv), 9), dtype=np.int8)
    for _ in range(4 * S + 20):
        nV, nU = f[Vv], f[Uu]
        ceil = (age_of(Vv, Uu) == S) & (age_of(nV, nU) == S)
        if ceil.any():
            ui = np.nonzero(ceil)[0]
            ud = np.stack([((Uu[ui] // pw[i]) % B) for i in range(9)], axis=1)
            nud = np.stack([((nU[ui] // pw[i]) % B) for i in range(9)], axis=1)
            dwell = (ud == 0) & (nud == 0)
            Qs = Q[ui]
            young = dwell & (Qs < S)
            hasd = dwell.any(axis=1)
            n_ceil += int(hasd.sum())
            ok += int((young.any(axis=1) & hasd).sum())
            bad = hasd & ~young.any(axis=1)
            viol += int(bad.sum())
            mq = np.where(dwell, Qs, 127).min(axis=1)
            for x in mq[hasd].tolist():
                minq[x] = minq.get(x, 0) + 1
            for j in np.nonzero(bad)[0]:
                viol_pairs.append((int(Vv[ui[j]]), int(Uu[ui[j]]), Qs[j].copy()))
        for i in range(9):
            Q[:, i] = np.where(((Uu // pw[i]) % B) == 0, Q[:, i] + 1, 0)
        Vv, Uu = nV, nU

    def digs(x):
        return [(int(x) // int(pw[k])) % B for k in range(9)]

    wit_any = wit_old = 0
    for v_, u_, Qv in viol_pairs:
        ud, nud = digs(u_), digs(int(f[u_]))
        dws = [i for i in range(9) if ud[i] == 0 and nud[i] == 0]
        found = None
        for k in range(1, len(dws) + 1):
            for su in combinations(dws, k):
                u2 = u_ + sum(S * int(pw[i]) for i in su)
                if int(f[int(u2)]) != int(f[u_]):
                    continue
                if age_of(np.array([v_]), np.array([int(u2)]))[0] == S - 1:
                    found = su
                    break
            if found:
                break
        if found is not None:
            wit_any += 1
            if all(Qv[i] >= S for i in found):
                wit_old += 1
    return dict(n_ceil=n_ceil, ok=ok, viol=viol, minq=minq,
                wit_any=wit_any, wit_old=wit_old, n_viol=len(viol_pairs))


def step4(cfg, ta, tp, NB):
    B = ta + tp + 1
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in NB[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % B for i, s in enumerate(cfg))


def maxq_4x4(ta, tp, NB, n_trials=400, T=80, seed=SEED_4X4):
    B = ta + tp + 1
    rng = np.random.default_rng(seed)
    glob = 0
    per_role = {"corner": 0, "edge": 0, "centre": 0}
    n_live = 0
    for _ in range(n_trials):
        cur = tuple(int(x) for x in rng.integers(0, B, size=16))
        Q = [0] * 16
        g = 0
        loc = {"corner": 0, "edge": 0, "centre": 0}
        alive = True
        for _t in range(T):
            if not any(cur):
                alive = False
                break
            Q = [Q[i] + 1 if cur[i] == 0 else 0 for i in range(16)]
            g = max(g, max(Q))
            for i in range(16):
                loc[ROLE4[i]] = max(loc[ROLE4[i]], Q[i])
            cur = step4(cur, ta, tp, NB)
        if alive:
            n_live += 1
            glob = max(glob, g)
            for k in per_role:
                per_role[k] = max(per_role[k], loc[k])
    return glob, per_role, n_live


def main():
    print("=== 1) the 2x2 closed form at 3x3 (exhaustive) ===")
    q3 = []
    for ta, tp in CELLS_3:
        S = ta + tp
        g, pr, nlive = maxq_3x3(ta, tp)
        form2 = min(-(-S // 2) + 1, tp + 2)
        q3.append(dict(ta=ta, tp=tp, maxQ=g, form2=form2,
                       corner=pr["corner"], edge=pr["edge"], centre=pr["centre"]))
        print(f"({ta},{tp}) S={S}: maxQ={g}  2x2-form={form2}  maxQ>S: {g > S}  "
              f"per-role {pr['corner']}/{pr['edge']}/{pr['centre']}  live={nlive:,}",
              flush=True)
        assert g == Q3_EXPECT[(ta, tp)], f"3x3 quiet run changed at ({ta},{tp})"
        assert g > form2, f"2x2 form must be FALSIFIED at 3x3 ({ta},{tp})"
        assert g > S, f"3x3 quiet runs must exceed S at ({ta},{tp})"
        exp = ROLE3_EXPECT[(ta, tp)]
        assert (pr["corner"], pr["edge"], pr["centre"]) == exp, \
            f"per-role maxima changed at ({ta},{tp})"
        assert pr["corner"] > pr["edge"] > pr["centre"], \
            f"degree ordering broke at ({ta},{tp})"

    print("\n=== 2) H1 at 3x3 (2,1), and are the violations witnessed? ===")
    h = h1_and_witnesses_3x3()
    print(f"ceiling holds with a dwelling u-cell: {h['n_ceil']:,}", flush=True)
    print(f"  H1 holds: {h['ok']:,}   VIOLATIONS: {h['viol']}")
    print(f"  min-quiet-run distribution at ceiling: {dict(sorted(h['minq'].items()))}")
    print(f"  violations witnessed: {h['wit_any']}/{h['n_viol']}   "
          f"witness uses only OLD cells: {h['wit_old']}/{h['n_viol']}")
    assert h["n_ceil"] == H1_3X3["n_ceil"], "3x3 ceiling count changed"
    assert h["viol"] == H1_3X3["viol"] > 0, "H1 must be FALSIFIED at 3x3"
    assert dict(sorted(h["minq"].items())) == H1_3X3["minq"], "min-Q census changed"
    assert h["wit_any"] == h["n_viol"], "covering lemma must survive"
    assert h["wit_old"] == h["n_viol"], "violations must be OLD-cell witnessed"

    print("\n=== 3) 4x4 sampled (seeded), open vs torus ===")
    q4 = []
    for name, NB in [("open", NB4O), ("torus", NB4T)]:
        for ta, tp in CELLS_4:
            S = ta + tp
            g, pr, nl = maxq_4x4(ta, tp, NB)
            q4.append(dict(lat=name, ta=ta, tp=tp, maxQ=g, n_live=nl,
                           corner=pr["corner"], edge=pr["edge"], centre=pr["centre"]))
            print(f"4x4 {name:>5} ({ta},{tp}) S={S}: maxQ={g} (>S: {g > S})  "
                  f"per-role {pr['corner']}/{pr['edge']}/{pr['centre']}  "
                  f"live orbits={nl}/400", flush=True)
            assert g == Q4_EXPECT[(name, ta, tp)], f"4x4 maxQ changed ({name},{ta},{tp})"
            assert g > S, f"4x4 quiet runs must exceed S ({name},{ta},{tp})"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "covering_lemma_quiet_runs_scaling.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"q3_{k}": np.array([r[k] for r in q3]) for k in q3[0]},
             h1_n_ceil=np.array([h["n_ceil"]]), h1_ok=np.array([h["ok"]]),
             h1_viol=np.array([h["viol"]]),
             h1_minq_keys=np.array(sorted(h["minq"])),
             h1_minq_vals=np.array([h["minq"][k] for k in sorted(h["minq"])]),
             h1_wit_any=np.array([h["wit_any"]]), h1_wit_old=np.array([h["wit_old"]]),
             q4_lat=np.array([r["lat"] for r in q4]),
             **{f"q4_{k}": np.array([r[k] for r in q4])
                for k in ("ta", "tp", "maxQ", "n_live", "corner", "edge", "centre")},
             seed_4x4=np.array([SEED_4X4]))
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| lattice | cell | S | max quiet run | 2×2 form | > S? | corner/edge/centre |")
    print("| :---: | :---: | ---: | ---: | ---: | :---: | :---: |")
    for r in q3:
        print(f"| 3×3 (exh.) | ({r['ta']}, {r['tp']}) | {r['ta']+r['tp']} | {r['maxQ']} | "
              f"{r['form2']} | yes | {r['corner']}/{r['edge']}/{r['centre']} |")
    for r in q4:
        S = r["ta"] + r["tp"]
        print(f"| 4×4 {r['lat']} (samp.) | ({r['ta']}, {r['tp']}) | {S} | {r['maxQ']} | "
              f"{min(-(-S//2)+1, r['tp']+2)} | yes | "
              f"{r['corner']}/{r['edge']}/{r['centre']} |")


if __name__ == "__main__":
    main()
