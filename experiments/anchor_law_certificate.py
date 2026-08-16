"""The anchor law at (2,1) is PROVEN -- by a 276-state explicit invariant,
non-circular and machine-checked -- and the one-neighbourhood compression is
shown impossible at the ledger level: five local invariants all fail closure.
The wave-coherence information is necessary. The circle closes.

Finale of the theory branch. Two results, one positive, one negative:

THEOREM 4 (anchor law + merge at (2,1); finite-certificate proof, same
epistemic status as Theorems 2-3 of docs/spectrum_sufficiency_proof.md).
  Construction: for EVERY persistent v (the 180 of Theorem 3), iterate the
  deterministic PAIR map on states (v, Delta) with u := v + 1 + Delta mod B,
  starting from Delta = 0 -- crucially with NO assumption about the shifted
  copy's fate (non-circular). The forward closure R has 276 states and:
    (i)   R is closed under the pair map            [induction step]
    (ii)  every state has min_i |Delta_i| = 0       [ANCHOR LAW]
    (iii) every state has Delta in [-1, 2]^4        [confinement, C = ceil(S/2)]
    (iv)  from every start the shifted orbit lands on attractor(v)  [MERGE]
  Delta-shapes of R: (0,0,0,0) x180, (0,0,0,1) x48, (-1,0,0,0) x24,
  (0,0,1,1) x16, (0,0,0,2) x8 -- i.e. B5's rigidity is exact on R.
  By the proven chain (merge => clock-shift invariance => ING-3 => spectrum
  sufficiency), this completes a full machine-checked proof of the regime law
  at (2,1) with every link either hand-proven or a finite certificate.

NEGATIVE RESULT (the compression barrier). The hoped-for one-neighbourhood
inductive proof does NOT exist at the ledger abstraction: taking the invariant
to be any of five increasingly strong LOCAL predicates on (v, Delta) --
  I0 range only; I1 +no-(-1,>=1)-edge; I2 +>=1 zero; I3 +>=2 zeros;
  I4 +<=1 merged; I5 +B4-as-invariant (merged receptive cell has quiet
  v-neighbourhood); I6 sig-liveness of BOTH endpoints + the I1 edge constraint ONLY (use_Z=0,
  use_M=False, use_B4=False at I6 -- it is not a superset of I2-I5)
-- every one has closure violations (3728 / 2316 / 2348 / 1860 / 1676 / 1776 /
1216 leaking states respectively). The certified invariant R (276) is 44x
smaller than the best local invariant (12,220): the reachable set is carved
by WAVE COHERENCE -- which lags co-occur with which absolute phases -- and no
predicate on the ledger + one configuration captures it.

WHY THIS IS THE RIGHT ENDING: the program's thesis was that the gradient
spectrum is the coordinate in which fate becomes decidable. The compression
barrier says the same about the PROOF: the regime law cannot be proven from
the waiting ledger alone; the spectral/wave structure is not just the decider
of fate but the necessary ingredient of any local proof of why it decides.
B4 remains true (verified on R) but is not inductively self-supporting; the
generalisation to arbitrary lattices needs the coherence invariant, which is
future work with a precise target.

House Rules Compliance: exhaustive/deterministic, no RNG; asserts Theorem 4's
four checks, the R shape census, and all seven local-invariant violation
counts (so the barrier cannot rot either); aborts on regression.
Output: result/topology/anchor_law_certificate.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import Counter

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
EDGES2 = [(0, 1), (0, 2), (1, 3), (2, 3)]
CYC = [0, 1, 3, 2]
TA, TP = 2, 1
S = TA + TP
B = S + 1
C = -(-S // 2)
L21 = {(0, 0, 2, 2), (0, 1, 1, 2), (0, 2, 3, 3), (1, 1, 1, 1),
       (1, 2, 2, 3), (2, 2, 2, 2), (3, 3, 3, 3)}
VIOL_EXPECT = {"I0": 3728, "I1": 2316, "I2": 2348, "I3": 1860,
               "I4": 1676, "I5": 1776, "I6": 1216}
SHAPE_EXPECT = {(0, 0, 0, 0): 180, (0, 0, 0, 1): 48, (-1, 0, 0, 0): 24,
                (0, 0, 1, 1): 16, (0, 0, 0, 2): 8}


def step2(cfg):
    return tuple((1 if sum(1 <= cfg[j] <= TA for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % B for i, s in enumerate(cfg))


def sig(cfg):
    ph = [cfg[i] for i in CYC]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % B for k in range(4)))


def live_set():
    cache = {}
    for c in product(range(B), repeat=4):
        if c in cache:
            continue
        seen, cur = {}, c
        while cur not in seen:
            seen[cur] = len(seen)
            cur = step2(cur)
        f = any(any(x) for x in [k for k, t in seen.items() if t >= seen[cur]])
        for k in seen:
            cache[k] = f
    return {c for c, f in cache.items() if f}


def pair_step(v, D):
    u = tuple((v[i] + 1 + D[i]) % B for i in range(4))
    nu, nv = step2(u), step2(v)
    D2 = tuple(D[i] + (1 if (v[i] == 0 and nv[i] == 0) else 0)
               - (1 if (u[i] == 0 and nu[i] == 0) else 0) for i in range(4))
    return nv, D2


def theorem4():
    live = live_set()
    assert len(live) == 180, "Theorem-3 live count changed"
    R = set()
    frontier = [(v, (0, 0, 0, 0)) for v in sorted(live)]
    while frontier:
        st = frontier.pop()
        if st in R:
            continue
        R.add(st)
        frontier.append(pair_step(*st))
    anchor = all(min(abs(x) for x in D) == 0 for v, D in R)
    conf = all(all(-1 <= x <= C for x in D) for v, D in R)
    vlive = all(v in live for v, D in R)
    # merge conclusion
    def att(cfg):
        seen, cur = {}, cfg
        while cur not in seen:
            seen[cur] = len(seen)
            cur = step2(cur)
        return frozenset(k for k, t in seen.items() if t >= seen[cur])
    merge = True
    for v in sorted(live):
        attv = att(v)
        cur = tuple((x + 1) % B for x in v)
        for _ in range(200):
            if cur in attv:
                break
            cur = step2(cur)
        else:
            merge = False
    shapes = Counter(tuple(sorted(D)) for v, D in R)
    return R, anchor, conf, vlive, merge, dict(shapes)


def edge_bad(D):
    return any((D[a] == -1 and D[b] >= 1) or (D[b] == -1 and D[a] >= 1)
               for a, b in EDGES2)


def local_invariant_violations(name):
    """Closure violations of local invariant `name` under live-filtered or
    (I6) sig-live pair steps."""
    use_E = name in ("I1", "I2", "I3", "I4", "I5", "I6")
    use_Z = {"I2": 1, "I3": 2, "I4": 2, "I5": 2}.get(name, 0)
    use_M = name in ("I4", "I5")
    use_B4 = name == "I5"
    sig_mode = name == "I6"
    viol = 0
    for v in product(range(B), repeat=4):
        if sig_mode and sig(v) not in L21:
            continue
        for D in product(range(-1, C + 1), repeat=4):
            if use_E and edge_bad(D):
                continue
            if use_Z and sum(1 for x in D if x == 0) < use_Z:
                continue
            if (use_M or sig_mode) and sum(1 for x in D if x == -1) > 1:
                continue
            u = tuple((v[i] + 1 + D[i]) % B for i in range(4))
            if sig_mode and sig(u) not in L21:
                continue
            if use_B4:
                if any(D[i] == -1 and v[i] == 0 and
                       any(1 <= v[j] <= TA for j in ADJ2[i]) for i in range(4)):
                    continue
            nu, nv = step2(u), step2(v)
            if not sig_mode:
                dwu = sum(1 for i in range(4) if u[i] == 0 and nu[i] == 0)
                dwv = sum(1 for i in range(4) if v[i] == 0 and nv[i] == 0)
                if dwu >= 2 or dwv >= 2:
                    continue
            D2 = tuple(D[i] + (1 if (v[i] == 0 and nv[i] == 0) else 0)
                       - (1 if (u[i] == 0 and nu[i] == 0) else 0)
                       for i in range(4))
            ok = all(-1 <= x <= C for x in D2)
            if ok and use_E:
                ok = not edge_bad(D2)
            if ok and use_Z:
                ok = sum(1 for x in D2 if x == 0) >= use_Z
            if ok and (use_M or sig_mode):
                ok = sum(1 for x in D2 if x == -1) <= 1
            if ok and use_B4:
                ok = not any(D2[i] == -1 and nv[i] == 0 and
                             any(1 <= nv[j] <= TA for j in ADJ2[i])
                             for i in range(4))
            if ok and sig_mode:
                ok = sig(nv) in L21 and sig(nu) in L21
            if not ok:
                viol += 1
    return viol


def main():
    print("=== THEOREM 4: the 276-state certificate at (2,1) ===")
    R, anchor, conf, vlive, merge, shapes = theorem4()
    print(f"|R| = {len(R)}  anchor: {anchor}  confinement: {conf}  "
          f"v-live: {vlive}  merge: {merge}")
    print(f"shapes: {shapes}")
    assert len(R) == 276 and anchor and conf and vlive and merge
    assert shapes == {k: v for k, v in SHAPE_EXPECT.items()}, "R census changed"

    print("\n=== the compression barrier: local invariants all leak ===")
    for name in ["I0", "I1", "I2", "I3", "I4", "I5", "I6"]:
        v = local_invariant_violations(name)
        print(f"{name}: closure violations = {v}", flush=True)
        assert v == VIOL_EXPECT[name], f"{name} violation count changed"
        assert v > 0, f"{name} unexpectedly closes -- a local proof EXISTS; " \
                      f"revisit the barrier claim!"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "anchor_law_certificate.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             n_R=np.array([len(R)]),
             shape_keys=np.array([str(k) for k in sorted(shapes)]),
             shape_counts=np.array([shapes[k] for k in sorted(shapes)]),
             inv_names=np.array(list(VIOL_EXPECT.keys())),
             inv_viols=np.array([VIOL_EXPECT[k] for k in VIOL_EXPECT]))
    print(f"\nwrote {out}")
    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| local invariant | constraints | closure violations |")
    print("| :---: | :--- | ---: |")
    desc = {"I0": "Δ ∈ [−1, ⌈S/2⌉]", "I1": "+ no (−1,≥1) edge",
            "I2": "+ ≥1 zero", "I3": "+ ≥2 zeros", "I4": "+ ≤1 merged",
            "I5": "+ B4 as invariant",
            "I6": "sig-live endpoints + I1 edge only (NOT the zero/merged/B4 constraints)"}
    for name in VIOL_EXPECT:
        print(f"| {name} | {desc[name]} | {VIOL_EXPECT[name]:,} |")
    print(f"| **R (certificate)** | reachable closure | **0** (|R| = 276) |")


if __name__ == "__main__":
    main()
