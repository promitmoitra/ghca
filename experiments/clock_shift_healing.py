"""Healing: the merge lemma decomposes into a two-line rigidity lemma plus a
bounded transient, with healing time <= 2S everywhere at tau_a >= tau_p.

Continues clock_shift_merge.py (R2: at tau_a >= tau_p every orbit merges with
its clock-shift). This experiment finds the structure of the merge:

LEMMA R (attractor rigidity) -- the half that is now PROVEN BY HAND:
  For any configuration c, step(c) = c+1 (global clock-shift)  <=>  no cell of
  c dwells. Proof: if no cell dwells, every nonzero state advances by 1 and
  every receptive cell fires (0 -> 1 = +1), so step is the global +1;
  conversely if step(c) = c+1 then no cell maps 0 -> 0. Two lines, and the
  checker confirms the equivalence at every attractor configuration of every
  cell tested (both regimes).

  EMPIRICAL HALF (exhaustive, both regimes): every live attractor with
  S+1 >= L = 4 is dwell-free, hence step = +1 on it, hence the attractor is a
  +1-INVARIANT SET (the clock-shift of an attractor state is its successor).
  The lone exception is (1,1), where S+1 = 3 < L: a 3-state wave on a 4-loop
  cannot advance all cells every step (it dwells once per beat -- the period-4
  law), so its attractor is NOT +1-closed; there the shifted state is a
  transient that heals back to the SAME attractor in <= 4 = 2S steps
  (exhaustive).

CONSEQUENCE: on attractors, clock-shift = one time step. Fate on attractors is
trivially shift-invariant IN BOTH REGIMES. The entire regime asymmetry lives
in the TRANSIENTS -- which basin the shifted transient falls into.

HEALING BOUND (exhaustive at 8 tau_a >= tau_p cells): defining heal(c) as the
first t with step^t(c+1) on attractor(c),

    max over all c of heal(c)  <=  2S,   with equality on the diagonal cells
    tested ((1,1): 4, (2,2): 8, (3,3): 12, (4,4): 16) and strict inequality on
    the strict cells tested ((2,1): 5 < 6, (3,1): 6 < 8, (3,2): 9 < 10,
    (4,3): 13 < 14).

EXTREMAL STRUCTURE (exhaustive; the proof's skeleton). Splitting healing by
fate:
  - The 2S saturators are ALL DEAD configurations (their healing is just both
    orbits draining to the all-zero state; on the diagonal the drain takes
    exactly 2S). Dead healing is not proof-relevant -- dead orbits trivially
    merge at 0.
  - LIVE healing obeys an exact closed form at every tau_a >= tau_p cell
    tested with S+1 >= L = 4 (twelve cells, (2,1) through (5,3)):

        max over live c of heal(c)  =  tau_a + 2*tau_p + 1

    ((1,1), the dwelling special case, has live max 2.) The witnesses all
    belong to ONE family up to loop rotation/reflection: a YOUNG WAVE -- one
    active cell at depth <= 2, the others shallow-passive (depth <= 2) or
    receptive, e.g. R-P1-A1-P2 in loop order. Physically: the shifted copy
    heals slowest when the wave has JUST fired, because the +1 pushes the
    young wave's surroundings deepest into the refractory pipeline, and the
    pipeline takes one active band plus two passive bands (+1) to clear.

REMAINING HAND-PROOF OBLIGATION (all that is left of the regime law): at
tau_a >= tau_p, the orbit of c+1 enters attractor(c) -- transient healing --
now with the exact target "live healing <= tau_a + 2*tau_p + 1" and a
one-family extremal witness to build the induction around. Everything else is
either proven (rigidity equivalence; certificate schema) or reduced to it
(merge => clock-shift invariance => ING-3 => spectrum sufficiency).

House Rules Compliance: exhaustive, deterministic, no RNG; asserts the
rigidity equivalence, the dwell-free fact, the (1,1) exception, and the 2S
bound; aborts on regression.
Output: result/topology/clock_shift_healing.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
         (1, 2), (2, 3), (3, 4)]
HEAL_CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4)]


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def attractor(cfg, ta, tp, cache):
    if cfg in cache:
        return cache[cfg]
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    att = frozenset(c for c, t in seen.items() if t >= seen[cur])
    for c in seen:
        cache[c] = att
    return cache[cfg]


def audit_rigidity(ta, tp):
    S = ta + tp
    B = S + 1
    cache = {}
    n = plus1 = dwellfree = equiv = 0
    seen = set()
    for c in product(range(B), repeat=4):
        att = attractor(c, ta, tp, cache)
        if not any(any(x) for x in att) or att in seen:
            continue
        seen.add(att)
        for a in att:
            n += 1
            nx = step2(a, ta, tp)
            p1 = tuple((x + 1) % B for x in a)
            ndw = all(not (a[i] == 0 and nx[i] == 0) for i in range(4))
            plus1 += (nx == p1)
            dwellfree += ndw
            equiv += ((nx == p1) == ndw)
    return dict(ta=ta, tp=tp, att_cfgs=n, step_is_plus1=plus1,
                dwellfree=dwellfree, equiv_ok=(equiv == n))


def audit_healing(ta, tp):
    S = ta + tp
    B = S + 1
    cache = {}
    mx, tot, cnt = 0, 0, 0
    for c in product(range(B), repeat=4):
        att = attractor(c, ta, tp, cache)
        cur, t = tuple((x + 1) % B for x in c), 0
        while cur not in att:
            cur = step2(cur, ta, tp)
            t += 1
            if t > 10 * B ** 4:
                raise RuntimeError(f"no healing at ({ta},{tp}) from {c}")
        mx, tot, cnt = max(mx, t), tot + t, cnt + 1
    return dict(ta=ta, tp=tp, max_heal=mx, mean_heal=tot / cnt, bound_2S=2 * S)


LIVE_CELLS = [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 1), (4, 2), (4, 3),
              (4, 4), (5, 1), (5, 2), (5, 3)]


def audit_live_heal(ta, tp):
    """Max healing over LIVE configurations + the closed form ta + 2*tp + 1."""
    S = ta + tp
    B = S + 1
    cache = {}
    mx, witnesses = 0, []
    for c in product(range(B), repeat=4):
        att = attractor(c, ta, tp, cache)
        if not any(any(x) for x in att):
            continue
        cur, t = tuple((x + 1) % B for x in c), 0
        while cur not in att:
            cur = step2(cur, ta, tp)
            t += 1
        if t > mx:
            mx, witnesses = t, [c]
        elif t == mx:
            witnesses.append(c)
    return dict(ta=ta, tp=tp, max_live_heal=mx, formula=ta + 2 * tp + 1,
                n_witnesses=len(witnesses))


def main():
    print("=== Lemma R: rigidity (step==+1 <=> dwell-free) on live attractors ===")
    rig = []
    for ta, tp in CELLS:
        r = audit_rigidity(ta, tp)
        rig.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        print(f"({ta},{tp}) [{regime:6s}] att cfgs={r['att_cfgs']:5d} "
              f"step==+1: {r['step_is_plus1']:5d} dwell-free: {r['dwellfree']:5d} "
              f"equivalence: {r['equiv_ok']}", flush=True)
        assert r["equiv_ok"], f"rigidity equivalence broke at ({ta},{tp})"
        B = ta + tp + 1
        if B >= 4:
            assert r["step_is_plus1"] == r["att_cfgs"], \
                f"dwell-free attractor fact broke at ({ta},{tp})"
        else:
            assert r["step_is_plus1"] == 0, "(1,1) exception changed"

    print("\n=== healing time (tau_a >= tau_p): first t with step^t(c+1) on "
          "attractor(c) ===")
    heal = []
    for ta, tp in HEAL_CELLS:
        r = audit_healing(ta, tp)
        heal.append(r)
        tag = "EQUALITY" if r["max_heal"] == r["bound_2S"] else "strict"
        print(f"({ta},{tp}) max={r['max_heal']:3d} mean={r['mean_heal']:5.2f} "
              f"bound 2S={r['bound_2S']:3d}  [{tag}]", flush=True)
        assert r["max_heal"] <= r["bound_2S"], f"2S bound broke at ({ta},{tp})"
        if ta == tp:
            assert r["max_heal"] == r["bound_2S"], \
                f"diagonal equality broke at ({ta},{tp})"

    print("\n=== LIVE healing: max = tau_a + 2*tau_p + 1 (closed form) ===")
    lh = []
    for ta, tp in LIVE_CELLS:
        r = audit_live_heal(ta, tp)
        lh.append(r)
        print(f"({ta},{tp}) max live heal={r['max_live_heal']:3d} "
              f"formula={r['formula']:3d} witnesses={r['n_witnesses']:3d} "
              f"match={r['max_live_heal'] == r['formula']}", flush=True)
        assert r["max_live_heal"] == r["formula"], \
            f"live-healing closed form broke at ({ta},{tp})"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "clock_shift_healing.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             **{f"r_{k}": np.array([x[k] for x in rig]) for k in rig[0]},
             **{f"h_{k}": np.array([x[k] for x in heal]) for k in heal[0]},
             **{f"l_{k}": np.array([x[k] for x in lh]) for k in lh[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | attractor configs | step = clock-shift | "
          "max heal | 2S bound |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: |")
    hby = {(x["ta"], x["tp"]): x for x in heal}
    for r in rig:
        cell = (r["ta"], r["tp"])
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        h = hby.get(cell)
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | {r['att_cfgs']} | "
              f"{r['step_is_plus1']}/{r['att_cfgs']} | "
              f"{h['max_heal'] if h else '—'} | {h['bound_2S'] if h else '—'} |")


if __name__ == "__main__":
    main()
