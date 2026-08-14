"""Discrete perturbation theory: the regime law is a damage-relaxation
transition, and the healing formula is the relaxation time of the +1 zero mode.

Answers "is there a discrete analogue of perturbative analysis?" concretely.
The formalism is DAMAGE SPREADING (Derrida-Pomeau; Bagnoli's CA Lyapunov
exponents): run a configuration c and its perturbed copy side by side and
track the difference field d_t = u_t - v_t (mod S+1, per cell). Here the
perturbation is not a random flip but the global clock-shift +1, which by
Lemma R (clock_shift_healing.py) is a SYMMETRY of dwell-free dynamics -- in
perturbative language, a ZERO MODE: a uniform d is invisible to the dynamics
wherever nothing dwells.

Life cycle of the damage, measured exhaustively:
  t = 0          d is uniform (= 1 everywhere) BY CONSTRUCTION.
  t_break        first time d becomes non-uniform. This happens exactly when a
                 dwell/wrap event fires in one copy and not the other (the
                 step-commutation analysis of clock_shift_merge.py): the
                 damage SCATTERS off a dwell.
  t_return       first time d is uniform again -- the damage has relaxed back
                 onto the zero mode; from there the orbits differ by a clock
                 shift, which on a dwell-free attractor is a time shift
                 (Lemma R), so they share the attractor.

FINDINGS (exhaustive; 10 cells):

1. At tau_a >= tau_p, EVERY pair that breaks returns, and
       max over live pairs of t_return  =  tau_a + 2*tau_p + 1
   -- the healing closed form IS the zero-mode relaxation time, exactly, at
   every cell tested ((2,1) through (4,4); asserted).
2. At tau_p > tau_a, the pairs that NEVER return are exactly the split-fate
   pairs of clock_shift_merge.py -- counts match one for one ((1,2): 32,
   (2,3): 80, (3,4): 144; asserted). Fate splits iff damage never relaxes.
3. So the regime boundary is a DAMAGE-HEALING TRANSITION in the
   damage-spreading sense: tau_a >= tau_p is the phase where the clock-shift
   zero mode always reabsorbs its scattered damage; tau_p > tau_a is the
   phase where damage can become permanent. The pipeline picture: a scattered
   dwell displacement needs one active band + two passive bands (+1) to
   flush, and tau_a >= tau_p is precisely when the wave's active support
   survives that flush.

Also verified: heal(c) (the entry of orbit(c+1) into attractor(c), the
quantity of clock_shift_healing.py) coincides with the plain transient length
of c+1 at every tau_a >= tau_p cell -- because the merge is total there, so
"entering attractor(c)" and "entering one's own attractor" are the same event.
The two experiments measure one relaxation.

House Rules Compliance: exhaustive, deterministic, no RNG; asserts findings
1-2 and the transient reduction; aborts on regression.
Output: result/topology/damage_relaxation.npz + printed doc table.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CELLS = [(2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
         (1, 2), (2, 3), (3, 4)]
# split-fate counts per tau_p > tau_a cell, from clock_shift_merge.npz
EXPECT_SPLIT = {(1, 2): 32, (2, 3): 80, (3, 4): 144}


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def att_and_entry(cfg, ta, tp, cache):
    if cfg in cache:
        return cache[cfg]
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    start = seen[cur]
    att = frozenset(c for c, t in seen.items() if t >= start)
    for c, t0 in seen.items():
        cache[c] = (att, max(0, start - t0))
    return cache[cfg]


def audit(ta, tp):
    S = ta + tp
    B = S + 1
    cache = {}
    cap = 10 * (ta + 2 * tp + 1) + 60
    n_break = never = 0
    mx_ret_live = 0
    mism_transient = 0
    mx_tr_live = 0
    for c in product(range(B), repeat=4):
        att_c, tr_c = att_and_entry(c, ta, tp, cache)
        alive = any(any(x) for x in att_c)
        c1 = tuple((x + 1) % B for x in c)
        att_1, tr_1 = att_and_entry(c1, ta, tp, cache)
        if alive:
            mx_tr_live = max(mx_tr_live, tr_c)
            if ta >= tp and att_1 != att_c:
                mism_transient += 1
        u, v = c1, c
        t_break = t_ret = None
        for t in range(cap):
            uni = len({(u[i] - v[i]) % B for i in range(4)}) == 1
            if t_break is None:
                if not uni:
                    t_break = t
            elif uni:
                t_ret = t
                break
            u, v = step2(u, ta, tp), step2(v, ta, tp)
        if t_break is not None:
            n_break += 1
            if t_ret is None:
                never += 1
            elif alive:
                mx_ret_live = max(mx_ret_live, t_ret)
    return dict(ta=ta, tp=tp, n_break=n_break, never_return=never,
                max_return_live=mx_ret_live, formula=ta + 2 * tp + 1,
                max_live_transient=mx_tr_live, merge_mismatch=mism_transient)


def main():
    print(f"{'(ta,tp)':>8} {'regime':>7} | {'breaks':>7} {'max ret (live)':>14} "
          f"{'formula':>8} | {'never-return':>12} {'expected split':>14}")
    recs = []
    for ta, tp in CELLS:
        r = audit(ta, tp)
        recs.append(r)
        regime = "ta>=tp" if ta >= tp else "ta<tp"
        exp = EXPECT_SPLIT.get((ta, tp), 0)
        print(f"{str((ta,tp)):>8} {regime:>7} | {r['n_break']:7d} "
              f"{r['max_return_live']:14d} {r['formula']:8d} | "
              f"{r['never_return']:12d} {exp:14d}", flush=True)

    for r in recs:
        cell = (r["ta"], r["tp"])
        if r["ta"] >= r["tp"]:
            assert r["never_return"] == 0, f"damage failed to relax at {cell}"
            assert r["max_return_live"] == r["formula"], \
                f"relaxation-time formula broke at {cell}"
            assert r["merge_mismatch"] == 0, f"merge totality broke at {cell}"
            assert r["max_live_transient"] == r["formula"], \
                f"transient reduction broke at {cell}"
        else:
            assert r["never_return"] == EXPECT_SPLIT[cell], \
                f"never-return != split-fate count at {cell}"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "damage_relaxation.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | pairs that scatter | max relaxation (live) | "
          "τa+2τp+1 | never relax | = split-fate pairs |")
    print("| :---: | :---: | ---: | ---: | ---: | ---: | :---: |")
    for r in recs:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        exp = EXPECT_SPLIT.get((r["ta"], r["tp"]))
        eq = "—" if exp is None else ("yes" if r["never_return"] == exp else "NO")
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | {r['n_break']} | "
              f"{r['max_return_live']} | {r['formula']} | "
              f"{r['never_return']} | {eq} |")


if __name__ == "__main__":
    main()
