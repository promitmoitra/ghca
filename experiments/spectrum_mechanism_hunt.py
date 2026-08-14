"""Why does the gradient spectrum decide persistence when tau_a >= tau_p?
Opening measurements: three mechanism candidates falsified, three exact laws found.

Theory branch companion to persistent_set_structure.md (Sections 3, 10). The
empirical fact to explain: persistence is a function of the sorted phase-gap
multiset ("gradient spectrum") exactly when tau_a >= tau_p (one-way; see (2,5)).
This script records the first round of mechanism candidates AGAINST the
exhaustive 2x2 oracle -- including the failures, which prune the theory space:

FALSIFIED (each was a plausible one-line mechanism):

F1. BRIDGING LEMMA: "on a live attractor with tau_a >= tau_p, a cell becoming
    receptive always has an active loop-neighbour (re-excitation is immediate,
    so relative lags suffice)." FALSE, and uncorrelated with the regime: over
    the 12 cells tested it fails at exactly ONE -- (1,1), the tightest cell,
    where 0 of 8 receptive events are bridged because the live wave has period
    4 > S+1 = 3, so cells DWELL at receptive with no active neighbour -- and
    holds everywhere else, INCLUDING every tau_p > tau_a cell where the
    spectrum is insufficient. Neither necessary nor sufficient for the regime.

F2. GLOBAL LUMPABILITY: "the spectrum partition is a congruence: sig(c1) =
    sig(c2) implies sig(step c1) = sig(step c2), so the dynamics descend to
    spectrum space." FALSE everywhere, both regimes -- e.g. (2,2): 11 of 14
    classes have multiple image spectra. Spectrum-sufficiency at tau_a >= tau_p
    is fate-constancy on classes WITHOUT the classes evolving coherently.

F3. LIVE-RESTRICTED LUMPABILITY: same, restricted to persistent configurations.
    FALSE everywhere except the trivial (1,1) (2 live classes).

EXACT LAWS (hold at every cell tested, BOTH regimes -- so none of them alone
can explain the regime boundary, but they are the invariant structure any
proof must run through):

L1. PERIOD LAW: on every attractor, for every cell i,
        T = k_i * S + d_i
    where T is the attractor period, k_i the number of firings of cell i per
    period, d_i its receptive-dwell steps per period, S = tau_a + tau_p the
    length of the non-receptive excursion. (Corollary: k_i*S + d_i is
    cell-independent.) Receptive dwell is a real degree of freedom of the wave
    -- the naive T = S+1 holds only when dwell is 0.

L2. SPECTRUM CONSTANCY ON ATTRACTORS: the gradient spectrum is constant along
    every attractor cycle (it is a genuine invariant of the asymptotic wave,
    even though it is NOT preserved on transients -- that is F2).

L3. DWELL IS SPECTRUM-DETERMINED: the sorted per-cell dwell profile of an
    attractor is a function of its (constant) spectrum, and even of the
    INITIAL configuration's spectrum, at every cell tested -- including
    tau_p > tau_a. So conditional on survival, the spectrum already determines
    the asymptotic wave's shape everywhere; what tau_p > tau_a breaks is only
    WHETHER the wave survives.

The sharpened open question these leave: on transients, the spectrum is not
preserved (F2) -- yet at tau_a >= tau_p the fate boundary in configuration
space is never crossed by two same-spectrum configurations. The mechanism must
therefore be a property of the TRANSIENT spectrum-trajectory bundle, not of a
static partition: candidate form, "at tau_a >= tau_p the set of spectra
reachable from a class is fate-homogeneous," which is weaker than lumpability
and is what the next experiment should test (reachable-spectrum-set as the
invariant).

House Rules Compliance: exhaustive, deterministic, no RNG; substrate/analysis
boundary as in the parent doc (all constructs analytic).
Output: result/topology/spectrum_mechanism_hunt.npz + printed doc table.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import product
from collections import defaultdict

ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
CELLS = [(1, 1), (2, 1), (2, 2), (3, 1), (3, 2), (3, 3), (4, 3), (4, 4),
         (1, 2), (2, 3), (2, 4), (3, 4)]


def step2(cfg, ta, tp):
    S = ta + tp
    return tuple((1 if sum(1 <= cfg[j] <= ta for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % (S + 1) for i, s in enumerate(cfg))


def sig(cfg, S):
    ph = [cfg[i] for i in CYC]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % (S + 1) for k in range(4)))


def orbit_full(cfg, ta, tp):
    seen, cur = {}, cfg
    while cur not in seen:
        seen[cur] = len(seen)
        cur = step2(cur, ta, tp)
    start = seen[cur]
    inv = {t: c for c, t in seen.items()}
    return [inv[t] for t in range(start)], [inv[t] for t in range(start, len(seen))]


def audit_cell(ta, tp):
    S = ta + tp
    r = dict(ta=ta, tp=tp)

    # --- F1: bridging lemma on live attractors -------------------------------
    seen, events, bridged = set(), 0, 0
    # --- L1/L2/L3 accumulators ------------------------------------------------
    law_ok, spec_const = True, True
    spec2dwell, init2dwell = defaultdict(set), defaultdict(set)
    for c in product(range(S + 1), repeat=4):
        tr, att = orbit_full(c, ta, tp)
        if not any(any(x) for x in att):
            continue
        key = tuple(sorted(att))
        if key not in seen:
            seen.add(key)
            T = len(att)
            for t, cfg in enumerate(att):
                nxt = att[(t + 1) % T]
                for i in range(4):
                    if cfg[i] == S:                     # becomes receptive next
                        events += 1
                        bridged += any(1 <= nxt[j] <= ta for j in ADJ2[i])
            for i in range(4):
                d_i = sum(1 for x in att if x[i] == 0)
                k_i = sum(1 for x in att if x[i] == 1)
                law_ok &= (T == k_i * S + d_i)
            specs = {sig(x, S) for x in att}
            spec_const &= (len(specs) == 1)
        T = len(att)
        prof = tuple(sorted(sum(1 for x in att if x[i] == 0) for i in range(4)))
        spec2dwell[frozenset(sig(x, S) for x in att)].add((T, prof))
        init2dwell[sig(c, S)].add((T, prof))
    r["n_attr"] = len(seen)
    r["bridge_events"] = events
    r["bridge_ok"] = bridged
    r["bridging_holds"] = (events == bridged)
    r["period_law"] = law_ok
    r["spec_const"] = spec_const
    r["dwell_fn_attr_spec"] = all(len(v) == 1 for v in spec2dwell.values())
    r["dwell_fn_init_spec"] = all(len(v) == 1 for v in init2dwell.values())

    # --- F2: global lumpability ----------------------------------------------
    img = defaultdict(set)
    for c in product(range(S + 1), repeat=4):
        img[sig(c, S)].add(sig(step2(c, ta, tp), S))
    r["n_classes"] = len(img)
    r["lump_viol"] = sum(1 for v in img.values() if len(v) > 1)

    # --- F3: live-restricted lumpability -------------------------------------
    fate = {}
    for c in product(range(S + 1), repeat=4):
        if c in fate:
            continue
        path, cur = [], c
        while cur not in fate and cur not in path:
            path.append(cur)
            cur = step2(cur, ta, tp)
        val = fate[cur] if cur in fate else any(
            any(x) for x in path[path.index(cur):])
        for p in path:
            fate[p] = val
    limg = defaultdict(set)
    for c, alive in fate.items():
        if alive:
            limg[sig(c, S)].add(sig(step2(c, ta, tp), S))
    r["live_classes"] = len(limg)
    r["live_lump_viol"] = sum(1 for v in limg.values() if len(v) > 1)
    return r


def main():
    recs = [audit_cell(ta, tp) for ta, tp in CELLS]
    for r in recs:
        regime = "ta>=tp" if r["ta"] >= r["tp"] else "ta<tp"
        print(f"({r['ta']},{r['tp']}) [{regime:6s}] attrs={r['n_attr']:3d} | "
              f"F1 bridging: {'holds' if r['bridging_holds'] else 'FAILS'} "
              f"({r['bridge_ok']}/{r['bridge_events']}) | "
              f"F2 lump viol: {r['lump_viol']}/{r['n_classes']} | "
              f"F3 live viol: {r['live_lump_viol']}/{r['live_classes']} | "
              f"L1: {r['period_law']} L2: {r['spec_const']} "
              f"L3: {r['dwell_fn_attr_spec']}/{r['dwell_fn_init_spec']}",
              flush=True)

    assert all(r["period_law"] for r in recs), "L1 violated"
    assert all(r["spec_const"] for r in recs), "L2 violated"
    assert all(r["dwell_fn_init_spec"] for r in recs), "L3 violated"
    assert not all(r["bridging_holds"] for r in recs), "F1 unexpectedly held"

    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "spectrum_mechanism_hunt.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out, **{k: np.array([r[k] for r in recs]) for k in recs[0]})
    print(f"\nwrote {out}")

    print("\n--- results-doc table (generated; do not hand-type) ---")
    print("| cell | regime | F1 bridging | F2 lump. viol. | F3 live viol. | "
          "L1 period law | L2 spec const | L3 dwell=f(spec) |")
    print("| :---: | :---: | :---: | ---: | ---: | :---: | :---: | :---: |")
    for r in recs:
        regime = "τa≥τp" if r["ta"] >= r["tp"] else "τa<τp"
        print(f"| ({r['ta']}, {r['tp']}) | {regime} | "
              f"{'holds' if r['bridging_holds'] else '**fails**'} | "
              f"{r['lump_viol']}/{r['n_classes']} | "
              f"{r['live_lump_viol']}/{r['live_classes']} | yes | yes | yes |")


if __name__ == "__main__":
    main()
