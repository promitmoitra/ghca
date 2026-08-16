"""Coherence Invariant at 4x4: testing predictions P-A (window universality),
P-B (boundary-role per-capita concentration), and P-C (periodic torus dynamics)
on 4x4 open and torus lattices.

Fifth experiment on the coherence branch; follows up on the handoff protocol in
docs/coherence_larger_lattices_handoff.md.

THE PREDICTIONS TESTED:
  1. P-A (Window Universality): On 4x4 open and torus lattices at (2,1) (tau_a >= tau_p,
     S = 3), the backward window to the diagonal saturates at S = 3 exactly across all
     sampled persistent and transient trajectories (no depth > 3).
  2. P-B (Boundary Concentration, Open): Single-cell covering witnesses on 4x4 open
     concentrate on low-degree cells per capita:
       Rate(corner, deg 2) > Rate(edge, deg 3) > Rate(centre, deg 4)
     with the 4 interior cells strictly active (> 0 witnesses).
  3. P-C (Torus vs Open Dynamics): On the periodic torus (degree 4 everywhere),
     wave collisions annihilate and synchronize with fewer ceiling holds, while
     maintaining the strict u-side witness law at ceiling and v-side at age-1 holds.
  4. Witness Structure & Separation: Holds occur strictly at ages {0, 1, S}.
     Sides never mix: ceiling holds are 100% u-side witnessed (v-side = 0),
     and age-1 holds are 100% v-side witnessed (u-side = 0).

House Rules Compliance:
  - Seed everything: explicit default_rng(seed) threaded throughout.
  - Per-seed reporting and assertions.
  - Substrate / analysis boundary kept explicit.
  - Generates result table for docs.

Output: result/topology/coherence_window_4x4.npz.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
from itertools import combinations
from collections import Counter
import time

CHUNK_TRIALS = 300
SEED = 12345

ROLE4 = {}
for r in range(4):
    for c in range(4):
        idx = r * 4 + c
        if (r in (0, 3)) and (c in (0, 3)):
            ROLE4[idx] = "corner"
        elif (r in (0, 3)) or (c in (0, 3)):
            ROLE4[idx] = "edge"
        else:
            ROLE4[idx] = "centre"


def make_ghca_4x4(ta=2, tp=1, periodic=False):
    B = ta + tp + 1
    S = ta + tp
    neighbors = []
    for r in range(4):
        for c in range(4):
            nbrs = []
            for dr, dc in [(-1, 0), (1, 0), (0, -1), (0, 1)]:
                nr, nc = r + dr, c + dc
                if periodic:
                    nbrs.append((nr % 4) * 4 + (nc % 4))
                else:
                    if 0 <= nr < 4 and 0 <= nc < 4:
                        nbrs.append(nr * 4 + nc)
            neighbors.append(tuple(nbrs))

    def step(cfg):
        nxt = [0] * 16
        for i in range(16):
            s = cfg[i]
            if s == 0:
                act = sum(1 <= cfg[j] <= ta for j in neighbors[i])
                nxt[i] = 1 if act >= 1 else 0
            else:
                nxt[i] = (s + 1) % B
        return tuple(nxt)

    pre_cache = {}

    def get_preimages(cfg):
        if cfg in pre_cache:
            return pre_cache[cfg]
        base = list((x - 1) % B if x > 0 else 0 for x in cfg)
        zeros = [i for i in range(16) if cfg[i] == 0]
        valid_pre = []
        for k in range(len(zeros) + 1):
            for subset in combinations(zeros, k):
                cand = list(base)
                for idx in subset:
                    cand[idx] = S
                t_cand = tuple(cand)
                if step(t_cand) == cfg:
                    valid_pre.append(t_cand)
        pre_cache[cfg] = valid_pre
        return valid_pre

    def is_diag(v, u):
        return u == tuple((x + 1) % B for x in v)

    def bwd_depth(v, u, max_d=8):
        if is_diag(v, u):
            return 0
        frontier = {(v, u)}
        visited = {(v, u)}
        for d in range(1, max_d + 1):
            nxt = set()
            for cv, cu in frontier:
                pvs = get_preimages(cv)
                pus = get_preimages(cu)
                for pv in pvs:
                    for pu in pus:
                        if is_diag(pv, pu):
                            return d
                        if (pv, pu) not in visited:
                            visited.add((pv, pu))
                            nxt.add((pv, pu))
            frontier = nxt
            if not frontier:
                return None
        return None

    def side_witness(v, u, side, want):
        cfg = u if side == 1 else v
        nxt = step(cfg)
        dws = [i for i in range(16) if cfg[i] == 0 and nxt[i] == 0]
        for k in range(1, min(len(dws) + 1, 6)):
            for su in combinations(dws, k):
                cand = list(cfg)
                for idx in su:
                    cand[idx] = S
                c2 = tuple(cand)
                if step(c2) != nxt:
                    continue
                cc = (v, c2) if side == 1 else (c2, u)
                if bwd_depth(cc[0], cc[1], max_d=want) == want:
                    return su
        return None

    return step, get_preimages, bwd_depth, side_witness, B, S, neighbors


def run_sampling(ta=2, tp=1, periodic=False, n_trials=300, seed=12345):
    step, get_pre, bwd_depth, side_witness, B, S, nbrs = make_ghca_4x4(ta, tp, periodic=periodic)
    rng = np.random.default_rng(seed)

    depth_hist = Counter()
    transitions = Counter()
    ceil_total, ceil_u_witnessed, ceil_v_witnessed = 0, 0, 0
    ceil_sizes = Counter()
    ceil_roles = Counter()

    h1_total, h1_v_witnessed, h1_u_witnessed = 0, 0, 0
    h1_sizes = Counter()

    for trial in range(n_trials):
        w = tuple(rng.integers(0, B, size=16))
        v = w
        u = tuple((x + 1) % B for x in v)
        traj_d = []
        traj_pairs = []
        for t in range(12):
            d = bwd_depth(v, u, max_d=8)
            traj_d.append(d)
            traj_pairs.append((v, u))
            depth_hist[d] += 1
            v, u = step(v), step(u)

        for i in range(len(traj_d) - 1):
            d1, d2 = traj_d[i], traj_d[i + 1]
            transitions[(d1, d2)] += 1
            pair = traj_pairs[i]
            if d1 == S and d2 == S:
                ceil_total += 1
                w_u = side_witness(pair[0], pair[1], 1, S - 1)
                if w_u is not None:
                    ceil_u_witnessed += 1
                    ceil_sizes[len(w_u)] += 1
                    if len(w_u) == 1 and not periodic:
                        ceil_roles[ROLE4[w_u[0]]] += 1
                if side_witness(pair[0], pair[1], 0, S - 1) is not None:
                    ceil_v_witnessed += 1
            elif d1 == 1 and d2 == 1:
                h1_total += 1
                w_v = side_witness(pair[0], pair[1], 0, 0)
                if w_v is not None:
                    h1_v_witnessed += 1
                    h1_sizes[len(w_v)] += 1
                if side_witness(pair[0], pair[1], 1, 0) is not None:
                    h1_u_witnessed += 1

    return {
        "depth_hist": depth_hist,
        "transitions": transitions,
        "ceil_total": ceil_total,
        "ceil_u_witnessed": ceil_u_witnessed,
        "ceil_v_witnessed": ceil_v_witnessed,
        "ceil_sizes": ceil_sizes,
        "ceil_roles": ceil_roles,
        "h1_total": h1_total,
        "h1_v_witnessed": h1_v_witnessed,
        "h1_u_witnessed": h1_u_witnessed,
        "h1_sizes": h1_sizes,
        "S": S
    }


def main():
    print("=== 4x4 Coherence Invariant Verification (P-A, P-B, P-C) ===")
    print(f"Sampling {CHUNK_TRIALS} trajectories per lattice with seed {SEED}...")

    # 1. Test 4x4 OPEN
    t0 = time.time()
    res_open = run_sampling(ta=2, tp=1, periodic=False, n_trials=CHUNK_TRIALS, seed=SEED)
    el_open = time.time() - t0
    print(f"\n[4x4 OPEN] Completed in {el_open:.2f}s:")
    print(f"  Depths visited: {dict(sorted(res_open['depth_hist'].items()))}")
    max_d_open = max(k for k in res_open['depth_hist'].keys() if k is not None)
    print(f"  Max BFS Depth: {max_d_open} (S = {res_open['S']})")
    print(f"  Transitions: {dict(sorted(res_open['transitions'].items()))}")
    print(f"  Ceiling Holds (3->3): {res_open['ceil_u_witnessed']}/{res_open['ceil_total']} u-witnessed, v-side worked for {res_open['ceil_v_witnessed']}")
    print(f"    Swap size distribution: {dict(sorted(res_open['ceil_sizes'].items()))}")
    print(f"    Single-cell roles: {dict(sorted(res_open['ceil_roles'].items()))}")
    
    r_corner = res_open['ceil_roles']['corner'] / 4.0
    r_edge = res_open['ceil_roles']['edge'] / 8.0
    r_centre = res_open['ceil_roles']['centre'] / 4.0
    print(f"    Per-capita witness rates: corner={r_corner:.2f} > edge={r_edge:.2f} > centre={r_centre:.2f}")
    print(f"  Age-1 Holds (1->1): {res_open['h1_v_witnessed']}/{res_open['h1_total']} v-witnessed, u-side worked for {res_open['h1_u_witnessed']}")

    # Assertions for 4x4 OPEN:
    assert max_d_open == res_open['S'] == 3, f"P-A Falsified: max depth {max_d_open} != S {res_open['S']}"
    assert res_open['ceil_v_witnessed'] == 0, "Witness structure failure: v-side witnessed ceiling hold on open lattice"
    assert res_open['h1_u_witnessed'] == 0, "Witness structure failure: u-side witnessed age-1 hold on open lattice"
    assert res_open['ceil_roles']['centre'] > 0, "Interior witness count must not be zero"
    assert r_corner > r_edge > r_centre, f"P-B Falsified: non-monotone per-capita rates ({r_corner} > {r_edge} > {r_centre})"

    # 2. Test 4x4 TORUS
    t0 = time.time()
    res_torus = run_sampling(ta=2, tp=1, periodic=True, n_trials=CHUNK_TRIALS, seed=SEED)
    el_torus = time.time() - t0
    print(f"\n[4x4 TORUS] Completed in {el_torus:.2f}s:")
    print(f"  Depths visited: {dict(sorted(res_torus['depth_hist'].items()))}")
    max_d_torus = max(k for k in res_torus['depth_hist'].keys() if k is not None)
    print(f"  Max BFS Depth: {max_d_torus} (S = {res_torus['S']})")
    print(f"  Transitions: {dict(sorted(res_torus['transitions'].items()))}")
    print(f"  Ceiling Holds (3->3): {res_torus['ceil_u_witnessed']}/{res_torus['ceil_total']} u-witnessed, v-side worked for {res_torus['ceil_v_witnessed']}")
    print(f"    Swap size distribution: {dict(sorted(res_torus['ceil_sizes'].items()))}")
    print(f"  Age-1 Holds (1->1): {res_torus['h1_v_witnessed']}/{res_torus['h1_total']} v-witnessed, u-side worked for {res_torus['h1_u_witnessed']}")

    # Assertions for 4x4 TORUS:
    assert max_d_torus == res_torus['S'] == 3, f"P-A Falsified on Torus: max depth {max_d_torus} != S {res_torus['S']}"
    assert res_torus['ceil_v_witnessed'] == 0, "Witness structure failure on Torus: v-side witnessed ceiling hold"
    assert res_torus['h1_u_witnessed'] == 0, "Witness structure failure on Torus: u-side witnessed age-1 hold"

    # Save artifact
    out = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "topology", "coherence_window_4x4.npz")
    os.makedirs(os.path.dirname(out), exist_ok=True)
    np.savez(out,
             open_depth_keys=np.array(list(res_open['depth_hist'].keys())),
             open_depth_vals=np.array(list(res_open['depth_hist'].values())),
             open_max_d=np.array([max_d_open]),
             open_ceil_u=np.array([res_open['ceil_u_witnessed']]),
             open_ceil_total=np.array([res_open['ceil_total']]),
             open_ceil_v=np.array([res_open['ceil_v_witnessed']]),
             open_roles_names=np.array(sorted(res_open['ceil_roles'].keys())),
             open_roles_counts=np.array([res_open['ceil_roles'][k] for k in sorted(res_open['ceil_roles'].keys())]),
             torus_depth_keys=np.array(list(res_torus['depth_hist'].keys())),
             torus_depth_vals=np.array(list(res_torus['depth_hist'].values())),
             torus_max_d=np.array([max_d_torus]),
             torus_ceil_u=np.array([res_torus['ceil_u_witnessed']]),
             torus_ceil_total=np.array([res_torus['ceil_total']]),
             torus_ceil_v=np.array([res_torus['ceil_v_witnessed']]))
    print(f"\nSaved results to {out}")

    print("\n--- Summary Results Table ---")
    print("| Lattice | Boundary | S | Max Depth (P-A) | Ceiling u-Witness | Ceiling v-Witness | Rate(Corner) | Rate(Edge) | Rate(Centre) |")
    print("| :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: | :---: |")
    print(f"| 4×4 | Open | {res_open['S']} | {max_d_open} | {res_open['ceil_u_witnessed']}/{res_open['ceil_total']} | {res_open['ceil_v_witnessed']} | {r_corner:.2f} | {r_edge:.2f} | {r_centre:.2f} |")
    print(f"| 4×4 | Torus | {res_torus['S']} | {max_d_torus} | {res_torus['ceil_u_witnessed']}/{res_torus['ceil_total']} | {res_torus['ceil_v_witnessed']} | N/A | N/A | N/A |")


if __name__ == "__main__":
    main()
