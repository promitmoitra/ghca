"""3a / P3b — seed & CI scale-up for the sigma-band / outcome-matrix headlines.

P1-P3 covered every dissociation-of-means headline (E-series, E8, C5 routing,
C6 necessity). What they explicitly deferred: C2/C3's achievable-BAND statistics
(range-across-policy / pooled-std, not a two-arm mean) and C4/C7's outcome
matrices + macro-sufficiency, because the band/matrix code has no outer seed
loop to bootstrap over. This is that per-seed refactor: each headline's own
reusable primitives (build/readouts/realize/decode -- unchanged, imported from
the original scripts) are called in a fresh outer seed loop here, producing a
per-seed ARRAY that feeds the same ghca_stats.bootstrap_ci machinery P1-P3 used.
Nothing in c2/c3/c4/c5/c7's own files is edited; their committed n=1..5 numbers
still reproduce bit-identically.

Usage: python experiments/stats_p3b.py c2 c3 c4 c5 c7
       python experiments/stats_p3b.py all
Optional: STATS_NC env var overrides the seed/replicate count (default 30).
"""

import os
import sys
import json

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

import ghca_stats as st  # noqa: E402

NC = int(os.environ.get("STATS_NC", "30"))
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def _save(key, groups, extra=None):
    rows = {lab: st.summarise(lab, np.asarray(a, float)) for lab, a in groups.items()}
    payload = {f"seed_{lab}": np.asarray(a, float) for lab, a in groups.items()}
    for k, v in (extra or {}).items():
        try:
            payload[k] = np.asarray(v, float)
        except (TypeError, ValueError):
            pass
    np.savez(os.path.join(OUT, f"{key}_stats.npz"), n=NC, **payload)
    with open(os.path.join(OUT, f"{key}_summary.json"), "w") as f:
        json.dump({"n": NC, "rows": rows}, f, indent=2, default=float)
    for lab in groups:
        print("    " + st.fmt_row(rows[lab]), flush=True)
    print(f"    wrote {key}_stats.npz + {key}_summary.json", flush=True)


# ---------------------------------------------------------------------------
# C2 -- do(W) achievable band (labeled-line vs collective), per substrate seed.
# achievable_band_t0 is pure math on the readout weight vector (no trial RNG),
# so the only per-seed cost is baseline_std's 800-step settle.
# ---------------------------------------------------------------------------

def run_c2():
    import c2_fat_handed as c2
    print(f"\n=== c2 achievable band (n={NC} substrate seeds) ===", flush=True)
    band_c = np.zeros(NC)
    band_l = np.zeros(NC)
    for s in range(NC):
        net, w_coll, w_lab = c2.build(seed=s)
        std_c, std_l = c2.baseline_std(net, w_coll, w_lab)
        lo_c, hi_c, _ = zip(*[c2.achievable_band_t0(w_coll, w) for w in c2.W_GRID])
        lo_l, hi_l, _ = zip(*[c2.achievable_band_t0(w_lab, w) for w in c2.W_GRID])
        band_c[s] = ((np.array(hi_c) - np.array(lo_c)) / std_c).mean()
        band_l[s] = ((np.array(hi_l) - np.array(lo_l)) / std_l).mean()
    groups = {"band_collective": band_c, "band_labeled": band_l}
    _save("c2_band", groups, {"ratio_labeled_over_collective": band_l / (band_c + 1e-9)})


# ---------------------------------------------------------------------------
# C3 -- do(theta) response at a fixed operating point (tau=14), per seed;
# the "do(theta) band" is the across-seed std over the within-seed trial std,
# the same sigma-unit convention as C2, computed once from the seed array.
# ---------------------------------------------------------------------------

def run_c3():
    import c3_do_theta as c3
    from ghca_causal import do_theta
    TAU_OP = 14
    print(f"\n=== c3 do(theta={TAU_OP}) response (n={NC} substrate seeds) ===", flush=True)
    Ec = np.zeros(NC)
    within_std = np.zeros(NC)
    for s in range(NC):
        net, w_coll, w_lab = c3.build(seed=s)
        rng = np.random.default_rng(1000 * s + TAU_OP)
        rc = []
        for _ in range(c3.N_TRIALS):
            net.phi[:] = 0
            net.seed_random(0.1, 0.1)
            do_theta(net, tau=TAU_OP, nodes=np.arange(c3.N_H))
            for _ in range(c3.T_PROP):
                net.step(None)
            c, _ = c3.readouts(net, w_coll, w_lab)
            rc.append(c)
        Ec[s] = np.mean(rc)
        within_std[s] = np.std(rc)
    band = float(Ec.std(ddof=1) / (within_std.mean() + 1e-9))
    print(f"    do(theta) band (across-seed / within-seed sigma) = {band:.4f}", flush=True)
    _save("c3_band", {"response_at_tau14": Ec}, {"within_seed_std": within_std,
                                                  "doTheta_band": np.array([band])})


# ---------------------------------------------------------------------------
# C4 -- A. outcome matrix (do(tau_gate) -> timing, do(g_route) -> identity) on
# a seeded E3 net (the original e3_net() hardcodes seed=0; this loops it).
# B. macro-sufficiency already takes `seed` -- just called at more seeds.
# ---------------------------------------------------------------------------

def run_c4():
    import c4_outcome_relativity as c4
    from ghca_net import Network
    import e3_timed_response as e3
    print(f"\n=== c4 outcome matrix + macro-sufficiency (n={NC} seeds) ===", flush=True)

    tau_grid = [8, 14, 20, 26]
    g_grid = [-0.4, -0.15, 0.15, 0.4]
    diag_tau = np.zeros(NC)     # do(tau) -> timing, normalized
    off_tau = np.zeros(NC)      # do(tau) -> identity, normalized
    diag_g = np.zeros(NC)       # do(g_route) -> identity, normalized
    off_g = np.zeros(NC)        # do(g_route) -> timing, normalized
    for s in range(NC):
        W, plastic, roles, theta, _ = e3.build(w_hm=0.15, w_gate=4.0, theta_m=5.0, seed=s)
        net = Network(W, act=e3.ACT, pas=e3.TAU0 - e3.ACT, theta=theta, p_s=0.0, seed=s)
        relay = np.concatenate([roles["hidden"]] + roles["motor"])
        net.tau[relay] = e3.ACT
        W0 = W.copy()
        lat_tau, id_tau = [], []
        for t in tau_grid:
            la, p0 = c4.run_e3(net, roles, W0, t, 0.0)
            lat_tau.append(la); id_tau.append(p0)
        lat_g, id_g = [], []
        for g in g_grid:
            la, p0 = c4.run_e3(net, roles, W0, e3.TAU0, g)
            lat_g.append(la); id_g.append(p0)
        eff = np.array([[np.ptp(lat_tau), np.ptp(id_tau)],
                        [np.ptp(lat_g), np.ptp(id_g)]])
        effn = eff / (eff.max(0, keepdims=True) + 1e-9)
        diag_tau[s], off_tau[s] = effn[0, 0], effn[0, 1]
        off_g[s], diag_g[s] = effn[1, 0], effn[1, 1]
        if (s + 1) % 10 == 0:
            print(f"    outcome matrix seed {s + 1}/{NC} done", flush=True)
    print("  A. outcome-relativity matrix (per seed):", flush=True)
    _save("c4_matrix", {"do_tau_timing": diag_tau, "do_tau_identity": off_tau,
                        "do_route_identity": diag_g, "do_route_timing": off_g})

    print("  B. macro-sufficiency:", flush=True)
    suff_c = np.zeros(NC)
    suff_l = np.zeros(NC)
    for s in range(NC):
        suff_c[s], suff_l[s], _ = c4.macro_sufficiency(seed=s)
        if (s + 1) % 10 == 0:
            print(f"    macro-sufficiency seed {s + 1}/{NC} done", flush=True)
    _save("c4_suff", {"suff_collective": suff_c, "suff_labeled": suff_l})


# ---------------------------------------------------------------------------
# C5 -- decode-by-policy sigma bands (center/tracked/global), replicated: same
# orchestration as measure_bands() but with a seed OFFSET per replicate (the
# original hardcodes base seed 0, so calling it twice is not independent).
# ---------------------------------------------------------------------------

def _measure_bands_rep(c5, base, n_trials=40):
    policies = [("centered", {}),
                ("displaced", {"d": 4}), ("displaced", {"d": 8}), ("displaced", {"d": 12}),
                ("pitch", {"pitch": 1.0}), ("pitch", {"pitch": 2.0}),
                ("noisy", {"jitter": 0.08})]
    net = c5.make_spiral()
    acc = {r: np.zeros(len(policies)) for r in c5.READERS}
    per_trial = {r: [np.zeros(n_trials) for _ in policies] for r in c5.READERS}
    for pi, (policy, kw) in enumerate(policies):
        for r in c5.READERS:
            for k in range(n_trials):
                rng = np.random.default_rng(base + 1000 * pi + k)
                chir = +1 if k % 2 == 0 else -1
                c5.realize(net, chir, policy, rng, **kw)
                g = c5.decode(net, r)
                per_trial[r][pi][k] = float(g == chir)
            acc[r][pi] = per_trial[r][pi].mean()
    bands = {}
    for r in c5.READERS:
        pooled = np.sqrt(np.mean([per_trial[r][pi].var() for pi in range(len(policies))]))
        bands[r] = (acc[r].max() - acc[r].min()) / (pooled + 1e-6)
    return bands


def run_c5():
    import c5_do_chirality as c5
    print(f"\n=== c5 decode-by-policy bands (n={NC} replicates) ===", flush=True)
    readers = list(c5.READERS)
    out = {f"band_{r}": np.zeros(NC) for r in readers}
    for rep in range(NC):
        bands = _measure_bands_rep(c5, base=rep * 10_000_000, n_trials=40)
        for r in readers:
            out[f"band_{r}"][rep] = bands[r]
        if (rep + 1) % 10 == 0:
            print(f"    replicate {rep + 1}/{NC} done", flush=True)
    _save("c5_bands", out)


# ---------------------------------------------------------------------------
# C7 -- outcome matrix + screening; main() already loops per-seed internally
# at n_seeds=3 -- this replays the same per-seed body (train_router / decode_
# context / outcomes_for_context, all unchanged) at NC seeds instead.
# ---------------------------------------------------------------------------

def run_c7():
    import c7_outcome_relativity as c7
    print(f"\n=== c7 outcome matrix + screening (n={NC} seeds) ===", flush=True)
    chi_rule = np.zeros(NC); chi_content = np.zeros(NC)
    route_rule = np.zeros(NC); route_content = np.zeros(NC)
    scr_pp = np.zeros(NC); scr_pm = np.zeros(NC)   # seed=+, inject +/-
    scr_mp = np.zeros(NC); scr_mm = np.zeros(NC)   # seed=-, inject +/-
    for s in range(NC):
        rnet, roles, spnet = c7.train_router(s)
        rng = np.random.default_rng(s + 31)

        gp = c7.decode_context(spnet, +1, rng); Op_rule, Op_con = c7.outcomes_for_context(rnet, roles, gp, rng)
        gm = c7.decode_context(spnet, -1, rng); Om_rule, Om_con = c7.outcomes_for_context(rnet, roles, gm, rng)
        d_chi_rule = abs(Op_rule - Om_rule); d_chi_con = abs(Op_con - Om_con)

        g = c7.decode_context(spnet, +1, rng)
        Oi_rule, Oi_con = c7.outcomes_for_context(rnet, roles, g, rng)
        Wsave = rnet.W.copy()
        rnet.W[rnet.plastic] = rnet.W[rnet.plastic].mean()
        Oe_rule, Oe_con = c7.outcomes_for_context(rnet, roles, g, rng)
        rnet.W[:] = Wsave
        d_g_rule = abs(Oi_rule - Oe_rule); d_g_con = abs(Oi_con - Oe_con)

        m = np.array([[d_chi_rule, d_chi_con], [d_g_rule, d_g_con]])
        mn = m / (m.max(0, keepdims=True) + 1e-9)
        chi_rule[s], chi_content[s] = mn[0, 0], mn[0, 1]
        route_content[s], route_rule[s] = mn[1, 1], mn[1, 0]

        for si, cj, dst in ((+1, +1, "pp"), (+1, -1, "pm"), (-1, +1, "mp"), (-1, -1, "mm")):
            c7.sp.seed_spiral(spnet, si, jitter=0.02, rng=rng)
            for _ in range(6):
                spnet.step(None)
            gg = c7.decode_context(spnet, cj, rng)
            o_rule, _ = c7.outcomes_for_context(rnet, roles, gg, rng, n=60)
            {"pp": scr_pp, "pm": scr_pm, "mp": scr_mp, "mm": scr_mm}[dst][s] = o_rule
        if (s + 1) % 5 == 0:
            print(f"    c7 seed {s + 1}/{NC} done", flush=True)
    print("  A. outcome matrix (normalized, per seed):", flush=True)
    _save("c7_matrix", {"do_chi_rule": chi_rule, "do_chi_content": chi_content,
                        "do_route_rule": route_rule, "do_route_content": route_content})
    print("  B. screening (O_rule by nucleation-seed x injected-chi):", flush=True)
    _save("c7_screen", {"seedplus_injplus": scr_pp, "seedplus_injminus": scr_pm,
                        "seedminus_injplus": scr_mp, "seedminus_injminus": scr_mm})


TARGETS = {"c2": run_c2, "c3": run_c3, "c4": run_c4, "c5": run_c5, "c7": run_c7}

if __name__ == "__main__":
    keys = sys.argv[1:] or ["all"]
    if keys == ["all"]:
        keys = list(TARGETS)
    for k in keys:
        TARGETS.get(k, lambda: print(f"unknown {k}"))()
