"""3a / P1 — seed & CI scale-up for the E-series headline dissociations.

Drives each experiment's committed per-seed run-functions at n=50 (publication
grade; see docs/stats_sweeps_plan.md), computes bootstrap 95% CIs, effect sizes,
and the automated bimodality flag from ghca_stats, and saves per-seed arrays +
summaries to result/stats/. It does NOT edit or re-run the experiments' own
committed n=5 numbers — it imports their run-functions and calls them fresh.

Usage:
    python experiments/stats_seed_scaleup.py e5           # one experiment
    python experiments/stats_seed_scaleup.py e1 e3f e5    # several
    python experiments/stats_seed_scaleup.py all
Optional: STATS_N env var overrides the seed count (default 50).
"""

import os
import sys
import json

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "50"))
OUTDIR = os.path.join(ROOT, "result", "stats")
os.makedirs(OUTDIR, exist_ok=True)


def _tail_mean(x, k=120):
    """Reduce a per-trial array (or scalar) to a headline scalar."""
    a = np.asarray(x, float).ravel()
    return float(a[-k:].mean()) if a.size else float("nan")


# ---- per-experiment adapters: return {label: per-seed 1-D array} -------------

def run_e1():
    import e1_conditioning as m
    groups = {L: np.array([_tail_mean(m.run_condition(L, s)[0], 120)
                           for s in range(N)]) for L in ("A", "B", "AB")}
    _log("E1", groups)
    return groups, {"effect_A_vs_B": st.effect_size(groups["A"], groups["B"])}


def run_e2():
    import e2_delayed_response as m
    # retention at the longest delay (D=200) is where the A/B dissociation is starkest
    full = {L: np.array([m.run_condition(L, s)[0] for s in range(N)])
            for L in ("A", "B", "AB")}
    groups = {L: full[L][:, -1] for L in full}          # D=200 column
    _log("E2 (retention @ D=200)", groups)
    extra = {f"{L}_full": full[L] for L in full}
    extra["delays"] = np.array(m.DELAYS_EVAL)
    extra["effect_B_vs_A"] = st.effect_size(groups["B"], groups["A"])
    return groups, extra


def run_e3f():
    import e3_timed_response as e3
    import e3_factored_credit as m
    cur = np.array([m.run_curriculum(s) for s in range(N)])       # (N, 2): ch, le
    fac = np.array([m.run_single("AB", "factored", s) for s in range(N)])
    ch_cur, le_cur = cur[:, 0], cur[:, 1]
    joint = ((ch_cur >= 0.75) & (le_cur <= e3.LAT_TOL)).astype(float)
    groups = {"curriculum_identity": ch_cur, "factored_identity": fac[:, 0]}
    k = int(joint.sum())
    p, lo, hi = st.wilson_ci(k, N)
    _log("E3 factored", groups)
    print(f"    joint-success rate = {k}/{N} = {p:.3f}  Wilson95=[{lo:.3f},{hi:.3f}]")
    return groups, {"curriculum_lat_err": le_cur, "joint_success": joint,
                    "joint_rate": p, "joint_ci": (lo, hi),
                    "TARGET_LAT": e3.TARGET_LAT, "LAT_TOL": e3.LAT_TOL}


def run_e3t():
    import e3_timed_response as m
    res = {L: np.array([m.run_condition(L, s)[:2] for s in range(N)])  # ch, lat
           for L in ("A", "B", "AB")}
    groups = {f"{L}_identity": res[L][:, 0] for L in res}
    _log("E3 timed (identity)", groups)
    extra = {f"{L}_latency": res[L][:, 1] for L in res}
    extra["effect_identity_A_vs_B"] = st.effect_size(groups["A_identity"],
                                                      groups["B_identity"])
    return groups, extra


def run_e5():
    import e5_executive as m
    nb, bl = 40, 25
    sw_i = np.array([_tail_mean(m.run_switching(s, False, nb, bl)[0], bl * 6) for s in range(N)])
    sw_a = np.array([_tail_mean(m.run_switching(s, True, nb, bl)[0], bl * 6) for s in range(N)])
    sr_i = np.array([_tail_mean(m.run_single_rule(s, False), 120) for s in range(N)])
    sr_a = np.array([_tail_mean(m.run_single_rule(s, True), 120) for s in range(N)])
    groups = {"switch_intact": sw_i, "switch_ablate": sw_a,
              "single_intact": sr_i, "single_ablate": sr_a}
    _log("E5", groups)
    return groups, {"effect_switch": st.effect_size(sw_i, sw_a),
                    "effect_single": st.effect_size(sr_i, sr_a)}


def run_e7():
    import e7_learning as m
    nb, bl = 30, 25
    sw_i = np.array([_tail_mean(m.run_switching(s, False, nb, bl)[0], bl * 4) for s in range(N)])
    sw_a = np.array([_tail_mean(m.run_switching(s, True, nb, bl)[0], bl * 4) for s in range(N)])
    sr_i = np.array([_tail_mean(m.run_single_rule(s, False, n_trials=600), 120) for s in range(N)])
    sr_a = np.array([_tail_mean(m.run_single_rule(s, True, n_trials=600), 120) for s in range(N)])
    groups = {"switch_intact": sw_i, "switch_ablate": sw_a,
              "single_intact": sr_i, "single_ablate": sr_a}
    _log("E7", groups)
    return groups, {"effect_switch": st.effect_size(sw_i, sw_a)}


def run_e9():
    import e9_emergent_conjunction as m
    nb, bl = 40, 25
    kinds = ["emergent", "wired", "frozen"]
    routing = {k: np.zeros(N) for k in kinds}
    sel_post = {k: np.zeros(N) for k in kinds}
    for s in range(N):
        for k in kinds:
            net, roles = m.make(s, kind=k)
            rng = np.random.default_rng(s + 3)
            m.probe_responses(net, roles, rng)
            if k == "emergent":
                m.selforg(net, roles, rng)
            R1 = m.probe_responses(net, roles, rng)
            sel_post[k][s] = m.conj_scores(R1).mean()
            acc, _ = m.run_switching(net, roles, s, nb, bl)
            routing[k][s] = _tail_mean(acc, bl * 6)
        print(f"    e9 seed {s + 1}/{N} done", flush=True)
    groups = {f"routing_{k}": routing[k] for k in kinds}
    groups.update({f"selectivity_{k}": sel_post[k] for k in kinds})
    _log("E9", groups)
    return groups, {"effect_emergent_vs_frozen":
                    st.effect_size(routing["emergent"], routing["frozen"])}


ADAPTERS = {"e1": run_e1, "e2": run_e2, "e3f": run_e3f, "e3t": run_e3t,
            "e5": run_e5, "e7": run_e7, "e9": run_e9}


def _log(tag, groups):
    print(f"  [{tag}] n={N}")
    for lab, arr in groups.items():
        print("    " + st.fmt_row(st.summarise(lab, arr)), flush=True)


def _save(key, groups, extra):
    rows = {lab: st.summarise(lab, arr) for lab, arr in groups.items()}
    payload = {f"seed_{lab}": np.asarray(arr, float) for lab, arr in groups.items()}
    for k, v in extra.items():
        try:
            payload[k] = np.asarray(v, float)
        except (TypeError, ValueError):
            pass
    np.savez(os.path.join(OUTDIR, f"{key}_stats.npz"), n=N, **payload)
    with open(os.path.join(OUTDIR, f"{key}_summary.json"), "w") as f:
        json.dump({"n": N, "rows": rows,
                   "extra": {k: (v if np.isscalar(v) else None)
                             for k, v in extra.items()}}, f, indent=2, default=float)
    print(f"    wrote {key}_stats.npz + {key}_summary.json", flush=True)


def main(keys):
    if keys == ["all"]:
        keys = list(ADAPTERS)
    for key in keys:
        if key not in ADAPTERS:
            print(f"  skip unknown '{key}'")
            continue
        print(f"\n=== {key} (n={N}) ===", flush=True)
        groups, extra = ADAPTERS[key]()
        _save(key, groups, extra)


if __name__ == "__main__":
    main(sys.argv[1:] or ["all"])
