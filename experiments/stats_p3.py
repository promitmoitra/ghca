"""3a / P3 — seed & CI scale-up for the E8 predictive series and the C-series.

E8 (E8.1-8.7) is the thinnest evidence in the repo (n=3, some deterministic
ridge panels) and is cheap to scale (linear readouts, no long sims): run its
headlines at n=50 with bootstrap CIs. The C-series headlines are heavier spiral
sims; C5 (behavioural fat-hand) and C6 (necessity) expose an `n_seeds` argument
and already return a SEM, so we call them at higher n and report a normal-approx
95% CI. C6 reuses E7's run_switching verbatim, so it is the E7 dissociation under
C6's windowing, not an independent result -- flagged as such.

Usage: python experiments/stats_p3.py e8         # fast (~minutes)
       python experiments/stats_p3.py c5 c6       # heavy spiral sims
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
NC = int(os.environ.get("STATS_NC", "30"))   # C-series seed count (heavier)
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
    np.savez(os.path.join(OUT, f"{key}_stats.npz"), n=N, **payload)
    with open(os.path.join(OUT, f"{key}_summary.json"), "w") as f:
        json.dump({"n_e8": N, "n_c": NC, "rows": rows}, f, indent=2, default=float)
    for lab in groups:
        print("    " + st.fmt_row(rows[lab]), flush=True)
    print(f"    wrote {key}_stats.npz", flush=True)


def run_e8():
    import e8_predictive as e8p
    import e8_nested as e8n
    import e8_conditional as e8c
    chance = 1.0 / e8p.M
    print(f"\n=== e8 (n={N}, chance={chance:.3f}) ===", flush=True)

    # E8.1/8.2 prediction accuracy by sequence type (periodic is deterministic)
    groups = {}
    per_tau = {}
    a_per, _, _, _, _, _ = e8p.eval_sequence(e8p.seq_periodic(6000, 4), tau=14)
    groups["periodic"] = np.full(N, a_per)  # deterministic
    for name, alpha in [("markov0.9", 0.9), ("markov0.5", 0.5), ("markov0.1", 0.1)]:
        groups[name] = np.array([e8p.eval_sequence(
            e8p.seq_markov(6000, alpha, np.random.default_rng(1000 + s)), tau=14)[0]
            for s in range(N)])
    groups["randomwalk"] = np.array([e8p.eval_sequence(
        e8p.seq_random_walk(6000, np.random.default_rng(2000 + s)), tau=14)[0]
        for s in range(N)])
    print("  E8.1/8.2 prediction:", flush=True)
    _save("e8_prediction", groups)

    # E8.3 integration window vs tau
    print("  E8.3 integration window vs tau:", flush=True)
    wgroups = {}
    for tau in [4, 8, 14, 20, 26]:
        wgroups[f"tau{tau}"] = np.array([e8p.integration_window(tau, seed=s)[0]
                                        for s in range(N)])
    _save("e8_window", wgroups)

    # E8.5 nested: intact (tau_ctx=70) vs ablate (6); nested vs within-context
    print("  E8.5 nested:", flush=True)
    ng = {
        "nested_intact": np.array([e8n.eval_nested(70, True, seed=s) for s in range(N)]),
        "nested_ablate": np.array([e8n.eval_nested(6, True, seed=s) for s in range(N)]),
        "within_intact": np.array([e8n.eval_nested(70, False, seed=s) for s in range(N)]),
    }
    _save("e8_nested", ng)

    # E8.7 conditional: grid+conj vs controls, at K=4
    print("  E8.7 conditional (K=4):", flush=True)
    cg = {
        "grid_conj": np.array([e8c.eval_condition("grid", 10, True, 4, seed=s) for s in range(N)]),
        "grid_noconj": np.array([e8c.eval_condition("grid", 10, False, 4, seed=s) for s in range(N)]),
        "trace_conj": np.array([e8c.eval_condition("trace", 26, True, 4, seed=s) for s in range(N)]),
    }
    _save("e8_conditional", cg)


def _ci_from_sem(mean, sem):
    return mean - 1.96 * sem, mean + 1.96 * sem


def run_c5():
    import c5_do_chirality as c5
    print(f"\n=== c5 routing fat-hand (n={NC}) ===", flush=True)
    res = c5.routing_confirmation(n_seeds=NC)   # {reader: (mean, sem)}
    out = {}
    for reader, (m, sem) in res.items():
        lo, hi = _ci_from_sem(m, sem)
        out[reader] = {"mean": float(m), "ci": [float(lo), float(hi)], "sem": float(sem)}
        print(f"  {reader:8s} {m:.3f} [{lo:.3f}, {hi:.3f}]  (normal-approx)", flush=True)
    with open(os.path.join(OUT, "c5_routing_stats.json"), "w") as f:
        json.dump({"n": NC, "readers": out}, f, indent=2, default=float)
    print("    wrote c5_routing_stats.json", flush=True)


def run_c6():
    import c6_do_theta_chi as c6
    print(f"\n=== c6 necessity (n={NC}; reuses E7 run_switching) ===", flush=True)
    nec = c6.necessity(n_seeds=NC)
    out = {}
    for k, (m, sem) in nec.items():
        lo, hi = _ci_from_sem(m, sem)
        out[k] = {"mean": float(m), "ci": [float(lo), float(hi)], "sem": float(sem)}
        print(f"  {k:16s} {m:.3f} [{lo:.3f}, {hi:.3f}]", flush=True)
    with open(os.path.join(OUT, "c6_necessity_stats.json"), "w") as f:
        json.dump({"n": NC, "note": "reuses e7.run_switching; the E7 dissociation "
                   "under C6 windowing, not independent", "groups": out}, f,
                  indent=2, default=float)
    print("    wrote c6_necessity_stats.json", flush=True)


TARGETS = {"e8": run_e8, "c5": run_c5, "c6": run_c6}

if __name__ == "__main__":
    keys = sys.argv[1:] or ["e8"]
    for k in keys:
        TARGETS.get(k, lambda: print(f"unknown {k}"))()
