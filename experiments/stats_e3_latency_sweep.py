"""3a / P2 — E3 composition across target latency (the "lucky resonance" test).

The audit's central E3 caveat: the composition headline rides a hand-picked
target latency whose required gate tau (latency ~ tau-2) happens to land in an
identity-learnable zone of the substrate. This sweeps TARGET_LAT and measures
joint-success rate (identity AND timing) at n=50 per point, overlaying the
identity-vs-fixed-tau resonance map. If joint-success tracks the resonance map,
the composition result is contingent on the operating point, exactly as flagged.

Monkeypatches e3_timed_response.TARGET_LAT (a module global that both the timed
and factored experiments read); modules are singletons so the patch propagates.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

import ghca_stats as st          # noqa: E402
import e3_timed_response as e3   # noqa: E402
import e3_factored_credit as m   # noqa: E402

N = int(os.environ.get("STATS_N", "50"))
LATS = [int(x) for x in os.environ.get("STATS_LATS", "10,12,14,16,18,20").split(",")]
REF_SEEDS = int(os.environ.get("STATS_REF_SEEDS", "10"))
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def main():
    rates, ident = {}, {}
    per_seed = {}
    for L in LATS:
        e3.TARGET_LAT = L
        ch = np.zeros(N)
        joint = np.zeros(N)
        for s in range(N):
            c, le = m.run_curriculum(s)
            ch[s] = c
            joint[s] = 1.0 if (c >= 0.75 and np.isfinite(le) and le <= e3.LAT_TOL) else 0.0
        k = int(joint.sum())
        p, lo, hi = st.wilson_ci(k, N)
        rates[L] = (k, p, lo, hi)
        ident[L] = st.summarise(f"L={L}", ch)
        per_seed[f"ident_L{L}"] = ch
        per_seed[f"joint_L{L}"] = joint
        print(f"  TARGET_LAT={L:2d} (~tau {L+2}): joint {k}/{N}={p:.3f} "
              f"[{lo:.3f},{hi:.3f}]  identity mean={ident[L]['mean']:.3f} "
              f"({ident[L]['shape'] if 'shape' in ident[L] else ''})", flush=True)

    # reference: identity-learnability vs fixed gate tau (higher-seed resonance map)
    taus = list(range(10, 25))
    e3.TARGET_LAT = 18  # restore default before the fixed-tau reference
    reson = m.resonance_map(np.array(taus), seeds=REF_SEEDS)
    print("  resonance map (identity vs fixed tau):",
          {t: round(float(a), 2) for t, a in zip(taus, reson)}, flush=True)

    np.savez(os.path.join(OUT, "e3_latency_sweep.npz"),
             lats=np.array(LATS), n=N,
             joint_rate=np.array([rates[L][1] for L in LATS]),
             joint_lo=np.array([rates[L][2] for L in LATS]),
             joint_hi=np.array([rates[L][3] for L in LATS]),
             ident_mean=np.array([ident[L]["mean"] for L in LATS]),
             taus=np.array(taus), resonance=reson, **per_seed)
    with open(os.path.join(OUT, "e3_latency_sweep.json"), "w") as f:
        json.dump({"n": N,
                   "points": {L: {"joint_k": rates[L][0], "joint_rate": rates[L][1],
                                  "joint_ci": [rates[L][2], rates[L][3]],
                                  "identity_mean": ident[L]["mean"]} for L in LATS},
                   "resonance": {t: float(a) for t, a in zip(taus, reson)}},
                  f, indent=2, default=float)
    print("  wrote e3_latency_sweep.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
