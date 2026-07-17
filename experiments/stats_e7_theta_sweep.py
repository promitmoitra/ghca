"""3a / P2 — E7 switching across the spiral-band threshold theta.

P1 found E7 switching softened at n=50 (0.86 -> 0.75) with a low-seed tail. This
asks whether the intact-vs-ablated switching dissociation is robust across the E0
spiral band, or an artifact of the single hand-chosen theta=4.0. Sweeps the
spiral substrate's THETA_LIVE (module constant read by make_spiral); ablation
keeps THETA_DEAD=5.0. As theta_live rises toward the death threshold the intact
spiral should degrade toward the ablated baseline -- that band edge is the point.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))

import ghca_stats as st         # noqa: E402
import e7_learning as e7        # noqa: E402

N = int(os.environ.get("STATS_N", "30"))
THETAS = [float(x) for x in os.environ.get("STATS_THETAS", "3.0,3.5,4.0,4.5").split(",")]
NB, BL = 30, 25
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def _tail(a, k=BL * 4):
    a = np.asarray(a, float).ravel()
    return float(a[-k:].mean()) if a.size else float("nan")


def main():
    per_seed = {}
    rows = []
    default = e7.THETA_LIVE
    for th in THETAS:
        e7.THETA_LIVE = th
        intact = np.array([_tail(e7.run_switching(s, False, NB, BL)[0]) for s in range(N)])
        ablate = np.array([_tail(e7.run_switching(s, True, NB, BL)[0]) for s in range(N)])
        per_seed[f"intact_th{th}"] = intact
        per_seed[f"ablate_th{th}"] = ablate
        ri = st.summarise(f"intact θ={th}", intact)
        gap = st.effect_size(intact, ablate)
        rows.append((th, ri, float(intact.mean()), float(ablate.mean()), gap))
        print(f"  θ={th}: intact {ri['mean']:.3f} [{ri['ci_lo']:.3f},{ri['ci_hi']:.3f}] "
              f"vs ablate {ablate.mean():.3f}  (Cohen d={gap:.2f})", flush=True)
    e7.THETA_LIVE = default

    np.savez(os.path.join(OUT, "e7_theta_sweep.npz"),
             thetas=np.array(THETAS), n=N,
             intact_mean=np.array([r[2] for r in rows]),
             ablate_mean=np.array([r[3] for r in rows]),
             intact_lo=np.array([r[1]["ci_lo"] for r in rows]),
             intact_hi=np.array([r[1]["ci_hi"] for r in rows]),
             cohen_d=np.array([r[4] for r in rows]), **per_seed)
    with open(os.path.join(OUT, "e7_theta_sweep.json"), "w") as f:
        json.dump({"n": N, "points": {r[0]: {"intact": r[2], "ablate": r[3],
                   "intact_ci": [r[1]["ci_lo"], r[1]["ci_hi"]], "cohen_d": r[4]}
                   for r in rows}}, f, indent=2, default=float)
    print("  wrote e7_theta_sweep.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
