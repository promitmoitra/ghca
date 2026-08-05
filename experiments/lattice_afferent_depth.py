"""How far into a recurrent medium does exogenous timing information penetrate?

`lattice_timescale_notes.md` showed that a **global** afferent timing pulse lets the
input-timing τ rule track a drive period exactly. That was the easy case, and it also
quietly departed from this repo's own convention: the layered experiments (E1–E6, 3d/3e)
carry `roles = {sensory, hidden, motor}`, in which the **sensory population *is* the
afferent channel** — a cell never has to be *told* which input is exogenous, because
afferent and lateral input arrive on different edges. That is the biologically ordinary
case (thalamic vs corticocortical input on distinct compartments).

So put the afferent where it belongs — a **localised sensory strip along one edge** —
and ask the question that then becomes interesting:

  Cells inside the strip have a clean exogenous timing signal. Cells deep in the medium
  have only lateral input, which is the self-confirming signal that fails. **How far from
  the sensory edge can τ still learn the right period?**

Geometry: periodic in y, **closed in x** (walls), so distance from the sensory edge is
simply the x coordinate and the gradient is 1-D and directly readable.

  * `x < W_S` — **sensory strip.** The transducer: τ is imposed at the drive period and
    the cells are force-fired on the beat, exactly as `layered_graph` gives sensory nodes
    fixed properties. Excluded from the τ statistics; it is not part of the learned medium.
  * `x >= W_S` — **hidden medium.** Plastic τ. Waves arrive from the strip.

Conditions (τ rule identical throughout; only the afferent's *extent* changes):
  local   — afferent sensed only inside the strip; the medium learns from lateral input
            (the honest localised case, and the penetration-depth measurement)
  global  — afferent sensed by every cell (the previously-run reference: should be flat)
  own     — the e10 self-firing rule (ratchet control)

Metric: |τ − P| as a function of distance from the sensory edge, plus a penetration
depth — the last distance at which the medium is still locked (|τ − P| < 1.5).
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("DEP_L", "96"))
W_S = int(os.environ.get("DEP_WS", "4"))          # sensory strip width
STEPS = int(os.environ.get("DEP_STEPS", "6000"))
N_SEEDS = int(os.environ.get("DEP_N", "5"))
PERIODS = [int(p) for p in os.environ.get("DEP_P", "6,12").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 34.0, 0.15
LOCK = 1.5                                        # |τ−P| below this counts as locked
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("DEP_TAG", "")


def nbr_count(active):
    """4-neighbour count: periodic in y, closed (walled) in x."""
    n = np.roll(active, 1, 0).astype(np.float32) + np.roll(active, -1, 0)
    xm = np.roll(active, 1, 1).astype(np.float32); xm[:, 0] = 0.0      # no wrap at x=0
    xp = np.roll(active, -1, 1).astype(np.float32); xp[:, -1] = 0.0    # no wrap at x=L-1
    return n + xm + xp


def run(mode, P, seed):
    rng = np.random.default_rng(seed)
    phi = np.zeros((L, L), np.int32)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)

    xs = np.arange(L)[None, :].repeat(L, 0)
    sensory = xs < W_S                             # the transducer strip
    medium = ~sensory
    tau[sensory] = P                               # imposed: sensory nodes keep the beat

    last_ev = np.full((L, L), -1.0, np.float32)
    last_fire = np.full((L, L), -1.0, np.float32)
    prev_lat = np.zeros((L, L), bool)
    prev_ev = np.zeros((L, L), bool)
    # timing-event threshold: a single planar wavefront presents ~1 active neighbour, so
    # requiring >=2 makes one passing wave register ONCE rather than several times.
    theta_ev = int(mode.split("-c")[1]) if "-c" in mode else int(THETA)

    for t in range(STEPS):
        active = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        nbr = nbr_count(active)

        beat = (t % P) == 0
        # the afferent timing signal: subthreshold, carries timing only.
        aff = np.zeros((L, L), bool)
        if beat:
            if mode.startswith("global"):
                aff[:] = True                      # every cell senses it (the easy case)
            else:
                aff[sensory] = True                # localised to the sensory strip
        # the strip is also the actual drive (suprathreshold) so waves exist
        drive = sensory & beat

        lateral = nbr >= THETA
        excite = rested & (lateral | drive)

        lat_edge = lateral & ~prev_lat
        prev_lat = lateral
        supra_ev = nbr >= theta_ev
        lat_ev = supra_ev & ~prev_ev
        prev_ev = supra_ev

        moving = phi >= 1
        phi[moving] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1

        if mode == "own":
            ev = (phi == 1) & medium
            prev = last_fire
        elif mode.endswith("-aff"):
            ev = aff & medium              # afferent-only rule
            prev = last_ev
        else:
            ev = (lat_ev | aff) & medium   # all-inputs rule (lat_ev == lat_edge if theta_ev==THETA)
            prev = last_ev

        seen = ev & (prev >= 0)
        if seen.any():
            dt = np.clip((t - prev[seen]).astype(np.float32), 0, 60.0)
            tau[seen] = np.clip(tau[seen] + ETA * (dt - tau[seen]), TMIN, TMAX)
        if mode == "own":
            last_fire[ev] = t
        else:
            last_ev[ev] = t

    err_by_x = np.abs(tau - P).mean(0)[W_S:]       # mean over y, medium columns only
    tau_by_x = tau.mean(0)[W_S:]
    return err_by_x, tau_by_x


def depth(err_by_x):
    """Last distance (in cells from the strip) at which the medium is still locked."""
    locked = err_by_x < LOCK
    if not locked[0]:
        return 0
    d = 0
    while d + 1 < len(locked) and locked[d + 1]:
        d += 1
    return d + 1


def main():
    out = {}
    for P in PERIODS:
        for mode in ("global-aff", "local", "local-c2", "local-c3", "local-c4", "own"):
            E = np.zeros((N_SEEDS, L - W_S)); T = np.zeros((N_SEEDS, L - W_S))
            for s in range(N_SEEDS):
                E[s], T[s] = run(mode, P, s)
            key = f"P{P}_{mode}"
            out[key] = {"err": E, "tau": T}
            e = E.mean(0)
            d = depth(e)
            print(f"P={P:2d} {mode:10s}: |τ−P| at x=0 {e[0]:5.2f} | x=8 {e[min(8,len(e)-1)]:5.2f} "
                  f"| x=24 {e[min(24,len(e)-1)]:5.2f} | far {e[-1]:5.2f} "
                  f"| locked to depth {d:3d} cells  (τ far = {T.mean(0)[-1]:.1f})", flush=True)
        print(flush=True)

    save = {"L": L, "W_S": W_S, "steps": STEPS, "n": N_SEEDS,
            "periods": np.array(PERIODS), "lock": LOCK}
    for k in out:
        save[f"{k}_err"] = out[k]["err"]; save[f"{k}_tau"] = out[k]["tau"]
    np.savez(os.path.join(OUT, f"lattice_afferent_depth{TAG}.npz"), **save)
    with open(os.path.join(OUT, f"lattice_afferent_depth{TAG}.json"), "w") as f:
        json.dump({"L": L, "W_S": W_S, "n": N_SEEDS, "periods": PERIODS, "lock": LOCK,
                   "rows": {k: {"err_profile": out[k]["err"].mean(0).round(3).tolist(),
                                "depth": depth(out[k]["err"].mean(0))} for k in out}},
                  f, indent=2)
    print(f"wrote lattice_afferent_depth{TAG}.{{npz,json}}", flush=True)


if __name__ == "__main__":
    main()
