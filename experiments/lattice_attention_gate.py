"""An attention strip, made of the same substrate, gating plasticity per column.

Where this comes from. `lattice_afferent_depth.py` established that exogenous timing
does **not** penetrate a recurrent medium: with the afferent confined to a sensory strip,
τ never locks anywhere (error 1.54 near the strip → 2.92 far away), because several
wavefronts arrive per drive cycle and a last-interval estimator cannot recover the
period. It also showed that **coincidence gating cannot fix it**: requiring ≥2 active
neighbours to count a timing event makes coincidences rare, so intervals become long and
τ ratchets to the ceiling (33.9 at θ_ev=2). So there is a bind — any amplitude threshold
low enough to sample every cycle admits spurious traffic, and any threshold high enough
to reject the traffic misses most cycles. **Amplitude cannot separate them; phase might.**

Hence a gate. This adds a **1-D attention strip of the same Greenberg–Hastings cells**,
running *orthogonal* to the sensory edge so that it has exactly one element per column,
i.e. one per distance-from-input. It is seeded at x=0 by the sensory layer on each beat
and propagates along x as a travelling pulse, so column x is gated once per drive period,
at a phase delayed in proportion to x.

This is E4's and E5's machinery pointed at plasticity: E4 showed selection is native to
the substrate (colliding waves annihilate — a winner-take-all), and E5 used a persistent
loop as an *option* gating spatial routing. Here the same idea gates **learning** rather
than routing — a three-factor-style rule (eligibility × gate).

Why it is not the global afferent in disguise: the gate is launched from the sensory
layer and carried by GH cells at finite speed, so it arrives at column x *later* than at
column 0. `gate_phase_by_x` is reported precisely so this is checkable — a propagating
reference shows arrival phase rising with x; a disguised global pulse would be flat.

Conditions:
  global-aff — global afferent, afferent-only rule (the trivial upper bound)
  local      — localised afferent, all-inputs rule (the failure to beat)
  att-gate   — timing event = the gate's **rising edge** (gate supplies the phase
               reference directly; the medium's own input is ignored)
  att-3f     — three-factor: the gate still supplies *when* (so intervals stay clean),
               but the update only lands on cells with recent lateral input — the gate
               times the plasticity, the medium's activity selects who gets it
  att-lat-wN — timing event = the **first** lateral edge inside each attention window,
               where the window is N steps wide (true gating: the cell still measures
               its own input, but only attended arrivals are timestamped — the closest
               thing to "attention determines learning")
  own        — the e10 self-firing rule (ratchet control)

Two things had to be got right, and both were already lessons from earlier runs:

  * **Edge, not level.** A gate that is *open* for N steps fires N times if you read the
    level, giving dt=1 and collapsing τ to the floor. Read its rising edge. (Same bug as
    the original lattice port, where level-detecting lateral input floored τ at 3.1.)
  * **One measurement per window.** Widening the window makes it likely that *some*
    lateral edge falls inside it — but several may, so latch the first and disarm. Too
    narrow a window and most periods contain no edge at all, and τ ratchets to the
    ceiling exactly as coincidence gating did (`lattice_afferent_depth.py`, θ_ev=2 → 33.9).

The attention strip's own τ is **hand-set** (just above the window width, always below
the drive period) so it recovers before the next beat. That is an honest hand-set choice,
not something learned; a plastic strip is a follow-up (it is the easy 1-D case, driven
from one end).
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("ATT_L", "96"))
W_S = int(os.environ.get("ATT_WS", "4"))
STEPS = int(os.environ.get("ATT_STEPS", "6000"))
N_SEEDS = int(os.environ.get("ATT_N", "3"))
PERIODS = [int(p) for p in os.environ.get("ATT_P", "6,12").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 34.0, 0.15
LOCK = 1.5
MODES = ("global-aff", "local", "att-gate", "att-3f",
         "att-lat-w1", "att-lat-w2", "att-lat-w4", "att-lat-w6", "own")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def nbr_count(active):
    """4-neighbour count: periodic in y, closed in x."""
    n = np.roll(active, 1, 0).astype(np.float32) + np.roll(active, -1, 0)
    xm = np.roll(active, 1, 1).astype(np.float32); xm[:, 0] = 0.0
    xp = np.roll(active, -1, 1).astype(np.float32); xp[:, -1] = 0.0
    return n + xm + xp


def run(mode, P, seed):
    rng = np.random.default_rng(seed)
    phi = np.zeros((L, L), np.int32)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)

    xs = np.arange(L)[None, :].repeat(L, 0)
    sensory = xs < W_S
    medium = ~sensory
    tau[sensory] = P

    # --- the attention strip: a 1-D GH chain, one cell per column -------------------
    # GW = attention window width (the strip's own active duration). tau_a sits just
    # above it and must stay below P so each column recovers before the next beat.
    # NB tau_a is set from GW alone — the drive period P never enters it, so the strip
    # is not covertly tuned to the thing the medium is supposed to learn. (For att-gate
    # and att-3f, GW=ACT=2 and tau_a=3 at every P tested.)
    GW = int(mode.split("-w")[1]) if "-w" in mode else int(ACT)
    tau_a = float(GW + 1)
    always_open = tau_a >= P                       # degenerate control: never gated
    phi_a = np.zeros(L, np.int32)
    gate_hits = np.zeros(L)                        # for the propagation check
    gate_phase_sum = np.zeros(L)
    prev_gate = np.zeros(L, bool)

    last_ev = np.full((L, L), -1.0, np.float32)
    last_fire = np.full((L, L), -1.0, np.float32)
    prev_lat = np.zeros((L, L), bool)
    armed = np.zeros((L, L), bool)                 # one measurement per window
    elig_t = np.full((L, L), -1000, np.int32)      # eligibility: last lateral arrival
    n_upd = np.zeros((L, L), np.int32)             # how often each cell actually learned

    for t in range(STEPS):
        active = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        nbr = nbr_count(active)

        beat = (t % P) == 0
        aff = np.zeros((L, L), bool)
        if beat:
            if mode == "global-aff":
                aff[:] = True
            else:
                aff[sensory] = True
        drive = sensory & beat

        # --- advance the attention chain: seeded at x=0 by the beat, travels along x --
        a_act = (phi_a >= 1) & (phi_a <= GW)
        a_rested = phi_a == 0
        a_in = np.zeros(L, bool)
        a_in[1:] = a_act[:-1]                     # cell x excited by cell x-1
        a_in[0] = beat                            # seeded from the sensory layer
        a_moving = phi_a >= 1
        phi_a[a_moving] += 1
        phi_a[phi_a > int(np.ceil(tau_a))] = 0
        phi_a[a_rested & a_in] = 1
        gate_open = (phi_a >= 1) & (phi_a <= GW)   # per-column gate, shape (L,)
        if always_open:
            gate_open = np.ones(L, bool)
        gate_edge = gate_open & ~prev_gate         # EDGE, not level
        prev_gate = gate_open
        gate_hits += gate_edge
        gate_phase_sum += gate_edge * (t % P)

        lateral = nbr >= THETA
        excite = rested & (lateral | drive)
        lat_edge = lateral & ~prev_lat
        prev_lat = lateral

        moving = phi >= 1
        phi[moving] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1

        gate2d = np.broadcast_to(gate_open[None, :], (L, L))
        # rearm each window; with an always-open gate every step is inside the window
        armed |= np.broadcast_to((gate_edge | always_open)[None, :], (L, L))
        if mode == "own":
            ev = (phi == 1) & medium
            prev = last_fire
        elif mode == "global-aff":
            ev = aff & medium
            prev = last_ev
        elif mode in ("att-gate", "att-3f"):
            # the gate's rising edge itself is the timing event
            ev = np.broadcast_to(gate_edge[None, :], (L, L)) & medium
            prev = last_ev
        elif mode.startswith("att-lat"):
            # own input, but only the FIRST arrival inside each attention window
            ev = lat_edge & gate2d & armed & medium
            prev = last_ev
        else:                                     # "local"
            ev = (lat_edge | aff) & medium
            prev = last_ev

        elig_t[lat_edge] = t
        seen = ev & (prev >= 0)
        if mode == "att-3f":
            # three-factor: the gate supplies *when* (timestamp taken every window,
            # so intervals stay clean), the cell's own recent input supplies *who*.
            seen = seen & ((t - elig_t) <= P)
        if seen.any():
            dt = np.clip((t - prev[seen]).astype(np.float32), 0, 60.0)
            tau[seen] = np.clip(tau[seen] + ETA * (dt - tau[seen]), TMIN, TMAX)
            n_upd[seen] += 1
        if mode == "own":
            last_fire[ev] = t
        else:
            last_ev[ev] = t
        armed &= ~ev                              # disarm until the next window

    err = np.abs(tau - P).mean(0)[W_S:]
    tv = tau.mean(0)[W_S:]
    phase = np.where(gate_hits > 0, gate_phase_sum / np.maximum(gate_hits, 1), np.nan)
    frac_upd = float((n_upd[medium] > 0).mean())    # selectivity of the plasticity
    return err, tv, phase[W_S:], frac_upd


def depth(err):
    if not (err[0] < LOCK):
        return 0
    d = 0
    while d + 1 < len(err) and err[d + 1] < LOCK:
        d += 1
    return d + 1


def main():
    out = {}
    for P in PERIODS:
        for mode in MODES:
            E = np.zeros((N_SEEDS, L - W_S)); T = np.zeros((N_SEEDS, L - W_S))
            PH = np.zeros((N_SEEDS, L - W_S)); FU = np.zeros(N_SEEDS)
            for s in range(N_SEEDS):
                E[s], T[s], PH[s], FU[s] = run(mode, P, s)
            out[f"P{P}_{mode}"] = {"err": E, "tau": T, "phase": PH, "frac_upd": FU}
            e = E.mean(0)
            print(f"P={P:2d} {mode:10s}: |τ−P| x=0 {e[0]:5.2f} | x=8 {e[min(8,len(e)-1)]:5.2f} "
                  f"| x=24 {e[min(24,len(e)-1)]:5.2f} | far {e[-1]:5.2f} "
                  f"| depth {depth(e):3d}  (τ far {T.mean(0)[-1]:.1f}, "
                  f"cells updated {FU.mean():.2f})", flush=True)
        # propagation check: gate arrival phase must RISE with x, not be flat
        ph = out[f"P{P}_att-gate"]["phase"].mean(0)
        print(f"          gate arrival phase (mod P) at x=0/8/24/far: "
              f"{ph[0]:.1f} / {ph[min(8,len(ph)-1)]:.1f} / {ph[min(24,len(ph)-1)]:.1f} / {ph[-1]:.1f}",
              flush=True)
        print(flush=True)

    save = {"L": L, "W_S": W_S, "steps": STEPS, "n": N_SEEDS,
            "periods": np.array(PERIODS), "lock": LOCK}
    for k in out:
        for f in ("err", "tau", "phase", "frac_upd"):
            save[f"{k}_{f}"] = out[k][f]
    np.savez(os.path.join(OUT, "lattice_attention_gate.npz"), **save)
    with open(os.path.join(OUT, "lattice_attention_gate.json"), "w") as f:
        json.dump({"L": L, "W_S": W_S, "n": N_SEEDS, "periods": PERIODS,
                   "rows": {k: {"err_profile": out[k]["err"].mean(0).round(3).tolist(),
                                "depth": depth(out[k]["err"].mean(0)),
                                "tau_far": float(out[k]["tau"].mean(0)[-1]),
                                "frac_cells_updated": float(out[k]["frac_upd"].mean()),
                                "gate_phase": out[k]["phase"].mean(0).round(2).tolist()}
                            for k in out}},
                  f, indent=2)
    print("wrote lattice_attention_gate.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
