"""Validation for the interactive demo — does the input-timing τ rule tile *space*?

3e.2 grew a fast/slow τ hierarchy in a non-spatial **pool** (every unit fanned in from
both drive sources; a population conscience split the pool ~50/50). This asks the
spatial question the interactive demo needs, and which the repo has not run:

  On a 2D excitable lattice driven by **sparse pacemaker sites** of two different
  periods, does the input-timing τ rule organise τ into spatial **domains** — fast τ
  near fast pacemakers, slow τ near slow ones — where the old own-firing rule cannot?

Why sparse pacemakers rather than global drive: driving every cell at two global
periods makes the lattice flash in unison, with no wave propagation and no spatial
structure. Scattering pacemakers means waves emanate from them, so the *local input
statistics* genuinely vary across space — which is exactly what a rule that learns
from input timing should be able to tile. It is also cheaper in a shader (no global
reduction, hence no population conscience needed: the drive geometry supplies the
symmetry breaking that the conscience had to supply in the pool).

Implementation is deliberately grid-native (`np.roll` neighbour sums, per-cell state
arrays) so it ports 1:1 to a fragment shader for the browser demo.

Rules compared, identical in every other respect:
  * **input** (new) — τ tracks the interval between *suprathreshold input events the
    cell senses*, registered regardless of the cell's own refractory state.
  * **own** (old, e10) — τ tracks the interval between the cell's *own* firings; this
    is the self-referential signal that ratchets τ upward once it overshoots.

Metrics: Sarle bimodality of the τ field; fraction of τ near each pacemaker period;
and — the spatial claim — Moran's I spatial autocorrelation of τ, plus mean τ as a
function of distance-to-nearest-fast-pacemaker.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
import ghca_stats as st  # noqa: E402

L = int(os.environ.get("DEMO_L", "96"))            # lattice side
STEPS = int(os.environ.get("DEMO_STEPS", "6000"))
N_SEEDS = int(os.environ.get("DEMO_N", "5"))
P_F, P_S = 6, 24                                   # fast / slow pacemaker periods
PACE_FRAC = 0.004                                  # each class, fraction of cells
ACT = 2                                            # global active duration (not learned)
THETA = float(os.environ.get("DEMO_THETA", "1"))   # active neighbours needed to excite
LAYOUT = os.environ.get("DEMO_LAYOUT", "mixed")    # "mixed" | "split" pacemaker geometry
TMIN, TMAX = 3.0, 34.0
ETA = 0.15
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def neighbour_active_count(active):
    """4-neighbour count on a torus — the exact operation a shader does."""
    return (np.roll(active, 1, 0) + np.roll(active, -1, 0)
            + np.roll(active, 1, 1) + np.roll(active, -1, 1)).astype(np.float32)


def run(rule, seed, collect_frames=False):
    rng = np.random.default_rng(seed)
    phi = np.zeros((L, L), dtype=np.int32)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)

    n_pace = max(1, int(PACE_FRAC * L * L))
    fast_mask = np.zeros((L, L), bool); slow_mask = np.zeros((L, L), bool)
    if LAYOUT == "split":
        # fast pacemakers confined to the left half, slow to the right: local input
        # statistics then vary across space by construction, and the question is
        # whether τ tiles to match (with a boundary), not whether symmetry breaks.
        for mask, lo, hi in ((fast_mask, 0, L // 2), (slow_mask, L // 2, L)):
            rr = rng.integers(0, L, n_pace)
            cc = rng.integers(lo, hi, n_pace)
            mask[rr, cc] = True
    else:
        idx = rng.choice(L * L, 2 * n_pace, replace=False)
        f = np.zeros(L * L, bool); f[idx[:n_pace]] = True
        s = np.zeros(L * L, bool); s[idx[n_pace:]] = True
        fast_mask, slow_mask = f.reshape(L, L), s.reshape(L, L)
    pace = fast_mask | slow_mask
    # pacemakers must be able to keep their own period: a cell whose τ exceeds its
    # drive period would still be refractory when the next pulse arrives.
    tau[fast_mask] = P_F
    tau[slow_mask] = P_S

    last_in = np.full((L, L), -1.0, dtype=np.float32)    # last input event time
    last_fire = np.full((L, L), -1.0, dtype=np.float32)  # last own firing time
    prev_supra = np.zeros((L, L), bool)                  # for rising-edge detection
    frames = []

    for t in range(STEPS):
        active = (phi >= 1) & (phi <= ACT)
        rested = (phi == 0)
        nbr = neighbour_active_count(active)

        # pacemaker sites are force-driven on their own period
        drive = np.zeros((L, L), bool)
        if t % P_F == 0:
            drive |= fast_mask
        if t % P_S == 0:
            drive |= slow_mask

        # An "input event" is the *onset* of suprathreshold presynaptic drive, sensed
        # regardless of whether this cell can fire (the fix for the e10 ratchet). Edge
        # rather than level detection is essential on a lattice: a passing wave holds
        # the input suprathreshold for `act` consecutive steps, and with dense local
        # coupling some neighbour is nearly always active — level detection therefore
        # measures ~1 step and destroys the rhythm. In the 3e.2 pool this distinction
        # was invisible because the drive sources fired as discrete pulses, so every
        # input already *was* an edge.
        supra = (nbr >= THETA) | drive
        input_event = supra & ~prev_supra
        prev_supra = supra
        excite = rested & supra

        # advance clocks (per-cell tau), wrap, then apply excitations
        moving = phi >= 1
        phi[moving] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1

        # ---- timescale plasticity (pacemakers excluded: their period is imposed) ----
        if rule == "input":
            ev = input_event & ~pace
            prev = last_in
        else:                                            # "own" — the e10 rule
            ev = (phi == 1) & ~pace                      # fired this step
            prev = last_fire
        seen = ev & (prev >= 0)
        if seen.any():
            dt = (t - prev[seen]).astype(np.float32)
            dt = np.clip(dt, 0, P_S * 2.0)               # ignore implausibly long gaps
            tau[seen] = np.clip(tau[seen] + ETA * (dt - tau[seen]), TMIN, TMAX)
        if rule == "input":
            last_in[ev] = t
        else:
            last_fire[ev] = t

        if collect_frames and t % 20 == 0:
            frames.append(tau.copy())

    return tau, fast_mask, slow_mask, frames


def morans_i(field):
    """Spatial autocorrelation of `field` on a torus (4-neighbour weights)."""
    x = field - field.mean()
    num = (x * (np.roll(x, 1, 0) + np.roll(x, -1, 0)
                + np.roll(x, 1, 1) + np.roll(x, -1, 1))).sum()
    den = (x * x).sum() * 4.0
    return float(num / den) if den > 0 else 0.0


def ring_distance_profile(tau, fast_mask, slow_mask, max_d=12):
    """Mean τ as a function of Chebyshev distance to the nearest fast / slow pacemaker."""
    def dist_to(mask):
        d = np.full((L, L), 1e9, np.float32)
        d[mask] = 0
        for _ in range(max_d):                      # BFS-ish dilation on a torus
            d = np.minimum.reduce([d,
                                   np.roll(d, 1, 0) + 1, np.roll(d, -1, 0) + 1,
                                   np.roll(d, 1, 1) + 1, np.roll(d, -1, 1) + 1])
        return d
    df, ds = dist_to(fast_mask), dist_to(slow_mask)
    closer_fast = df < ds
    return (float(tau[closer_fast].mean()) if closer_fast.any() else float("nan"),
            float(tau[~closer_fast].mean()) if (~closer_fast).any() else float("nan"))


def main():
    out = {}
    for rule in ("input", "own"):
        rows = {"bc": [], "near_f": [], "near_s": [], "moran": [],
                "tau_near_fast": [], "tau_near_slow": [], "mean": []}
        tau0 = None
        for s in range(N_SEEDS):
            tau, fm, sm, _ = run(rule, s)
            if s == 0:
                tau0 = tau.copy()
            flat = tau.ravel()
            rows["bc"].append(st.bimodality(flat)["bc"])
            rows["near_f"].append(float((np.abs(flat - P_F) <= 2).mean()))
            rows["near_s"].append(float((np.abs(flat - P_S) <= 3).mean()))
            rows["moran"].append(morans_i(tau))
            nf, ns = ring_distance_profile(tau, fm, sm)
            rows["tau_near_fast"].append(nf); rows["tau_near_slow"].append(ns)
            rows["mean"].append(float(flat.mean()))
        out[rule] = {k: np.array(v) for k, v in rows.items()}
        out[rule]["tau0"] = tau0
        m = out[rule]
        print(f"=== rule={rule} (L={L}, {STEPS} steps, n={N_SEEDS}) ===", flush=True)
        print(f"  mean τ            {m['mean'].mean():.1f}", flush=True)
        print(f"  near P_f={P_F}       {m['near_f'].mean():.2f}", flush=True)
        print(f"  near P_s={P_S}      {m['near_s'].mean():.2f}", flush=True)
        print(f"  bimodality BC     {m['bc'].mean():.2f}", flush=True)
        print(f"  Moran's I (space) {m['moran'].mean():.3f}", flush=True)
        print(f"  τ near fast pace  {m['tau_near_fast'].mean():.1f}"
              f"   vs near slow pace  {m['tau_near_slow'].mean():.1f}"
              f"   (Δ={m['tau_near_slow'].mean()-m['tau_near_fast'].mean():+.1f})", flush=True)

    np.savez(os.path.join(OUT, "lattice_timescale_demo.npz"), L=L, steps=STEPS,
             P_F=P_F, P_S=P_S, n=N_SEEDS,
             **{f"{r}_{k}": out[r][k] for r in out for k in out[r]})
    with open(os.path.join(OUT, "lattice_timescale_demo.json"), "w") as f:
        json.dump({"L": L, "steps": STEPS, "P_F": P_F, "P_S": P_S, "n": N_SEEDS,
                   "rows": {r: {k: float(np.mean(v)) for k, v in out[r].items()
                                if k != "tau0"} for r in out}}, f, indent=2)
    print("\nwrote lattice_timescale_demo.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
