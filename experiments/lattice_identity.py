"""Identity as a plastic variable: can a depth gradient emerge instead of being wired?

The whole layered arc still hand-draws its anatomy — which cells are wave-carriers, which are
coincidence detectors, which read out. The proposal is to make *identity* itself a plastic cell
property and let structure organise from use, with the discipline agreed in advance: **one rate,
one context dependence**, and a falsifiable prediction committed before running.

Identity here is the firing threshold **θ**. This is not an arbitrary choice — it is the one cell
property this repo has already tied to function. Negative 1 established that θ=1 propagates planar
waves and θ≥2 blocks them, so θ *is* the propagator↔coincidence-detector axis. A low-θ cell
carries waves; a high-θ cell fires only on convergent input (an A/M-style reader). Letting θ move
is letting a cell change what kind of cell it is.

The transition rule is homeostatic intrinsic plasticity — the minimal one-rate/one-context law:

    ā_i ← (1−λ) ā_i + λ · [cell i active this step]          (its own recent activity)
    θ_i ← clip( θ_i + r · (ā_i − a*),  1,  θ_max )            (one rate r, one context ā_i)

`a*` is a single activity set-point; `ā_i` is the one context signal (the cell's own activity, so
no extra spatial-averaging parameter). The rule is **self-limiting by construction**: too active →
θ up → fires less; too quiet → θ down → fires more. That is the conserved-flow property argued for
over birth–death — there is no runaway direction, because the feedback opposes the deviation.

The sweep is over r, and it doubles as the timescale-separation test promised earlier. τ also
learns concurrently (the input-timing rule tracking the drive period P), so identity churn and
timing learning share the substrate: what counts as a timing edge depends on θ, so if θ moves fast
it moves the ground under τ. Sweeping r therefore measures both whether identity self-organises and
whether it must be slow relative to learning.

Init θ = 1 everywhere (all propagators). Drive is a periodic patch at the left edge, so activity is
high near it and decays with distance x — a stable geometry for identity to organise against.

PREDICTIONS, committed before running (house rule: state them, then let the data judge):

  1. Differentiation. θ rises where drive is high (near the patch) and stays at the floor where it
     is low (far), so a **spatial gradient emerges** — corr(θ, x) negative — with no wiring. If θ
     stays flat at every r, homeostatic identity does not self-organise and the idea is wrong.
  2. Fixed point vs speed. The gradient's *magnitude* is set by geometry (where ā = a*), not by r;
     r sets how fast it is reached. So small r under-converges within the run, medium r reaches a
     stable gradient, large r overshoots and **oscillates** (θ churn stays high late in the run).
     Expect an inverted-U in stability.
  3. Timescale separation. τ-learning stays intact for small/medium r and **degrades above some
     r\*** where θ churns faster than τ can exploit a stable substrate. If |τ−P| is flat across all
     r, identity plasticity is free and the "must be slow" worry was wrong — a real possible
     outcome, recorded as such.
  4. Function. Wave propagation to the far half is preserved in the window and may destabilise at
     very high r (the sheet's excitability map never settles).
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("ID_L", "64"))
WS = int(os.environ.get("ID_WS", "3"))
STEPS = int(os.environ.get("ID_STEPS", "6000"))
N_SEEDS = int(os.environ.get("ID_N", "3"))
P = int(os.environ.get("ID_P", "10"))               # drive period
ACT = 2.0
TMIN, TMAX, ETA = 3.0, 40.0, 0.15                   # τ range + input-timing learning rate
THETA_MIN, THETA_MAX = 1.0, 4.0                     # the identity axis
A_STAR = float(os.environ.get("ID_ASTAR", "0.08"))  # activity set-point (the one context constant)
LAM = float(os.environ.get("ID_LAM", "0.02"))       # activity EMA rate
THETA_INIT = float(os.environ.get("ID_INIT", "1.0"))  # starting identity (floor=propagator)
R_LIST = [float(x) for x in os.environ.get("ID_R", "0,3,10,30,100,300").split(",")]
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("ID_TAG", "")


def nbr4(a):
    n = np.zeros(a.shape, np.float32)
    n[1:, :] += a[:-1, :]
    n[:-1, :] += a[1:, :]
    n[:, 1:] += a[:, :-1]
    n[:, :-1] += a[:, 1:]
    return n


def run(r, seed):
    rng = np.random.default_rng(seed)
    phi = np.zeros((L, L), np.int32)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    theta = np.full((L, L), THETA_INIT, np.float32)       # starting identity
    abar = np.zeros((L, L), np.float32)                   # EMA of own activity
    prev_supra = np.zeros((L, L), bool)
    last_edge = np.full((L, L), -1.0, np.float32)

    xs = np.arange(L)[None, :].repeat(L, 0).astype(np.float32)
    patch = xs < WS
    theta_hist = []                                       # snapshots for late-churn measurement

    for t in range(STEPS):
        active = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        nbr = nbr4(active)
        drive = patch & ((t % P) == 0)

        supra = nbr >= theta                              # per-cell threshold == identity
        excite = rested & (supra | drive)

        # input-timing τ rule: learn from the rising edge of suprathreshold input
        edge = supra & ~prev_supra
        prev_supra = supra
        seen = edge & (last_edge >= 0)
        if seen.any():
            dt = np.clip(t - last_edge[seen], 0, 60.0)
            tau[seen] = np.clip(tau[seen] + ETA * (dt - tau[seen]), TMIN, TMAX)
        last_edge[edge] = t

        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1

        # homeostatic identity: one rate r, one context ā (own activity) vs set-point a*
        abar = (1 - LAM) * abar + LAM * active.astype(np.float32)
        if r > 0:
            theta = np.clip(theta + r * LAM * (abar - A_STAR), THETA_MIN, THETA_MAX)

        if t % 200 == 0 and t > STEPS // 2:
            theta_hist.append(theta.copy())

    # --- metrics ---------------------------------------------------------------------
    x_of = np.arange(L)[None, :].repeat(L, 0).astype(np.float32)
    th = theta.reshape(-1)
    xx = x_of.reshape(-1)
    corr = float(np.corrcoef(th, xx)[0, 1]) if th.std() > 1e-6 else 0.0
    churn = (float(np.mean([np.abs(theta_hist[i + 1] - theta_hist[i]).mean()
                            for i in range(len(theta_hist) - 1)]))
             if len(theta_hist) > 1 else float("nan"))
    reach = float(active_reach(phi, tau, theta, patch, rng))
    # τ error over cells that actually track the drive (near the patch, where a rhythm exists)
    near = (x_of < L // 3) & ~patch
    tau_err = float(np.abs(tau[near] - P).mean())
    return {
        "theta_mean": float(theta.mean()), "theta_std": float(theta.std()),
        "corr_theta_x": corr, "churn_late": churn,
        "act_mean": float(abar.mean()), "diff_frac": float((theta > THETA_MIN + 0.3).mean()),
        "tau_err": tau_err, "reach": reach,
    }


def active_reach(phi, tau, theta, patch, rng):
    """After learning, does a fresh drive wave still reach the far half of the sheet?"""
    ph = np.zeros((L, L), np.int32)
    prev = np.zeros((L, L), bool)
    xs = np.arange(L)[None, :].repeat(L, 0)
    far = xs >= (2 * L // 3)
    far_hits = np.zeros((L, L), bool)
    for t in range(P * 6):
        act = (ph >= 1) & (ph <= ACT)
        drive = patch & ((t % P) == 0)
        excite = (ph == 0) & ((nbr4(act) >= theta) | drive)
        ph[ph >= 1] += 1
        ph[ph > np.ceil(tau).astype(np.int32)] = 0
        ph[excite] = 1
        far_hits |= act & far
    return far_hits.mean() / max(far.mean(), 1e-9)


def main():
    out = {}
    print(f"identity sweep: L={L} P={P} steps={STEPS} a*={A_STAR} λ={LAM} n={N_SEEDS}", flush=True)
    print(f"init θ={THETA_MIN} (all propagators), range [{THETA_MIN},{THETA_MAX}]\n", flush=True)
    for r in R_LIST:
        rs = [run(r, s) for s in range(N_SEEDS)]
        agg = {k: float(np.nanmean([x[k] for x in rs])) for k in rs[0]}
        out[f"r{r:g}"] = {k: [x[k] for x in rs] for k in rs[0]}
        print(f"r={r:6g}: θ̄ {agg['theta_mean']:4.2f} ±{agg['theta_std']:4.2f}"
              f" | diff {agg['diff_frac']:4.2f}"
              f" | corr(θ,x) {agg['corr_theta_x']:+5.2f}"
              f" | churn {agg['churn_late']:.4f}"
              f" | act {agg['act_mean']:.3f} (a*={A_STAR})"
              f" | |τ−P| {agg['tau_err']:4.1f}"
              f" | reach {agg['reach']:.2f}", flush=True)

    with open(os.path.join(OUT, f"lattice_identity{TAG}.json"), "w") as f:
        json.dump({"L": L, "P": P, "steps": STEPS, "a_star": A_STAR, "lam": LAM,
                   "n": N_SEEDS, "theta_range": [THETA_MIN, THETA_MAX],
                   "r_list": R_LIST, "rows": out}, f, indent=2)
    print(f"\nwrote lattice_identity{TAG}.json", flush=True)


if __name__ == "__main__":
    main()
