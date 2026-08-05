"""Foundations for a sensorimotor loop: what IS the action, physically?

`lattice_layers.py` reported a "motor readout" — H cells recovering synchronously at the
predicted reward time (91% of recovery events within ±3 steps). Read carefully, that is not
an output. A Greenberg-Hastings cell returning to rest emits nothing; recovery is only
visible to an observer counting it. The scalar `rec_hist` was a statistic, not a population,
and calling it motor was generous.

What the synchronous recovery actually *is* is a synchronous window in which the medium
becomes **transmissive**. Refractoriness gates propagation — that is the one thing excitable
media do natively — so a learned τ field is a learned schedule of *when input can get
through*. That reframes the action primitive:

    not  "emit at the expected time"
    but  "be transmissive at the expected time"

which makes the motor layer definable without inventing new machinery: M is a population
downstream of H that fires when H transmits to it. The medium acts by passing or blocking.

This experiment tests the reframing, and lays the interface the environment will need:

  1. Train H exactly as in the winning layered condition (stimulus patch at t=0, diffuse
     value signal at t=D, τ ← τ + η((t_now − t_fire) − τ), ungated).
  2. Then, on separate probe trials, re-present the stimulus at a swept time t_probe and
     measure how much of it gets through — H's evoked firing, and the spikes of a real
     motor sheet M (an L×L sheet of the same cells, 3x3 convergent input from H, θ_M > 1,
     no plasticity of its own).
  3. Sweep t_probe across and past D.

**Prediction, and note it is a step rather than a peak.** Before D the trained cells are
still refractory, so the probe wave is blocked and M is silent. At D they recover together
and the probe propagates, so M fires. After D they simply stay rested — nothing re-refracts
them — so transmission stays high. The signature is therefore a **transmission edge whose
position is the learned interval**, not a band-pass tuning curve. Getting a genuine peak
would need tonic drive to re-refract the sheet, which is a later question, not this one.

**Why this matters for bootstrapping.** Roadmap item 13 worried that contingent reward
cannot get started: if reward requires a well-timed action and the action is initially
untimed, reward never arrives and nothing learns. If the action is transmission, that worry
may be unfounded — an untrained sheet has τ scattered across its whole range, so at *any*
probe time some fraction of cells is excitable and transmission is partial rather than zero.
Partial transmission is a graded signal, which is exactly what shaping needs. The measurement
here is whether the untrained transmission curve is genuinely intermediate (bootstrapping is
free) or flat-and-saturated (it is not, and shaping has to be designed).

Conditions: trained / untrained (random τ) / unpaired (τ trained on random reward times).
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("SM_L", "64"))
WS = int(os.environ.get("SM_WS", "3"))
PATCH = int(os.environ.get("SM_PATCH", "4"))
ITI = int(os.environ.get("SM_ITI", "260"))
N_TRIALS = int(os.environ.get("SM_TRIALS", "40"))
N_SEEDS = int(os.environ.get("SM_N", "3"))
D_LIST = [int(d) for d in os.environ.get("SM_D", "30,50,70").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_M = float(os.environ.get("SM_THETAM", "3"))   # convergence needed to fire a motor cell
RESP_W = int(os.environ.get("SM_RESPW", "12"))      # window after the probe to count response
MODES = ("trained", "untrained", "unpaired")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def nbr4(a):
    n = np.zeros(a.shape, np.float32)
    n[1:, :] += a[:-1, :]
    n[:-1, :] += a[1:, :]
    n[:, 1:] += a[:, :-1]
    n[:, :-1] += a[:, 1:]
    return n


def box3(a):
    f = a.astype(np.float32)
    p = np.zeros((a.shape[0] + 2, a.shape[1] + 2), np.float32)
    p[1:-1, 1:-1] = f
    return (p[:-2, :-2] + p[:-2, 1:-1] + p[:-2, 2:] +
            p[1:-1, :-2] + p[1:-1, 1:-1] + p[1:-1, 2:] +
            p[2:, :-2] + p[2:, 1:-1] + p[2:, 2:])


def patch_mask():
    y0 = L // 2
    p = np.zeros((L, L), bool)
    p[y0 - PATCH:y0 + PATCH + 1, :WS] = True
    return p


def train(mode, D, seed):
    """Return a τ field trained under `mode` (the winning ungated layered condition)."""
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    if mode == "untrained":
        return tau
    pa = patch_mask()
    for _ in range(N_TRIALS):
        d_rew = int(rng.integers(TMIN + 2, ITI // 2)) if mode == "unpaired" else D
        phi = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            excite = (phi == 0) & ((nbr4(act) >= THETA) | (pa if t == 0 else False))
            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t
            if t == d_rew:                       # diffuse value signal: reward now
                sel = t_fire >= 0
                if sel.any():
                    dt = (t - t_fire[sel]).astype(np.float32)
                    tau[sel] = np.clip(tau[sel] + ETA * (dt - tau[sel]), TMIN, TMAX)
    return tau


def probe(tau, t_probe):
    """The probe's *causal* contribution: run the trial with and without it and difference.

    Counting raw activity after t_probe does not work — the conditioning wave launched at t=0
    is still crossing the sheet, and its front advances into never-fired cells regardless of τ,
    so it swamps the probe response and is identical across conditions. The dynamics are
    deterministic given τ, so the with-minus-without difference is exactly what the probe caused.
    """
    a = _run(tau, t_probe, True)
    b = _run(tau, t_probe, False)
    return a[0] - b[0], a[1] - b[1]


def _run(tau, t_probe, with_probe):
    pa = patch_mask()
    phi = np.zeros((L, L), np.int32)
    phi_m = np.zeros((L, L), np.int32)
    t_fire = np.full((L, L), -1, np.int32)
    evoked, m_spikes = 0.0, 0.0
    t_end = min(ITI, t_probe + RESP_W + 1)
    for t in range(t_end):
        act = (phi >= 1) & (phi <= ACT)
        drive = np.zeros((L, L), bool)
        if t == 0:
            drive |= pa                          # the conditioning stimulus
        if with_probe and t == t_probe:
            drive |= pa                          # the probe, at the swept time
        excite = (phi == 0) & ((nbr4(act) >= THETA) | drive)
        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        newly = excite & (t_fire < 0)
        t_fire[newly] = t

        # motor sheet: real cells, convergent input from H, no plasticity of its own
        act_m = (phi_m >= 1) & (phi_m <= ACT)
        exc_m = (phi_m == 0) & (box3(act) >= THETA_M)
        phi_m[phi_m >= 1] += 1
        phi_m[phi_m > int(ACT + 1)] = 0
        phi_m[exc_m] = 1

        if t >= t_probe:                         # only count what the PROBE evoked
            evoked += float(excite.sum())
            m_spikes += float(exc_m.sum())
    return evoked / (L * L), m_spikes / (L * L)


def main():
    out = {}
    for D in D_LIST:
        offs = [-24, -16, -10, -6, -3, 0, 3, 6, 10, 16, 24, 40]
        ts = [D + o for o in offs if 1 <= D + o < ITI - RESP_W]
        for mode in MODES:
            ev = np.zeros((N_SEEDS, len(ts)))
            ms = np.zeros((N_SEEDS, len(ts)))
            for s in range(N_SEEDS):
                tau = train(mode, D, s)
                for i, tp in enumerate(ts):
                    ev[s, i], ms[s, i] = probe(tau, tp)
            out[f"D{D}_{mode}"] = {"t": ts, "evoked": ev.tolist(), "motor": ms.tolist()}
            e, m = ev.mean(0), ms.mean(0)
            # where does transmission cross halfway between its floor and ceiling?
            lo, hi = m.min(), m.max()
            half = lo + 0.5 * (hi - lo)
            cross = float("nan")
            for i in range(1, len(m)):
                if m[i - 1] < half <= m[i]:
                    f = (half - m[i - 1]) / max(m[i] - m[i - 1], 1e-9)
                    cross = ts[i - 1] + f * (ts[i] - ts[i - 1])
                    break
            out[f"D{D}_{mode}"]["cross"] = cross
            print(f"D={D:2d} {mode:9s}: motor spikes/cell by probe time")
            print("           t " + " ".join(f"{t:5d}" for t in ts))
            print("       motor " + " ".join(f"{v:5.2f}" for v in m))
            print("      evoked " + " ".join(f"{v:5.2f}" for v in e))
            print(f"      half-max crossing {cross:6.1f}  (reward time {D})"
                  f"   floor {lo:.2f} ceiling {hi:.2f}", flush=True)
        print(flush=True)

    print("does the transmission edge sit at the learned interval?", flush=True)
    for mode in MODES:
        c = [out[f"D{D}_{mode}"]["cross"] for D in D_LIST]
        print(f"  {mode:9s}: crossing " + " / ".join(f"{v:6.1f}" for v in c)
              + "   vs D " + " / ".join(f"{D:6d}" for D in D_LIST), flush=True)

    with open(os.path.join(OUT, "lattice_sensorimotor.json"), "w") as f:
        json.dump({"L": L, "trials": N_TRIALS, "n": N_SEEDS, "delays": D_LIST,
                   "theta_m": THETA_M, "resp_w": RESP_W, "rows": out}, f, indent=2)
    print("wrote lattice_sensorimotor.json", flush=True)


if __name__ == "__main__":
    main()
