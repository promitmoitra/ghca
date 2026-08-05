"""Is avoidance sign-symmetric with approach?

`lattice_sensorimotor.py` established the action primitive: the medium is transmissive when
its cells have recovered, so a learned τ field is a schedule of when input gets through, and
the transmission edge sits at the learned interval (31.5 / 51.5 / 71.5 for D = 30 / 50 / 70).
That is *approach*: let the expected event through.

Avoidance asks for the opposite — **block** the event at D. The question is whether the same
machinery does it with a sign flip, or whether approach and avoidance are structurally
different problems.

There is a reason to expect the latter, and it is about the substrate rather than the rule.
A cell is refractory over `[t_fire, t_fire + τ]` and rested afterwards — refractoriness is a
**prefix** of the time following a firing, never a window in the middle of it. So a cell can
express "pass after T" by setting τ ≈ T − t_fire, but it has no way at all to express "block
near T and pass on either side". If that is right then no choice of learning rule repairs it,
and avoidance can only ever be implemented as *approach's complement* — block everything up
to T — rather than as selective avoidance.

The decisive measurement is therefore not whether transmission at D goes down. It is whether
any rule produces a **dip**: transmission low at D while still high at D−δ *and* D+δ. A
shifted step is not avoidance, it is postponed approach.

Four rules, identical in every other respect:

  approach       τ ← τ + η(dt − τ)                     the established rule
  avoid-flip     τ ← τ − η(dt − τ)                     the naive sign flip
  avoid-ratchet  on failure only (τ ≤ dt):              error-driven: push τ past the
                 τ ← τ + η((dt + M) − τ)               interval only when the event got through
  avoid-margin   τ ← τ + η((dt + M) − τ)               aim to be refractory by a margin M

Predictions recorded before running:
  * approach — a clean step up at D; τ ≈ interval; little mass at either bound.
  * avoid-flip — **bistable, not avoidant.** If τ < dt then (dt − τ) > 0 and the flipped rule
    *decreases* τ, driving the cell further below the interval, i.e. recovered and transmissive
    at D — the wrong answer. If τ > dt it increases τ, which is the right answer. So the rule
    should split the population to both bounds depending on initialisation, and roughly half of
    it should be actively wrong.
  * avoid-ratchet — τ → ceiling; blocks everything; no timing information survives, which is
    the e10 ratchet reappearing as the "correct" answer to a degenerate objective.
  * avoid-margin — works, but only by moving the step later by M. Needs an extra designed
    constant that approach does not, and still blocks everything before D + M.
  * **and no rule produces a dip**, because of the prefix argument above.
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("AV2_L", "64"))
WS = int(os.environ.get("AV2_WS", "3"))
PATCH = int(os.environ.get("AV2_PATCH", "4"))
ITI = int(os.environ.get("AV2_ITI", "260"))
N_TRIALS = int(os.environ.get("AV2_TRIALS", "40"))
N_SEEDS = int(os.environ.get("AV2_N", "3"))
D_LIST = [int(d) for d in os.environ.get("AV2_D", "30,50").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
MARG = float(os.environ.get("AV2_MARG", "10"))       # the avoidance margin (a designed constant)
THETA_M = 3.0
RESP_W = int(os.environ.get("AV2_RESPW", "12"))
DELTA = int(os.environ.get("AV2_DELTA", "10"))       # ±δ for the dip test
MODES = ("approach", "avoid-flip", "avoid-ratchet", "avoid-margin")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("AV2_TAG", "")


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
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    pa = patch_mask()
    for _ in range(N_TRIALS):
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
            if t == D:
                sel = t_fire >= 0
                if not sel.any():
                    continue
                dt = (t - t_fire[sel]).astype(np.float32)
                cur = tau[sel]
                if mode == "approach":
                    new = cur + ETA * (dt - cur)
                elif mode == "avoid-flip":
                    new = cur - ETA * (dt - cur)
                elif mode == "avoid-margin":
                    new = cur + ETA * ((dt + MARG) - cur)
                else:                                  # avoid-ratchet: only on failure
                    new = cur.copy()
                    failed = cur <= dt                 # recovered by the event, so it passed
                    new[failed] = (cur[failed]
                                   + ETA * ((dt[failed] + MARG) - cur[failed]))
                tau[sel] = np.clip(new, TMIN, TMAX)
    return tau


def _run(tau, t_probe, with_probe):
    pa = patch_mask()
    phi = np.zeros((L, L), np.int32)
    phi_m = np.zeros((L, L), np.int32)
    m_spikes = 0.0
    for t in range(min(ITI, t_probe + RESP_W + 1)):
        act = (phi >= 1) & (phi <= ACT)
        drive = np.zeros((L, L), bool)
        if t == 0:
            drive |= pa
        if with_probe and t == t_probe:
            drive |= pa
        excite = (phi == 0) & ((nbr4(act) >= THETA) | drive)
        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        exc_m = (phi_m == 0) & (box3(act) >= THETA_M)
        phi_m[phi_m >= 1] += 1
        phi_m[phi_m > int(ACT + 1)] = 0
        phi_m[exc_m] = 1
        if t >= t_probe:
            m_spikes += float(exc_m.sum())
    return m_spikes / (L * L)


def transmission(tau, t_probe):
    """Causal contribution of the probe (see lattice_sensorimotor.py for why differencing)."""
    return _run(tau, t_probe, True) - _run(tau, t_probe, False)


def main():
    out = {}
    for D in D_LIST:
        ts = sorted({max(1, D + o) for o in
                     (-30, -20, -DELTA, -4, 0, 4, DELTA, 20, 30, 45, 60)})
        ts = [t for t in ts if t < ITI - RESP_W]
        for mode in MODES:
            T = np.zeros((N_SEEDS, len(ts)))
            floor_f, ceil_f, taus = [], [], []
            for s in range(N_SEEDS):
                tau = train(mode, D, s)
                for i, tp in enumerate(ts):
                    T[s, i] = transmission(tau, tp)
                floor_f.append(float((tau <= TMIN + 1e-3).mean()))
                ceil_f.append(float((tau >= TMAX - 1e-3).mean()))
                taus.append(float(tau.mean()))
            m = T.mean(0)
            iD = ts.index(D) if D in ts else int(np.argmin([abs(t - D) for t in ts]))
            iL = int(np.argmin([abs(t - (D - DELTA)) for t in ts]))
            iR = int(np.argmin([abs(t - (D + DELTA)) for t in ts]))
            # A dip requires transmission low at D and high on BOTH sides, so the
            # denominator must be the MINIMUM of the flanks. Using their mean scores any
            # rising step as a dip -- it rated plain approach at +1.00.
            flank = min(m[iL], m[iR])
            dip = float(1.0 - m[iD] / flank) if flank > 1e-3 else float("nan")
            # monotonicity: the largest DECREASE anywhere along the curve. The prefix
            # argument says this must be ~0 for any tau field under transient dynamics.
            mono = float(max([m[i] - m[i + 1] for i in range(len(m) - 1)] + [0.0]))
            out[f"D{D}_{mode}"] = {"t": ts, "trans": T.tolist(), "dip": dip,
                                  "floor": floor_f, "ceil": ceil_f, "tau": taus,
                                  "mono_viol": mono}
            print(f"D={D:2d} {mode:14s}: transmission by probe time")
            print("            t " + " ".join(f"{t:5d}" for t in ts))
            print("        trans " + " ".join(f"{v:5.2f}" for v in m))
            print(f"        at D {m[iD]:5.2f} | flanks D±{DELTA} {m[iL]:5.2f}/{m[iR]:5.2f}"
                  f" | DIP {dip:+6.2f} | max drop {mono:+5.2f}"
                  f" | τ {np.mean(taus):5.1f} floor {np.mean(floor_f):.2f}"
                  f" ceil {np.mean(ceil_f):.2f}", flush=True)
        print(flush=True)

    print("Does any rule produce a DIP (block at D, pass on both sides)?", flush=True)
    print("  (dip = nan means the flanks were not both high, so no dip is even detectable)",
          flush=True)
    for mode in MODES:
        d = [out[f"D{D}_{mode}"]["dip"] for D in D_LIST]
        mv = [out[f"D{D}_{mode}"]["mono_viol"] for D in D_LIST]
        fl = [float(np.mean(out[f"D{D}_{mode}"]["floor"])) for D in D_LIST]
        ce = [float(np.mean(out[f"D{D}_{mode}"]["ceil"])) for D in D_LIST]
        print(f"  {mode:14s}: dip " + " / ".join(f"{v:+5.2f}" for v in d)
              + "   τ at floor " + " / ".join(f"{v:.2f}" for v in fl)
              + "   at ceiling " + " / ".join(f"{v:.2f}" for v in ce)
              + "   max drop " + " / ".join(f"{v:+.2f}" for v in mv), flush=True)

    with open(os.path.join(OUT, f"lattice_avoidance{TAG}.json"), "w") as f:
        json.dump({"L": L, "trials": N_TRIALS, "n": N_SEEDS, "delays": D_LIST,
                   "margin": MARG, "delta": DELTA, "rows": out}, f, indent=2)
    print(f"wrote lattice_avoidance{TAG}.json", flush=True)


if __name__ == "__main__":
    main()
