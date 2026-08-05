"""Tonic drive: one structural change under two stalled threads.

Two separate dead ends in this arc trace to the same assumption.

  * `lattice_sensorimotor.py` — the transmission curve is a monotone **step**, not a band-pass
    peak. The medium becomes transmissive at the learned interval and then *stays* transmissive,
    so it expresses "pass after T" rather than "pass at T".
  * `lattice_avoidance.py` — selective avoidance is impossible, and the argument is a theorem:
    for a dip you need a cell rested at T−δ and refractory at T, but a rested cell **stays**
    rested because nothing re-fires it. Hence transmission is monotone non-decreasing in probe
    time (measured max decrease across all rules and delays: +0.00).

Both rest on **single firing**: each cell fires once per trial and then sits rested forever.
Weak tonic drive breaks exactly that. A rested cell that receives a little background excitation
fires again, re-enters refractoriness, and becomes periodically — rather than monotonically —
excitable. If the assumption is what binds, one structural change should unlock both a band-pass
transmission peak and a genuine avoidance dip.

The risk is equally specific, and it is this arc's oldest pathology. A self-re-firing medium is
one where each cell's input arrives at interval ≈ its own τ, which is precisely the
**self-confirming fixed point** that sank the lateral-input rule (Negative 2): τ locks and stops
tracking the drive because the rule measures its own output. So tonic drive may buy band-pass
timing at the cost of the timing being learnable at all. That trade is the real content here, and
it is measured directly rather than assumed away.

Design. Identical to the sensorimotor setup — stimulus patch at t=0, diffuse value signal at
t=D, τ ← τ + η((t_now − t_fire) − τ) — plus a tonic term. On each step every rested cell receives
background excitation with probability `q`. Sweeping q from 0 (the established regime) upward is
the one knob. Learning uses the FIRST firing per trial as its eligibility timestamp, so the τ rule
is unchanged by tonic drive; only the dynamics change.

Measurements, each aimed at one claim:
  * `peakiness`  — band-pass: (transmission at D) / (transmission at D+40). >1 means a peak; ~1
                   a step; the step regime gives ~1 by construction.
  * `max_drop`   — the monotonicity violation. Under single firing this is provably 0. Any
                   reliably positive value means tonic drive broke the theorem's assumption.
  * `dip`        — with the avoid-margin rule: transmission at D relative to the MINIMUM of both
                   flanks (using the mean scores any rising step as a dip — the instrument error
                   caught in `lattice_avoidance.py`).
  * `tau_err`    — |τ − true interval|, the cost side of the trade. If this degrades with q, the
                   self-confirming fixed point is reasserting itself.
  * `reentrant`  — fraction of steps with activity late in the trial, long after the stimulus
                   wave should have died. Diagnoses runaway self-sustaining activity.

PREDICTIONS, committed before running:
  1. `max_drop` rises above 0 for q > 0 — the monotonicity theorem's assumption is broken, as its
     own statement said it would be. This is close to a sanity check; if it fails, the tonic term
     is not doing anything and the implementation is wrong.
  2. `peakiness` rises above 1 over some intermediate q range. A cell that re-fires becomes
     periodically excitable, so transmission should acquire structure rather than a single edge.
  3. `tau_err` degrades as q grows, because re-firing at interval ≈ τ is the self-confirming
     signal. **I expect a narrow or empty window** where 2 holds and 3 has not yet broken — the
     honest prior after four failures of this shape, and if the window is empty that is the result.
  4. `reentrant` grows with q; above some q the sheet never quiets and trials stop being separable.
  5. The avoidance dip becomes non-nan for some q > 0. If it never does, selective avoidance is
     unreachable even with re-firing, and the structural claim is stronger than first stated.
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("TO_L", "64"))
WS = int(os.environ.get("TO_WS", "3"))
PATCH = int(os.environ.get("TO_PATCH", "4"))
ITI = int(os.environ.get("TO_ITI", "260"))
N_TRIALS = int(os.environ.get("TO_TRIALS", "40"))
N_SEEDS = int(os.environ.get("TO_N", "5"))
D = int(os.environ.get("TO_D", "50"))
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_M = 3.0
RESP_W = int(os.environ.get("TO_RESPW", "12"))
DELTA = int(os.environ.get("TO_DELTA", "10"))
MARG = float(os.environ.get("TO_MARG", "10"))
Q_LIST = [float(x) for x in os.environ.get(
    "TO_Q", "0,1e-7,3e-7,1e-6,3e-6,1e-5,3e-5,1e-4").split(",")]
RULES = ("approach", "avoid-margin")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("TO_TAG", "")


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


def train(rule, q, seed):
    """Train τ with tonic drive present. rng threaded explicitly (house rule)."""
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    pa = patch_mask()
    reentrant = []
    for _ in range(N_TRIALS):
        phi = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        late_act = 0
        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            rested = phi == 0
            tonic = (rng.random((L, L)) < q) if q > 0 else False
            excite = rested & ((nbr4(act) >= THETA) | (pa if t == 0 else False) | tonic)
            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t              # FIRST firing only: the τ rule is unchanged
            if t > ITI - 60:
                late_act += int(act.sum() > 0)
            if t == D:
                sel = t_fire >= 0
                if not sel.any():
                    continue
                dt = (t - t_fire[sel]).astype(np.float32)
                cur = tau[sel]
                if rule == "approach":
                    new = cur + ETA * (dt - cur)
                else:                       # avoid-margin
                    new = cur + ETA * ((dt + MARG) - cur)
                tau[sel] = np.clip(new, TMIN, TMAX)
        reentrant.append(late_act / 60.0)
    return tau, float(np.mean(reentrant))


def _run(tau, q, t_probe, with_probe, seed, want_sat=False):
    """One probe trial. Same rng stream for both arms so the difference is causal."""
    rng = np.random.default_rng(10_000 + seed)
    pa = patch_mask()
    phi = np.zeros((L, L), np.int32)
    phi_m = np.zeros((L, L), np.int32)
    m = 0.0
    sat_acc = []
    for t in range(min(ITI, t_probe + RESP_W + 1)):
        act = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        tonic = (rng.random((L, L)) < q) if q > 0 else False
        drive = np.zeros((L, L), bool)
        if t == 0:
            drive |= pa
        if with_probe and t == t_probe:
            drive |= pa
        excite = rested & ((nbr4(act) >= THETA) | drive | tonic)
        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        exc_m = (phi_m == 0) & (box3(act) >= THETA_M)
        phi_m[phi_m >= 1] += 1
        phi_m[phi_m > int(ACT + 1)] = 0
        phi_m[exc_m] = 1
        if t >= t_probe:
            m += float(exc_m.sum())
            sat_acc.append(float(act.mean()))
    if want_sat:
        return m / (L * L), float(np.mean(sat_acc)) if sat_acc else 0.0
    return m / (L * L)


def alignment(tau, q, seed, offset=0.0):
    """|(t_rew - t_fire) - tau| measured on a probe trial, as in the rest of the arc.

    The earlier version of this metric compared tau against an analytic guess (D - x), which is
    wrong for any cell off the patch row since t_fire ~ x + |dy|. It therefore could not see
    learning quality at all. This measures the real interval each cell experienced.
    """
    rng = np.random.default_rng(20_000 + seed)
    pa = patch_mask()
    phi = np.zeros((L, L), np.int32)
    t_fire = np.full((L, L), -1, np.int32)
    for t in range(D + 1):
        act = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        tonic = (rng.random((L, L)) < q) if q > 0 else False
        excite = rested & ((nbr4(act) >= THETA) | (pa if t == 0 else False) | tonic)
        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        newly = excite & (t_fire < 0)
        t_fire[newly] = t
        if t == D:
            sel = t_fire >= 0
            if not sel.any():
                return float("nan"), 0.0
            dt = (D - t_fire[sel]).astype(np.float32) + offset
            e = np.abs(dt - tau[sel])
            return float(e.mean()), float((e <= 2.0).mean())
    return float("nan"), 0.0


def transmission(tau, q, t_probe, seed):
    """Probe's causal contribution, plus baseline activity so 0.0 is not ambiguous.

    A saturated sheet also yields transmission ~0 -- not because the probe was blocked but
    because every cell is already firing and there is nothing left to recruit. Returning the
    no-probe arm's mean activity distinguishes the two.
    """
    with_p = _run(tau, q, t_probe, True, seed)
    without_p, sat = _run(tau, q, t_probe, False, seed, want_sat=True)
    return with_p - without_p, sat


def main():
    ts = sorted({max(1, D + o) for o in (-30, -20, -DELTA, -4, 0, 4, DELTA, 20, 40, 60)})
    ts = [t for t in ts if t < ITI - RESP_W]
    iD = ts.index(D)
    iL = int(np.argmin([abs(t - (D - DELTA)) for t in ts]))
    iR = int(np.argmin([abs(t - (D + DELTA)) for t in ts]))
    iFar = int(np.argmin([abs(t - (D + 40)) for t in ts]))
    out = {}
    print(f"tonic sweep: L={L} D={D} trials={N_TRIALS} n={N_SEEDS}", flush=True)
    print(f"probe times: {ts}\n", flush=True)

    for rule in RULES:
        for q in Q_LIST:
            T = np.zeros((N_SEEDS, len(ts)))
            S = np.zeros((N_SEEDS, len(ts)))
            errs, fracs, reents = [], [], []
            for s in range(N_SEEDS):
                tau, reent = train(rule, q, s)
                for i, tp in enumerate(ts):
                    T[s, i], S[s, i] = transmission(tau, q, tp, s)
                # cost side, measured not guessed: does τ still encode the interval each
                # cell actually experienced?
                e, fr = alignment(tau, q, s,
                                  MARG if rule == "avoid-margin" else 0.0)
                errs.append(e)
                fracs.append(fr)
                reents.append(reent)
            mt = T.mean(0)
            sd = T.std(0)
            peak = float(mt[iD] / mt[iFar]) if mt[iFar] > 1e-3 else float("nan")
            drop = float(max([mt[i] - mt[i + 1] for i in range(len(mt) - 1)] + [0.0]))
            flank = min(mt[iL], mt[iR])
            dip = float(1.0 - mt[iD] / flank) if flank > 1e-3 else float("nan")
            satm = float(S.mean())
            out[f"{rule}_q{q:g}"] = {"t": ts, "trans": T.tolist(), "peakiness": peak,
                                     "max_drop": drop, "dip": dip, "sat": satm,
                                     "tau_err": errs, "aligned_frac": fracs,
                                     "reentrant": reents}
            print(f"{rule:12s} q={q:<7g} " + " ".join(f"{v:5.2f}" for v in mt), flush=True)
            print(f"{'':12s}   ±sd     " + " ".join(f"{v:5.2f}" for v in sd), flush=True)
            print(f"{'':12s}   peakiness {peak:5.2f} | max_drop {drop:+5.2f}"
                  f" | dip {dip:+6.2f}"
                  f" | τ_err {np.nanmean(errs):5.1f}±{np.nanstd(errs):4.1f}"
                  f" (within±2 {np.nanmean(fracs):.2f})"
                  f" | reentrant {np.mean(reents):.2f}"
                  f" | SAT {satm:.3f}", flush=True)
        print(flush=True)

    print("=== the trade: does band-pass timing cost learnability? ===", flush=True)
    for rule in RULES:
        pk = [out[f"{rule}_q{q:g}"]["peakiness"] for q in Q_LIST]
        dr = [out[f"{rule}_q{q:g}"]["max_drop"] for q in Q_LIST]
        te = [float(np.nanmean(out[f"{rule}_q{q:g}"]["tau_err"])) for q in Q_LIST]
        wf = [float(np.nanmean(out[f"{rule}_q{q:g}"]["aligned_frac"])) for q in Q_LIST]
        st = [out[f"{rule}_q{q:g}"]["sat"] for q in Q_LIST]
        dp = [out[f"{rule}_q{q:g}"]["dip"] for q in Q_LIST]
        rn = [float(np.mean(out[f"{rule}_q{q:g}"]["reentrant"])) for q in Q_LIST]
        print(f"  {rule:12s} q      " + " ".join(f"{q:>6g}" for q in Q_LIST), flush=True)
        print(f"  {'':12s} peaky  " + " ".join(f"{v:6.2f}" for v in pk), flush=True)
        print(f"  {'':12s} drop   " + " ".join(f"{v:+6.2f}" for v in dr), flush=True)
        print(f"  {'':12s} τ_err  " + " ".join(f"{v:6.1f}" for v in te), flush=True)
        print(f"  {'':12s} w±2    " + " ".join(f"{v:6.2f}" for v in wf), flush=True)
        print(f"  {'':12s} sat    " + " ".join(f"{v:6.3f}" for v in st), flush=True)
        print(f"  {'':12s} dip    " + " ".join(f"{v:6.2f}" for v in dp), flush=True)
        print(f"  {'':12s} reent  " + " ".join(f"{v:6.2f}" for v in rn), flush=True)

    with open(os.path.join(OUT, f"lattice_tonic{TAG}.json"), "w") as f:
        json.dump({"L": L, "D": D, "trials": N_TRIALS, "n": N_SEEDS, "q_list": Q_LIST,
                   "delta": DELTA, "margin": MARG, "rows": out}, f, indent=2)
    print("\nwrote lattice_tonic.json", flush=True)


if __name__ == "__main__":
    main()
