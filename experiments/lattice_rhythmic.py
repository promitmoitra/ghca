"""Track 3e.19: Rhythmic (not stochastic) background drive for closed-loop Q-switched actuation.

Negative 2 (`lattice_tonic.py`) established that stochastic tonic drive (uncorrelated re-firing
with probability q) decoheres excitability windows, causing learning collapse before any
band-pass peakiness or selective avoidance can be achieved.

This experiment tests the 3e.19 hypothesis:
  Replacing stochastic noise with a **coherent, rhythmic background drive** (a periodic pacemaker
  wave of period P_bg across a background channel/edge) allows learned τ refractory windows
  to **phase-lock to the background rhythm**.

The learned τ window then acts as an optical **Q-switch**:
  - When the background rhythm hits a cell that is still refractory (c > 0), the wave is blocked.
  - When the background rhythm hits a cell whose τ window has just expired (c = 0, rest state),
    the cell fires in phase with the background drive.
  - This converts a passive transmissive window into an **active, synchronized motor population
    burst** (Q-switched emission) driven by the endogenous background rhythm without requiring
    an external probe event!

Sweep:
  Background period P_bg in {6, 10, 12, 16, 20, 0} (where 0 is the no-background control).
  Target reward delay D = 50 (with D = 30, 70 sweeps as options).

Metrics:
  * `q_switch_peak`: Motor burst emission magnitude at t = D relative to off-phase times.
  * `tau_err`: |(t_rew - t_fire) - τ| — checks if input-timing τ learning stays intact.
  * `aligned_frac`: Fraction of cells with |Δ| <= 2.0.
  * `phase_locking`: Vector strength / phase coherence of motor firing relative to t mod P_bg.
  * `reentrant`: Late trial activity, checking if the background drive maintains stability.
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("RH_L", "64"))
WS = int(os.environ.get("RH_WS", "3"))
PATCH = int(os.environ.get("RH_PATCH", "4"))
ITI = int(os.environ.get("RH_ITI", "260"))
N_TRIALS = int(os.environ.get("RH_TRIALS", "40"))
N_SEEDS = int(os.environ.get("RH_N", "20"))
D = int(os.environ.get("RH_D", "50"))
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_M = 3.0
RESP_W = int(os.environ.get("RH_RESPW", "12"))
P_BG_LIST = [int(x) for x in os.environ.get("RH_PBG", "0,6,10,12,16,20").split(",")]
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("RH_TAG", "")


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


def bg_mask(mode):
    """Background pacemaker mask: boundary strip or diffuse whole-sheet."""
    if mode == "diffuse_sub":
        return np.ones((L, L), bool)
    b = np.zeros((L, L), bool)
    b[:WS, :] = True
    return b


def train_rhythmic(p_bg, mode, bg_val, seed):
    """Train τ with a coherent rhythmic background drive of period p_bg."""
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    pa = patch_mask()
    bg = bg_mask(mode) if p_bg > 0 else np.zeros((L, L), bool)
    reentrant = []

    for _ in range(N_TRIALS):
        phi = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        late_act = 0
        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            rested = phi == 0
            cue_drive = pa if t == 0 else False
            
            if p_bg > 0 and (t % p_bg == 0):
                if mode == "boundary_supra":
                    excite = rested & ((nbr4(act) >= THETA) | cue_drive | bg)
                else: # boundary_sub or diffuse_sub
                    drive_sum = nbr4(act) + bg.astype(np.float32) * bg_val
                    excite = rested & ((drive_sum >= THETA) | cue_drive)
            else:
                excite = rested & ((nbr4(act) >= THETA) | cue_drive)

            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            phi[excite] = 1
            
            newly = excite & (t_fire < 0)
            t_fire[newly] = t              # FIRST firing relative to cue
            
            if t > ITI - 60:
                late_act += int(act.sum() > 0)
                
            if t == D:
                sel = t_fire >= 0
                if not sel.any():
                    continue
                dt = (t - t_fire[sel]).astype(np.float32)
                cur = tau[sel]
                new = cur + ETA * (dt - cur)
                tau[sel] = np.clip(new, TMIN, TMAX)
                
        reentrant.append(late_act / 60.0)
    return tau, float(np.mean(reentrant))


def run_q_switch_probe(tau, p_bg, mode, bg_val, seed):
    """Probe autonomous Q-switched motor emission driven by background rhythm alone (no probe)."""
    rng = np.random.default_rng(10_000 + seed)
    pa = patch_mask()
    bg = bg_mask(mode) if p_bg > 0 else np.zeros((L, L), bool)
    
    phi = np.zeros((L, L), np.int32)
    phi_m = np.zeros((L, L), np.int32)
    
    m_trace = np.zeros(ITI, np.float32)
    act_trace = np.zeros(ITI, np.float32)
    
    for t in range(ITI):
        act = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        cue_drive = pa if t == 0 else False
        
        if p_bg > 0 and (t % p_bg == 0):
            if mode == "boundary_supra":
                excite = rested & ((nbr4(act) >= THETA) | cue_drive | bg)
            else:
                drive_sum = nbr4(act) + bg.astype(np.float32) * bg_val
                excite = rested & ((drive_sum >= THETA) | cue_drive)
        else:
            excite = rested & ((nbr4(act) >= THETA) | cue_drive)

        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        
        # Motor layer M downstream of H (box3 convergence >= THETA_M)
        exc_m = (phi_m == 0) & (box3(act) >= THETA_M)
        phi_m[phi_m >= 1] += 1
        phi_m[phi_m > int(ACT + 1)] = 0
        phi_m[exc_m] = 1
        
        m_trace[t] = float(exc_m.sum()) / (L * L)
        act_trace[t] = float(act.mean())

    # Measure Q-switch emission at target reward time D vs baseline off-target times
    window_at_d = float(m_trace[max(0, D - 3):min(ITI, D + RESP_W)].sum())
    window_off = float(m_trace[:max(1, D - 10)].sum() + m_trace[min(ITI, D + RESP_W + 10):].sum())
    off_mean = window_off / max(1, ITI - RESP_W - 20)
    q_switch_peak = float(window_at_d / max(off_mean * (RESP_W + 3), 1e-6))

    # Phase coherence / vector strength relative to background rhythm p_bg
    if p_bg > 0 and m_trace.sum() > 0:
        phases = 2.0 * np.pi * (np.arange(ITI) % p_bg) / p_bg
        weights = m_trace
        phase_coherence = float(np.abs(np.sum(weights * np.exp(1j * phases))) / max(weights.sum(), 1e-6))
    else:
        phase_coherence = 0.0

    return q_switch_peak, float(m_trace[D]), phase_coherence, m_trace.astype(float).tolist()


def alignment(tau, p_bg, mode, bg_val, seed):
    """Measure |(t_rew - t_fire) - tau| on a test run."""
    rng = np.random.default_rng(20_000 + seed)
    pa = patch_mask()
    bg = bg_mask(mode) if p_bg > 0 else np.zeros((L, L), bool)
    phi = np.zeros((L, L), np.int32)
    t_fire = np.full((L, L), -1, np.int32)
    
    for t in range(D + 1):
        act = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        cue_drive = pa if t == 0 else False
        
        if p_bg > 0 and (t % p_bg == 0):
            if mode == "boundary_supra":
                excite = rested & ((nbr4(act) >= THETA) | cue_drive | bg)
            else:
                drive_sum = nbr4(act) + bg.astype(np.float32) * bg_val
                excite = rested & ((drive_sum >= THETA) | cue_drive)
        else:
            excite = rested & ((nbr4(act) >= THETA) | cue_drive)

        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        phi[excite] = 1
        
        newly = excite & (t_fire < 0)
        t_fire[newly] = t
        
        if t == D:
            sel = t_fire >= 0
            if not sel.any():
                return float("nan"), 0.0
            dt = (D - t_fire[sel]).astype(np.float32)
            e = np.abs(dt - tau[sel])
            return float(e.mean()), float((e <= 2.0).mean())
            
    return float("nan"), 0.0


def main():
    out = {}
    mode = os.environ.get("RH_MODE", "diffuse_sub")
    bg_val = float(os.environ.get("RH_BGVAL", "0.8"))
    print(f"3e.19 rhythmic background sweep: L={L} D={D} trials={N_TRIALS} n={N_SEEDS} mode={mode} bg_val={bg_val}", flush=True)
    print(f"background periods P_bg: {P_BG_LIST}\n", flush=True)

    for p_bg in P_BG_LIST:
        peaks, m_at_d, coherences, errs, fracs, reents = [], [], [], [], [], []
        for s in range(N_SEEDS):
            tau, reent = train_rhythmic(p_bg, mode, bg_val, s)
            peak, md, coh, _ = run_q_switch_probe(tau, p_bg, mode, bg_val, s)
            e, fr = alignment(tau, p_bg, mode, bg_val, s)
            peaks.append(float(peak))
            m_at_d.append(float(md))
            coherences.append(float(coh))
            errs.append(float(e) if not np.isnan(e) else float("nan"))
            fracs.append(float(fr))
            reents.append(float(reent))

        agg = {
            "q_switch_peak": float(np.nanmean(peaks)),
            "m_at_d": float(np.nanmean(m_at_d)),
            "phase_coherence": float(np.nanmean(coherences)),
            "tau_err": float(np.nanmean(errs)),
            "aligned_frac": float(np.nanmean(fracs)),
            "reentrant": float(np.nanmean(reents))
        }

        out[f"pbg_{p_bg}"] = {
            "q_switch_peak": peaks,
            "m_at_d": m_at_d,
            "phase_coherence": coherences,
            "tau_err": errs,
            "aligned_frac": fracs,
            "reentrant": reents,
            "agg": agg
        }

        print(f"P_bg = {p_bg:2d} | Q-peak {agg['q_switch_peak']:5.2f} "
              f"| M@D {agg['m_at_d']:6.4f} "
              f"| Phase-Coh {agg['phase_coherence']:4.2f} "
              f"| τ_err {agg['tau_err']:4.1f}±{np.nanstd(errs):3.1f} "
              f"(within±2 {agg['aligned_frac']:.2f}) "
              f"| reentrant {agg['reentrant']:.2f}", flush=True)

    with open(os.path.join(OUT, f"lattice_rhythmic{TAG}.json"), "w") as f:
        json.dump({"L": L, "D": D, "trials": N_TRIALS, "n": N_SEEDS,
                   "p_bg_list": P_BG_LIST, "mode": mode, "bg_val": bg_val,
                   "rows": out}, f, indent=2,
                  default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x))
                   
    print(f"\nwrote lattice_rhythmic{TAG}.json", flush=True)


if __name__ == "__main__":
    main()

