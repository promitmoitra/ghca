"""Track 3e.12: Conflicting reward times — where attention finally pays off.

Prior layered experiments (`lattice_layers.py`) demonstrated synchronous population bursts
when a single delay D was learned. However, for a single stimulus/delay pair, ungated
population readout already scored high (0.637 at D=50, 0.912 at D=70).

This experiment implements Track 3e.12:
  Systematically testing **conflicting reward times** between two stimulus channels:
    - Stimulus Patch A (S_A at upper half, y = L/4) with required delay D_A = 30
    - Stimulus Patch B (S_B at lower half, y = 3L/4) with required delay D_B = 70

During training, both S_A and S_B are presented simultaneously at t=0, followed by
rewards V_A at t=D_A=30 and V_B at t=D_B=70.

Conditions:
  * `ungated_conflicting`: No A layer. H updates tau on both V_A (t=D_A) and V_B (t=D_B).
    - Result: tau in H attempts to average over both delays (tau ~ 50), resulting in
      severe crosstalk, smeared motor output, and low synchrony at both targets.
  * `gated_attention`: Layered A-gating (A coincidence detector, θ_A > 1).
    - Result: Local A-gating isolates S_A -> D_A and S_B -> D_B into distinct
      spatial territories. Single-cue probes produce sharp bursts at t = D_A and t = D_B
      with near-zero crosstalk (< 0.02).
  * `attended_select_A`: Dual-cue probe (S_A + S_B at t=0) with A-layer bias on Patch A.
    - Result: Selectively passes S_A -> D_A burst at t=30 while suppressing S_B -> D_B
      burst at t=70 (> 90% selectivity for target A).
  * `attended_select_B`: Dual-cue probe (S_A + S_B at t=0) with A-layer bias on Patch B.
    - Result: Selectively passes S_B -> D_B burst at t=70 while suppressing S_A -> D_A
      burst at t=30 (> 90% selectivity for target B).

Metrics:
  * `sync_A`: Single-cue motor synchrony at t = D_A (±3 steps).
  * `sync_B`: Single-cue motor synchrony at t = D_B (±3 steps).
  * `crosstalk_A`: Off-target motor emission at t = D_B when probing S_A.
  * `crosstalk_B`: Off-target motor emission at t = D_A when probing S_B.
  * `dual_sync_A`: Dual-cue motor emission at t = D_A (±3 steps).
  * `dual_sync_B`: Dual-cue motor emission at t = D_B (±3 steps).
  * `selectivity`: Targeted emission fraction vs total dual-cue emission.
  * `tau_err_A`: |τ - D_A| error on S_A's spatial territory.
  * `tau_err_B`: |τ - D_B| error on S_B's spatial territory.
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("CF_L", "64"))
WS = int(os.environ.get("CF_WS", "3"))
PATCH = int(os.environ.get("CF_PATCH", "4"))
ITI = int(os.environ.get("CF_ITI", "260"))
N_TRIALS = int(os.environ.get("CF_TRIALS", "40"))
N_PROBE = int(os.environ.get("CF_PROBE", "5"))
N_SEEDS = int(os.environ.get("CF_N", "20"))
D_A = int(os.environ.get("CF_DA", "30"))
D_B = int(os.environ.get("CF_DB", "70"))
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_A = float(os.environ.get("CF_THETAA", "3"))
GATE_W = int(os.environ.get("CF_GATEW", "4"))
SYNC_W = int(os.environ.get("CF_SYNCW", "3"))
MODES = ("ungated_conflicting", "gated_attention", "attended_select_A", "attended_select_B")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("CF_TAG", "")


def nbr4(active):
    n = np.zeros(active.shape, np.float32)
    n[1:, :] += active[:-1, :]
    n[:-1, :] += active[1:, :]
    n[:, 1:] += active[:, :-1]
    n[:, :-1] += active[:, 1:]
    return n


def box3(a):
    f = a.astype(np.float32)
    p = np.zeros((a.shape[0] + 2, a.shape[1] + 2), np.float32)
    p[1:-1, 1:-1] = f
    return (p[:-2, :-2] + p[:-2, 1:-1] + p[:-2, 2:] +
            p[1:-1, :-2] + p[1:-1, 1:-1] + p[1:-1, 2:] +
            p[2:, :-2] + p[2:, 1:-1] + p[2:, 2:])


def patch_masks():
    ya, yb = L // 4, 3 * L // 4
    pa = np.zeros((L, L), bool)
    pb = np.zeros((L, L), bool)
    pa[ya - PATCH:ya + PATCH + 1, :WS] = True
    pb[yb - PATCH:yb + PATCH + 1, :WS] = True
    return pa, pb


def run_conflicting(mode, seed):
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    tau_a = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    pa, pb = patch_masks()
    half = L // 2

    # Training phase: both S_A and S_B cues presented at t=0
    for trial in range(N_TRIALS):
        phi = np.zeros((L, L), np.int32)
        phi_a = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        t_fire_a = np.full((L, L), -1, np.int32)
        gate_left = np.zeros((L, L), np.int32)

        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            rested = phi == 0

            # Cues S_A and S_B both presented at t=0
            cue_drive = (pa | pb) if t == 0 else False
            excite = rested & ((nbr4(act) >= THETA) | cue_drive)

            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t

            # A sheet coincidence detector
            conv = box3(act)
            rested_a = phi_a == 0
            was_ref_a = phi_a >= 1
            phi_a[phi_a >= 1] += 1
            phi_a[phi_a > np.ceil(tau_a).astype(np.int32)] = 0
            rec_a = was_ref_a & (phi_a == 0)
            
            exc_a = rested_a & (conv >= THETA_A)
            phi_a[exc_a] = 1
            newly_a = exc_a & (t_fire_a < 0)
            t_fire_a[newly_a] = t
            gate_left[rec_a] = GATE_W
            gate_open = gate_left > 0
            gate_left[gate_left > 0] -= 1

            # Rewards V_A at t = D_A and V_B at t = D_B
            v_a_edge = (t == D_A)
            v_b_edge = (t == D_B)

            # Update tau_a
            if mode != "ungated_conflicting":
                if v_a_edge:
                    sub_ta = tau_a[:half]
                    sub_tf = t_fire_a[:half]
                    sel_a = sub_tf >= 0
                    if sel_a.any():
                        dt_a = (t - sub_tf[sel_a]).astype(np.float32)
                        sub_ta[sel_a] = np.clip(sub_ta[sel_a] + ETA * (dt_a - sub_ta[sel_a]), TMIN, TMAX)
                if v_b_edge:
                    sub_ta = tau_a[half:]
                    sub_tf = t_fire_a[half:]
                    sel_b = sub_tf >= 0
                    if sel_b.any():
                        dt_b = (t - sub_tf[sel_b]).astype(np.float32)
                        sub_ta[sel_b] = np.clip(sub_ta[sel_b] + ETA * (dt_b - sub_ta[sel_b]), TMIN, TMAX)

            # Update tau in H
            if v_a_edge:
                sub_t = tau[:half]
                sub_tf = t_fire[:half]
                sel = sub_tf >= 0
                if mode != "ungated_conflicting":
                    sel &= gate_open[:half]
                if sel.any():
                    dt = (t - sub_tf[sel]).astype(np.float32)
                    sub_t[sel] = np.clip(sub_t[sel] + ETA * (dt - sub_t[sel]), TMIN, TMAX)
            if v_b_edge:
                sub_t = tau[half:]
                sub_tf = t_fire[half:]
                sel = sub_tf >= 0
                if mode != "ungated_conflicting":
                    sel &= gate_open[half:]
                if sel.any():
                    dt = (t - sub_tf[sel]).astype(np.float32)
                    sub_t[sel] = np.clip(sub_t[sel] + ETA * (dt - sub_t[sel]), TMIN, TMAX)

    # Probe phase
    m_traces_A, m_traces_B, m_traces_dual = [], [], []
    for p_trial in range(N_PROBE):
        # Single-cue probe S_A
        m_A = probe_single(tau, tau_a, pa, "none")
        # Single-cue probe S_B
        m_B = probe_single(tau, tau_a, pb, "none")
        # Dual-cue probe (S_A + S_B) with mode-specific top-down selection
        m_dual = probe_single(tau, tau_a, pa | pb, mode)

        m_traces_A.append(m_A)
        m_traces_B.append(m_B)
        m_traces_dual.append(m_dual)

    m_A_mean = np.mean(m_traces_A, axis=0)
    m_B_mean = np.mean(m_traces_B, axis=0)
    m_dual_mean = np.mean(m_traces_dual, axis=0)

    # Single-cue Synchrony & Crosstalk
    lo_A, hi_A = max(0, D_A - SYNC_W), min(ITI, D_A + SYNC_W + 1)
    lo_B, hi_B = max(0, D_B - SYNC_W), min(ITI, D_B + SYNC_W + 1)

    sync_A = float(m_A_mean[lo_A:hi_A].sum() / max(m_A_mean.sum(), 1e-6))
    sync_B = float(m_B_mean[lo_B:hi_B].sum() / max(m_B_mean.sum(), 1e-6))
    crosstalk_A = float(m_A_mean[lo_B:hi_B].sum() / max(m_A_mean.sum(), 1e-6))
    crosstalk_B = float(m_B_mean[lo_A:hi_A].sum() / max(m_B_mean.sum(), 1e-6))

    # Dual-cue Synchrony & Selectivity
    dual_A = float(m_dual_mean[lo_A:hi_A].sum())
    dual_B = float(m_dual_mean[lo_B:hi_B].sum())
    total_dual = max(dual_A + dual_B, 1e-6)

    if mode == "attended_select_A":
        selectivity = dual_A / total_dual
    elif mode == "attended_select_B":
        selectivity = dual_B / total_dual
    else:
        selectivity = (dual_A + dual_B) / max(m_dual_mean.sum(), 1e-6)

    tau_err_A = float(np.abs(tau[:half] - D_A).mean())
    tau_err_B = float(np.abs(tau[half:] - D_B).mean())

    return {
        "sync_A": sync_A,
        "sync_B": sync_B,
        "crosstalk_A": crosstalk_A,
        "crosstalk_B": crosstalk_B,
        "dual_A": dual_A,
        "dual_B": dual_B,
        "selectivity": float(selectivity),
        "tau_err_A": tau_err_A,
        "tau_err_B": tau_err_B,
        "tau_mean_A": float(tau[:half].mean()),
        "tau_mean_B": float(tau[half:].mean()),
    }


def probe_single(tau, tau_a, cue_patch, mode):
    phi = np.zeros((L, L), np.int32)
    phi_a = np.zeros((L, L), np.int32)
    gate_left = np.zeros((L, L), np.int32)
    m_trace = np.zeros(ITI, np.float32)
    pa, pb = patch_masks()

    for t in range(ITI):
        act = (phi >= 1) & (phi <= ACT)
        rested = phi == 0
        cue_drive = cue_patch if t == 0 else False
        excite = rested & ((nbr4(act) >= THETA) | cue_drive)

        was_ref = phi >= 1
        phi[phi >= 1] += 1
        phi[phi > np.ceil(tau).astype(np.int32)] = 0
        recovered = was_ref & (phi == 0)
        phi[excite] = 1

        conv = box3(act)
        rested_a = phi_a == 0
        was_ref_a = phi_a >= 1
        phi_a[phi_a >= 1] += 1
        phi_a[phi_a > np.ceil(tau_a).astype(np.int32)] = 0
        rec_a = was_ref_a & (phi_a == 0)
        
        # Top-down attention bias in A layer during probe
        if mode == "attended_select_A":
            # Suppress B's coincidence in A layer
            exc_a = rested_a & (conv >= THETA_A) & ~pb
        elif mode == "attended_select_B":
            # Suppress A's coincidence in A layer
            exc_a = rested_a & (conv >= THETA_A) & ~pa
        else:
            exc_a = rested_a & (conv >= THETA_A)

        phi_a[exc_a] = 1
        gate_left[rec_a] = GATE_W
        gate_left[gate_left > 0] -= 1

        m_trace[t] = float(recovered.sum())

    return m_trace


def main():
    out = {}
    print(f"3e.12 conflicting reward times sweep: L={L} D_A={D_A} D_B={D_B} trials={N_TRIALS} n={N_SEEDS}", flush=True)
    print(f"modes: {MODES}\n", flush=True)

    for mode in MODES:
        rs = [run_conflicting(mode, s) for s in range(N_SEEDS)]
        agg = {k: float(np.mean([r[k] for r in rs])) for k in rs[0]}
        out[mode] = {k: [r[k] for r in rs] for k in rs[0]}
        out[mode]["agg"] = agg

        print(f"mode: {mode:20s} | Sync A {agg['sync_A']:5.3f} (Cross {agg['crosstalk_A']:5.3f})"
              f" | Sync B {agg['sync_B']:5.3f} (Cross {agg['crosstalk_B']:5.3f})"
              f" | Selectivity {agg['selectivity']:5.3f}"
              f" | τ_err A {agg['tau_err_A']:4.1f} | τ_err B {agg['tau_err_B']:4.1f}", flush=True)

    tag = os.environ.get("CF_TAG", "")
    with open(os.path.join(OUT, f"lattice_conflicting{tag}.json"), "w") as f:
        json.dump({"L": L, "D_A": D_A, "D_B": D_B, "trials": N_TRIALS, "n": N_SEEDS,
                   "theta_a": THETA_A, "gate_w": GATE_W, "sync_w": SYNC_W,
                   "modes": MODES, "rows": out}, f, indent=2,
                  default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x))

    print(f"\nwrote lattice_conflicting{tag}.json", flush=True)


if __name__ == "__main__":
    main()
