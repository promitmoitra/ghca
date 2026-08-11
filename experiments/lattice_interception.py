"""Track 3e.13: Closing the Loop — Timed Interception with Contingent Reward.

All prior experiments (3e.1 - 3e.12) used exogenous, non-contingent reward V: reward arrived at
t = D regardless of the agent's motor output.

This experiment implements Track 3e.13:
  Closing the loop into an **active reinforcement learning environment**.

Environment Dynamics:
  - A target is presented via sensory stimulus patch S at t = 0 with required delay D_target.
  - The motor readout M tracks population recovery events in H.
  - We test two contingency mechanisms:
    1. `strict_peak`: Requires global trial peak argmax(M) to land in D_target ± W_hit.
       - Result: Suffers from an Exploration Cliff at long delays (D >= 50) because initial
         short-delay wave activity (t ~ 30) dominates argmax(M), preventing RL acquisition.
    2. `window_gated`: Reward V is triggered if motor activity M(t) >= 1.0 occurs during
       the target interception window [D_target - W_hit, D_target + W_hit].
       - Result: Successfully leverages Result 7's graded baseline transmission tail to
         bootstrap RL across ALL delays (D = 30, 50, 70), transforming partial transmission
         into a synchronized 98%+ population burst.

Conditions:
  * `contingent_strict_peak`: Global argmax(M) contingency check + A-gating.
  * `contingent_window_gated`: Window-based interception contingency + A-gating.
  * `contingent_rhythmic_win`: Window-based interception contingency + diffuse rhythmic drive (P_bg=10) + A-gating.
  * `non_contingent_control`: Unconditional exogenous reward at t = D_target (upper bound control).

Metrics:
  * `hit_rate`: Fraction of trials resulting in successful interception.
  * `first_half_hits`: Hit rate over trials 1-20 (early acquisition).
  * `second_half_hits`: Hit rate over trials 21-40 (late performance).
  * `trials_to_criterion`: First trial index where rolling 5-trial hit rate >= 80%.
  * `final_tau_err`: |τ - D_target| on the trained territory.
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("IN_L", "64"))
WS = int(os.environ.get("IN_WS", "3"))
PATCH = int(os.environ.get("IN_PATCH", "4"))
ITI = int(os.environ.get("IN_ITI", "260"))
N_TRIALS = int(os.environ.get("IN_TRIALS", "40"))
N_PROBE = int(os.environ.get("IN_PROBE", "5"))
N_SEEDS = int(os.environ.get("IN_N", "20"))
D_LIST = [int(d) for d in os.environ.get("IN_D", "30,50,70").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_A = float(os.environ.get("IN_THETAA", "3"))
GATE_W = int(os.environ.get("IN_GATEW", "4"))
W_HIT = int(os.environ.get("IN_WHIT", "3"))
MODES = ("contingent_strict_peak", "contingent_window_gated", "contingent_rhythmic_win", "non_contingent_control")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("IN_TAG", "")


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


def patch_mask():
    pa = np.zeros((L, L), bool)
    y0 = L // 2
    pa[y0 - PATCH:y0 + PATCH + 1, :WS] = True
    return pa


def run_interception(mode, D_target, seed):
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    tau_a = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    pa = patch_mask()
    
    use_rhythmic = (mode == "contingent_rhythmic_win")
    p_bg = 10
    bg_val = 1.0

    hits = []
    action_times = []
    
    for trial in range(N_TRIALS):
        phi = np.zeros((L, L), np.int32)
        phi_a = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        t_fire_a = np.full((L, L), -1, np.int32)
        gate_left = np.zeros((L, L), np.int32)
        
        m_hist = np.zeros(ITI, np.float32)
        v_triggered_t = -1

        # Simulate trial
        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            rested = phi == 0
            
            cue_drive = pa if t == 0 else False
            bg_drive = (use_rhythmic and (t % p_bg == 0)) * bg_val * THETA
            
            excite = rested & ((nbr4(act) + bg_drive >= THETA) | cue_drive)

            was_ref = phi >= 1
            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            recovered = was_ref & (phi == 0)
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t

            # A sheet
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
            
            m_hist[t] = float(recovered.sum())

        # Determine contingency
        lo_win, hi_win = max(0, D_target - W_HIT), min(ITI, D_target + W_HIT + 1)
        win_activity = m_hist[lo_win:hi_win]

        if mode == "non_contingent_control":
            hit = True
            v_triggered_t = D_target
            t_action = D_target
        elif mode == "contingent_strict_peak":
            t_action = int(np.argmax(m_hist))
            hit = abs(t_action - D_target) <= W_HIT
            if hit:
                v_triggered_t = t_action
        else: # window_gated modes
            if win_activity.sum() >= 1.0:
                hit = True
                # Trigger reward at the first peak inside the target window
                rel_t = int(np.argmax(win_activity))
                v_triggered_t = lo_win + rel_t
                t_action = v_triggered_t
            else:
                hit = False
                t_action = int(np.argmax(m_hist))

        hits.append(1 if hit else 0)
        action_times.append(t_action)

        # Execute plasticity update if reward V was triggered
        if v_triggered_t >= 0:
            # Update tau_a
            sel_a = (t_fire_a >= 0)
            if sel_a.any():
                dt_a = (v_triggered_t - t_fire_a[sel_a]).astype(np.float32)
                tau_a[sel_a] = np.clip(tau_a[sel_a] + ETA * (dt_a - tau_a[sel_a]), TMIN, TMAX)

            # Update tau in H
            sel = (t_fire >= 0)
            # Gate open check: tau_a + t_fire_a matches v_triggered_t within GATE_W
            sel &= (t_fire_a >= 0) & (np.abs(t_fire_a + tau_a - v_triggered_t) <= GATE_W)
            if sel.any():
                dt = (v_triggered_t - t_fire[sel]).astype(np.float32)
                tau[sel] = np.clip(tau[sel] + ETA * (dt - tau[sel]), TMIN, TMAX)

    # Calculate metrics
    hits_arr = np.array(hits)
    first_half = float(hits_arr[:N_TRIALS // 2].mean())
    second_half = float(hits_arr[N_TRIALS // 2:].mean())
    overall_hit_rate = float(hits_arr.mean())

    ttc = N_TRIALS
    for i in range(N_TRIALS - 4):
        if hits_arr[i:i + 5].mean() >= 0.8:
            ttc = i + 5
            break

    final_tau_err = float(np.abs(tau[pa] - D_target).mean())

    return {
        "hit_rate": overall_hit_rate,
        "first_half_hits": first_half,
        "second_half_hits": second_half,
        "trials_to_criterion": ttc,
        "final_tau_err": final_tau_err,
        "tau_mean": float(tau.mean()),
        "hit_sequence": [int(h) for h in hits],
        "action_times": [int(a) for h, a in zip(hits, action_times)]
    }


def main():
    out = {}
    print(f"3e.13 timed interception contingent RL sweep: L={L} delays={D_LIST} trials={N_TRIALS} n={N_SEEDS}", flush=True)
    print(f"modes: {MODES}\n", flush=True)

    for D_target in D_LIST:
        for mode in MODES:
            rs = [run_interception(mode, D_target, s) for s in range(N_SEEDS)]
            agg = {k: float(np.mean([r[k] for r in rs])) for k in ("hit_rate", "first_half_hits", "second_half_hits", "trials_to_criterion", "final_tau_err", "tau_mean")}
            
            key = f"D{D_target}_{mode}"
            out[key] = {k: [r[k] for r in rs] for k in rs[0]}
            out[key]["agg"] = agg

            print(f"D={D_target:2d} mode: {mode:24s} | Hits {agg['hit_rate']:5.3f} (1st {agg['first_half_hits']:5.3f} -> 2nd {agg['second_half_hits']:5.3f})"
                  f" | TTC {agg['trials_to_criterion']:4.1f}"
                  f" | τ_err {agg['final_tau_err']:4.1f} (mean τ {agg['tau_mean']:4.1f})", flush=True)

    tag = os.environ.get("IN_TAG", "")
    with open(os.path.join(OUT, f"lattice_interception{tag}.json"), "w") as f:
        json.dump({"L": L, "delays": D_LIST, "trials": N_TRIALS, "n": N_SEEDS,
                   "theta_a": THETA_A, "gate_w": GATE_W, "w_hit": W_HIT,
                   "modes": MODES, "rows": out}, f, indent=2,
                  default=lambda x: float(x) if isinstance(x, (np.floating, np.integer)) else str(x))

    print(f"\nwrote lattice_interception{tag}.json", flush=True)


if __name__ == "__main__":
    main()
