"""Layers instead of edges: sheets of the same cells, stacked and modulating each other.

The four-edge work (`lattice_reward_edges.py`, `lattice_attention_value.py`) got content and
placement to be learned rather than designed, but left three things standing:

  1. **the anatomy was hand-built** — a 1-D delay line laid out by hand, whose speed was the
     designed constant the value chain then had to remove;
  2. **attention could only select a column**, because a chain has one cell per column;
  3. **the value signal could only travel one direction**, so delays downstream of the write
     site were never trained at all — the mechanism's reach was visible in its own parameters.

All three are artifacts of using a 1-D edge for a job the cortex does with a sheet. So: replace
both chains with **2-D layers of the same Greenberg-Hastings cells**, coupled vertically.

    V  value / neuromodulator sheet — carries "reward now", diffuse or propagating
    A  attention sheet   — fires on convergent H activity, plastic τ_A; its RECOVERY gates H
    H  hidden sheet      — plastic τ, receives the stimulus patch, waves propagate
    M  motor readout     — H's recovery events, counted per step

Two design choices are forced by results already in this repo rather than chosen:

  * **Vertical coupling is modulatory, not excitatory.** Negative 1 (the fastest pacemaker
    entrains the whole medium) and the two failed cross-frequency designs in 3e.2b both say
    that coupling rhythmic populations excitatorily produces rigid entrainment; only
    *excitability/plasticity modulation* worked there. `layered-exc` is the control that
    reproduces the failure on purpose.
  * **The attention sheet is a coincidence detector**, driven by the count of active H cells
    in a 3x3 receptive field with θ_A > 1. Negative 1 found θ=2 blocks planar-wave propagation
    in a 4-neighbourhood — a wall there, the useful property here: A cannot sustain waves of
    its own, so it lights up only where H delivers convergent drive, and never floods.

**The learning hierarchy has no bootstrapping problem.** V teaches A *ungated* (A learns when
reward is due relative to its own firing), and A's recovery then gates V's teaching of H. So
the gate does not have to be right before it can learn to be right.

**The motor readout is the point of the whole arc.** If every H cell's τ equals its own
stimulus-to-reward interval, then cells that fired at different times — a travelling wave,
staggered across the sheet — all become excitable again at the *same* moment, the predicted
reward time. A spatially distributed, temporally staggered response converges into a
synchronous population event. That is measurable (`sync`: the fraction of recovery events
landing within ±3 steps of the true reward time) and it is the first thing in this arc that
the medium *does* with what it has learned, rather than merely representing.

Note on what is still not emergent. The remaining constants (θ_A, the receptive field size,
the gate window) are **cell properties applied homogeneously**, not a topology drawn by hand.
That is a real improvement on a delay line, and it is all it is. Reward is still exogenous and
uncontingent — nothing the medium does changes whether reward arrives, so the value signal
remains informational rather than appetitive, and this is still not reinforcement learning.

Conditions:
  layered        — modulatory A gate; V diffuse (a neuromodulator arriving everywhere at once)
  layered-vwave  — same, but V propagates as a wave from the reward edge. Does value need to
                   travel, or is diffuse enough?
  layered-exc    — A projects *excitatorily* into H instead of gating it: the entrainment
                   control, predicted to flood
  layered-nogate — no A layer; every H cell updates on reward. The ungated baseline.
  unpaired       — reward time redrawn each trial: same signal count, no contingency
  untrained      — no plasticity at all, to measure M's baseline synchrony
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("LY_L", "64"))
WS = int(os.environ.get("LY_WS", "3"))
PATCH = int(os.environ.get("LY_PATCH", "4"))
ITI = int(os.environ.get("LY_ITI", "260"))
N_TRIALS = int(os.environ.get("LY_TRIALS", "40"))
N_PROBE = int(os.environ.get("LY_PROBE", "3"))
N_SEEDS = int(os.environ.get("LY_N", "3"))
D_LIST = [int(d) for d in os.environ.get("LY_D", "30,50,70").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
THETA_A = float(os.environ.get("LY_THETAA", "3"))   # convergence needed to fire an A cell
GATE_W = int(os.environ.get("LY_GATEW", "4"))       # steps the gate stays open after recovery
SYNC_W = int(os.environ.get("LY_SYNCW", "3"))       # ± window for counting a recovery as synchronous
MODES = ("layered", "layered-vwave", "layered-exc", "layered-nogate",
         "distract-g2", "distract-g4", "distract-g8", "distract-g16",
         "distract-nogate", "unpaired", "untrained")
JIT_MAX = int(os.environ.get("LY_JIT", "40"))       # distractor onset jitter
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("LY_TAG", "")


def nbr4(active):
    """4-neighbour count on a sheet with walls on all four sides."""
    n = np.zeros(active.shape, np.float32)
    n[1:, :] += active[:-1, :]
    n[:-1, :] += active[1:, :]
    n[:, 1:] += active[:, :-1]
    n[:, :-1] += active[:, 1:]
    return n


def box3(a):
    """Count over a 3x3 receptive field — the vertical projection H -> A."""
    f = a.astype(np.float32)
    p = np.zeros((a.shape[0] + 2, a.shape[1] + 2), np.float32)
    p[1:-1, 1:-1] = f
    return (p[:-2, :-2] + p[:-2, 1:-1] + p[:-2, 2:] +
            p[1:-1, :-2] + p[1:-1, 1:-1] + p[1:-1, 2:] +
            p[2:, :-2] + p[2:, 1:-1] + p[2:, 2:])


def run(mode, D, seed):
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)      # H timescales
    tau_a = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)    # A timescales
    learn = mode != "untrained"
    # unpaired is ungated so it is a clean contingency control against layered-nogate
    gated = not (mode.endswith("nogate") or mode == "unpaired")
    distract = mode.startswith("distract")
    # the gate window trades coverage against selectivity; "-gN" sweeps it
    gw = int(mode.split("-g")[1]) if "-g" in mode else GATE_W
    # one patch whose onset is fixed (so its stimulus-to-reward interval is predictable) and,
    # in the distractor conditions, a second patch whose onset jitters (so its interval is not).
    # Without a distractor there is nothing for attention to exclude and gating can only subtract.
    pa = np.zeros((L, L), bool)
    pb = np.zeros((L, L), bool)
    ya, yb = L // 4, 3 * L // 4
    pa[ya - PATCH:ya + PATCH + 1, :WS] = True
    pb[yb - PATCH:yb + PATCH + 1, :WS] = True
    if not distract:
        pa[:] = False
        y0 = L // 2
        pa[y0 - PATCH:y0 + PATCH + 1, :WS] = True
        pb[:] = False

    aligns, afracs, fw, syncs, peaks, aflood = [], [], [], [], [], []
    sync_a, sync_b, wr_a, wr_b = [], [], [], []

    for trial in range(N_TRIALS + N_PROBE):
        probe = trial >= N_TRIALS
        d_rew = int(rng.integers(TMIN + 2, ITI // 2)) if (mode == "unpaired" and not probe) else D
        jit = int(rng.integers(0, JIT_MAX + 1)) if distract else 0

        phi = np.zeros((L, L), np.int32)
        phi_a = np.zeros((L, L), np.int32)
        v_phi = np.zeros((L, L), np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        t_fire_a = np.full((L, L), -1, np.int32)
        gate_left = np.zeros((L, L), np.int32)      # steps of gate still open
        n_upd = np.zeros((L, L), np.int32)
        rec_hist = np.zeros(ITI, np.float64)        # the motor readout: recoveries per step
        # split by half so the distractor conditions can show *where* learning was suppressed:
        # the predictable patch sits in the upper half, the jittered one in the lower half
        hist_pred = np.zeros(ITI, np.float64)
        hist_dist = np.zeros(ITI, np.float64)
        half = L // 2
        a_act_sum, v_prev = 0.0, np.zeros((L, L), bool)

        for t in range(ITI):
            act = (phi >= 1) & (phi <= ACT)
            act_a = (phi_a >= 1) & (phi_a <= ACT)
            a_act_sum += float(act_a.mean())

            # --- H sheet: stimulus patch, lateral waves, optional excitatory input from A ---
            rested = phi == 0
            drive = np.zeros((L, L), bool)
            if t == 0:
                drive |= pa
            if distract and t == jit:
                drive |= pb
            exc = (nbr4(act) >= THETA) | drive
            if mode == "layered-exc":
                # the entrainment control: A delivers the SAME signal at the SAME time as the
                # modulatory version, but as excitation into H rather than as a plasticity gate
                exc = exc | (gate_left > 0)
            excite = rested & exc
            was_ref = phi >= 1
            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            recovered = was_ref & (phi == 0)         # recovery edge == the motor event
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t
            rec = recovered & (t_fire >= 0)
            rec_hist[t] = float(rec.sum())
            hist_pred[t] = float(rec[:half].sum())
            hist_dist[t] = float(rec[half:].sum())

            # --- A sheet: coincidence detector on H, cannot sustain its own waves ----------
            conv = box3(act)
            rested_a = phi_a == 0
            was_ref_a = phi_a >= 1
            phi_a[phi_a >= 1] += 1
            phi_a[phi_a > np.ceil(tau_a).astype(np.int32)] = 0
            rec_a = was_ref_a & (phi_a == 0)         # A's recovery == the gate opening
            exc_a = rested_a & (conv >= THETA_A)
            phi_a[exc_a] = 1
            newly_a = exc_a & (t_fire_a < 0)
            t_fire_a[newly_a] = t
            gate_left[rec_a] = gw
            gate_open = gate_left > 0
            gate_left[gate_left > 0] -= 1

            # --- V sheet: "reward now", diffuse or propagating ----------------------------
            if mode == "layered-vwave":
                v_in = np.zeros((L, L), bool)
                v_in[:, -1] = (t == d_rew)
                v_in |= nbr4((v_phi >= 1) & (v_phi <= ACT)) >= 1
                rested_v = v_phi == 0
                v_phi[v_phi >= 1] += 1
                v_phi[v_phi > int(ACT + 1)] = 0
                v_phi[rested_v & v_in] = 1
                v_here = (v_phi >= 1) & (v_phi <= ACT)
            else:
                v_here = np.full((L, L), t == d_rew)   # diffuse: everywhere at once
            v_edge = v_here & ~v_prev
            v_prev = v_here

            if v_edge.any():
                # V teaches A ungated: A learns when reward is due relative to its own firing
                selA = v_edge & (t_fire_a >= 0)
                if learn and not probe and selA.any():
                    dta = (t - t_fire_a[selA]).astype(np.float32)
                    tau_a[selA] = np.clip(tau_a[selA] + ETA * (dta - tau_a[selA]), TMIN, TMAX)
                # A's gate then decides which H cells V may teach
                sel = v_edge & (t_fire >= 0)
                if gated:
                    sel &= gate_open
                if sel.any():
                    dt = (t - t_fire[sel]).astype(np.float32)
                    if probe:
                        e = np.abs(dt - tau[sel])
                        aligns.append(float(e.mean()))
                        afracs.append(float((e <= 2.0).mean()))
                    elif learn:
                        tau[sel] = np.clip(tau[sel] + ETA * (dt - tau[sel]), TMIN, TMAX)
                        n_upd[sel] += 1

        lo, hi = max(0, D - SYNC_W), min(ITI, D + SYNC_W + 1)
        frac_in = lambda h: (float(h[lo:hi].sum() / h.sum()) if h.sum() > 0 else np.nan)
        if probe:
            syncs.append(frac_in(rec_hist))
            sync_a.append(frac_in(hist_pred))
            sync_b.append(frac_in(hist_dist))
            peaks.append(float(np.argmax(rec_hist)) if rec_hist.sum() > 0 else np.nan)
        else:
            fw.append(float((n_upd > 0).mean()))
            wr_a.append(float((n_upd[:half] > 0).mean()))
            wr_b.append(float((n_upd[half:] > 0).mean()))
        aflood.append(a_act_sum / ITI)

    return {
        "align": float(np.nanmean(aligns)) if aligns else float("nan"),
        "aligned_frac": float(np.nanmean(afracs)) if afracs else float("nan"),
        "frac_written": float(np.nanmean(fw)) if fw else float("nan"),
        "sync": float(np.nanmean(syncs)) if syncs else float("nan"),
        "peak": float(np.nanmean(peaks)) if peaks else float("nan"),
        "a_active": float(np.nanmean(aflood)),
        "sync_pred": float(np.nanmean(sync_a)) if sync_a else float("nan"),
        "sync_distract": float(np.nanmean(sync_b)) if sync_b else float("nan"),
        "written_pred": float(np.nanmean(wr_a)) if wr_a else float("nan"),
        "written_distract": float(np.nanmean(wr_b)) if wr_b else float("nan"),
        "tau_mean": float(tau.mean()),
        "tau_a_mean": float(tau_a.mean()),
    }


def main():
    out = {}
    for D in D_LIST:
        for mode in MODES:
            rs = [run(mode, D, s) for s in range(N_SEEDS)]
            agg = {k: [r[k] for r in rs] for k in rs[0]}
            out[f"D{D}_{mode}"] = agg
            m = lambda k: float(np.nanmean(agg[k]))
            print(f"D={D:2d} {mode:15s}: probe |Δ| {m('align'):6.2f}"
                  f" within ±2 {m('aligned_frac'):5.2f}"
                  f" | H written {m('frac_written'):5.3f}"
                  f" | A active {m('a_active'):5.3f}"
                  f" | motor peak {m('peak'):6.1f} (true {D})"
                  f" sync {m('sync'):5.3f}"
                  f" | τ {m('tau_mean'):5.1f} τ_A {m('tau_a_mean'):5.1f}", flush=True)
        print(flush=True)

    print("motor synchrony — does the learned τ field time a population burst to reward?",
          flush=True)
    for mode in MODES:
        s = [float(np.nanmean(out[f"D{D}_{mode}"]["sync"])) for D in D_LIST]
        pk = [float(np.nanmean(out[f"D{D}_{mode}"]["peak"])) for D in D_LIST]
        err = [abs(p - D) for p, D in zip(pk, D_LIST)]
        print(f"  {mode:15s}: sync " + " / ".join(f"{v:5.3f}" for v in s)
              + "   |peak − D| " + " / ".join(f"{v:5.1f}" for v in err), flush=True)

    print("\nwith a distractor — does gating suppress the unpredictable half?", flush=True)
    for mode in [m for m in MODES if m.startswith("distract")]:
        for D in D_LIST:
            r = out[f"D{D}_{mode}"]
            g = lambda k: float(np.nanmean(r[k]))
            print(f"  {mode:16s} D={D:2d}: written pred {g('written_pred'):5.3f}"
                  f" vs distractor {g('written_distract'):5.3f}"
                  f" | sync pred {g('sync_pred'):5.3f}"
                  f" vs distractor {g('sync_distract'):5.3f}"
                  f" | overall sync {g('sync'):5.3f}", flush=True)

    with open(os.path.join(OUT, f"lattice_layers{TAG}.json"), "w") as f:
        json.dump({"L": L, "trials": N_TRIALS, "n": N_SEEDS, "delays": D_LIST,
                   "theta_a": THETA_A, "gate_w": GATE_W, "sync_w": SYNC_W,
                   "tau_range": [TMIN, TMAX], "rows": out}, f, indent=2)
    print(f"wrote lattice_layers{TAG}.json", flush=True)


if __name__ == "__main__":
    main()
