"""The fourth edge: a value signal that teaches the attention chain where to gate.

`lattice_reward_edges.py` closed three quarters of the four-edge proposal. Its gap was
specific and measurable: the attention chain's speed K — one delay-line cell per medium
column, so the gate reaches column x at time Kx — is a **designer's constant**. At K=1 the
gate meets reward exactly where the stimulus wave and the reward arrive together, the
interval there is ~0, τ pins to the floor, and the write carries no interval information at
all. Hand-setting K=3 fixes it. That hand-setting is the last designed thing in the loop.

So give the attention chain its own value input, on the remaining edge:

      ┌──── attention chain: per-column delays d[x], travels → ──────────────┐
   stimulus                     plastic medium                        (reward at t=D)
      └──── value chain: travels ←, seeded where the write happened ────────┘

**Why the value signal has to be an edge, not a scalar.** The gate's arrival time at column
x is the *cumulative* delay of every chain cell upstream of it. A write that lands in a bad
place is therefore not the fault of the cell that fired — it is the fault of the delays
behind it. A global scalar cannot express that. A pulse launched backward from the write
site and updating each delay it passes assigns credit exactly where it belongs, and it is
built from the same cells as everything else. The geometry does the credit assignment.

**What the value signal is.** When a gated write happens at column x*, the cells written
either land inside τ's usable range or pin against a bound. Saturation is the substrate
failing to represent the interval it was shown, and it is locally detectable:

    sat = (fraction of written cells at τ_min) − (fraction at τ_max)

sat > 0 means the interval was too short — the gate met reward too far to the right, so the
chain must be *slower* to move the meeting point left. sat < 0 means the opposite. A
backward pulse from x* carries sat, and each chain cell it passes does
`d[j] ← clip(d[j] + η_a · sat, 1, d_max)`.

Being honest about what kind of signal this is: it is **homeostatic, not appetitive**. It is
derived from the reward-gated update rather than from reward magnitude, and it says "do not
saturate your own representation", not "get more reward". There is no action and no policy
here, so calling it reinforcement learning would be wrong. What it does test is whether the
last hand-set constant in the loop can be removed.

The risk this run exists to measure. `lattice_timescale_notes.md` predicted that closing a
recursive gate↔medium loop would reproduce the self-confirming fixed point that sank the
lateral-input rule — the gate shaping the learning that shapes the gate. There is a reason
to think it might not here: the signal is a saturation *boundary*, not a rhythm, so it
cannot confirm itself the way a period can. That is the hypothesis, and d could equally
well run away or oscillate.

Conditions:
  plastic          — the full loop, d starts at 1 (the degenerate value) and must move
  plastic-unpaired — same loop, reward delay redrawn each trial. No contingency, so no
                     stable place to put the gate: d should not settle anywhere useful.
  fixed-k1         — d pinned at 1: the degenerate baseline to escape
  fixed-k3         — d pinned at 3: the hand-set value the loop is trying to discover
  plastic-two      — delays D and D+16 on alternating trials. One chain, two intervals:
                     does a single learned placement serve both, or compromise?
"""

import os
import sys
import json
import warnings
import numpy as np

warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("AV_L", "64"))
WS = int(os.environ.get("AV_WS", "3"))
PATCH = int(os.environ.get("AV_PATCH", "4"))
ITI = int(os.environ.get("AV_ITI", "260"))
N_TRIALS = int(os.environ.get("AV_TRIALS", "50"))
N_PROBE = int(os.environ.get("AV_PROBE", "3"))
N_SEEDS = int(os.environ.get("AV_N", "3"))
D_LIST = [int(d) for d in os.environ.get("AV_D", "8,20,32").split(",")]
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
ETA_A = float(os.environ.get("AV_ETAA", "0.6"))      # value-driven delay learning rate
DMAX = float(os.environ.get("AV_DMAX", "12"))        # ceiling on a per-column delay
MODES = ("plastic", "plastic-unpaired", "fixed-k1", "fixed-k3", "plastic-two")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def nbr_count(active):
    n = np.zeros(active.shape, np.float32)
    n[1:, :] += active[:-1, :]
    n[:-1, :] += active[1:, :]
    n[:, 1:] += active[:, :-1]
    n[:, :-1] += active[:, 1:]
    return n


def step_chain(phi, inp, act):
    """Advance a 1-D GH chain in place; return which cells are active after the step."""
    rested = phi == 0
    phi[phi >= 1] += 1
    phi[phi > int(np.ceil(act + 1.0))] = 0
    phi[rested & inp] = 1
    return (phi >= 1) & (phi <= act)


def run(mode, D, seed):
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    plastic_d = mode.startswith("plastic")
    d = np.ones(L, np.float64) * (3.0 if mode == "fixed-k3" else 1.0)

    ys0 = L // 2
    patch = np.zeros((L, L), bool)
    patch[ys0 - PATCH:ys0 + PATCH + 1, :WS] = True

    traj = []                                  # (trial, write column, sat, mean d)
    n_written = 0
    # arm 0 is delay D; arm 1 is the second delay of the interleaved condition. Probing
    # both is the whole point of plastic-two — scoring only arm 0 would look like success
    # while saying nothing about whether one chain can hold two intervals.
    pa = {0: [], 1: []}
    pf = {0: [], 1: []}
    wcol = {0: [], 1: []}

    for trial in range(N_TRIALS + N_PROBE):
        probe = trial >= N_TRIALS
        two = mode == "plastic-two"
        arm = 1 if (two and trial % 2 == 1) else 0
        if mode == "plastic-unpaired" and not probe:
            d_rew = int(rng.integers(4, 2 * max(D_LIST)))
        else:
            d_rew = D + (16 if arm else 0)

        # build the attention delay line from the current per-column delays
        K = np.maximum(1, np.rint(d)).astype(np.int64)
        cum = np.concatenate([[0], np.cumsum(K[:-1])])          # gate cell index per column
        GW = float(max(ACT, K.max() + 1))                       # pulse width; see reward-edges
        att_phi = np.zeros(int(cum[-1]) + 2, np.int32)

        phi = np.zeros((L, L), np.int32)
        rew_phi = np.zeros(L, np.int32)
        val_phi = np.zeros(L, np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        prev_rew = np.zeros(L, bool)
        prev_val = np.zeros(L, bool)
        sat = 0.0
        wrote_at = -1

        for t in range(ITI):
            active = (phi >= 1) & (phi <= ACT)
            rested = phi == 0
            drive = patch if t == 0 else False
            excite = rested & ((nbr_count(active) >= THETA) | drive)
            phi[phi >= 1] += 1
            phi[phi > np.ceil(tau).astype(np.int32)] = 0
            phi[excite] = 1
            newly = excite & (t_fire < 0)
            t_fire[newly] = t

            # --- attention chain (delays d), reward chain, value chain ----------------
            a_in = np.zeros(att_phi.size, bool)
            a_in[1:] = ((att_phi[:-1] >= 1) & (att_phi[:-1] <= GW))
            a_in[0] = (t == 0)
            att_all = step_chain(att_phi, a_in, GW)
            att_gate = att_all[cum]

            r_in = np.zeros(L, bool)
            r_in[:-1] = ((rew_phi[1:] >= 1) & (rew_phi[1:] <= ACT))
            r_in[-1] = (t == d_rew)
            rew_act = step_chain(rew_phi, r_in, ACT)
            r_edge = rew_act & ~prev_rew
            prev_rew = rew_act

            # --- the gated write, and the saturation it reveals ------------------------
            cols = r_edge & att_gate
            if cols.any():
                sel = np.zeros((L, L), bool)
                sel[:, cols] = True
                sel &= t_fire >= 0
                if sel.any():
                    dt = (t - t_fire[sel]).astype(np.float32)
                    if probe:
                        e = np.abs(dt - tau[sel])
                        pa[arm].append(float(e.mean()))
                        pf[arm].append(float((e <= 2.0).mean()))
                        wcol[arm].append(float(np.flatnonzero(cols)[0]))
                    else:
                        tau[sel] = np.clip(tau[sel] + ETA * (dt - tau[sel]), TMIN, TMAX)
                        lo = float((tau[sel] <= TMIN + 1e-3).mean())
                        hi = float((tau[sel] >= TMAX - 1e-3).mean())
                        sat = lo - hi
                        wrote_at = int(np.flatnonzero(cols)[0])
                        val_phi[wrote_at] = 1          # launch the value pulse from here

            # --- the value chain: travels ← from the write, retuning delays as it goes --
            v_in = np.zeros(L, bool)
            v_in[:-1] = ((val_phi[1:] >= 1) & (val_phi[1:] <= ACT))
            val_act = step_chain(val_phi, v_in, ACT)
            v_edge = val_act & ~prev_val
            prev_val = val_act
            if plastic_d and not probe and v_edge.any() and sat != 0.0:
                d[v_edge] = np.clip(d[v_edge] + ETA_A * sat, 1.0, DMAX)

        if not probe:
            traj.append((trial, wrote_at, sat, float(d.mean())))
            n_written += int(wrote_at >= 0)

    tr = np.array([(w, s, dm) for _, w, s, dm in traj], float)
    late = tr[-10:] if len(tr) >= 10 else tr
    wl = float(np.nanmean(np.where(late[:, 0] >= 0, late[:, 0], np.nan)))
    # the value pulse only ever travels left from the write site, so the delays it can
    # reach are those upstream of it. Report that region separately from the whole chain.
    cut = int(wl) if np.isfinite(wl) and wl >= 1 else 0
    return {
        "d_final_mean": float(d.mean()),
        "d_upstream": float(d[:cut].mean()) if cut else float("nan"),
        "d_downstream": float(d[cut:].mean()) if cut < L else float("nan"),
        "d_prof": [float(d[int(f * (L - 1))]) for f in (0, .25, .5, .75, 1.0)],
        "write_col_late": wl,
        "write_col_first": float(tr[0, 0]),
        "sat_late": float(np.nanmean(late[:, 1])),
        "sat_first": float(tr[0, 1]),
        "d_traj": [float(x) for x in tr[:, 2]],
        "w_traj": [float(x) for x in tr[:, 0]],
        "trials_written": n_written,
        "align": float(np.mean(pa[0])) if pa[0] else float("nan"),
        "aligned_frac": float(np.mean(pf[0])) if pf[0] else float("nan"),
        "align_b": float(np.mean(pa[1])) if pa[1] else float("nan"),
        "aligned_frac_b": float(np.mean(pf[1])) if pf[1] else float("nan"),
        "wcol_a": float(np.mean(wcol[0])) if wcol[0] else float("nan"),
        "wcol_b": float(np.mean(wcol[1])) if wcol[1] else float("nan"),
    }


def main():
    out = {}
    for D in D_LIST:
        for mode in MODES:
            rs = [run(mode, D, s) for s in range(N_SEEDS)]
            agg = {k: [r[k] for r in rs] for k in rs[0]}
            out[f"D{D}_{mode}"] = agg
            m = lambda k: float(np.nanmean(agg[k]))
            print(f"D={D:2d} {mode:17s}: d̄ {m('d_final_mean'):5.2f}"
                  f" (up {m('d_upstream'):4.2f} / down {m('d_downstream'):4.2f})"
                  f" | write col {m('write_col_first'):5.1f} → {m('write_col_late'):5.1f}"
                  f" | sat {m('sat_first'):+5.2f} → {m('sat_late'):+5.2f}"
                  f" | probe |Δ| {m('align'):6.2f} within ±2 {m('aligned_frac'):5.2f}"
                  f" | wrote {m('trials_written'):4.1f}/{N_TRIALS}", flush=True)
        print(flush=True)

    print(f"learned gate placement vs reward delay (D = {D_LIST}):", flush=True)
    for mode in ("plastic", "fixed-k1", "fixed-k3"):
        w = [float(np.nanmean(out[f"D{D}_{mode}"]["write_col_late"])) for D in D_LIST]
        dd = [float(np.nanmean(out[f"D{D}_{mode}"]["d_final_mean"])) for D in D_LIST]
        print(f"  {mode:17s}: write col " + " / ".join(f"{v:5.1f}" for v in w)
              + "   d̄ " + " / ".join(f"{v:4.2f}" for v in dd), flush=True)

    print("\ntwo interleaved delays — is each arm represented, or compromised?", flush=True)
    for D in D_LIST:
        a = out[f"D{D}_plastic-two"]
        print(f"  D={D:2d} vs {D+16:2d}: arm A |Δ| {np.nanmean(a['align']):5.2f}"
              f" (within ±2 {np.nanmean(a['aligned_frac']):4.2f}) at col"
              f" {np.nanmean(a['wcol_a']):5.1f}"
              f"   |   arm B |Δ| {np.nanmean(a['align_b']):5.2f}"
              f" (within ±2 {np.nanmean(a['aligned_frac_b']):4.2f}) at col"
              f" {np.nanmean(a['wcol_b']):5.1f}", flush=True)

    with open(os.path.join(OUT, "lattice_attention_value.json"), "w") as f:
        json.dump({"L": L, "trials": N_TRIALS, "n": N_SEEDS, "delays": D_LIST,
                   "eta_a": ETA_A, "dmax": DMAX, "tau_range": [TMIN, TMAX],
                   "rows": out}, f, indent=2)
    print("wrote lattice_attention_value.json", flush=True)


if __name__ == "__main__":
    main()
