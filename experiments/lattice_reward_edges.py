"""Reward as a fourth edge: can the substrate learn a stimulus-reward interval?

Where this comes from. `lattice_attention_gate.py` established that a 1-D attention chain
made of the same Greenberg-Hastings cells can carry a timing reference to arbitrary depth
into a recurrent medium, where lateral input alone cannot. But that result is thin in one
specific way: the chain's τ is hand-set, and all the medium does is *copy a clock*. The
content of the representation is designed, not discovered.

This asks whether **reward** can supply content instead. Geometry (all four edges used):

      ┌──────── attention chain (1-D, travels →, seeded by stimulus onset) ────────┐
      │                                                                           │
   stimulus                          plastic medium                            (reward
    patch                          (τ learned per cell)                         onset)
      │                                                                           │
      └──────── reward chain (1-D, travels ←, seeded at the right edge at t=D) ────┘

Both chains are **labelled lines by construction** — separate rows, exactly as this repo's
`roles = {sensory, hidden, motor}` convention licenses. That matters, because the obvious
cheaper design fails: if reward were a suprathreshold wave *inside* the medium, then "input
arriving from my right" would not mean "reward". A rightward-travelling wave excites cell x
from the left at t, and one step later cell x-1 sees x active on its right — so every cell
in a passing wave measures an interval of ~1 and floors its τ. Direction is not pathway
identity, the same way amplitude was not (coincidence gating) and phase was not (gated
own-input). A dedicated chain sidesteps it.

The rule is three-factor, and each factor is a different physical thing:
  *who*   — eligibility: the cell's own firing time this trial, set by the stimulus wave
            reaching it, so the medium's 2-D geometry selects who can learn
  *when*  — the attention chain, optionally gating which columns may update
  *what*  — the reward chain: on reward arrival, τ ← τ + η((t_now − t_fire) − τ)

So τ converges to *this cell's own stimulus-to-reward interval*. That is not a clock copy:
the interval differs cell by cell because it depends on how long the stimulus wave took to
arrive, which is wave kinematics, not a designed value.

The kinematics also give a coincidence for free. The attention pulse is at column t at
time t; reward reaches column x at t = D + (Lx-1-x). They coincide where
x* = (D + Lx - 1)/2 — so **attention-gating the update writes the interval into a single
column whose position encodes the delay**, while leaving the gate off writes a distributed
map. Both are reported.

Conditions:
  paired      — reward at a fixed delay D after stimulus; no attention gate (distributed)
  paired-att  — same, but only columns where the attention gate is open may update (sparse)
  unpaired    — reward delay redrawn at random every trial, same rate, no contingency.
                THE control that matters: the reward line still delivers timed input, so
                anything that survives here is about extra input rather than the association.
  no-reward   — reward line never fires; τ keeps its random initialisation (floor baseline)

Test (deliberately not the training objective, which τ satisfies by construction). After
training, run probe trials at the canonical delay and measure how well each cell's τ still
predicts its own stimulus-to-reward interval: alignment = |(t_rew − t_fire) − τ|. A
paired-trained medium should be aligned; an unpaired-trained one should not, even though it
received exactly as many reward events.

Note on τ range: intervals here run to ~Lx + D, well past the [3, 34] used elsewhere in the
lattice work, so the range is widened to [3, 90]. Nothing else changes.
"""

import os
import sys
import json
import warnings
import numpy as np

# the no-reward condition deliberately produces empty slices; nan is the intended answer
warnings.filterwarnings("ignore", category=RuntimeWarning)

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)

L = int(os.environ.get("REW_L", "64"))            # square sheet, walls on all four sides
WS = int(os.environ.get("REW_WS", "3"))           # stimulus patch width (in x)
PATCH = int(os.environ.get("REW_PATCH", "4"))     # stimulus patch half-height (in y)
ITI = int(os.environ.get("REW_ITI", "220"))       # trial length; waves must die out
N_TRIALS = int(os.environ.get("REW_TRIALS", "30"))
N_PROBE = int(os.environ.get("REW_PROBE", "3"))
N_SEEDS = int(os.environ.get("REW_N", "3"))
D_LIST = [int(d) for d in os.environ.get("REW_D", "8,20,32").split(",")]
JIT = int(os.environ.get("REW_JIT", "0"))         # per-trial jitter of the stimulus patch (y)
ACT, THETA = 2.0, 1.0
TMIN, TMAX, ETA = 3.0, 90.0, 0.4
TAU_CH = 3.0                                      # both chains' own τ (hand-set, as before)
# "-kN" slows the attention chain to one cell per N columns of medium — a physically longer
# delay line, N GH cells per column, not a fudge factor on its speed.
MODES = ("paired", "paired-att", "paired-att-k3", "paired-att-k6", "unpaired", "no-reward")
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)
TAG = os.environ.get("REW_TAG", "")


def nbr_count(active):
    """4-neighbour count on a sheet with walls on all four sides."""
    n = np.zeros(active.shape, np.float32)
    n[1:, :] += active[:-1, :]
    n[:-1, :] += active[1:, :]
    n[:, 1:] += active[:, :-1]
    n[:, :-1] += active[:, 1:]
    return n


def step_chain(phi, inp, act=ACT):
    """Advance a 1-D GH chain in place; return which cells are active after the step.

    `act` is the pulse width. It has to be at least K+1 for a K-slowed attention chain:
    the gate and the reward pulse approach at a relative speed of (K+1)/K columns per step,
    so with a 2-step pulse they cross *between* columns and never overlap at any single one
    unless (D+L-1) mod (K+1) <= 1. Widening the pulse to K+1 guarantees a meeting column.
    Same lesson as the attention-window sweep in `lattice_attention_gate.py`.
    """
    rested = phi == 0
    phi[phi >= 1] += 1
    phi[phi > int(np.ceil(act + 1.0))] = 0
    phi[rested & inp] = 1
    return (phi >= 1) & (phi <= act)


def run(mode, D, seed):
    rng = np.random.default_rng(seed)
    tau = rng.uniform(TMIN, TMAX, (L, L)).astype(np.float32)
    tau0 = tau.copy()
    K = int(mode.split("-k")[1]) if "-k" in mode else 1      # attention cells per column
    GW = float(max(ACT, K + 1))                              # gate pulse width (see step_chain)
    n_upd = np.zeros((L, L), np.int32)
    upd_cols = np.zeros(L, np.int64)
    align, aligned_frac, probe_n = [], [], 0

    for trial in range(N_TRIALS + N_PROBE):
        probe = trial >= N_TRIALS
        if mode == "unpaired" and not probe:
            d = int(rng.integers(4, 2 * max(D_LIST)))   # same rate, no contingency
        else:
            d = D

        # the stimulus patch, jittered in y per trial when JIT > 0
        ys = L // 2 + (int(rng.integers(-JIT, JIT + 1)) if JIT else 0)
        patch = np.zeros((L, L), bool)
        patch[max(0, ys - PATCH):ys + PATCH + 1, :WS] = True

        phi = np.zeros((L, L), np.int32)
        att_phi = np.zeros(K * L, np.int32)
        rew_phi = np.zeros(L, np.int32)
        t_fire = np.full((L, L), -1, np.int32)
        prev_rew = np.zeros(L, bool)

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

            # --- the two chains, both made of the same cells -------------------------
            a_in = np.zeros(K * L, bool)
            a_in[1:] = ((att_phi[:-1] >= 1) & (att_phi[:-1] <= GW))
            a_in[0] = (t == 0)                       # seeded by stimulus onset
            att_gate = step_chain(att_phi, a_in, GW)[::K]   # one gate per medium column

            r_in = np.zeros(L, bool)
            r_in[:-1] = ((rew_phi[1:] >= 1) & (rew_phi[1:] <= ACT))
            # no-reward: the line is silent during training, but fires on probe trials so
            # that an untrained τ can still be scored against the true intervals
            r_in[-1] = (t == d) and (mode != "no-reward" or probe)
            rew_act = step_chain(rew_phi, r_in)

            # --- reward arrival is the learning event ---------------------------------
            r_edge = rew_act & ~prev_rew
            prev_rew = rew_act
            cols = r_edge & att_gate if "att" in mode else r_edge
            if cols.any():
                sel = np.zeros((L, L), bool)
                sel[:, cols] = True
                sel &= t_fire >= 0                   # eligibility: I fired this trial
                if sel.any():
                    dt = (t - t_fire[sel]).astype(np.float32)
                    if probe:
                        e = np.abs(dt - tau[sel])
                        align.append(float(e.mean()))
                        aligned_frac.append(float((e <= 2.0).mean()))
                    else:
                        tau[sel] = np.clip(tau[sel] + ETA * (dt - tau[sel]), TMIN, TMAX)
                        n_upd[sel] += 1
                        upd_cols[cols] += 1
        probe_n += int(probe)

    touched = n_upd > 0
    if not touched.any():                       # no-reward: score the whole sheet instead
        touched = np.ones((L, L), bool)
    # spatial structure: does τ vary along y at fixed x? (if not, the 2-D medium is idle)
    ystd = []
    for x in range(L):
        col = tau[touched[:, x], x]
        if col.size > 3:
            ystd.append(col.std())
    cols_used = np.flatnonzero(upd_cols > 0)
    return {
        "align": float(np.mean(align)) if align else float("nan"),
        "aligned_frac": float(np.mean(aligned_frac)) if aligned_frac else float("nan"),
        "frac_touched": float(touched.mean()),
        "tau_touched": float(tau[touched].mean()) if touched.any() else float("nan"),
        "tau_init": float(tau0[touched].mean()) if touched.any() else float("nan"),
        "y_std": float(np.mean(ystd)) if ystd else float("nan"),
        "col_lo": int(cols_used.min()) if cols_used.size else -1,
        "col_hi": int(cols_used.max()) if cols_used.size else -1,
        "col_mid": float(cols_used.mean()) if cols_used.size else float("nan"),
    }


def main():
    out = {}
    for D in D_LIST:
        for mode in MODES:
            rs = [run(mode, D, s) for s in range(N_SEEDS)]
            agg = {k: [r[k] for r in rs] for k in rs[0]}
            out[f"D{D}_{mode}"] = agg
            m = lambda k: float(np.nanmean(agg[k]))
            sd = lambda k: float(np.nanstd(agg[k]))
            print(f"D={D:2d} {mode:11s}: probe |Δ| {m('align'):6.2f} ±{sd('align'):4.2f}"
                  f" | within ±2 {m('aligned_frac'):5.2f}"
                  f" | cells touched {m('frac_touched'):5.3f}"
                  f" | τ {m('tau_touched'):5.1f} (init {m('tau_init'):5.1f})"
                  f" | σ_y {m('y_std'):5.2f}"
                  f" | cols {int(m('col_lo')):2d}–{int(m('col_hi')):2d}"
                  f" mid {m('col_mid'):5.1f}", flush=True)
        print(flush=True)

    # Where does an attention-gated write land? The gate reaches column x at t=Kx and reward
    # at t=D+L-1-x, so they meet at x* = (D+L-1)/(K+1) — position slope 1/(K+1) per unit delay.
    print(f"attention-gated write column vs delay (D = {D_LIST}):", flush=True)
    ds = np.array(D_LIST, float)
    for mode in [m for m in MODES if "att" in m]:
        K = int(mode.split("-k")[1]) if "-k" in mode else 1
        mids = np.array([np.nanmean(out[f"D{D}_{mode}"]["col_mid"]) for D in D_LIST])
        got = " / ".join(f"{v:5.1f}" for v in mids)
        pr = " / ".join(f"{(D + L - 1)/(K + 1):5.1f}" for D in D_LIST)
        line = f"  {mode:14s} (K={K}): measured {got}   predicted {pr}"
        if np.isfinite(mids).all() and len(ds) > 1:
            line += (f"   slope {np.polyfit(ds, mids, 1)[0]:+.3f}"
                     f" (predicted {1/(K + 1):+.3f})")
        print(line, flush=True)

    with open(os.path.join(OUT, f"lattice_reward_edges{TAG}.json"), "w") as f:
        json.dump({"L": L, "WS": WS, "patch": PATCH, "iti": ITI, "trials": N_TRIALS,
                   "probe": N_PROBE, "n": N_SEEDS, "delays": D_LIST,
                   "tau_range": [TMIN, TMAX], "eta": ETA, "rows": out}, f, indent=2)
    print(f"wrote lattice_reward_edges{TAG}.json", flush=True)


if __name__ == "__main__":
    main()
