"""Track 3d (emergent arm) — a *learned* timescale basis, grown not hand-set.

The 3d wired arm showed a hand-set graded-τ basis buys continual-learning capacity
on temporal (delay-keyed) tasks that a homogeneous basis cannot represent. This arm
asks the harder question — the one Track 4a is blocked on: can that τ spread be
*grown* by a local, reward-free plasticity rule instead of wired by hand?

`e10_notes.md` diagnosed why the existing Line B rule can't: `τ ← τ + η(interfire − τ)`
reads a node's *own* inter-fire interval, which is corrupted once τ overshoots the
period (it reads a multiple), so τ only ratchets upward. The fix it proposed — and
that this implements — is a **bidirectional, input-timing-driven** rule: a node nudges
its τ toward the *delays it is externally probed at*, not toward its own firing
interval. Combined with E9's competitive grouping (k-WTA + DeSieno conscience) so
nodes specialise to different delays, the τ distribution self-organises to *tile* the
delay range. The teaching signal (the probe delay Δ) is external and the update is
bidirectional, so the ratchet never arises.

Self-organisation (reward-free, label-free): repeatedly present a stimulus, wait a
delay Δ drawn from the bins, probe. Only nodes that actually *respond* to the probe
(rested, τ < Δ — an activity gate, so the rule is substrate-driven) are eligible; a
k-WTA + conscience competition picks winners, which nudge their τ toward Δ. From a
near-homogeneous start the τ distribution spreads to cover the delay range.

Then freeze τ and run the 3d continual-learning sweep exactly as the wired arm, over
three bases — homogeneous / **emergent** (grown here) / wired (hand-set) — under both
head modes. Prediction (mirroring the E9 bridge): emergent ≈ wired ≫ homogeneous.
Only this arm makes the *dynamics* (not just a readout) carry the task structure, so
only it retires the plastic-dynamics caveat.

τ is grown once per seed and cached (self-organisation depends only on the seed, not
on T or head mode or the reward task), so the emergent basis costs one selforg pass
per seed plus the same continual-learning sweep as the other bases.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import continual_temporal_saturation as B  # noqa: E402  (the wired harness)
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "20"))
T_LIST = [int(t) for t in os.environ.get("STATS_TLIST", "2,3,4,5,6").split(",")]
N_PRES = int(os.environ.get("STATS_NPRES", "1200"))     # selforg presentations
ETA_TAU = 0.20
WIN_FRAC = 0.25
CONSC_C = 6.0
CONSC_RATE = 0.05
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)

_orig_make = B.make
EMERGENT_TAU = {}                                        # seed -> grown hidden-τ vector


def _probe_response(net, roles, delay):
    """Run one stimulus→delay→probe and return per-hidden probe-evoked activity."""
    hid = roles["hidden"]
    net.phi[:] = 0; net.t = 0
    for _ in range(B.SETTLE):
        net.step(None)
    drive = net.sensory_drive(0)
    for _ in range(B.CUE):
        net.step(drive)
    for _ in range(delay):
        net.step(None)
    acc = np.zeros(len(hid))
    for k in range(B.PROBE + B.WWIN):
        net.step(drive if k < B.PROBE else None)
        acc += net.active_mask()[hid].astype(float)
    return acc / (B.PROBE + B.WWIN)


def grow_emergent_tau(seed):
    """Self-organise the hidden-τ distribution from a near-homogeneous start with the
    input-timing-driven competitive rule. Returns the grown τ (hidden nodes)."""
    net = _orig_make(seed, "homog")                      # all hidden τ = TAU_HOMOG
    hid = net.roles["hidden"]; nH = len(hid)
    rng = np.random.default_rng(seed + 3)
    net.tau[hid] = np.clip(B.TAU_HOMOG + 0.5 * rng.standard_normal(nH), 2, 40)  # break ties
    p_win = np.full(nH, WIN_FRAC)
    n_win = max(1, int(WIN_FRAC * nH))
    for _ in range(N_PRES):
        delay = B.DELAYS[int(rng.integers(B.K_DELAY))]
        resp = _probe_response(net, net.roles, delay)     # activity gate (substrate-driven)
        elig = resp > 0
        score = resp - CONSC_C * (p_win - WIN_FRAC)
        score[~elig] = -1e9
        win = np.argsort(-score)[:n_win]
        win = win[elig[win]]
        won = np.zeros(nH, bool); won[win] = True
        p_win += CONSC_RATE * (won.astype(float) - p_win)
        # winners nudge τ toward the *external* probe delay (bidirectional; no ratchet)
        net.tau[hid[win]] = np.clip(net.tau[hid[win]] + ETA_TAU * (delay - net.tau[hid[win]]),
                                    2, 40)
    return net.tau[hid].copy()


def make_with_emergent(seed, basis):
    if basis == "emergent":
        net = _orig_make(seed, "homog")
        net.tau[net.roles["hidden"]] = EMERGENT_TAU[seed]
        net.tau_base = net.tau.copy()
        return net
    return _orig_make(seed, basis)


B.make = make_with_emergent                              # so B.run_sequence sees "emergent"

BASES = ["homog", "emergent", "graded"]
HEADS = ["shared", "per-task"]


def main():
    # 1) grow (and cache) the emergent τ once per seed
    tau_hist = {b: [] for b in BASES}
    for s in range(N):
        EMERGENT_TAU[s] = grow_emergent_tau(s)
    print(f"grew emergent τ for {N} seeds", flush=True)
    # record τ distributions (seed 0) for the figure
    for b in BASES:
        net = B.make(0, b)
        tau_hist[b] = net.tau[net.roles["hidden"]].copy()

    # 2) continual-learning sweep, three bases × two head modes
    out = {f"{b}+{h}": {"avg": [], "bwt": []} for b in BASES for h in HEADS}
    for basis in BASES:
        for head in HEADS:
            arm = f"{basis}+{head}"
            for T in T_LIST:
                avg = np.zeros(N); bwt = np.zeros(N)
                for s in range(N):
                    R, jt = B.run_sequence(s, basis, head, T)
                    avg[s], bwt[s] = B.cl_metrics(R, jt)
                out[arm]["avg"].append(avg); out[arm]["bwt"].append(bwt)
                am, alo, ahi = st.bootstrap_ci(avg)
                bm, blo, bhi = st.bootstrap_ci(bwt)
                print(f"{arm:18s} T={T}: avg={am:.3f}[{alo:.3f},{ahi:.3f}] "
                      f"bwt={bm:+.3f}[{blo:+.3f},{bhi:+.3f}]", flush=True)

    save = {"n": N, "chance": B.CHANCE, "n_h": B.N_H, "k_delay": B.K_DELAY,
            "delays": np.array(B.DELAYS), "T_list": np.array(T_LIST)}
    for arm in out:
        save[f"{arm}_avg"] = np.array(out[arm]["avg"])
        save[f"{arm}_bwt"] = np.array(out[arm]["bwt"])
    for b in BASES:
        save[f"tau_{b}"] = np.array(tau_hist[b])
    np.savez(os.path.join(OUT, "continual_temporal_emergent.npz"), **save)
    summ = {"n": N, "chance": B.CHANCE, "n_h": B.N_H, "delays": B.DELAYS,
            "tau_std": {b: float(np.std(tau_hist[b])) for b in BASES}, "rows": {}}
    for arm in out:
        for i, T in enumerate(T_LIST):
            am, alo, ahi = st.bootstrap_ci(out[arm]["avg"][i])
            bm, blo, bhi = st.bootstrap_ci(out[arm]["bwt"][i])
            summ["rows"][f"{arm}:T{T}"] = {"avg": [am, alo, ahi], "bwt": [bm, blo, bhi]}
    with open(os.path.join(OUT, "continual_temporal_emergent.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote continual_temporal_emergent.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
