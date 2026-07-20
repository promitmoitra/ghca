"""Track 3e.2 / closes 4a — an emergent fast/slow timescale hierarchy.

Track 4a asked: under a two-rhythm drive, do per-node timescales `τ` self-organise
into a fast/slow **hierarchy** (a bimodal `τ` distribution, clusters at the two drive
periods)? `docs/e10_notes.md` found the existing Line B rule **cannot**: `τ ← τ +
η(interfire − τ)` reads a node's *own* inter-fire interval, which once `τ` overshoots
the period reports a *multiple* of it, ratcheting `τ` upward to the ceiling — no fast
population forms. That was the block.

The fix (validated in 3d-emergent for delays, applied here to periods): update `τ`
from the interval between **input arrivals** the node senses — presynaptic drive,
registered *regardless of the node's own refractory state* — not its own firing. The
overshoot that corrupted the old signal is invisible to this one, so `τ` locks to the
true drive period and can move *down* as well as up.

Two channel-assignment modes:
  * **wired** — first half of the pool fed only the fast source, second half only the
    slow (the `e10_diagnostics.py` diagnostic-3 setup). Isolates the *rule*.
  * **emergent** — every node fans in from *both* sources; a competitive Hebbian rule
    with a population **conscience** (the fix for e10 diagnostic-2, where the more
    frequent fast rhythm otherwise swamps the slow) splits the pool ~50/50 and each
    group's `τ` locks to its channel — a fully self-organised hierarchy.

Rules compared: **input** (new, input-interval) vs **own** (old, e10 self-referential).
Drives: **two** (fast+slow) vs single-rhythm **fast** / **slow** controls (should stay
unimodal). Metric: Sarle bimodality coefficient (BC > 5/9) + fraction of `τ` near each
period. Cross-frequency coupling (theta–gamma nesting) needs inter-population coupling
this pool does not have — deferred, and stated as such.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ghca_net import Network  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "20"))
NPOOL = int(os.environ.get("STATS_NPOOL", "120"))
STEPS = int(os.environ.get("STATS_STEPS", "10000"))
P_F, P_S = 6, 24                 # fast / slow drive periods
ACT = 2
TMIN, TMAX = 3, 34
ETA_TAU, ETA_W, CONSC = 0.15, 0.03, 0.6
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def run(rule, channel, drive, seed):
    """rule in {input,own}; channel in {wired,emergent}; drive in {two,fast,slow}."""
    rng = np.random.default_rng(seed)
    M = NPOOL + 2
    fast, slow = NPOOL, NPOOL + 1
    W = np.zeros((M, M))
    if channel == "wired":
        W[:NPOOL // 2, fast] = 1.5
        W[NPOOL // 2:NPOOL, slow] = 1.5
    else:                                    # emergent: fan in from both, plastic
        W[:NPOOL, fast] = 1.2 * (1 + 0.2 * rng.standard_normal(NPOOL))
        W[:NPOOL, slow] = 1.2 * (1 + 0.2 * rng.standard_normal(NPOOL))
    net = Network(W, act=ACT, pas=8, theta=np.full(M, 1.0), p_s=0.0, seed=seed)
    net.tau = net.tau.astype(float)
    net.tau[:NPOOL] = rng.uniform(TMIN, TMAX, NPOOL)
    net.tau[fast], net.tau[slow] = P_F, P_S

    last_fire = np.full(NPOOL, -1.0)
    last_in = np.full(NPOOL, -1.0)
    drive_fast = drive in ("two", "fast")
    drive_slow = drive in ("two", "slow")

    for t in range(STEPS):
        d = np.zeros(M, bool)
        if drive_fast and t % P_F == 0:
            d[fast] = True
        if drive_slow and t % P_S == 0:
            d[slow] = True
        prev = net.phi.copy()
        net.step(d)
        ff, sf = bool(d[fast]), bool(d[slow])

        if rule == "own":
            for i in np.where((net.phi[:NPOOL] == 1) & (prev[:NPOOL] == 0))[0]:
                if last_fire[i] >= 0:
                    net.tau[i] = np.clip(net.tau[i] + ETA_TAU * ((t - last_fire[i]) - net.tau[i]),
                                         TMIN, TMAX)
                last_fire[i] = t
            continue
        if not (ff or sf):
            continue

        # input rule: commit each node to a channel (conscience-balanced for emergent),
        # then track the inter-arrival interval of that committed channel only.
        if channel == "emergent":
            wf, ws = net.W[:NPOOL, fast], net.W[:NPOOL, slow]
            p_fast = float((wf > ws).mean())
            dom_fast = (wf - CONSC * (p_fast - 0.5)) >= (ws - CONSC * ((1 - p_fast) - 0.5))
        else:
            dom_fast = np.arange(NPOOL) < NPOOL // 2
        for i in range(NPOOL):
            domf = bool(dom_fast[i])
            if channel == "emergent":
                if domf:
                    net.W[i, fast] = min(net.W[i, fast] + ETA_W, 2.5)
                    net.W[i, slow] = max(net.W[i, slow] - ETA_W, 0.2)
                else:
                    net.W[i, slow] = min(net.W[i, slow] + ETA_W, 2.5)
                    net.W[i, fast] = max(net.W[i, fast] - ETA_W, 0.2)
            fired_dom = (domf and ff) or ((not domf) and sf)
            if fired_dom:
                if last_in[i] >= 0 and (t - last_in[i]) <= P_S * 1.5:
                    net.tau[i] = np.clip(net.tau[i] + ETA_TAU * ((t - last_in[i]) - net.tau[i]),
                                         TMIN, TMAX)
                last_in[i] = t

    tau = net.tau[:NPOOL].copy()
    dom = (net.W[:NPOOL, fast] > net.W[:NPOOL, slow]) if channel == "emergent" \
        else (np.arange(NPOOL) < NPOOL // 2)
    return tau, dom


def summarise(tau):
    bc = st.bimodality(tau)
    return {"bc": bc["bc"], "bimodal": bool(bc["bimodal"]),
            "near_f": float((np.abs(tau - P_F) <= 2).mean()),
            "near_s": float((np.abs(tau - P_S) <= 3).mean()),
            "mean": float(tau.mean()), "sd": float(tau.std())}


CONDS = [
    ("input", "wired", "two"),      # rule works with channels given (repairs diag-3)
    ("input", "emergent", "two"),   # fully self-organised hierarchy (the 4a goal)
    ("own", "wired", "two"),        # old rule (reproduces e10 ratchet failure)
    ("own", "emergent", "two"),
    ("input", "emergent", "fast"),  # single-rhythm controls: should stay unimodal
    ("input", "emergent", "slow"),
]


def main():
    out = {}
    for rule, channel, drive in CONDS:
        key = f"{rule}+{channel}+{drive}"
        bcs = np.zeros(N); nf = np.zeros(N); ns = np.zeros(N)
        taus = None
        for s in range(N):
            tau, dom = run(rule, channel, drive, s)
            m = summarise(tau)
            bcs[s], nf[s], ns[s] = m["bc"], m["near_f"], m["near_s"]
            if s == 0:
                taus = tau
        out[key] = {"bc": bcs, "near_f": nf, "near_s": ns, "tau0": taus}
        bm, blo, bhi = st.bootstrap_ci(bcs)
        print(f"{key:24s}: BC={bm:.2f}[{blo:.2f},{bhi:.2f}] "
              f"nearP_f={nf.mean():.2f} nearP_s={ns.mean():.2f} "
              f"bimodal={(bcs > 5/9).mean():.0%} of seeds", flush=True)

    save = {"n": N, "P_F": P_F, "P_S": P_S, "npool": NPOOL}
    for k in out:
        save[f"{k}_bc"] = out[k]["bc"]
        save[f"{k}_near_f"] = out[k]["near_f"]
        save[f"{k}_near_s"] = out[k]["near_s"]
        save[f"{k}_tau0"] = out[k]["tau0"]
    np.savez(os.path.join(OUT, "timescale_hierarchy.npz"), **save)
    summ = {"n": N, "P_F": P_F, "P_S": P_S, "rows": {}}
    for k in out:
        bm, blo, bhi = st.bootstrap_ci(out[k]["bc"])
        summ["rows"][k] = {"bc": [bm, blo, bhi], "near_f": float(out[k]["near_f"].mean()),
                           "near_s": float(out[k]["near_s"].mean()),
                           "frac_bimodal": float((out[k]["bc"] > 5/9).mean())}
    with open(os.path.join(OUT, "timescale_hierarchy.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote timescale_hierarchy.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
