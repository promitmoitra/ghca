"""Track 3e.1 — re-tiling under a shifting delay distribution.

The capability the *wired* timescale basis provably lacks: adaptation. 3d-emergent
grew a τ basis to tile a fixed delay range; here the delay distribution *shifts*
between regimes, and we ask two things a hand-set basis cannot answer:

  1. **Plasticity** — does the emergent basis *reallocate* τ to tile the new regime?
  2. **Stability** — does adapting to the new regime *destroy* the old regime's
     tiling? i.e. does the plastic *representation* have its own stability–plasticity
     frontier — the 3c interference story one level down, in the dynamics rather than
     the readout?

Protocol: from a near-homogeneous start, self-organise τ (the 3d-emergent
input-timing-driven competitive rule, reward-free) through a sequence of delay
regimes SHORT → LONG → SHORT. After each stage, measure delay-decodability on *every*
regime's delay bins → an accuracy matrix M[after_stage, eval_regime]. Compared against
a wired-static basis (τ hand-set to tile SHORT, never updated) — the no-adaptation
reference.

Read-outs: plasticity = the current-regime diagonal (does it adapt?); representation
backward transfer = decode on SHORT after adapting to LONG, minus decode on SHORT
right after the first SHORT stage (negative ⇒ the representation forgot). The τ
histogram after each stage shows the migration directly.
"""

import os
import sys
import json
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
import continual_temporal_saturation as B      # make(), constants  # noqa: E402
import ghca_stats as st  # noqa: E402

N = int(os.environ.get("STATS_N", "20"))
N_PRES = int(os.environ.get("STATS_NPRES", "1000"))     # selforg presentations per stage
ETA_TAU, WIN_FRAC, CONSC_C, CONSC_RATE = 0.20, 0.25, 6.0, 0.05
SHORT = [2, 5, 8, 11]
LONG = [18, 24, 30, 36]
REGIMES = {"SHORT": SHORT, "LONG": LONG}
STAGES = ["SHORT", "LONG", "SHORT"]                      # the shifting schedule
OUT = os.path.join(ROOT, "result", "stats")
os.makedirs(OUT, exist_ok=True)


def probe_response(net, roles, delay):
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


def selforg(net, roles, delays, rng, n_pres=N_PRES):
    """Input-timing-driven competitive τ rule (3d-emergent), on `delays`."""
    hid = roles["hidden"]; nH = len(hid)
    p_win = np.full(nH, WIN_FRAC); n_win = max(1, int(WIN_FRAC * nH))
    for _ in range(n_pres):
        d = delays[int(rng.integers(len(delays)))]
        resp = probe_response(net, roles, d)
        elig = resp > 0
        score = resp - CONSC_C * (p_win - WIN_FRAC)
        score[~elig] = -1e9
        win = np.argsort(-score)[:n_win]; win = win[elig[win]]
        won = np.zeros(nH, bool); won[win] = True
        p_win += CONSC_RATE * (won.astype(float) - p_win)
        net.tau[hid[win]] = np.clip(net.tau[hid[win]] + ETA_TAU * (d - net.tau[hid[win]]),
                                    2, 40)


def decode(net, roles, delays, rng, reps=25):
    """Nearest-centroid delay-bin decodability of the probe-evoked pattern."""
    X, y = [], []
    for di, d in enumerate(delays):
        for _ in range(reps):
            X.append(probe_response(net, roles, d)); y.append(di)
    X = np.array(X); y = np.array(y)
    cents = np.array([X[y == c].mean(0) for c in range(len(delays))])
    pred = np.argmin(((X[:, None] - cents[None]) ** 2).sum(-1), 1)
    return float((pred == y).mean())


def run_seed(seed, adaptive):
    net = B.make(seed, "homog")
    hid = net.roles["hidden"]
    rng = np.random.default_rng(seed + 5)
    net.tau[hid] = np.clip(B.TAU_HOMOG + 0.5 * rng.standard_normal(len(hid)), 2, 40)
    if not adaptive:                                    # wired-static: tile SHORT, freeze
        net.tau[hid] = np.linspace(SHORT[0], SHORT[-1], len(hid))
    net.tau_base = net.tau.copy()
    reg_names = list(REGIMES)
    M = np.zeros((len(STAGES), len(reg_names)))         # [after_stage, eval_regime]
    tau_snaps = []
    srng = np.random.default_rng(seed + 11)
    for si, stage in enumerate(STAGES):
        if adaptive:
            selforg(net, net.roles, REGIMES[stage], srng)
        tau_snaps.append(net.tau[hid].copy())
        for ri, rn in enumerate(reg_names):
            M[si, ri] = decode(net, net.roles, REGIMES[rn], np.random.default_rng(seed + 20 + ri))
    return M, np.array(tau_snaps)


def main():
    reg_names = list(REGIMES)
    out = {}
    for adaptive, name in [(True, "emergent"), (False, "wired-static")]:
        Ms = np.zeros((N, len(STAGES), len(reg_names)))
        taus = None
        for s in range(N):
            M, ts = run_seed(s, adaptive)
            Ms[s] = M
            if s == 0:
                taus = ts
        out[name] = {"M": Ms, "tau_snaps": taus}
        Mm = Ms.mean(0)
        print(f"=== {name} (n={N}) — decode M[after_stage, eval_regime] ===", flush=True)
        print("            " + "  ".join(f"{r:>6s}" for r in reg_names), flush=True)
        for si, stage in enumerate(STAGES):
            print(f"  after {stage:5s} " + "  ".join(f"{Mm[si, ri]:6.2f}"
                  for ri in range(len(reg_names))), flush=True)
        # representation backward transfer: SHORT decode after LONG minus after first SHORT
        si_first = 0; si_after_long = 1; short_idx = reg_names.index("SHORT")
        bwt = Ms[:, si_after_long, short_idx] - Ms[:, si_first, short_idx]
        m, lo, hi = st.bootstrap_ci(bwt)
        print(f"  representation backward transfer (SHORT: afterLONG − afterSHORT) = "
              f"{m:+.3f} [{lo:+.3f}, {hi:+.3f}]", flush=True)

    save = {"n": N, "stages": np.array(STAGES, dtype=object), "regimes": reg_names,
            "short": np.array(SHORT), "long": np.array(LONG)}
    for name in out:
        save[f"{name}_M"] = out[name]["M"]
        save[f"{name}_tau"] = out[name]["tau_snaps"]
    np.savez(os.path.join(OUT, "continual_temporal_retile.npz"), **save)
    summ = {"n": N, "stages": STAGES, "regimes": reg_names, "rows": {}}
    for name in out:
        Mm = out[name]["M"].mean(0)
        summ["rows"][name] = {f"after_{STAGES[si]}_{si}": {reg_names[ri]: float(Mm[si, ri])
                              for ri in range(len(reg_names))} for si in range(len(STAGES))}
    with open(os.path.join(OUT, "continual_temporal_retile.json"), "w") as f:
        json.dump(summ, f, indent=2, default=float)
    print("\nwrote continual_temporal_retile.{npz,json}", flush=True)


if __name__ == "__main__":
    main()
