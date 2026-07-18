"""Track 3b — does the E0 substrate story generalise off the lattice?

E0 characterised the self-sustaining band and the period~tau law on `lattice2d`
only; the design doc's `smallworld` default was never used. This re-runs the two
topology-agnostic E0 observables on `lattice2d` (reference), `smallworld`, and
`rgg`, matched on node count and (approximately) mean degree:

  1. Self-sustaining band: A_ss (late active fraction, p_s=0) vs theta -- is there
     a live excitable band, and where?
  2. period ~ tau: at a persistent operating point, does the dominant global
     period track the local timescale tau (E0's `period = 1.00*tau + 0.95`)?

Random topologies (smallworld, rgg) are run over several graph seeds with CIs
(ghca_stats); lattice2d is deterministic. Spiral-core structure is lattice-only
(needs 2-D geometry) and is out of scope here -- this asks whether the *excitable
dynamics* generalise, not the 2-D spiral.
"""

import os
import sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ghca_net import Network, lattice2d, smallworld, rgg, dominant_period  # noqa: E402
import ghca_stats as st  # noqa: E402

OUT = os.path.join(ROOT, "result", "e0_topo")
FIG = os.path.join(ROOT, "docs", "figures")
os.makedirs(OUT, exist_ok=True)
os.makedirs(FIG, exist_ok=True)

L = 40
N = L * L                      # 1600 nodes, matched across topologies
ACT, PAS = 6, 8                # tau = 14 (E0 operating point)
NSEED = int(os.environ.get("STATS_N", "8"))   # graph seeds for random topologies


def build(topo, seed=0):
    if topo == "lattice2d":
        return lattice2d(L, r=2, periodic=True)
    if topo == "smallworld":
        return smallworld(N, k=12, beta=0.1, seed=seed)   # k~lattice r=2 degree
    if topo == "rgg":
        return rgg(N, radius=0.05, seed=seed)             # mean degree ~12
    raise ValueError(topo)


def mean_degree(W):
    return float((W > 0).sum(1).mean())


def a_ss(W, theta, act=ACT, pas=PAS, seed=0, T=600, burn=400):
    net = Network(W, act=act, pas=pas, theta=float(theta), p_s=0.0, seed=seed)
    net.seed_random(0.30, 0.30)
    out = net.run(T, record=False)
    return float(out["A"][burn:].mean())


PERIOD_THETA = 1.0   # E0 measures the period~tau law in the loop-dominated regime


def period_at(W, tau, theta=PERIOD_THETA, act=2, seed=0, T=600, burn=350):
    # match E0.period_vs_tau exactly: act=2 fixed, pas = tau - act, theta=1
    net = Network(W, act=act, pas=int(tau) - act, theta=float(theta), p_s=0.0, seed=seed)
    net.seed_random(0.30, 0.30)
    out = net.run(T, record=False)
    return dominant_period(out["A"][burn:])


def band_sweep(topo, thetas):
    seeds = [0] if topo == "lattice2d" else list(range(NSEED))
    grid = np.zeros((len(seeds), len(thetas)))
    for si, s in enumerate(seeds):
        W = build(topo, s)
        for ti, th in enumerate(thetas):
            grid[si, ti] = a_ss(W, th, seed=100 + s)
    return grid  # (seeds, thetas)


def period_sweep(topo, taus):
    """period~tau in the loop-dominated regime (theta=1), matching E0."""
    seeds = [0] if topo == "lattice2d" else list(range(NSEED))
    grid = np.full((len(seeds), len(taus)), np.nan)
    for si, s in enumerate(seeds):
        W = build(topo, s)
        for ki, tau in enumerate(taus):
            grid[si, ki] = period_at(W, tau, seed=200 + s)
    return grid


def main():
    thetas = np.arange(1, 14)
    taus = np.array([4, 6, 8, 10, 12, 14, 16, 18, 20])
    topos = ["lattice2d", "smallworld", "rgg"]

    degs = {t: mean_degree(build(t, 0)) for t in topos}
    print("mean degree:", {t: round(d, 1) for t, d in degs.items()}, flush=True)

    bands = {}
    band_theta = {}
    for t in topos:
        g = band_sweep(t, thetas)
        bands[t] = g
        m = g.mean(0)
        # pick a persistent operating point: highest theta with A_ss in a live band
        live = np.where((m > 0.02) & (m < 0.45))[0]
        band_theta[t] = int(thetas[live[-1]]) if len(live) else int(thetas[np.argmax(m)])
        print(f"[{t}] deg={degs[t]:.1f} A_ss vs theta: "
              + " ".join(f"{th}:{v:.2f}" for th, v in zip(thetas, m))
              + f"  -> op theta={band_theta[t]}", flush=True)

    periods = {}
    fits = {}
    for t in topos:
        g = period_sweep(t, taus)
        periods[t] = g
        m = np.nanmean(g, 0)
        fin = np.isfinite(m)
        if fin.sum() >= 2:
            sl, ic = np.polyfit(taus[fin], m[fin], 1)
            r = np.corrcoef(taus[fin], m[fin])[0, 1]
        else:
            sl = ic = r = np.nan
        fits[t] = (sl, ic, r)
        print(f"[{t}] period~tau (theta={PERIOD_THETA:.0f}): slope={sl:.3f} "
              f"intercept={ic:.3f} r={r:.4f} ({fin.sum()}/{len(taus)} oscillating)",
              flush=True)

    np.savez(os.path.join(OUT, "e0_topo.npz"),
             thetas=thetas, taus=taus, degrees=np.array([degs[t] for t in topos]),
             op_theta=np.array([band_theta[t] for t in topos]),
             **{f"band_{t}": bands[t] for t in topos},
             **{f"period_{t}": periods[t] for t in topos},
             **{f"fit_{t}": np.array(fits[t]) for t in topos})
    print("wrote e0_topo.npz", flush=True)


if __name__ == "__main__":
    main()
