"""Topology — the GGH winding number as the sustain criterion (refines the
cycle-space capacity bound).

`topology_cycle_capacity.py` bounds how many reentrant loops a topology admits,
using E2's *length gate* (`tau < L`) to decide which cycles can carry a
circulating wave. That gate is an approximation. The exact criterion is older
and comes from the combinatorial companion to the GH model itself:

  Greenberg, Greene & Hastings, "A combinatorial problem arising in the study of
  reaction-diffusion equations", SIAM J. Alg. Disc. Meth. 1(1), 1980.

GGH characterise which initial conditions lead to periodic patterns via the
**winding number of a continuous cycle**: walk once around a closed path summing
the phase increments; if the total is a nonzero multiple of 2*pi, that cycle
carries a rotating wave. The repo already applies exactly this construction to
2x2 lattice plaquettes in `ghca_net.count_phase_singularities`. This experiment
generalises it to an *arbitrary cycle in the graph*, which is what the capacity
packing needs (it works on general cycles, not lattice plaquettes).

What this establishes
---------------------
1. **The invariant behaves as GGH predict.** On a directed ring carrying one
   rotating pulse the winding is a nonzero constant; on a dead ring it is 0.
2. **Winding is an exact sustain predictor where the length gate is not.** On a
   (L, tau) grid of directed rings, steady-state winding predicts persistence
   45/45, while the length gate scores 40/45. Every length-gate miss is the case
   `tau == L` -- precisely the death boundary that `e2_results.md` flags as
   marginal and that the length-gate packing excludes conservatively via a
   strict `>`. Winding resolves that boundary from first principles instead of
   sidestepping it.
3. **Measurement timing matters, and is a real caveat.** Winding must be read at
   steady state. Measured too early (t=5) it scores 39/45 -- WORSE than the
   length gate -- because a slow pulse has not yet wrapped enough of the ring to
   register a full turn. This is a property of the measurement, not of the
   invariant, and it is why the criterion is applied here to settled dynamics.

Substrate / analysis boundary
-----------------------------
Unlike `topology_cycle_capacity.py` (pure graph analysis), the winding number is
a **dynamical** quantity: it is read off a phase configuration, so it requires
running the substrate. That is the cost of the added exactness. It therefore
cannot replace the length gate inside a purely topological pre-dynamical bound;
it *calibrates* that bound, and tells us the length gate's only systematic error
is at `tau == L`.

Outputs
-------
docs/figures/topology_winding_capacity.png
result/topology/winding_capacity.npz
"""

import os
import sys

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
from ghca_net import Network  # noqa: E402

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "topology")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

LS = [8, 12, 16, 24, 32]
TAUS = [4, 6, 8, 12, 16, 20, 24, 28, 32]
T_RUN = 300
SS = slice(250, 300)          # steady-state window for the winding read
EARLY = 5                     # deliberately-too-early read, for the caveat
SEED = 0


def cycle_winding(phi, cycle, tau):
    """GGH winding number of a continuous cycle: turns of the phase per circuit.

    Sums wrapped phase increments along the closed vertex path `cycle` and
    divides by 2*pi. Nonzero => the cycle carries a rotating wave. This is the
    same construction `ghca_net.count_phase_singularities` applies to 2x2
    plaquettes, lifted to an arbitrary cycle.
    """
    ang = 2 * np.pi * np.asarray(phi, dtype=float) / tau
    wrap = lambda d: (d + np.pi) % (2 * np.pi) - np.pi          # noqa: E731
    total = sum(wrap(ang[cycle[(i + 1) % len(cycle)]] - ang[cycle[i]])
                for i in range(len(cycle)))
    return int(round(total / (2 * np.pi)))


def directed_ring(L):
    """Directed ring, matching e2_delayed_response.two_ring().

    A bare pulse on a BIDIRECTIONAL ring splits and the two fronts annihilate
    head-on, so reentry tests must use a directed ring.
    """
    W = np.zeros((L, L))
    for i in range(L):
        W[i, (i - 1) % L] = 1.0
    return W


def probe(L, tau, seed=SEED):
    """Run one seeded pulse on a directed ring; return the three verdicts."""
    net = Network(directed_ring(L), act=2, pas=int(tau) - 2, theta=1.0,
                  p_s=0.0, seed=seed)
    net.phi[:] = 0
    net.phi[0] = 1
    out = net.run(T_RUN, record=True)
    alive = bool(out["A"][-30:].mean() > 0)
    cyc = list(range(L))
    ws = [cycle_winding(out["phi"][t], cyc, tau) for t in range(*SS.indices(T_RUN))]
    w_ss = max(set(ws), key=ws.count)                 # modal steady-state winding
    w_early = cycle_winding(out["phi"][EARLY], cyc, tau)
    return dict(alive=alive, w_ss=w_ss, w_early=w_early, pred_len=bool(tau < L))


def main():
    print("GGH winding number vs the E2 length gate, on directed rings")
    print(f"{'L':>4}{'tau':>5}{'len gate':>10}{'winding':>9}{'actual':>9}")
    rows = []
    for L in LS:
        for tau in TAUS:
            r = probe(L, tau)
            rows.append((L, tau, r["pred_len"], r["w_ss"] != 0,
                         r["w_early"] != 0, r["alive"]))
            flag = "" if r["pred_len"] == (r["w_ss"] != 0) == r["alive"] else \
                   ("  len WRONG" if (r["w_ss"] != 0) == r["alive"] else "  wind WRONG")
            print(f"{L:>4}{tau:>5}{str(r['pred_len']):>10}"
                  f"{str(r['w_ss'] != 0):>9}{str(r['alive']):>9}{flag}")

    n = len(rows)
    acc_len = sum(1 for r in rows if r[2] == r[5]) / n
    acc_ss = sum(1 for r in rows if r[3] == r[5]) / n
    acc_early = sum(1 for r in rows if r[4] == r[5]) / n
    miss_len = [(r[0], r[1]) for r in rows if r[2] != r[5]]
    miss_ss = [(r[0], r[1]) for r in rows if r[3] != r[5]]

    print(f"\n  length gate         : {acc_len:.3f}  ({n - len(miss_len)}/{n})"
          f"   misses {miss_len}")
    print(f"  winding, steady     : {acc_ss:.3f}  ({n - len(miss_ss)}/{n})"
          f"   misses {miss_ss}")
    print(f"  winding, read early : {acc_early:.3f}  (caveat: too-early reads"
          f" underperform the length gate)")
    print(f"  every length-gate miss is tau == L: "
          f"{all(L == t for L, t in miss_len)}")

    plot(rows, acc_len, acc_ss, acc_early)

    np.savez(os.path.join(DATADIR, "winding_capacity.npz"),
             Ls=np.array([r[0] for r in rows]),
             taus=np.array([r[1] for r in rows]),
             pred_len=np.array([r[2] for r in rows]),
             pred_wind_ss=np.array([r[3] for r in rows]),
             pred_wind_early=np.array([r[4] for r in rows]),
             alive=np.array([r[5] for r in rows]),
             acc_len=acc_len, acc_wind_ss=acc_ss, acc_wind_early=acc_early)
    print("\nwrote", os.path.join(DATADIR, "winding_capacity.npz"))


def plot(rows, acc_len, acc_ss, acc_early):
    Ls, Ts = LS, TAUS
    grid_a = np.zeros((len(Ls), len(Ts)))
    grid_l = np.zeros_like(grid_a)
    grid_w = np.zeros_like(grid_a)
    for L, tau, pl, pw, _, a in rows:
        i, j = Ls.index(L), Ts.index(tau)
        grid_a[i, j], grid_l[i, j], grid_w[i, j] = a, pl, pw

    fig, axes = plt.subplots(1, 3, figsize=(10.6, 3.2))
    for ax, g, title in zip(
            axes, [grid_a, grid_l, grid_w],
            ["Observed: does the pulse persist?",
             f"E2 length gate  τ<L   ({acc_len:.0%})",
             f"GGH winding ≠ 0   ({acc_ss:.0%})"]):
        ax.imshow(g, cmap="Blues", vmin=0, vmax=1, aspect="auto")
        ax.set_xticks(range(len(Ts)))
        ax.set_xticklabels(Ts, fontsize=6)
        ax.set_yticks(range(len(Ls)))
        ax.set_yticklabels(Ls, fontsize=6)
        ax.set_xlabel(r"local timescale $\tau$", fontsize=7)
        ax.set_title(title, loc="left", fontsize=7.5)
        # mark the cells where the panel's prediction differs from observation
        if g is not grid_a:
            for i in range(len(Ls)):
                for j in range(len(Ts)):
                    if g[i, j] != grid_a[i, j]:
                        ax.plot(j, i, "x", color="#dc2626", ms=7, mew=1.8)
    axes[0].set_ylabel("ring length $L$", fontsize=7)
    fig.suptitle("The GGH winding number resolves the τ = L boundary the length "
                 "gate misses (red × = wrong cell)",
                 fontsize=8.5, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    path = os.path.join(FIGDIR, "topology_winding_capacity.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


if __name__ == "__main__":
    main()
