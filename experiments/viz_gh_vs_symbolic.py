"""Animation: Greenberg-Hastings dynamics vs its symbolic shadow, side by side.

Visual companion to docs/spectrum_sufficiency_proof.md at (tau_a,tau_p)=(2,1).
Two runs, three synchronized views each, one frame per time step:

  left    the 2x2 medium itself (white=receptive, reds=active, blue=passive);
          a dwelling cell (receptive now, still receptive next step) is
          outlined orange, and the step before death carries the
          "DOUBLE DWELL -> death" flag -- the P4 mechanism, visible.
  middle  the space-time raster in loop order (0,1,3,2) with a time cursor:
          a live run shows the diagonal banding of a circulating wave, a dead
          run frays and flatlines.
  right   the full 10-symbol automaton of (2,1) (green = uniformly-live
          classes, black = the zero sink), current symbol haloed, realized
          symbol trail in red. The live run performs exactly one class
          exchange (0022 <-> 1223, the two-class cycle of Theorem 3) and then
          parks on a constant signature (law L2); the dead run wanders grey
          symbols into 0000.

Run selection is deterministic: the live run is the lexicographically first
persistent configuration with signature (0,0,2,2) whose trajectory visits both
exchange classes; the dead run is the first dying configuration with signature
(0,0,1,3), nonzero, whose activity survives >= 5 steps (so the transient is
visible). The dead run triple-dwells at t=4 (asserted).

House Rules Compliance: deterministic, no RNG; the substrate/analysis boundary
is the point of the figure -- left/middle are the dynamics, right is the
analysis construct.
Output: result/topology/gh_vs_symbolic.gif + gh_vs_symbolic_death_frame.png.
"""
import sys, os
sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.animation import FuncAnimation, PillowWriter
from matplotlib.patches import Circle
from itertools import product
from collections import defaultdict

TA, TP = 2, 1
S = TA + TP
B = S + 1
T = 13
ADJ2 = {0: (1, 2), 1: (0, 3), 2: (0, 3), 3: (1, 2)}
CYC = [0, 1, 3, 2]
STATE_COL = {0: "#ffffff", 1: "#ff6b57", 2: "#c21807", 3: "#4477cc"}
POS = {0: (0, 1), 1: (1, 1), 2: (0, 0), 3: (1, 0)}


def step2(cfg):
    return tuple((1 if sum(1 <= cfg[j] <= TA for j in ADJ2[i]) >= 1 else 0)
                 if s == 0 else (s + 1) % B for i, s in enumerate(cfg))


def sig(cfg):
    ph = [cfg[i] for i in CYC]
    return tuple(sorted((ph[(k + 1) % 4] - ph[k]) % B for k in range(4)))


def dwellers(c):
    n = step2(c)
    return [i for i in range(4) if c[i] == 0 and n[i] == 0]


def build():
    fate = {}
    for c in product(range(B), repeat=4):
        if c in fate:
            continue
        path, cur = [], c
        while cur not in fate and cur not in path:
            path.append(cur)
            cur = step2(cur)
        val = fate[cur] if cur in fate else any(
            any(x) for x in path[path.index(cur):])
        for p in path:
            fate[p] = val
    cls_fate, edges = defaultdict(set), defaultdict(set)
    for c, f in fate.items():
        s = sig(c)
        cls_fate[s].add(f)
        edges[s].add(sig(step2(c)))
    return fate, cls_fate, edges


def traj(c0):
    out = [c0]
    for _ in range(T):
        out.append(step2(out[-1]))
    return out


def pick_runs(fate):
    live0 = next(c for c in sorted(fate)
                 if fate[c] and sig(c) == (0, 0, 2, 2)
                 and len({sig(x) for x in traj(c)}) > 1)
    dead0 = next(c for c in sorted(fate)
                 if not fate[c] and sig(c) == (0, 0, 1, 3) and any(c)
                 and sum(1 for x in traj(c) if any(x)) >= 5)
    return live0, dead0


def main():
    fate, cls_fate, edges = build()
    live0, dead0 = pick_runs(fate)
    tr = {"live": traj(live0), "dead": traj(dead0)}
    dd = [t for t in range(T) if len(dwellers(tr["dead"][t])) >= 2]
    assert dd and dd[0] == 4 and tr["dead"][4] == (3, 0, 0, 0), \
        "dead-run anatomy changed -- update the docstring"
    print(f"live start {live0} | dead start {dead0} | dead triple-dwell t={dd[0]}")

    symbols = sorted(cls_fate)
    ang = np.linspace(np.pi/2, np.pi/2 - 2*np.pi, len(symbols), endpoint=False)
    NPOS = {s: (np.cos(a), np.sin(a)) for s, a in zip(symbols, ang)}

    def scol(s):
        f = cls_fate[s]
        return ("#2e8b57" if f == {True}
                else "#222222" if s == (0, 0, 0, 0) else "#aaaaaa")

    fig, axes = plt.subplots(2, 3, figsize=(11.6, 7.6),
                             gridspec_kw=dict(width_ratios=[1, 1.35, 1.5]))
    fig.suptitle("Greenberg\u2013Hastings vs its symbolic shadow \u2014 "
                 "(\u03c4a,\u03c4p)=(2,1), 2\u00d72 core", fontsize=13, y=0.985)

    def draw_frame(t):
        for r, key in enumerate(["live", "dead"]):
            cfgs = tr[key]
            c = cfgs[min(t, T)]
            dw = dwellers(c)
            ax = axes[r, 0]
            ax.clear()
            for i in range(4):
                x, y = POS[i]
                ax.add_patch(plt.Rectangle((x, y), 0.92, 0.92,
                                           fc=STATE_COL[c[i]], ec="k", lw=1.2))
                ax.text(x + 0.46, y + 0.46, str(c[i]), ha="center", va="center",
                        fontsize=15, color="w" if c[i] in (2, 3) else "k")
                if i in dw:
                    ax.add_patch(plt.Rectangle((x - 0.04, y - 0.04), 1.0, 1.0,
                                               fill=False, ec="#ff9900", lw=3.0))
            if len(dw) >= 2:
                ax.text(0.96, -0.28, "DOUBLE DWELL \u2192 death", ha="center",
                        fontsize=10, color="#cc4400", fontweight="bold")
            elif dw:
                ax.text(0.96, -0.28, "one dweller (survivable)", ha="center",
                        fontsize=9, color="#ff9900")
            ax.set_xlim(-0.25, 2.2)
            ax.set_ylim(-0.45, 2.15)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(f"{key.upper()} run \u2014 medium, t={min(t, T)}",
                         fontsize=10)

            ax = axes[r, 1]
            ax.clear()
            M = np.array([[cfgs[u][i] for u in range(T + 1)] for i in CYC])
            img = np.zeros(M.shape + (3,))
            for v, col in STATE_COL.items():
                img[M == v] = matplotlib.colors.to_rgb(col)
            ax.imshow(img, aspect="auto", interpolation="nearest")
            ax.axvline(min(t, T), color="k", lw=2)
            ax.set_yticks(range(4))
            ax.set_yticklabels([f"cell {i}" for i in CYC], fontsize=8)
            ax.set_xlabel("t", fontsize=8)
            ax.tick_params(labelsize=8)
            ax.set_title("space\u2013time (loop order 0,1,3,2)", fontsize=10)

            ax = axes[r, 2]
            ax.clear()
            for s in symbols:
                for u in edges[s]:
                    x0, y0 = NPOS[s]
                    x1, y1 = NPOS[u]
                    if s == u:
                        ax.add_patch(Circle((x0*1.16, y0*1.16), 0.09,
                                            fill=False, ec="#cccccc", lw=0.8))
                    else:
                        ax.annotate("", xy=(x1*0.92, y1*0.92),
                                    xytext=(x0*0.92, y0*0.92),
                                    arrowprops=dict(arrowstyle="-|>",
                                                    color="#cccccc", lw=0.8,
                                                    shrinkA=10, shrinkB=10))
            trail = [sig(x) for x in cfgs[:min(t, T) + 1]]
            for a, b in zip(trail, trail[1:]):
                if a != b:
                    x0, y0 = NPOS[a]
                    x1, y1 = NPOS[b]
                    ax.annotate("", xy=(x1*0.92, y1*0.92),
                                xytext=(x0*0.92, y0*0.92),
                                arrowprops=dict(arrowstyle="-|>",
                                                color="#cc4400", lw=2.0,
                                                shrinkA=10, shrinkB=10))
            for s in symbols:
                x, y = NPOS[s]
                cur = (s == trail[-1])
                ax.add_patch(Circle((x, y), 0.13 if cur else 0.10, fc=scol(s),
                                    ec="#ff9900" if cur else "k",
                                    lw=3 if cur else 0.8, zorder=5))
                ax.text(x*1.38, y*1.38, "".join(map(str, s)), ha="center",
                        va="center", fontsize=7.5,
                        color=scol(s) if scol(s) != "#aaaaaa" else "#777777",
                        fontweight="bold" if cls_fate[s] == {True} else "normal")
            ax.set_xlim(-1.75, 1.75)
            ax.set_ylim(-1.75, 1.75)
            ax.set_aspect("equal")
            ax.axis("off")
            ax.set_title(f"symbol automaton \u2014 "
                         f"sig={''.join(map(str, trail[-1]))}", fontsize=10)
        fig.subplots_adjust(top=0.90, bottom=0.05, hspace=0.35, wspace=0.25)

    outdir = os.path.join(os.path.dirname(os.path.dirname(
        os.path.abspath(__file__))), "result", "topology")
    os.makedirs(outdir, exist_ok=True)
    anim = FuncAnimation(fig, draw_frame, frames=T + 2, interval=900)
    gif = os.path.join(outdir, "gh_vs_symbolic.gif")
    anim.save(gif, writer=PillowWriter(fps=1.2), dpi=105)
    draw_frame(dd[0])
    fig.savefig(os.path.join(outdir, "gh_vs_symbolic_death_frame.png"), dpi=105)
    plt.close(fig)
    sz = os.path.getsize(gif) / 1e6
    print(f"wrote {gif} ({sz:.2f} MB)")
    assert sz < 1.0, "gif exceeds the repo ~1MB guardrail"


if __name__ == "__main__":
    main()
