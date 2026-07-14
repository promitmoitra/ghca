"""Visualisation / animation for the Greenberg-Hastings network substrate.

The network experiments (E-series, C-series) all run on ``ghca_net.Network`` but
only ever emit *static* summary plots (learning curves, bar charts, phase
diagrams). The thesis, though, is about *dynamics* -- "the medium generates its
own repertoire of patterns." This module renders that: given a rollout of node
phases over time it animates the excitable activity -- travelling waves on a
chain, rotating pulses on a ring, spirals on a lattice -- and writes a GIF you
can actually watch.

A *rollout* is a ``(T, N)`` integer array of node phases ``phi`` over time --
exactly what ``Network.run(T, record=True)["phi"]`` returns. Given a rollout, a
node layout, and the substrate's ``(act, tau)``, :func:`animate` writes a GIF.

State encoding (from ``ghca_net``): ``phi == 0`` rested; ``1..act`` active;
``act+1..tau`` refractory (counting down).

Example
-------
    from ghca_net import Network, lattice2d
    from ghca_net_viz import animate
    net = Network(lattice2d(40, r=1), act=2, pas=8, theta=1.0, p_s=1e-3, seed=1)
    net.seed_random()
    roll = net.run(120, record=True)["phi"]
    animate(roll, net.act, net.tau, layout="grid", L=40, out="spiral.gif")
"""

import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib import animation


# ---------------------------------------------------------------------------
# State -> colour
# ---------------------------------------------------------------------------

# rested is a light neutral; active is a bright red; refractory fades from a
# deep red down toward near-black as the phase counts back to rest -- so a wave
# reads as a bright crest with a darkening tail, the way excitable media look.
_RESTED = np.array([0.90, 0.90, 0.90, 1.0])
_ACTIVE = np.array([0.83, 0.09, 0.11, 1.0])
_REFR_HEAD = np.array([0.45, 0.00, 0.05, 1.0])
_REFR_TAIL = np.array([0.12, 0.12, 0.17, 1.0])


def state_colors(phi_row, act, tau):
    """Map a phase row ``(N,)`` to an ``(N, 4)`` RGBA array by GH state."""
    phi_row = np.asarray(phi_row)
    tau = np.broadcast_to(np.asarray(tau), phi_row.shape)
    rgba = np.tile(_RESTED, (len(phi_row), 1))
    active = (phi_row >= 1) & (phi_row <= act)
    refr = phi_row > act
    rgba[active] = _ACTIVE
    if refr.any():
        span = np.maximum(tau[refr] - act, 1)
        frac = ((phi_row[refr] - act) / span)[:, None]      # 0 head -> 1 tail
        rgba[refr] = (1 - frac) * _REFR_HEAD + frac * _REFR_TAIL
    return rgba


# ---------------------------------------------------------------------------
# Layouts: map N nodes to 2D positions (scatter modes) or a grid (imshow).
# ---------------------------------------------------------------------------

def layout_positions(layout, N, L=None, pos=None):
    """Return ``(pos, mode)``. mode is 'grid' or 'scatter'.

    layout: 'line' | 'ring' | 'grid' | 'free' (uses `pos`).
    """
    if layout == "grid":
        if L is None:
            L = int(round(np.sqrt(N)))
        return None, "grid"
    if layout == "line":
        p = np.column_stack([np.arange(N), np.zeros(N)])
    elif layout == "ring":
        ang = 2 * np.pi * np.arange(N) / N
        p = np.column_stack([np.cos(ang), np.sin(ang)])
    elif layout == "free":
        p = np.asarray(pos, float)
    else:
        raise ValueError(f"unknown layout {layout!r}")
    return p, "scatter"


# ---------------------------------------------------------------------------
# The animator
# ---------------------------------------------------------------------------

def animate(rollout, act, tau, layout="line", L=None, pos=None, out="anim.gif",
            fps=12, stride=1, title=None, captions=None, colors_rollout=None,
            marker_size=None, figsize=None, dpi=90, annotations=None,
            grid_shape=None):
    """Animate a phase rollout and write a GIF.

    Parameters
    ----------
    rollout : (T, N) int array   -- phases over time (Network.run(...)['phi']).
    act : int                    -- active duration.
    tau : int or (N,) array      -- per-node cycle length (net.tau).
    layout : str                 -- 'line' | 'ring' | 'grid' | 'free'.
    L : int                      -- grid side (layout='grid').
    pos : (N,2) array            -- explicit positions (layout='free').
    out : str                    -- output GIF path.
    fps : int                    -- frames per second.
    stride : int                 -- keep every `stride`-th frame (thins long runs).
    title : str                  -- static suptitle.
    captions : list[str] or None -- per-frame subtitle text (narration); indexed
                                    by the *original* time step.
    colors_rollout : (T,N,4) or None -- override per-frame colours (e.g. to tint
                                    two competing waves); default colours-by-state.
    annotations : list[(x, y, text)] or None -- static labels drawn in data
                                    coordinates (scatter layouts), e.g. to title
                                    side-by-side panels.
    grid_shape : (rows, cols) or None -- for layout='grid', reshape each frame to
                                    this instead of the default square (L, L). Use
                                    to render a non-square field or two lattices
                                    side by side (cols = 2*L + gap).
    """
    rollout = np.asarray(rollout)
    T, N = rollout.shape
    frames = list(range(0, T, stride))
    p, mode = layout_positions(layout, N, L=L, pos=pos)

    def frame_colors(t):
        if colors_rollout is not None:
            return np.asarray(colors_rollout[t])
        return state_colors(rollout[t], act, tau)

    if mode == "grid":
        if grid_shape is None:
            if L is None:
                L = int(round(np.sqrt(N)))
            grid_shape = (L, L)
        rows, cols = grid_shape
        fig, ax = plt.subplots(figsize=figsize or (5.2, 5.6))
        img = ax.imshow(frame_colors(0).reshape(rows, cols, 4), interpolation="nearest")
        ax.set_xticks([]); ax.set_yticks([])
    else:
        span_x = np.ptp(p[:, 0]) or 1.0
        span_y = np.ptp(p[:, 1]) or 1.0
        if figsize is None:
            if layout == "line":
                figsize = (min(14, 2 + N * 0.22), 2.4)
            elif layout == "ring":
                figsize = (5.4, 5.8)
            else:
                figsize = (6.0, 6.2)
        if marker_size is None:
            marker_size = 380 if layout == "line" else (260 if layout == "ring" else 90)
        fig, ax = plt.subplots(figsize=figsize)
        scat = ax.scatter(p[:, 0], p[:, 1], s=marker_size, c=frame_colors(0),
                          edgecolors="none")
        pad_x = 0.06 * span_x + 0.5
        pad_y = 0.06 * span_y + 0.5
        ax.set_xlim(p[:, 0].min() - pad_x, p[:, 0].max() + pad_x)
        ax.set_ylim(p[:, 1].min() - pad_y, p[:, 1].max() + pad_y)
        ax.set_aspect("equal")
        ax.axis("off")

    for (ax_x, ax_y, txt) in (annotations or []):
        ax.text(ax_x, ax_y, txt, ha="center", va="center", fontsize=10,
                fontweight="bold")

    sub = ax.set_title("", fontsize=10)
    if title:
        fig.suptitle(title, fontsize=12, y=0.99)

    def update(t):
        c = frame_colors(t)
        if mode == "grid":
            img.set_data(c.reshape(grid_shape[0], grid_shape[1], 4))
            artists = [img]
        else:
            scat.set_facecolors(c)
            artists = [scat]
        cap = captions[t] if captions is not None and t < len(captions) else f"t = {t}"
        sub.set_text(cap)
        artists.append(sub)
        return artists

    ani = animation.FuncAnimation(fig, update, frames=frames, blit=False)
    fig.tight_layout()
    ani.save(out, writer=animation.PillowWriter(fps=fps))
    plt.close(fig)
    return out


if __name__ == "__main__":
    # demo: a spiral-forming lattice at the E0 organised operating point.
    import os
    from ghca_net import Network, lattice2d
    ROOT = os.path.dirname(os.path.abspath(__file__))
    net = Network(lattice2d(44, r=2), act=6, pas=8, theta=4.0, p_s=5e-3, seed=3)
    net.seed_random(0.15, 0.15)
    roll = net.run(150, record=True)["phi"]
    out = os.path.join(ROOT, "docs", "figures", "demo_lattice.gif")
    animate(roll[40:], net.act, net.tau, layout="grid", L=44, out=out,
            fps=12, stride=2, dpi=80,
            title="GH substrate: self-organised spiral waves (E0 point)")
    print("wrote", out)
