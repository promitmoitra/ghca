"""Animations for the three capacity mechanisms on a 2-D lattice.

Illustrates the claims in [`lattice_capacities.md`](../docs/lattice_capacities.md):
each cognitive capacity uses a *different substrate primitive*, and each one is
visible in the dynamics.

  1. E2 memory     -- a plane wave circulating a torus row. The loop transit
                      length is L, so it sustains iff tau < L. Two panels at the
                      same L, one either side of the boundary: the sustaining
                      loop keeps going, the over-slow one dies.
  2. E4 attention  -- two wavefronts launched from opposite edges collide and
                      annihilate in each other's refractory wake. The bias moves
                      the annihilation locus, which decides who owns the centre.
                      No inhibitory units exist anywhere in the substrate.
  3. E5 control    -- the context ring as a held option: alive vs dead ring,
                      same stimulus, different downstream gating.

Outputs
-------
docs/figures/lattice_e2_memory.gif
docs/figures/lattice_e4_attention.gif
docs/figures/lattice_e5_option.gif
"""

import os
import sys

import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network, lattice2d  # noqa: E402
from ghca_net_viz import animate  # noqa: E402

FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

GAP = 2                      # blank columns between side-by-side lattices


def _plane_wave(net, L, tau):
    """Seed a full-row wavefront with a refractory tail behind it."""
    net.phi[:] = 0
    g = net.phi.reshape(L, L)
    g[0, :] = 1
    for k in range(1, int(tau)):
        g[(-k) % L, :] = min(k + 1, int(tau))
    net.phi = g.reshape(-1)
    return net


def side_by_side(roll_a, roll_b, L, fill=0):
    """Stitch two (T, L*L) rollouts into one (T, L*(2L+GAP)) grid field."""
    T = min(len(roll_a), len(roll_b))
    cols = 2 * L + GAP
    out = np.full((T, L * cols), fill, dtype=np.int64)
    for t in range(T):
        a = np.asarray(roll_a[t]).reshape(L, L)
        b = np.asarray(roll_b[t]).reshape(L, L)
        f = np.full((L, cols), fill, dtype=np.int64)
        f[:, :L] = a
        f[:, L + GAP:] = b
        out[t] = f.reshape(-1)
    return out, cols


# ---------------------------------------------------------------------------
def anim_e2_memory(L=16, tau_live=10, tau_dead=22, T=90, seed=0):
    """Reentry on a torus row: tau < L sustains, tau >= L dies."""
    rolls = []
    for tau in (tau_live, tau_dead):
        net = Network(lattice2d(L, r=1, periodic=True), act=2, pas=int(tau) - 2,
                      theta=1.0, p_s=0.0, seed=seed)
        _plane_wave(net, L, tau)
        rolls.append(net.run(T, record=True)["phi"])
    field, cols = side_by_side(rolls[0], rolls[1], L)
    caps = [f"t={t}   left: tau={tau_live} < L={L} (loop holds)   |   "
            f"right: tau={tau_dead} > L={L} (dies)" for t in range(len(field))]
    out = os.path.join(FIGDIR, "lattice_e2_memory.gif")
    animate(field, act=2, tau=max(tau_live, tau_dead), layout="grid",
            grid_shape=(L, cols), out=out, fps=9, stride=2,
            title="E2 memory on a 2-D torus: a circulating plane wave, "
                  "gated by tau < L",
            captions=caps, marker_size=26, figsize=(9.2, 4.2), dpi=85)
    return out


def anim_e4_attention(L=31, bias=3, T=None, seed=1):
    """Two edge-launched wavefronts collide; bias shifts the annihilation locus."""
    T = T or int(1.6 * L)
    rolls = []
    for b in (0, bias):
        net = Network(lattice2d(L, r=1, periodic=False), act=6, pas=6,
                      theta=1.0, p_s=0.0, seed=seed)
        net.phi[:] = 0
        N = L * L
        left = np.arange(L) * L + 0
        right = np.arange(L) * L + (L - 1)
        s0, s1 = max(0, -b), max(0, b)
        frames = []
        for t in range(T):
            drive = np.zeros(N, bool)
            if t == s0:
                drive[left] = True
            if t == s1:
                drive[right] = True
            net.step(drive if drive.any() else None)
            frames.append(net.phi.copy())
        rolls.append(np.array(frames))
    field, cols = side_by_side(rolls[0], rolls[1], L)
    caps = [f"t={t}   left: no bias (collision at centre)   |   "
            f"right: bias={bias} to the LEFT stream (locus pushed right)"
            for t in range(len(field))]
    out = os.path.join(FIGDIR, "lattice_e4_attention.gif")
    animate(field, act=6, tau=12, layout="grid", grid_shape=(L, cols),
            out=out, fps=10, stride=1,
            title="E4 attention on a 2-D arena: competition BY wave annihilation "
                  "(no inhibitory units)",
            captions=caps, marker_size=16, figsize=(9.2, 4.2), dpi=85)
    return out


def anim_e5_option(L=16, tau_live=10, tau_dead=20, T=80, seed=0):
    """The context ring as a held option: alive vs dead, on a torus row."""
    rolls = []
    for tau in (tau_live, tau_dead):
        net = Network(lattice2d(L, r=1, periodic=True), act=2, pas=int(tau) - 2,
                      theta=1.0, p_s=0.0, seed=seed)
        net.phi[:] = 0
        g = net.phi.reshape(L, L)
        g[0, :L // 2] = 1                        # a partial front = the cue
        for k in range(1, int(tau)):
            g[(-k) % L, :L // 2] = min(k + 1, int(tau))
        net.phi = g.reshape(-1)
        rolls.append(net.run(T, record=True)["phi"])
    field, cols = side_by_side(rolls[0], rolls[1], L)
    caps = [f"t={t}   left: option HELD (tau={tau_live} < L) -> gate open   |   "
            f"right: option LAPSED (tau={tau_dead} > L) -> gate shut"
            for t in range(len(field))]
    out = os.path.join(FIGDIR, "lattice_e5_option.gif")
    animate(field, act=2, tau=max(tau_live, tau_dead), layout="grid",
            grid_shape=(L, cols), out=out, fps=8, stride=2,
            title="E5 executive control: the context loop is a HELD OPTION "
                  "(E2's primitive, reused)",
            captions=caps, marker_size=26, figsize=(9.2, 4.2), dpi=85)
    return out


def main():
    for fn in (anim_e2_memory, anim_e4_attention, anim_e5_option):
        p = fn()
        print(f"wrote {p}  ({os.path.getsize(p)/1e6:.2f} MB)")


if __name__ == "__main__":
    main()
