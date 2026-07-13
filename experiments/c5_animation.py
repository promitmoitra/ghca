"""C5 mechanism animation -- do(chirality) is fat-handed at a fixed locus.

Renders the C5 result (docs/c5_results.md) as a narrated GIF: a single spiral of
TRUE chirality +1 (CCW), nucleated with its core DISPLACED from the lattice
centre (the "disp-8" realization from the C5 table, where the fixed-locus reader
scores 0.10 -- essentially wrong). Two probe loops are overlaid on the same
frame:

  - a FIXED probe at the lattice centre (blue outline) -- the naive "read the
    wave at a fixed cortical locus" reader E7's router uses.
  - a TRACKED probe that follows the actual core (orange outline) -- the
    topology-aware reader.

Each probe's outline is colour-coded green/red by whether its local winding
readout currently agrees with the TRUE chirality. The punchline: the true
rotation never changes, but the fixed reader loses it (drifts red) once the
core has moved away from the lattice centre, while the tracked reader stays
locked on (green) throughout -- fat-handedness is a property of the (variable,
reader) pair, not the chirality itself.

Uses the real C5/E7 machinery (`e7_spiral_option.local_winding` /
`signed_charge`, `c5_do_chirality.realize`), so the readouts are the genuine
substrate computation.

Output
------
docs/figures/c5_fixed_vs_tracked.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
from ghca_net_viz import state_colors, animate
import e7_spiral_option as sp
import c5_do_chirality as c5

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

L = sp.L
TAU = sp.TAU
C = L // 2
D = 14                  # core displacement (beyond the C5 "disp-12" failing case)
CHIR = +1               # true chirality: CCW
PROBE_RAD = 6           # matches e7_spiral_option.local_winding default
T_RUN = 90

# fixed-locus reader: green/crimson: tracked reader: teal/orange (so identity is
# readable by hue even where the two rings overlap, correctness by which pair)
_FIXED_OK, _FIXED_BAD = np.array([0.16, 0.70, 0.20]), np.array([0.80, 0.08, 0.10])
_TRACK_OK, _TRACK_BAD = np.array([0.05, 0.55, 0.55]), np.array([0.90, 0.45, 0.05])


def _ring_pixels(cx, cy, rad):
    pts = []
    for x in range(cx - rad, cx + rad + 1):
        pts.append((x, cy - rad)); pts.append((x, cy + rad))
    for y in range(cy - rad, cy + rad + 1):
        pts.append((cx - rad, y)); pts.append((cx + rad, y))
    return [(min(max(x, 0), L - 1), min(max(y, 0), L - 1)) for x, y in pts]


def main():
    net = c5.make_spiral()
    # nucleate once, off-centre by D at a fixed angle (reproducible for the GIF)
    rng = np.random.default_rng(0)
    theta = 0.0                                    # displace straight along +x
    cx0 = int(round(C + D * np.cos(theta)))
    cy0 = int(round(C + D * np.sin(theta)))
    c5.realize(net, CHIR, "displaced", rng, d=D)
    # realize() picks a random angle internally; force the exact reproducible core
    xs = np.arange(L)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    ang = np.arctan2(Y - cy0, X - cx0)
    frac = ((CHIR * ang) / (2 * np.pi)) % 1.0
    net.phi = np.clip(np.floor(frac * (TAU + 1)).astype(np.int64), 0, TAU).reshape(-1).copy()

    phi_roll = np.zeros((T_RUN, L * L), np.int64)
    fixed_ok, tracked_ok = [], []
    for t in range(T_RUN):
        net.step(None)
        phi_roll[t] = net.phi
        fixed_sign = np.sign(sp.local_winding(net.phi, C, C, rad=PROBE_RAD))
        _, _, _, core = sp.signed_charge(net.phi)
        tcx, tcy = core
        if np.isfinite(tcx):
            tracked_sign = np.sign(sp.local_winding(net.phi, int(round(tcx)), int(round(tcy)), rad=PROBE_RAD))
        else:
            tracked_sign = 0
        fixed_ok.append(bool(fixed_sign == CHIR))
        tracked_ok.append(bool(tracked_sign == CHIR))

    # build colour rollout: base GH state colouring + two probe-ring overlays
    C_roll = np.zeros((T_RUN, L, L, 4))
    for t in range(T_RUN):
        base = state_colors(phi_roll[t], sp.ACT, TAU).reshape(L, L, 4)
        C_roll[t] = base
        for (px, py) in _ring_pixels(C, C, PROBE_RAD):
            C_roll[t, px, py, :3] = _FIXED_OK if fixed_ok[t] else _FIXED_BAD
            C_roll[t, px, py, 3] = 1.0
        _, _, _, core = sp.signed_charge(phi_roll[t].reshape(-1))
        tcx, tcy = core
        if np.isfinite(tcx):
            for (px, py) in _ring_pixels(int(round(tcx)), int(round(tcy)), PROBE_RAD):
                col = _TRACK_OK if tracked_ok[t] else _TRACK_BAD
                C_roll[t, px, py, :3] = col       # tracked ring drawn on top, undiluted
                C_roll[t, px, py, 3] = 1.0

    captions = []
    for t in range(T_RUN):
        f = "OK" if fixed_ok[t] else "WRONG"
        tr = "OK" if tracked_ok[t] else "WRONG"
        captions.append(f"t={t}  true chirality=CCW  |  fixed-locus reader: {f}  |  "
                         f"tracked reader: {tr}")

    dummy = np.zeros((T_RUN, L * L), np.int64)
    out = os.path.join(FIGDIR, "c5_fixed_vs_tracked.gif")
    animate(dummy, sp.ACT, TAU, layout="grid", L=L, out=out,
            fps=10, stride=2, dpi=95, colors_rollout=C_roll.reshape(T_RUN, L * L, 4),
            figsize=(7.2, 7.8),
            title="C5: do(chirality) is fat-handed at a fixed locus, not when tracked",
            captions=captions)
    print(f"fixed reader correct: {np.mean(fixed_ok):.2f}   "
          f"tracked reader correct: {np.mean(tracked_ok):.2f}")
    print("wrote", out)


if __name__ == "__main__":
    main()
