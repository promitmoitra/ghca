"""E7 mechanism animation -- rotation direction as the rule (2D spiral option).

Renders the E7 mechanism (docs/e7_results.md) as a narrated GIF: two genuine 2D
spiral waves on the E0 organised-band lattice (no-flux boundaries), seeded with
opposite handedness. The left core rotates counter-clockwise (CCW, net charge +1)
and the right clockwise (CW, net charge -1). In E7 Phase B these two handednesses
*are* the two task rules (CCW -> identity, CW -> reversal), so the picture is
literally "the rule is which way the spiral turns" -- the 2D analogue of E5's
option loop, on the substrate's native medium.

Uses the real E7 nucleation (`e7_spiral_option.make` / `seed_spiral`), so the
spirals and their persistence are the genuine substrate behaviour.

Output
------
docs/figures/e7_spiral_rule.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
from ghca_net_viz import state_colors, animate
import e7_spiral_option as e7

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

L = e7.L
GAP = 3                      # white columns between the two lattices
T_RUN = 140
WARMUP = 20                  # drop the nucleation transient


def spiral_rollout(chirality, seed=0):
    net = e7.make(seed=seed)
    e7.seed_spiral(net, chirality)
    roll = np.zeros((T_RUN, L * L), np.int64)
    for t in range(T_RUN):
        net.step(None)
        roll[t] = net.phi
    return roll[WARMUP:]


def main():
    ccw = spiral_rollout(+1)          # net charge +1 -> rule 0 (identity)
    cw = spiral_rollout(-1)           # net charge -1 -> rule 1 (reversal)
    T = ccw.shape[0]
    tau = e7.TAU

    # build a combined (L, 2L+GAP) colour field per frame: CCW | gap | CW
    cols = 2 * L + GAP
    white = np.array([1.0, 1.0, 1.0, 1.0])
    C = np.zeros((T, L, cols, 4))
    for t in range(T):
        left = state_colors(ccw[t], e7.ACT, tau).reshape(L, L, 4)
        right = state_colors(cw[t], e7.ACT, tau).reshape(L, L, 4)
        C[t, :, :L] = left
        C[t, :, L:L + GAP] = white
        C[t, :, L + GAP:] = right
    C = C.reshape(T, L * cols, 4)                       # (T, N, 4) colour rollout

    captions = [f"t={t + WARMUP}    left: CCW (+1) = identity rule    |"
                f"    right: CW (-1) = reversal rule" for t in range(T)]

    dummy = np.zeros((T, L * cols), np.int64)           # colours come from colors_rollout
    out = os.path.join(FIGDIR, "e7_spiral_rule.gif")
    animate(dummy, e7.ACT, tau, layout="grid", grid_shape=(L, cols), out=out,
            fps=12, stride=3, dpi=90, colors_rollout=C, figsize=(9.5, 3.6),
            title="E7: rotation direction is the rule (a 2D spiral option)",
            captions=captions)
    print(f"frames={T} (warmup {WARMUP} dropped); grid {L}x{cols}")
    print("wrote", out)


if __name__ == "__main__":
    main()
