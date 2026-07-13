"""C6 mechanism animation -- the persistent core is necessary, not sufficient.

Renders C6 Result B (docs/c6_results.md) as a narrated GIF: two spirals nucleated
identically (same true chirality, CCW = rule 0 = identity), differing only in the
lattice threshold theta -- the do(theta) ablation. Left is INTACT (theta=4.0, the
E0 organised band): the core persists and the rule stays readable for the whole
block. Right is ABLATED (theta=5.0): the core dies almost immediately after
nucleation and the rule becomes unreadable for the rest of the block -- exactly
why C6 finds switching collapses (0.85 -> 0.52) while the re-nucleated single-rule
control is spared (0.90 vs 0.89): the ablation removes persistence, not the
spiral's ability to encode the rule at all.

This is the 2D-spiral analogue of the E2 (ring memory) and E5 (ring persistence)
animations: same true content, only the *timescale/threshold parameter* differs,
and that parameter alone determines whether the context is available later. Uses
the real C6/E7 machinery (`e7_learning.make_spiral`, `e7_spiral_option.seed_spiral`
/ `local_winding`) -- the readouts are genuine substrate computation.

Output
------
docs/figures/c6_necessity.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
from ghca_net_viz import state_colors, animate
import e7_spiral_option as sp
import e7_learning as e7

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

L = sp.L
GAP = 3
CHIR = +1                # true chirality: CCW = rule 0 = identity
T_RUN = 110               # roughly one block (block_len=25 trials x 4 steps/trial)
C = L // 2


def rollout(ablate, seed=0):
    net = e7.make_spiral(seed, ablate=ablate)
    rng = np.random.default_rng(seed)
    sp.seed_spiral(net, CHIR, jitter=0.02, rng=rng)
    phi = np.zeros((T_RUN, L * L), np.int64)
    rule_ok = []
    for t in range(T_RUN):
        net.step(None)
        phi[t] = net.phi
        w = sp.local_winding(net.phi, C, C)
        g = 1 if w > 0.5 else (-1 if w < -0.5 else 0)
        rule_ok.append(bool(g == CHIR))
    return phi, rule_ok, net.tau


def main():
    intact, intact_ok, tau = rollout(ablate=False)
    ablated, ablated_ok, _ = rollout(ablate=True)

    cols = 2 * L + GAP
    white = np.array([1.0, 1.0, 1.0, 1.0])
    Croll = np.zeros((T_RUN, L, cols, 4))
    for t in range(T_RUN):
        left = state_colors(intact[t], sp.ACT, tau).reshape(L, L, 4)
        right = state_colors(ablated[t], sp.ACT, tau).reshape(L, L, 4)
        Croll[t, :, :L] = left
        Croll[t, :, L:L + GAP] = white
        Croll[t, :, L + GAP:] = right
    Croll = Croll.reshape(T_RUN, L * cols, 4)

    captions = []
    for t in range(T_RUN):
        li = "rule HELD" if intact_ok[t] else "rule lost"
        la = "rule HELD" if ablated_ok[t] else "rule LOST"
        captions.append(f"t={t}   left intact: {li}   |   right do(theta) ablated: {la}")

    dummy = np.zeros((T_RUN, L * cols), np.int64)
    out = os.path.join(FIGDIR, "c6_necessity.gif")
    animate(dummy, sp.ACT, tau, layout="grid", grid_shape=(L, cols), out=out,
            fps=9, stride=2, dpi=85, colors_rollout=Croll, figsize=(7.6, 4.2),
            title="C6: do(theta) sets whether the core -- and the rule -- persists",
            captions=captions)
    death_t = next((i for i in range(len(ablated_ok)) if not any(ablated_ok[i:i + 8])), None)
    print(f"intact holds the rule for all {T_RUN} steps: {all(intact_ok)}")
    print(f"ablated rule lost by step: {death_t}")
    print("wrote", out)


if __name__ == "__main__":
    main()
