"""E5 mechanism animation -- the slow loop as an option holding the rule.

Renders the core of the E5 executive-control mechanism (docs/e5_results.md): each
of the two rules owns a slow context ring. At a block start the active rule's
ring is ignited and *persists* (tau < L, the E2 mechanism), supplying a standing
context signal for the whole block; the other ring stays dark. At the next block
the roles flip. So the "rule in force" is literally which loop is rotating -- a
slow reentrant loop acting as an option that gates the fast routing.

Uses the real E5 substrate (`e5_executive.make` / `ignite_block`); the rings
rotate autonomously once ignited (no learning needed to show the gating), so this
is the genuine substrate behaviour.

Output
------
docs/figures/e5_options.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
from ghca_net_viz import animate
import e5_executive as e5

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

N_BLOCKS = 4
STEPS_PER_BLOCK = 26


def main():
    net, roles = e5.make(0, ablate=False)
    L = e5.L_RING
    r0, r1 = roles["rings"][0], roles["rings"][1]
    phi0, phi1, rule_at = [], [], []
    for b in range(N_BLOCKS):
        rule = b % 2
        e5.ignite_block(net, roles, rule)
        for _ in range(STEPS_PER_BLOCK):
            net.step(None)
            phi0.append(net.phi[r0].copy())
            phi1.append(net.phi[r1].copy())
            rule_at.append(rule)
    phi0 = np.array(phi0); phi1 = np.array(phi1)
    T = len(rule_at)

    # two rings side by side ('free' layout)
    ang = 2 * np.pi * np.arange(L) / L
    unit = np.column_stack([np.cos(ang), np.sin(ang)])
    pos = np.vstack([unit + [-1.7, 0.0], unit + [1.7, 0.0]])
    rollout = np.concatenate([phi0, phi1], axis=1)                 # (T, 2L)
    tau_vec = np.full(2 * L, e5.TAU_SLOW)

    captions = []
    for t in range(T):
        b = t // STEPS_PER_BLOCK
        rule = rule_at[t]
        holder = "LEFT ring" if rule == 0 else "RIGHT ring"
        captions.append(f"E5  block {b}   rule {rule} in force   ->   {holder} "
                        f"holds it (the other loop is extinguished)")

    annotations = [(-1.7, 1.5, "Rule 0 loop"), (1.7, 1.5, "Rule 1 loop")]

    out = os.path.join(FIGDIR, "e5_options.gif")
    animate(rollout, e5.ACT, tau_vec, layout="free", pos=pos, out=out,
            fps=12, stride=2, marker_size=120, figsize=(9, 5.2), dpi=85,
            title="E5: the rule in force = which slow loop is rotating (an option)",
            captions=captions, annotations=annotations)
    print(f"blocks={N_BLOCKS}, steps/block={STEPS_PER_BLOCK}, frames={T}")
    print("wrote", out)


if __name__ == "__main__":
    main()
