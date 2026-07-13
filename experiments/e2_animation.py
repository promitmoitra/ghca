"""E2 mechanism animation -- working memory as a tau-tuned reentrant loop.

Renders the E2 mechanism (docs/e2_results.md) as a narrated GIF: two identical
stimulus-specific directed rings are ignited the same way, but with different
local timescales tau. When tau < L (loop transit time) the rotating pulse
sustains indefinitely -- the stimulus is *held* across the delay (memory). When
tau >= L the pulse dies within ~L steps -- the memory is lost. tau is the
memory-duration knob; this is the causal lever behind the E2 learning result
(Line B tunes tau below L to hold memory).

Both rings use the *real* E2 substrate (`e2_delayed_response.two_ring`) and the
real ignition (drive the first IGNITE nodes for CUE steps), so the death at
tau>=L is the genuine substrate behaviour, not a hand-tuned cartoon.

Output
------
docs/figures/e2_ring_memory.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), "experiments"))
from ghca_net import Network
from ghca_net_viz import animate
import e2_delayed_response as e2

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

TAU_HOLD = 20      # < L = 24 : loop sustains  (memory held)
TAU_DIE = 28       # >= L     : loop dies ~L steps after the cue (memory lost)
T_POST = 130       # steps after the cue to animate


def ring_rollout(tau, x=0, seed=0):
    """Ignite ring x on the real E2 substrate and record ring x's phases."""
    W, plastic, roles, theta, _ = e2.two_ring(seed=seed)
    net = Network(W, act=e2.ACT, pas=int(tau - e2.ACT), theta=theta, p_s=0.0, seed=seed)
    net.phi[:] = 0
    for _ in range(e2.CUE):                              # real ignition
        d = np.zeros(net.N, bool); d[roles["sensory"][x]] = True
        net.step(d)
    ring_nodes = roles["hidden"][x * e2.L:(x + 1) * e2.L]
    roll = np.zeros((T_POST, e2.L), np.int64)
    for t in range(T_POST):
        net.step(None)
        roll[t] = net.phi[ring_nodes]
    return roll


def main():
    hold = ring_rollout(TAU_HOLD)
    die = ring_rollout(TAU_DIE)
    L = e2.L

    # layout: two rings side by side on the 'free' layout (2L nodes)
    ang = 2 * np.pi * np.arange(L) / L
    unit = np.column_stack([np.cos(ang), np.sin(ang)])
    left = unit + np.array([-1.7, 0.0])
    right = unit + np.array([1.7, 0.0])
    pos = np.vstack([left, right])
    rollout = np.concatenate([hold, die], axis=1)        # (T, 2L)
    tau_vec = np.concatenate([np.full(L, TAU_HOLD), np.full(L, TAU_DIE)])

    # narration: report the delay and each ring's status (alive => still holding)
    captions = []
    die_alive = (die >= 1).any(axis=1)
    death_t = int(np.argmax(~die_alive)) if (~die_alive).any() else T_POST
    for t in range(T_POST):
        r_state = "holding" if die_alive[t] else "FORGOTTEN"
        msg = f"E2  delay D={t}   left tau={TAU_HOLD}<L: holding   |   right tau={TAU_DIE}≥L: {r_state}"
        captions.append(msg)

    annotations = [(-1.7, 1.45, f"tau={TAU_HOLD} < L  (memory holds)"),
                   (1.7, 1.45, f"tau={TAU_DIE} ≥ L  (memory decays)")]

    out = os.path.join(FIGDIR, "e2_ring_memory.gif")
    animate(rollout, e2.ACT, tau_vec, layout="free", pos=pos, out=out,
            fps=12, stride=2, marker_size=90, figsize=(9, 5.2), dpi=85,
            title="E2: working memory = a tau-tuned reentrant loop (ring L=24)",
            captions=captions, annotations=annotations)
    print(f"right ring (tau={TAU_DIE}) dies at delay D={death_t}; "
          f"left ring (tau={TAU_HOLD}) sustains to D={T_POST}")
    print("wrote", out)


if __name__ == "__main__":
    main()
