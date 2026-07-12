"""E4 mechanism animation -- selective attention as wave annihilation.

Renders the E4 mechanism (docs/e4_results.md) as a narrated GIF: two waves
ignite at opposite ends of an excitable chain, travel inward, collide, and
annihilate in each other's refractory wake. A top-down *bias* gives the attended
stream a head start, so it captures the centre node -- attention as a biased
winner-take-all, with no inhibitory machinery (only excitation + refractoriness).

This is the animated counterpart to panel C of e4_attention.py (the static
space-time raster). Two-colour tinting tracks which stream each node belongs to,
so you watch the blue and orange fronts meet and cancel.

Output
------
docs/figures/e4_annihilation.gif
"""

import os
import sys
import numpy as np

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network
from ghca_net_viz import animate
from e4_attention import chain, L, CENTER, ACT, PAS, T_RUN

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
os.makedirs(FIGDIR, exist_ok=True)

# two-stream palette (colour-blind-friendly blue / orange), + neutral red for
# an as-yet-unlabelled active node and light grey for rested.
_BLUE = np.array([0.13, 0.40, 0.72, 1.0])       # stream 0 (left, attended)
_ORANGE = np.array([0.90, 0.52, 0.10, 1.0])     # stream 1 (right)
_RED = np.array([0.83, 0.09, 0.11, 1.0])        # active, unlabelled
_GREY = np.array([0.90, 0.90, 0.90, 1.0])       # rested
_STREAM_RGB = {0: _BLUE[:3], 1: _ORANGE[:3]}


def record_trial(bias, seed=0):
    """Re-run one deterministic E4 trial, recording phi and the wave label of
    every node at every step (label propagates from the igniting stream)."""
    net = Network(chain(L), act=ACT, pas=PAS, theta=1.0, p_s=0.0, seed=seed)
    net.phi[:] = 0
    s0, s1 = max(0, -bias), max(0, bias)         # no jitter: clean mechanism shot
    label = np.full(L, -1)
    first = np.full(L, -1)
    phi_roll = np.zeros((T_RUN, L), np.int64)
    lab_roll = np.zeros((T_RUN, L), np.int64)
    collision_t = None
    for t in range(T_RUN):
        prev_active = net.active_mask()
        drive = np.zeros(L, bool)
        if t == s0:
            drive[0] = True
        if t == s1:
            drive[L - 1] = True
        net.step(drive if drive.any() else None)
        am = net.active_mask()
        for i in np.where(am & (first < 0))[0]:
            if i == 0 and t == s0:
                label[i] = 0
            elif i == L - 1 and t == s1:
                label[i] = 1
            else:
                nb = [j for j in (i - 1, i + 1)
                      if 0 <= j < L and prev_active[j] and label[j] >= 0]
                if len(nb) == 1:
                    label[i] = label[nb[0]]
                elif len(nb) == 2 and collision_t is None:
                    collision_t = t
            first[i] = t
        phi_roll[t] = net.phi
        lab_roll[t] = label
    winner = int(label[CENTER]) if label[CENTER] >= 0 else -1
    return phi_roll, lab_roll, winner, collision_t


def colors_from_labels(phi_roll, lab_roll):
    """Build a (T, N, 4) colour rollout: active nodes tinted by stream label,
    refractory nodes a darkened tint of the stream that passed, rested grey."""
    T, N = phi_roll.shape
    C = np.tile(_GREY, (T, N, 1))
    for t in range(T):
        phi = phi_roll[t]
        lab = lab_roll[t]
        active = (phi >= 1) & (phi <= ACT)
        refr = phi > ACT
        for i in np.where(active)[0]:
            C[t, i, :3] = _STREAM_RGB.get(lab[i], _RED[:3]); C[t, i, 3] = 1.0
        for i in np.where(refr)[0]:
            base = _STREAM_RGB.get(lab[i], _RED[:3])
            frac = (phi[i] - ACT) / max(net_tau - ACT, 1)     # 0 head -> 1 tail
            C[t, i, :3] = (1 - frac) * base + frac * np.array([0.13, 0.13, 0.17])
            C[t, i, 3] = 1.0
    return C


def main():
    bias = 6                                     # left stream attended (head start)
    phi_roll, lab_roll, winner, coll = record_trial(bias, seed=0)
    global net_tau
    net_tau = ACT + PAS
    C = colors_from_labels(phi_roll, lab_roll)

    # trim to the informative window (a few steps past annihilation)
    end = min(T_RUN, (coll or T_RUN // 2) + 10)
    phi_roll, C = phi_roll[:end], C[:end]

    s0, s1 = max(0, -bias), max(0, bias)
    win_name = {0: "LEFT (attended)", 1: "right", -1: "tie"}[winner]
    captions = []
    for t in range(end):
        if t < max(s0, s1):
            msg = "two streams about to ignite the chain ends"
        elif coll is not None and t < coll:
            msg = "waves travel inward (blue=attended left, orange=right)"
        elif coll is not None and t <= coll + 2:
            msg = "collision: waves annihilate in each other's refractory wake"
        else:
            msg = f"annihilation left of centre -> {win_name} stream captures centre"
        captions.append(f"E4  bias=+{bias}   t={t}   {msg}")

    out = os.path.join(FIGDIR, "e4_annihilation.gif")
    animate(phi_roll, ACT, net_tau, layout="line", out=out, fps=5, stride=1,
            title="E4: attention = biased winner-take-all by wave annihilation",
            colors_rollout=C, captions=captions, marker_size=300)
    print(f"winner={win_name}  collision_t={coll}  frames={end}")
    print("wrote", out)


if __name__ == "__main__":
    main()
