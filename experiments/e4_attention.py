"""
E4 - Selective attention as biased winner-take-all by wave annihilation.

Two stimulus streams are presented simultaneously (conflict); each would drive a
different action. A top-down bias favours one stream. The claim (design doc E4):
attention here is biased winner-take-all among competing waves, and it needs NO
added inhibitory machinery -- the substrate's own refractory annihilation does
the competition.

Mechanism: a shared 1D excitable chain is the arena. Stream 0 ignites the left
end, stream 1 the right end; the two waves travel inward, collide, and annihilate
in each other's refractory wake. Whichever wave captures the CENTRE node wins
(its action fires). A top-down bias (a small timing advantage to the attended
stream) shifts the collision locus, hence the winner. Against sensory noise
(trial-to-trial jitter in the ignition times) the winner is probabilistic near
zero bias -> a psychometric curve whose width is set by the noise, exactly the
signal-detection picture of attention as a bias on a noisy competition.

We show: (A) the psychometric curve P(attended wins) vs bias, with accuracy -> 1
for modest bias; (B) the collision (annihilation) locus is linear in the bias --
the competition is resolved spatially by where the waves annihilate; (C) a
space-time raster of one trial showing the two waves meeting and annihilating.
No node in the substrate is inhibitory (only excitatory weights + refractoriness).

Outputs
-------
docs/figures/e4_attention.png
result/e4/e4_data.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e4")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

L = 41
CENTER = L // 2
ACT, PAS = 1, 3
T_RUN = 60


def chain(L, w=1.0):
    W = np.zeros((L, L))
    for i in range(L):
        if i > 0:
            W[i, i - 1] = w
        if i < L - 1:
            W[i, i + 1] = w
    return W


def trial(bias, jitter=0.0, seed=0, record=False):
    """bias>0 gives stream 0 a timing advantage (stream 1 ignites later).
    `jitter` = std of Gaussian sensory noise on each stream's ignition time.
    Propagation is deterministic (p_s=0); the only noise is the ignition times.
    Returns (winner, annihilation_locus, raster or None).

    Wave identity is tracked by label propagation: a firing node inherits the
    label of the active neighbour that triggered it, so we read which stream
    captured the centre."""
    net = Network(chain(L), act=ACT, pas=PAS, theta=1.0, p_s=0.0, seed=seed)
    net.phi[:] = 0
    rng = np.random.default_rng(seed)
    j0, j1 = (rng.standard_normal(2) * jitter).round().astype(int)
    s0 = max(0, -bias) + max(0, int(j0))
    s1 = max(0, bias) + max(0, int(j1))
    label = np.full(L, -1)                # -1 none, 0 stream0, 1 stream1
    first = np.full(L, -1)
    raster = np.zeros((T_RUN, L), np.int8) if record else None
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
                nb = [j for j in (i - 1, i + 1) if 0 <= j < L and prev_active[j] and label[j] >= 0]
                if len(nb) == 1:
                    label[i] = label[nb[0]]
                elif len(nb) == 2:                    # collision: symmetric tie-break
                    label[i] = label[nb[int(rng.integers(2))]]
            first[i] = t
        if record:
            raster[t] = net.phi
    winner = int(label[CENTER]) if label[CENTER] >= 0 else -1
    xc = int(np.argmax(first))            # last node to fire = annihilation locus
    return winner, xc, raster


def main():
    biases = np.arange(-10, 11, 1)
    jitters = [1.0, 2.5, 5.0]              # sensory-noise levels (ignition jitter std)
    n_trials = 200

    # A. psychometric: P(stream 0 wins) vs bias, per sensory-noise level
    psy = {j: np.zeros(len(biases)) for j in jitters}
    xc_mean = np.zeros(len(biases))
    for j in jitters:
        for i, b in enumerate(biases):
            wins = np.array([trial(b, jitter=j, seed=1000 * i + k)[0] for k in range(n_trials)])
            valid = wins >= 0
            psy[j][i] = np.mean(wins[valid] == 0) if valid.any() else np.nan
            if j == jitters[1]:
                xc_mean[i] = np.mean([trial(b, jitter=j, seed=1000 * i + k)[1]
                                      for k in range(n_trials)])

    i4p, i4m = np.where(biases == 4)[0][0], np.where(biases == -4)[0][0]
    acc4 = np.nanmean([psy[jitters[1]][i4p], 1 - psy[jitters[1]][i4m]])
    sl, ic = np.polyfit(biases, xc_mean, 1)
    print(f"attention accuracy at |bias|=4 (jitter={jitters[1]}): {acc4:.2f}")
    print(f"annihilation locus vs bias: slope={sl:.2f} (expect ~0.5), "
          f"intercept={ic:.1f} (expect ~centre={CENTER})")

    # C. one clean raster (bias favouring stream 0, no jitter) for the mechanism
    _, _, raster = trial(6, jitter=0.0, seed=0, record=True)

    fig, axes = plt.subplots(1, 3, figsize=(16, 4.6))
    for j in jitters:
        axes[0].plot(biases, psy[j], "o-", label=f"sensory noise σ={j}")
    axes[0].axhline(0.5, ls=":", color="k", alpha=0.4)
    axes[0].axvline(0, ls=":", color="k", alpha=0.4)
    axes[0].set_xlabel("top-down bias  (>0 favours stream 0)")
    axes[0].set_ylabel("P(stream 0 wins)")
    axes[0].set_title("A. Psychometric: bias selects the winner\n(width set by sensory noise)")
    axes[0].legend(fontsize=8)

    axes[1].plot(biases, xc_mean, "o-", color="purple")
    axes[1].axhline(CENTER, ls=":", color="k", alpha=0.4, label=f"centre={CENTER}")
    axes[1].set_xlabel("top-down bias")
    axes[1].set_ylabel("collision (annihilation) locus")
    axes[1].set_title("B. The annihilation locus tracks the bias\n(competition resolved spatially)")
    axes[1].legend(fontsize=8)

    im = axes[2].imshow((raster > 0).T, aspect="auto", cmap="hot", origin="lower")
    axes[2].set_xlabel("time step")
    axes[2].set_ylabel("chain position")
    axes[2].axhline(CENTER, ls=":", color="cyan", alpha=0.7)
    axes[2].set_title("C. Two waves collide & annihilate\n(bias +6: stream 0 captures centre)")

    fig.suptitle("E4: selective attention = biased winner-take-all by wave "
                 "annihilation (no inhibitory machinery)")
    fig.tight_layout()
    p_ = os.path.join(FIGDIR, "e4_attention.png")
    fig.savefig(p_, dpi=110)
    print("wrote", p_)

    np.savez(os.path.join(DATADIR, "e4_data.npz"),
             biases=biases, jitters=np.array(jitters),
             psy=np.array([psy[j] for j in jitters]),
             xc_mean=xc_mean, acc4=acc4, slope=sl, intercept=ic)
    print("wrote", os.path.join(DATADIR, "e4_data.npz"))


if __name__ == "__main__":
    main()
