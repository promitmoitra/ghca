"""
E2 (information addendum) - working memory as a low information-destruction rate.

Casts E2's delayed-response result in the information-theoretic language of
discrete-diffusion learning (Casado Noguerales et al., 2026): the irreducible
cost of a noising process is the rate at which it destroys information about the
clean data, -d/dt I(Z0; Zt). Here the "clean data" is the stimulus identity X and
the "noised state" is the network state at delay D. A reentrant loop is then a
noising process whose information-destruction rate is tuned by the local
timescale tau: tau < L sustains the loop (near-zero destruction = memory), tau
>= L kills it (fast destruction = forgetting).

We measure I(X ; readout of the ring state) in bits vs delay D, for several fixed
tau, with no learning (mechanism), reusing the two-ring substrate of E2. The
readout is which ring carries more activity over a full loop period (phase-
insensitive); I(X;Y) is computed exactly from the 2x(<=3) confusion counts.

Connections: this is the temporal-information-loss companion to C4's
coarse-graining information loss (I(B;W) vs I(B;S)); together they place E2 and
the C-series in one information currency.

Outputs
-------
docs/figures/e2_information.png
result/e2/e2_information.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, ROOT)
sys.path.insert(0, os.path.join(ROOT, "experiments"))
from ghca_net import Network
import e2_delayed_response as e2

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e2")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

L = e2.L                 # ring length (loop transit time)
ACT = e2.ACT
TAUS = [16, 20, 22, 24, 28]
DELAYS = [0, 5, 10, 20, 40, 80, 160]
N_TRIALS = 80


def mutual_information_bits(X, Y):
    """Exact I(X;Y) in bits from samples of two discrete variables."""
    xs = np.unique(X); ys = np.unique(Y)
    n = len(X)
    I = 0.0
    for xv in xs:
        px = np.mean(X == xv)
        for yv in ys:
            py = np.mean(Y == yv)
            pxy = np.mean((X == xv) & (Y == yv))
            if pxy > 0:
                I += pxy * np.log2(pxy / (px * py))
    return I


def retained_info(tau, delay, seed0=0):
    """I(stimulus X ; ring-readout Y) at a given delay, for fixed uniform tau."""
    Xs, Ys = [], []
    for tr in range(N_TRIALS):
        Wt, _, roles, theta, _ = e2.two_ring(seed=tr + seed0)
        net = Network(Wt, act=ACT, pas=int(tau - ACT), theta=theta, p_s=0.0, seed=tr)
        x = tr % 2
        net.phi[:] = 0
        for _ in range(e2.CUE):                       # ignite ring x
            d = np.zeros(net.N, bool); d[roles["sensory"][x]] = True
            net.step(d)
        for _ in range(delay):                        # stimulus-free delay (p_s=0)
            net.step(None)
        ring = [roles["hidden"][:L], roles["hidden"][L:2 * L]]
        sc = np.zeros(2)
        for _ in range(L):                            # integrate over a loop period
            net.step(None)
            am = net.active_mask()
            sc += [am[ring[0]].sum(), am[ring[1]].sum()]
        y = int(np.argmax(sc)) if sc.sum() > 0 else 2  # 2 = "no response"
        Xs.append(x); Ys.append(y)
    return mutual_information_bits(np.array(Xs), np.array(Ys))


def main():
    I = np.zeros((len(TAUS), len(DELAYS)))
    for i, tau in enumerate(TAUS):
        for j, d in enumerate(DELAYS):
            I[i, j] = retained_info(tau, d)
        print(f"tau={tau:2d} (L={L}): I(X;state) vs delay = "
              + " ".join(f"{v:.2f}" for v in I[i]))

    # information-destruction summary: retained info at a long delay, and a
    # memory half-life (first delay where I drops below 0.5 bit)
    D = np.array(DELAYS)
    long_idx = np.argmin(np.abs(D - 80))
    retained_long = I[:, long_idx]
    half_life = []
    for i in range(len(TAUS)):
        below = np.where(I[i] < 0.5)[0]
        half_life.append(D[below[0]] if len(below) else D[-1] + 1)

    print("\nretained info at D=80 (bits):",
          {t: round(r, 2) for t, r in zip(TAUS, retained_long)})
    print("memory half-life (delay to <0.5 bit):",
          {t: int(h) for t, h in zip(TAUS, half_life)})

    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    cmap = plt.cm.viridis(np.linspace(0, 1, len(TAUS)))
    for i, tau in enumerate(TAUS):
        axes[0].plot(DELAYS, I[i], "o-", color=cmap[i],
                     label=f"tau={tau}" + (" (<L)" if tau < L else " (>=L)"))
    axes[0].axhline(1.0, ls=":", color="k", alpha=0.4)
    axes[0].set_xlabel("delay D (steps)")
    axes[0].set_ylabel("I(stimulus ; network state)  [bits]")
    axes[0].set_title("Information retained about the stimulus vs delay\n"
                      "memory = near-zero information-destruction rate")
    axes[0].legend(fontsize=8)

    axes[1].plot(TAUS, retained_long, "o-", color="crimson")
    axes[1].axvline(L, ls="--", color="k", alpha=0.5, label=f"tau = L = {L}")
    axes[1].set_xlabel("local timescale tau")
    axes[1].set_ylabel("info retained at D=80 (bits)")
    axes[1].set_title("tau controls the information-destruction rate\n"
                      "(sharp memory transition at tau = L)")
    axes[1].legend()
    fig.suptitle("E2 information addendum: the reentrant loop is a noising "
                 "process whose forgetting rate is set by tau")
    fig.tight_layout()
    p = os.path.join(FIGDIR, "e2_information.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e2_information.npz"),
             taus=np.array(TAUS), delays=D, I=I,
             retained_long=retained_long, half_life=np.array(half_life), L=L)
    print("wrote", os.path.join(DATADIR, "e2_information.npz"))


if __name__ == "__main__":
    main()
