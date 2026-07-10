"""
E0 - Substrate characterisation (no learning).

Locates the self-sustaining "spiral band" of the Greenberg-Hastings network
substrate and checks that the dominant loop period tracks the local timescale
tau (validating Line B's control variable). See docs/learning_experiments.md,
experiment E0.

Outputs
-------
docs/figures/e0_range_death.png     : self-sustained activity vs theta for
                                      several neighbourhood ranges (r=1 dies)
docs/figures/e0_phase_diagram.png   : (theta x p_s) heatmaps of active fraction,
                                      coherence, dominant period, spiral cores
docs/figures/e0_period_vs_tau.png   : dominant period vs tau (band validation)
result/e0/e0_data.npz               : raw swept arrays
docs/e0_results.md                  : written separately after inspecting output

Operating point (found here): r=2, act=6, pas=8 (tau=14). Threshold theta is
the regime knob: theta<=3 turbulent (many cores), theta~4 few organised spiral
cores, theta>=5 death. Range r=1 (von Neumann) fixates and cannot self-sustain.
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network, lattice2d, smallworld, dominant_period, \
    count_phase_singularities

FIGDIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                      "docs", "figures")
DATADIR = os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))),
                       "result", "e0")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)


def characterise(net, T=400, burn=200, grid_L=None):
    """Run one network and return steady-state observables."""
    net.seed_random(frac_active=0.05, frac_refractory=0.10)
    out = net.run(T, record=True)
    A_ss = out["A"][burn:].mean()
    R_ss = out["R"][burn:].mean()
    period = dominant_period(out["A"][burn:])
    cores = np.nan
    if grid_L is not None:
        # average spiral-core count over the last few frames
        cs = [count_phase_singularities(out["phi"][t].reshape(grid_L, grid_L),
                                        net.act, int(net.tau[0]))
              for t in range(T - 5, T)]
        cores = float(np.mean(cs))
    return A_ss, R_ss, period, cores


def range_death(L=40, act=6, pas=8, ranges=(1, 2, 3), thetas=None, seed=7):
    """Self-sustained activity (p_s=0) vs theta for several ranges.

    Demonstrates the E0 discriminator: range-1 (von Neumann) fixates; larger
    ranges have a live self-sustaining band. Returns {r: (thetas, A_late)}.
    """
    if thetas is None:
        thetas = np.arange(1, 9)
    res = {}
    for r in ranges:
        W = lattice2d(L, r=r, periodic=True)
        A_late = np.zeros(len(thetas))
        for k, th in enumerate(thetas):
            net = Network(W, act=act, pas=pas, theta=float(th), p_s=0.0, seed=seed)
            net.seed_random(0.30, 0.30)
            out = net.run(600, record=False)
            A_late[k] = out["A"][400:].mean()
        res[r] = (thetas, A_late)
        print(f"  range r={r} done")
    return res


def phase_diagram(L=40, act=6, pas=8, r=2, thetas=None, p_ss=None, seed=1):
    """Sweep (theta, p_s) on a torus; return observable grids."""
    if thetas is None:
        thetas = np.arange(1, 9).astype(float)
    if p_ss is None:
        p_ss = np.logspace(-4, -1, 8)
    W = lattice2d(L, r=r, periodic=True)
    A = np.zeros((len(thetas), len(p_ss)))
    R = np.zeros_like(A)
    P = np.zeros_like(A)
    C = np.zeros_like(A)
    for i, th in enumerate(thetas):
        for j, ps in enumerate(p_ss):
            net = Network(W, act=act, pas=pas, theta=th, p_s=ps, seed=seed)
            A[i, j], R[i, j], P[i, j], C[i, j] = characterise(net, grid_L=L)
        print(f"  theta={th:.1f} done")
    return thetas, p_ss, A, R, P, C


def period_vs_tau(L=40, act=2, taus=None, theta=1.0, p_s=0.0, r=2, seed=2):
    """Sweep tau (via pas) in the self-sustaining turbulent regime and measure
    the dominant global period. In this regime the local cycle dominates the
    global oscillation, so the period should track tau (validating Line B's
    control variable)."""
    if taus is None:
        taus = np.arange(4, 22, 2)
    W = lattice2d(L, r=r, periodic=True)
    periods = np.zeros(len(taus))
    for k, tau in enumerate(taus):
        pas = int(tau) - act
        net = Network(W, act=act, pas=pas, theta=theta, p_s=p_s, seed=seed)
        net.seed_random(0.30, 0.30)
        out = net.run(600, record=False)
        periods[k] = dominant_period(out["A"][350:])
    return taus, periods


def plot_phase_diagram(thetas, p_ss, A, R, P, C):
    fig, axes = plt.subplots(2, 2, figsize=(11, 9))
    ext = [np.log10(p_ss[0]), np.log10(p_ss[-1]), thetas[0], thetas[-1]]
    fields = [(A, "active fraction  A_ss", "viridis"),
              (R, "phase coherence  R_ss", "magma"),
              (np.where(np.isfinite(P), P, np.nan), "dominant period", "cividis"),
              (C, "spiral-core count", "inferno")]
    for ax, (F, title, cm) in zip(axes.ravel(), fields):
        im = ax.imshow(F, origin="lower", aspect="auto", extent=ext, cmap=cm)
        ax.set_xlabel(r"$\log_{10} p_s$")
        ax.set_ylabel(r"threshold $\theta$")
        ax.set_title(title)
        fig.colorbar(im, ax=ax)
    fig.suptitle("E0 phase diagram - GH torus (L=40, r=2, a=6, tau=14)", y=0.99)
    fig.tight_layout()
    path = os.path.join(FIGDIR, "e0_phase_diagram.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


def plot_range_death(res):
    fig, ax = plt.subplots(figsize=(6, 5))
    for r, (thetas, A_late) in sorted(res.items()):
        ax.plot(thetas, A_late, "o-", label=f"range r={r}")
    ax.set_xlabel(r"threshold $\theta$")
    ax.set_ylabel("self-sustained active fraction (p_s=0)")
    ax.set_title("E0 - range-1 fixates; larger ranges self-sustain")
    ax.legend()
    fig.tight_layout()
    path = os.path.join(FIGDIR, "e0_range_death.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


def plot_period_vs_tau(taus, periods):
    fig, ax = plt.subplots(figsize=(6, 5))
    finite = np.isfinite(periods)
    ax.plot(taus, taus, "k--", alpha=0.5, label=r"period $=\tau$")
    ax.plot(taus[finite], periods[finite], "o-", color="crimson",
            label="measured dominant period")
    ax.set_xlabel(r"local timescale $\tau = a + p$")
    ax.set_ylabel("dominant period (steps)")
    ax.set_title("E0 - dominant loop period tracks tau")
    ax.legend()
    fig.tight_layout()
    path = os.path.join(FIGDIR, "e0_period_vs_tau.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


def band_summary(thetas, p_ss, A, R, P):
    """Simple rule for the 'useful band': intermediate activity + coherence
    + a finite dominant period. Print the qualifying (theta, p_s) cells."""
    band = (A > 0.03) & (A < 0.35) & (R > 0.35) & np.isfinite(P)
    print("\n-- candidate self-sustaining band cells (theta, p_s) --")
    for i, th in enumerate(thetas):
        for j, ps in enumerate(p_ss):
            if band[i, j]:
                print(f"   theta={th:.1f}  p_s={ps:.1e}  "
                      f"A={A[i,j]:.3f} R={R[i,j]:.3f} period={P[i,j]:.1f}")
    return band


def main():
    print("E0: range-death test (self-sustaining vs range) ...")
    res = range_death()
    print("E0: sweeping (theta, p_s) phase diagram on torus ...")
    thetas, p_ss, A, R, P, C = phase_diagram()
    print("E0: sweeping period vs tau ...")
    taus, periods = period_vs_tau()

    plot_range_death(res)
    plot_phase_diagram(thetas, p_ss, A, R, P, C)
    plot_period_vs_tau(taus, periods)
    band = band_summary(thetas, p_ss, A, R, P)

    np.savez(os.path.join(DATADIR, "e0_data.npz"),
             thetas=thetas, p_ss=p_ss, A=A, R=R, P=P, C=C,
             taus=taus, periods=periods, band=band,
             **{f"range_{r}_A": v[1] for r, v in res.items()},
             range_thetas=next(iter(res.values()))[0])
    print("wrote", os.path.join(DATADIR, "e0_data.npz"))


if __name__ == "__main__":
    main()
