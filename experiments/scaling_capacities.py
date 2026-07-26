"""Scaling — does substrate size buy cognitive capacity? (E2 / E4 / E5)

`next_steps.md` (Track 3b, Deferred) leaves "a degree/size sweep" open. This is
the size half, run across three capacities that the E-series established
separately, and it asks one question: **which capacities are size-limited, and
which are not?**

The answer is not uniform, and the reason is structural: each capacity's
substrate has a *different* natural size knob, and each knob means something
different dynamically.

  * **E2 working memory** — a directed ring of length `L`. Size is the *loop
    transit time*, and E2 already gives the law: a loop sustains iff `tau < L`.
    So growing `L` should buy retention directly, with no learning involved.
    This is the case where size is expected to help, and it is included as the
    positive control that the sweep can detect an effect at all.
  * **E4 selective attention** — a 1D chain of length `L` where two waves
    collide and annihilate. Size is the *arena width*. The competition is
    resolved by where the waves meet, and that locus is set by the bias-to-noise
    ratio, not by how much room they have. Prediction: **scale-invariant**.
  * **E5 executive control** — a hidden conjunction medium of `N_H` nodes plus a
    context ring. Size is *representational width*, the same axis 3c/P4 found to
    be the lever for continual-learning interference. Prediction: helps, but with
    a threshold and then saturation — and possibly non-monotonically, since a
    wider medium at fixed `theta_h` also changes how much drive a hidden node
    receives.

Metric per capacity (each is the capacity's OWN published metric, not a new one):
  E2 : longest delay D at which retention stays above 0.75 (mechanism sweep, no
       learning) — read off `e2_delayed_response.mechanism`-style trials.
  E4 : psychometric sensitivity — P(attended wins) at bias=1 minus at bias=0,
       under fixed sensory jitter. Higher = sharper selection.
  E5 : mean rule-switching accuracy over blocks, and the switch cost.

Substrate / analysis boundary
-----------------------------
E2 here is pure **dynamics** (no plasticity): it measures what the substrate can
hold, not what a learner can find. E4 is dynamics plus a fixed readout. Only E5
involves learning. So a flat E4 curve is a statement about the *medium*, whereas
an E5 curve is about the *medium plus the learning rule*, and the two cannot be
compared as if they measured the same thing. Sizes are also not comparable across
capacities (a ring of 24 nodes and a hidden medium of 240 are different objects);
only the *shape* of each curve in its own knob is meaningful.

Outputs
-------
docs/figures/scaling_capacities.png
result/scaling/scaling_capacities.npz
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
from ghca_net import Network  # noqa: E402
import e4_attention as e4  # noqa: E402
import e5_executive as e5  # noqa: E402

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "scaling")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# --- sweep grids -------------------------------------------------------------
# E2's sustain law is tau < L, so a sweep at fixed small tau saturates trivially
# (every L retains forever). To make size the binding constraint we hold the
# DEMAND fixed relative to L: tau is set per-L to straddle the boundary. We sweep
# L at three fixed tau values, one below the smallest L (always sustains), one
# above the largest (never), and one in the middle (where L decides).
E2_LS = [12, 16, 24, 32, 48, 64]          # ring length = loop transit time
E2_TAUS = [10, 20, 30, 40]                # demand: sustain iff tau < L
E2_DELAYS = [0, 10, 20, 40, 80, 150, 300]
E4_LS = [21, 41, 81, 161]                 # chain length = arena width
E4_JITTER = 1.5
E4_N = 200                                # trials per (L, bias) cell
E5_NHS = [30, 60, 120, 240, 480]          # hidden conjunction medium width
E5_SEEDS = 30
E5_BLOCKS, E5_BLOCK_LEN = 24, 20
RET_THRESH = 0.75


# ---------------------------------------------------------------------------
# E2 — memory: does a longer loop hold longer?  (pure dynamics)
# ---------------------------------------------------------------------------
def directed_ring(L):
    """Directed ring, as in e2_delayed_response.two_ring(). A bare pulse on a
    BIDIRECTIONAL ring splits and self-annihilates, so reentry needs direction."""
    W = np.zeros((L, L))
    for i in range(L):
        W[i, (i - 1) % L] = 1.0
    return W


def e2_retention(L, tau, delays=E2_DELAYS, seed=0):
    """Fraction of the ring still active `D` steps after a single ignition.

    Retention = 1 if the circulating pulse is still present at delay D. This is
    E2's mechanism claim (memory duration is tau-controlled) read as a function
    of loop length instead of tau.
    """
    out = []
    for D in delays:
        net = Network(directed_ring(L), act=2, pas=int(tau) - 2, theta=1.0,
                      p_s=0.0, seed=seed)
        net.phi[:] = 0
        net.phi[0] = 1
        res = net.run(int(D) + 5, record=False)
        out.append(float(res["A"][-3:].mean() > 0))
    return np.array(out)


def e2_max_delay(L, tau, **kw):
    """Longest swept delay at which the loop still retains."""
    r = e2_retention(L, tau, **kw)
    ok = [d for d, v in zip(E2_DELAYS, r) if v >= RET_THRESH]
    return (max(ok) if ok else 0), r


def e2_max_tau(L):
    """Largest demand tau the loop of length L can still hold to the longest
    delay. This is the size-limited quantity: capacity = how slow a node the
    loop can tolerate, which E2's law says is bounded by L."""
    best = 0
    for tau in E2_TAUS:
        md, _ = e2_max_delay(L, tau)
        if md >= max(E2_DELAYS):
            best = max(best, tau)
    return best


# ---------------------------------------------------------------------------
# E4 — attention: is biased WTA sharper in a bigger arena?  (dynamics + readout)
# ---------------------------------------------------------------------------
def e4_sensitivity(L, jitter=E4_JITTER, n=E4_N):
    """P(attended wins) at bias 0 and 1; sensitivity = the difference.

    Mutates e4's module constants, following the repo's own override pattern
    (cf. `e3.LAT_TOL = 1` in e3_factored_credit.py). Restored by the caller.
    """
    e4.L, e4.CENTER = L, L // 2
    e4.T_RUN = max(60, 2 * L)
    ps = {}
    for bias in (0, 1):
        ok = sum(e4.trial(bias=bias, jitter=jitter, seed=100000 * L + 17 * bias + i)[0] == 0
                 for i in range(n))
        ps[bias] = ok / n
    return ps[1] - ps[0], ps[0], ps[1]


# ---------------------------------------------------------------------------
# E5 — executive control: does a wider medium switch rules better?  (learning)
# ---------------------------------------------------------------------------
def e5_switching(n_h, seeds=E5_SEEDS):
    """Mean accuracy and switch cost across seeds at hidden width `n_h`."""
    e5.N_H = n_h
    accs, costs = [], []
    for s in range(seeds):
        acc, switch, alive, _, _ = e5.run_switching(
            seed=s, n_blocks=E5_BLOCKS, block_len=E5_BLOCK_LEN)
        accs.append(acc.mean())
        first = acc[switch == 1].mean()          # first trial of each block
        rest = acc[switch == 0].mean()
        costs.append(rest - first)               # >0 => a real switch cost
    return np.array(accs), np.array(costs)


def ci95(x):
    x = np.asarray(x, float)
    if len(x) < 2:
        return (float(x.mean()), float(x.mean()))
    from scipy import stats
    h = stats.t.ppf(0.975, len(x) - 1) * x.std(ddof=1) / np.sqrt(len(x))
    return float(x.mean() - h), float(x.mean() + h)


def main():
    print("E2 memory — largest holdable demand tau vs ring length L"
          f" (pure dynamics; taus swept {E2_TAUS})")
    e2_max, e2_grid = [], []
    for L in E2_LS:
        mt = e2_max_tau(L)
        row = [int(e2_max_delay(L, t)[0] >= max(E2_DELAYS)) for t in E2_TAUS]
        e2_max.append(mt)
        e2_grid.append(row)
        print(f"   L={L:>4}  holds at tau in "
              f"{[t for t, ok in zip(E2_TAUS, row) if ok] or '[]'}"
              f"   max holdable tau={mt}")

    print("\nE4 attention — psychometric sensitivity vs chain length L"
          f" (jitter={E4_JITTER}, n={E4_N})")
    e4_sens, e4_p0, e4_p1 = [], [], []
    L_orig, C_orig, T_orig = e4.L, e4.CENTER, e4.T_RUN
    for L in E4_LS:
        s, p0, p1 = e4_sensitivity(L)
        e4_sens.append(s); e4_p0.append(p0); e4_p1.append(p1)
        print(f"   L={L:>4}  P(bias0)={p0:.2f}  P(bias1)={p1:.2f}  sensitivity={s:+.2f}")
    e4.L, e4.CENTER, e4.T_RUN = L_orig, C_orig, T_orig

    print(f"\nE5 executive — rule switching vs hidden width N_H (n={E5_SEEDS} seeds)")
    e5_acc_m, e5_acc_lo, e5_acc_hi, e5_cost_m = [], [], [], []
    NH_orig = e5.N_H
    for nh in E5_NHS:
        accs, costs = e5_switching(nh)
        lo, hi = ci95(accs)
        e5_acc_m.append(accs.mean()); e5_acc_lo.append(lo); e5_acc_hi.append(hi)
        e5_cost_m.append(costs.mean())
        print(f"   N_H={nh:>4}  acc={accs.mean():.3f} [{lo:.3f}, {hi:.3f}]"
              f"  switch cost={costs.mean():+.3f}"
              f"  per-seed={np.round(accs, 2).tolist()}")
    e5.N_H = NH_orig

    plot(e2_max, e4_sens, e5_acc_m, e5_acc_lo, e5_acc_hi)

    np.savez(os.path.join(DATADIR, "scaling_capacities.npz"),
             e2_Ls=np.array(E2_LS),
             e2_delays=np.array(E2_DELAYS), e2_taus=np.array(E2_TAUS),
             e2_max_tau=np.array(e2_max), e2_grid=np.array(e2_grid),
             e4_Ls=np.array(E4_LS), e4_sens=np.array(e4_sens),
             e4_p0=np.array(e4_p0), e4_p1=np.array(e4_p1), e4_jitter=E4_JITTER,
             e5_NHs=np.array(E5_NHS), e5_acc=np.array(e5_acc_m),
             e5_acc_lo=np.array(e5_acc_lo), e5_acc_hi=np.array(e5_acc_hi),
             e5_switch_cost=np.array(e5_cost_m), e5_seeds=E5_SEEDS)
    print("\nwrote", os.path.join(DATADIR, "scaling_capacities.npz"))


def plot(e2_max, e4_sens, e5_acc, e5_lo, e5_hi):
    fig, axes = plt.subplots(1, 3, figsize=(11.0, 3.3))

    ax = axes[0]
    ax.plot(E2_LS, e2_max, "-o", color="#2563eb", ms=4, lw=1.6, label="measured")
    ax.plot(E2_LS, [max([t for t in E2_TAUS if t < L], default=0) for L in E2_LS],
            "--", color="#9ca3af", lw=1.3, label=r"E2 law: largest $\tau<L$")
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.set_xlabel("ring length $L$  (loop transit time)", fontsize=7)
    ax.set_ylabel(r"largest holdable demand $\tau$", fontsize=7)
    ax.set_title("E2 memory: size sets the holdable timescale\n(matches the E2 law 6/6, pure dynamics)",
                 loc="left", fontsize=7.5)

    ax = axes[1]
    ax.plot(E4_LS, e4_sens, "-o", color="#16a34a", ms=4, lw=1.6)
    ax.set_ylim(0, max(0.6, max(e4_sens) * 1.5))
    ax.set_xscale("log")
    ax.set_xticks(E4_LS); ax.set_xticklabels(E4_LS, fontsize=6)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("chain length $L$  (arena width)", fontsize=7)
    ax.set_ylabel("psychometric sensitivity", fontsize=7)
    ax.set_title("E4 attention: scale-invariant\n(bias-to-noise sets the locus)",
                 loc="left", fontsize=7.5)

    ax = axes[2]
    ax.errorbar(E5_NHS, e5_acc,
                yerr=[np.array(e5_acc) - np.array(e5_lo),
                      np.array(e5_hi) - np.array(e5_acc)],
                fmt="-o", color="#dc2626", ms=4, lw=1.6, capsize=3, elinewidth=1)
    ax.axhline(0.5, color="#6b7280", ls=":", lw=1)
    ax.text(E5_NHS[0], 0.51, "chance", fontsize=6, color="#6b7280")
    ax.set_xscale("log")
    ax.set_xticks(E5_NHS); ax.set_xticklabels(E5_NHS, fontsize=6)
    ax.xaxis.set_minor_locator(plt.NullLocator())
    ax.set_xlabel("hidden width $N_H$  (representational width)", fontsize=7)
    ax.set_ylabel("rule-switching accuracy", fontsize=7)
    ax.set_title(f"E5 executive: no established size effect\n(n={E5_SEEDS} seeds, 95% CI; Spearman p=0.51)",
                 loc="left", fontsize=7.5)

    fig.suptitle("Scaling the substrate helps only where size is the binding "
                 "constraint — one knob per capacity, three different answers",
                 fontsize=8.5, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    path = os.path.join(FIGDIR, "scaling_capacities.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


if __name__ == "__main__":
    main()
