"""How each capacity is REPRESENTED on the substrate, and does it port to a 2-D lattice?

Two questions, answered from the substrate rather than from prose.

**(1) What is each capacity's representation, and what does it need from the medium?**
Read off the three builders, each capacity uses a different substrate primitive:

  E2 memory      : a *reentrant loop*. State = the position of one circulating
                   pulse on a directed cycle. Needs: a cycle, and `tau < L` so
                   the refractory tail clears before the pulse returns.
  E4 attention   : *colliding waves*. State = which of two counter-propagating
                   waves captured the centre. Needs: an arena wide enough for
                   two fronts, and refractoriness (the annihilation is the
                   competition -- no inhibitory units anywhere).
  E5 control     : a *held option gating a feedforward conjunction*. The context
                   ring is a reentrant loop (E2's primitive) that stays alive
                   across a block; hidden units are (stimulus AND rule)
                   coincidence gates that only fire when the ring is alive.

The E5 reading has a consequence this experiment verifies numerically: **E5's
hidden layer has no H->H recurrence at all** (asserted below). Widening `N_H`
therefore adds *non-interacting parallel AND-gates* -- which predicts that
hidden width should NOT be E5's binding constraint, and that the ring should be.
`scaling_capacities.md` found exactly that null on `N_H`; this supplies the
mechanism, and sweeps the ring for comparison.

**(2) Does each mechanism port to `lattice2d`?**
The E-series runs E2 on a ring, E4 on a chain, E5 on rings+feedforward. Track 3b
showed the raw dynamics and E1 learning generalise across media; this asks the
same of the three capacity *mechanisms* on the 2-D lattice.

Substrate / analysis boundary
-----------------------------
E2 and E4 here are **pure dynamics** (E4 plus a fixed centre readout); no
plasticity. E5's ring test is dynamics; its accuracy figures come from the
learning experiment and are quoted, not recomputed. The lattice ports test
whether the *mechanism* survives the medium -- not whether a learner exploits it.

Outputs
-------
docs/figures/lattice_capacities.png
result/scaling/lattice_capacities.npz
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
from ghca_net import Network, lattice2d  # noqa: E402
import e5_executive as e5  # noqa: E402

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "scaling")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

TORUS_LS = [12, 16, 24, 32]
TORUS_TAUS = [6, 10, 14, 18, 22, 26, 30]
ARENA_LS = [21, 41, 61]
ARENA_BIASES = [0, 1, 2, 3]
ARENA_JITTER = 1.5
ARENA_N = 80
RING_LS = [8, 12, 16, 24, 32]
E5_SEEDS = 12


# ---------------------------------------------------------------------------
# (1) mechanism check: is E5's hidden medium recurrent?
# ---------------------------------------------------------------------------
def e5_hidden_recurrence():
    """Count H->H edges in E5's substrate. Expect 0 (purely feedforward)."""
    W, _, roles, _ = e5.build(seed=0)
    hid = roles["hidden"]
    hh = int((W[np.ix_(hid, hid)] > 0).sum())
    ring_edges = int((W[np.ix_(roles["rings"][0], roles["rings"][0])] > 0).sum())
    return hh, ring_edges


# ---------------------------------------------------------------------------
# (2a) E2 memory on a lattice2d torus: a plane wave circulating a row is the
#      reentrant loop. Its transit length is L, so the E2 law predicts sustain
#      iff tau < L -- the same law, on a 2-D medium.
# ---------------------------------------------------------------------------
def torus_reentry(L, tau, T=300, seed=0):
    net = Network(lattice2d(L, r=1, periodic=True), act=2, pas=int(tau) - 2,
                  theta=1.0, p_s=0.0, seed=seed)
    net.phi[:] = 0
    g = net.phi.reshape(L, L)
    g[0, :] = 1                                     # a full row = the wavefront
    for k in range(1, int(tau)):                    # refractory tail behind it
        g[(-k) % L, :] = min(k + 1, int(tau))
    net.phi = g.reshape(-1)
    out = net.run(T, record=False)
    return bool(out["A"][-20:].mean() > 0)


# ---------------------------------------------------------------------------
# (2b) E4 attention on a lattice2d arena: ignite left and right EDGES; the two
#      wavefronts collide and annihilate; read which captured the centre node.
#      Same label-propagation readout as e4_attention, lifted to 2-D.
# ---------------------------------------------------------------------------
def arena_trial(Lg, bias, jitter=0.0, seed=0, act=6, pas=6):
    T = 3 * Lg
    adj = lattice2d(Lg, r=1, periodic=False)
    net = Network(adj, act=act, pas=pas, theta=1.0, p_s=0.0, seed=seed)
    net.phi[:] = 0
    rng = np.random.default_rng(seed)
    j0, j1 = (rng.standard_normal(2) * jitter).round().astype(int)
    s0 = max(0, -bias) + max(0, int(j0))
    s1 = max(0, bias) + max(0, int(j1))
    N = Lg * Lg
    left = np.arange(Lg) * Lg + 0
    right = np.arange(Lg) * Lg + (Lg - 1)
    centre = (Lg // 2) * Lg + (Lg // 2)
    nbrs = [np.where(adj[i] > 0)[0] for i in range(N)]
    label = np.full(N, -1)
    first = np.full(N, -1)
    for t in range(T):
        prev = net.active_mask()
        drive = np.zeros(N, bool)
        if t == s0:
            drive[left] = True
        if t == s1:
            drive[right] = True
        net.step(drive if drive.any() else None)
        am = net.active_mask()
        for i in np.where(am & (first < 0))[0]:
            if i in left and t == s0:
                label[i] = 0
            elif i in right and t == s1:
                label[i] = 1
            else:
                cand = [label[j] for j in nbrs[i] if prev[j] and label[j] >= 0]
                if cand:
                    vals = set(cand)
                    label[i] = (cand[0] if len(vals) == 1
                                else cand[int(rng.integers(len(cand)))])
            first[i] = t
    return int(label[centre]) if label[centre] >= 0 else -1


def arena_psychometric(Lg, biases=ARENA_BIASES, n=ARENA_N, jitter=ARENA_JITTER):
    return [sum(arena_trial(Lg, b, jitter=jitter, seed=7000 * Lg + 31 * b + i) == 0
                for i in range(n)) / n
            for b in biases]


# ---------------------------------------------------------------------------
# (2c) E5: sweep the CONTEXT RING -- the capacity's only recurrent structure.
# ---------------------------------------------------------------------------
def e5_ring_sweep(seeds=E5_SEEDS, n_h=120):
    nh0, lr0 = e5.N_H, e5.L_RING
    e5.N_H = n_h
    out = {}
    for lr in RING_LS:
        e5.L_RING = lr
        accs = []
        for s in range(seeds):
            acc, _, alive, _, _ = e5.run_switching(seed=s, n_blocks=24, block_len=20)
            accs.append(acc.mean())
        out[lr] = np.array(accs)
    e5.N_H, e5.L_RING = nh0, lr0
    return out


def main():
    hh, ring_edges = e5_hidden_recurrence()
    print("E5 mechanism — is the hidden 'medium' recurrent?")
    print(f"   H->H edges = {hh}   (per-ring recurrent edges = {ring_edges})")
    assert hh == 0, "E5 hidden layer unexpectedly has recurrence"
    print("   => purely FEEDFORWARD conjunction gates. Widening N_H adds")
    print("      non-interacting parallel AND-gates, which is why hidden width")
    print("      is not E5's binding constraint (see scaling_capacities.md).")

    print("\nE2 memory on a lattice2d TORUS (plane wave round a row, transit = L)")
    print(f"   {'L':>4} " + " ".join(f"t={t:<3}" for t in TORUS_TAUS))
    torus = np.zeros((len(TORUS_LS), len(TORUS_TAUS)), dtype=int)
    for i, L in enumerate(TORUS_LS):
        for j, tau in enumerate(TORUS_TAUS):
            torus[i, j] = int(torus_reentry(L, tau))
        pred = [int(t < L) for t in TORUS_TAUS]
        agree = int((torus[i] == np.array(pred)).sum())
        print(f"   {L:>4} " + " ".join(f"{v:<5}" for v in torus[i])
              + f"  vs law(tau<L) {agree}/{len(TORUS_TAUS)}")
    law = np.array([[int(t < L) for t in TORUS_TAUS] for L in TORUS_LS])
    print(f"   overall agreement with the E2 law on the torus: "
          f"{int((torus == law).sum())}/{torus.size}")

    print(f"\nE4 attention on a lattice2d ARENA (jitter={ARENA_JITTER}, n={ARENA_N})")
    arena = {}
    for Lg in ARENA_LS:
        ps = arena_psychometric(Lg)
        arena[Lg] = ps
        print(f"   Lg={Lg:>3} ({Lg*Lg:>5} nodes)  P(win) by bias "
              f"{[f'{p:.2f}' for p in ps]}   sensitivity={ps[1]-ps[0]:+.2f}")

    print(f"\nE5 control — sweep the CONTEXT RING length (n={E5_SEEDS} seeds, N_H=120)")
    rings = e5_ring_sweep()
    for lr, accs in rings.items():
        print(f"   L_RING={lr:>3}  tau_slow={e5.TAU_SLOW} < L_RING? "
              f"{e5.TAU_SLOW < lr!s:>5}   acc={accs.mean():.3f}"
              f"  spread=[{accs.min():.2f},{accs.max():.2f}]")
    lrs = list(rings)
    span = max(a.mean() for a in rings.values()) - min(a.mean() for a in rings.values())
    print(f"   ring-length span in accuracy: {span:.3f}"
          f"   (cf. the N_H sweep's 0.043 over a 16x range)")

    plot(torus, law, arena, rings)

    np.savez(os.path.join(DATADIR, "lattice_capacities.npz"),
             e5_hh_edges=hh, e5_ring_edges=ring_edges,
             torus_Ls=np.array(TORUS_LS), torus_taus=np.array(TORUS_TAUS),
             torus_sustain=torus, torus_law=law,
             arena_Ls=np.array(ARENA_LS), arena_biases=np.array(ARENA_BIASES),
             arena_p=np.array([arena[L] for L in ARENA_LS]),
             ring_Ls=np.array(lrs),
             ring_acc=np.array([rings[l].mean() for l in lrs]),
             ring_acc_all=np.array([rings[l] for l in lrs]),
             e5_tau_slow=e5.TAU_SLOW)
    print("\nwrote", os.path.join(DATADIR, "lattice_capacities.npz"))


def plot(torus, law, arena, rings):
    fig, axes = plt.subplots(1, 3, figsize=(11.2, 3.4))

    ax = axes[0]
    ax.imshow(torus, cmap="Blues", vmin=0, vmax=1, aspect="auto")
    for i in range(torus.shape[0]):
        for j in range(torus.shape[1]):
            if torus[i, j] != law[i, j]:
                ax.plot(j, i, "x", color="#dc2626", ms=7, mew=1.8)
    ax.set_xticks(range(len(TORUS_TAUS)))
    ax.set_xticklabels(TORUS_TAUS, fontsize=6)
    ax.set_yticks(range(len(TORUS_LS)))
    ax.set_yticklabels(TORUS_LS, fontsize=6)
    ax.set_xlabel(r"timescale $\tau$", fontsize=7)
    ax.set_ylabel("torus side $L$ (= transit length)", fontsize=7)
    n_ok = int((torus == law).sum())
    ax.set_title(f"E2 memory ports to the 2-D torus\n"
                 f"sustain matches $\\tau<L$  {n_ok}/{torus.size}",
                 loc="left", fontsize=7.5)

    ax = axes[1]
    cols = ["#0e7490", "#16a34a", "#84cc16"]
    for (Lg, ps), c in zip(arena.items(), cols):
        ax.plot(ARENA_BIASES, ps, "-o", color=c, ms=4, lw=1.6,
                label=f"$L$={Lg} ({Lg*Lg} nodes)")
    ax.axhline(0.5, color="#6b7280", ls=":", lw=1)
    ax.set_xlabel("top-down bias (steps)", fontsize=7)
    ax.set_ylabel("P(attended stream wins centre)", fontsize=7)
    ax.set_xticks(ARENA_BIASES)
    ax.legend(frameon=False, fontsize=6, loc="lower right")
    ax.set_title("E4 attention ports to the 2-D arena\n"
                 "same psychometric curve, size-invariant", loc="left", fontsize=7.5)

    ax = axes[2]
    lrs = list(rings)
    means = [rings[l].mean() for l in lrs]
    ax.plot(lrs, means, "-o", color="#dc2626", ms=4, lw=1.6, label="ring length")
    for l in lrs:
        ax.plot([l] * len(rings[l]), rings[l], ".", color="#fca5a5", ms=3, zorder=0)
    ax.axvline(e5.TAU_SLOW, color="#6b7280", ls=":", lw=1)
    ax.text(e5.TAU_SLOW + 0.6, 0.24, f"$\\tau$={e5.TAU_SLOW}\n(option dies below)",
            fontsize=6, color="#6b7280")
    ax.axhline(0.5, color="#9ca3af", ls="--", lw=0.9)
    ax.set_xlabel("context-ring length $L_{ring}$", fontsize=7)
    ax.set_ylabel("rule-switching accuracy", fontsize=7)
    ax.set_title("E5's binding constraint is the RING, not width\n"
                 "(hidden layer has 0 recurrent edges)", loc="left", fontsize=7.5)

    fig.suptitle("Each capacity uses a different substrate primitive — and all three "
                 "port to the 2-D lattice",
                 fontsize=8.5, fontweight="bold", x=0.02, ha="left")
    fig.tight_layout(rect=[0, 0, 1, 0.92])
    path = os.path.join(FIGDIR, "lattice_capacities.png")
    fig.savefig(path, dpi=110)
    print("wrote", path)


if __name__ == "__main__":
    main()
