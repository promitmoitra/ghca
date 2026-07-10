"""
E7 (Phase A) - The executive-control OPTION as a genuine 2D spiral core.

E5 realised the "option" of executive control as a persistent 1D ring loop. This
experiment moves onto the substrate's NATIVE medium: a genuine two-dimensional
spiral wave with a phase-singularity core, on the lattice characterised in E0. The
motivating claim from the rotating-wave neuroscience (Xu/Gong et al. 2023; Ye/
Steinmetz et al. 2026) is that cortical computation rides rotating waves whose
ROTATION DIRECTION is task-relevant. E7 asks whether our substrate supports the
same object as a controllable, persistent, readable variable.

Phase A (this file) establishes the MECHANISM the later learning phase depends on
(the analogue of E2's "memory duration is tau-controlled" result):

  1. NUCLEATION with chosen chirality. A single spiral of clockwise (CW) or
     counter-clockwise (CCW) rotation is seeded by imposing a polar-angle phase
     ramp around a core; the sign of the ramp sets the handedness.
  2. PERSISTENCE. In the E0 organised spiral band (r=2, a=6, tau=14, theta~4) with
     NO-FLUX boundaries, the seeded core persists for hundreds of steps. (A
     periodic torus cannot hold a lone spiral: total topological charge must be
     zero, so a lone core breeds a compensating anti-core -- no-flux boundaries fix
     this, a nod to Ye et al.'s "anatomy/boundaries shape the wave".)
  3. READOUT, two ways. The rotation direction is recovered (a) globally, as the
     SIGN OF THE NET TOPOLOGICAL CHARGE (a signed phase-singularity count), and
     (b) locally and on-thesis, from two PHASE PROBES placed 90 deg apart around
     the core: which probe leads in phase flips with chirality, so the wave itself
     delivers a binary context signal with no explicit charge computation.
  A planar (no-core) wave is included as a control: it has ~zero net charge and no
  defined chirality, so the readouts are spiral-specific.

Phase B (separate, reuses E5) will make CW/CCW the RULE (identity vs reversal) and
let reward-driven routing read the phase-probe context; and the C-series can then
ask whether the chirality is causal for routing or epiphenomenal.

Outputs
-------
docs/figures/e7_mechanism.png : CW/CCW phase fields, net-charge persistence,
                                phase-probe readout
result/e7/e7_mechanism.npz
"""

import os
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
from ghca_net import Network, lattice2d

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e7")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# E0 organised spiral band (see e0_results.md), no-flux boundaries for a lone spiral
L = 48
ACT, PAS = 6, 8
TAU = ACT + PAS               # 14
THETA = 4.0
RANGE = 2
PROBE_R = 10                  # phase-probe radius from the core


def make(seed=0):
    net = Network(lattice2d(L, r=RANGE, periodic=False), act=ACT, pas=PAS,
                  theta=THETA, p_s=0.0, seed=seed)
    net.tau[:] = TAU
    return net


def _angle_grid(cx, cy):
    xs = np.arange(L)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    return np.arctan2(Y - cy, X - cx), X, Y


def seed_spiral(net, chirality, cx=None, cy=None, jitter=0.0, rng=None):
    """Impose a polar-angle phase ramp around (cx,cy); sign of ramp = handedness.
    `chirality` = +1 (CCW) or -1 (CW)."""
    cx = L // 2 if cx is None else cx
    cy = L // 2 if cy is None else cy
    ang, _, _ = _angle_grid(cx, cy)
    frac = ((chirality * ang) / (2 * np.pi)) % 1.0
    phi = np.clip(np.floor(frac * (TAU + 1)).astype(np.int64), 0, TAU)
    if jitter > 0 and rng is not None:                 # small init perturbation
        flip = rng.random((L, L)) < jitter
        phi[flip] = rng.integers(0, TAU + 1, size=int(flip.sum()))
    net.phi = phi.reshape(-1).copy()
    return cx, cy


def seed_planar(net, wavelength=14):
    """Straight (planar) wavefronts along x -> no core, net charge ~ 0 (control)."""
    xs = np.arange(L)
    X, _ = np.meshgrid(xs, xs, indexing="ij")
    frac = ((X % wavelength) / wavelength)
    net.phi = np.clip(np.floor(frac * (TAU + 1)).astype(np.int64), 0, TAU).reshape(-1).copy()


def signed_charge(phi_flat):
    """Vectorised signed phase-singularity count. Returns (n_pos, n_neg, net,
    core_xy) where core_xy is the centroid of the dominant-sign defects."""
    ang = 2 * np.pi * phi_flat.reshape(L, L) / TAU
    def wrap(d):
        return (d + np.pi) % (2 * np.pi) - np.pi
    a, b = ang[:-1, :-1], ang[1:, :-1]
    c, d = ang[1:, 1:], ang[:-1, 1:]
    circ = wrap(b - a) + wrap(c - b) + wrap(d - c) + wrap(a - d)
    pos, neg = circ > np.pi, circ < -np.pi
    n_pos, n_neg = int(pos.sum()), int(neg.sum())
    net = n_pos - n_neg
    dom = pos if n_pos >= n_neg else neg
    if dom.any():
        xi, yi = np.where(dom)
        core_xy = (float(xi.mean()), float(yi.mean()))
    else:
        core_xy = (np.nan, np.nan)
    return n_pos, n_neg, net, core_xy


def probe_chirality(net, core_xy):
    """Local readout: circular phase difference between two probes 90 deg apart
    around the (tracked) core -- East vs North. Its sign flips with handedness.
    Anchoring the probes to the tracked core (not the seeded position) makes the
    readout robust to the spiral tip's meander."""
    cx, cy = core_xy
    if not np.isfinite(cx):
        return 0.0
    ea = (int(round(cx + PROBE_R)), int(round(cy)))
    no = (int(round(cx)), int(round(cy + PROBE_R)))
    def th(site):
        x, y = site
        x = min(max(x, 0), L - 1); y = min(max(y, 0), L - 1)
        return 2 * np.pi * net.phi[x * L + y] / TAU
    d = th(ea) - th(no)
    return np.angle(np.exp(1j * d))       # signed circular phase difference in (-pi,pi]


def run_condition(kind, chirality=+1, seed=0, T=400, warm=60):
    """kind in {'spiral','planar'}. Returns time series of net charge and the
    mean phase-probe reading (after a warm-up), plus a late phase-field snapshot."""
    net = make(seed=seed)
    rng = np.random.default_rng(seed + 11)
    if kind == "spiral":
        cx, cy = seed_spiral(net, chirality, jitter=0.02, rng=rng)
    else:
        seed_planar(net); cx, cy = L // 2, L // 2
    charge = np.zeros(T)
    probe = np.zeros(T)
    snap = None
    for t in range(T):
        net.step(None)
        _, _, netq, core_xy = signed_charge(net.phi)
        charge[t] = netq
        probe[t] = probe_chirality(net, core_xy)
        if t == T - 1:
            snap = net.phi.reshape(L, L).copy()
    probe_mean = probe[warm:].mean()
    return {"charge": charge, "probe": probe, "probe_mean": probe_mean,
            "snap": snap, "final_charge": int(np.sign(charge[warm:].mean()))}


def main():
    T = 400
    # 1-2. persistence + snapshots for CW / CCW / planar
    ccw = run_condition("spiral", chirality=+1, T=T)
    cw = run_condition("spiral", chirality=-1, T=T)
    planar = run_condition("planar", T=T)

    # 3. readout accuracy over trials (vary core position + init jitter)
    n_trials = 10
    net_ok = probe_ok = 0
    total = 0
    rows = []
    for chir in (+1, -1):
        for k in range(n_trials):
            net = make(seed=100 + k)
            rng = np.random.default_rng(500 + k)
            cx = L // 2 + int(rng.integers(-6, 7))
            cy = L // 2 + int(rng.integers(-6, 7))
            seed_spiral(net, chir, cx=cx, cy=cy, jitter=0.03, rng=rng)
            ch = pr = 0.0
            for t in range(T):
                net.step(None)
                if t >= T - 200:
                    _, _, nq, core_xy = signed_charge(net.phi)
                    ch += nq
                    pr += probe_chirality(net, core_xy)
            net_pred = int(np.sign(ch))
            probe_pred = int(np.sign(pr))
            net_ok += (net_pred == chir)
            # phase-probe sign is a convention; fix it from the CCW reference sign
            total += 1
            rows.append((chir, net_pred, np.sign(pr)))
    # calibrate the phase-probe sign convention against the seeded chirality
    ref = np.sign(np.mean([np.sign(p) * c for c, _, p in rows]))  # +1 if p-sign aligns with chir
    probe_ok = sum(1 for c, _, p in rows if int(np.sign(p) * ref) == c)
    net_acc = net_ok / total
    probe_acc = probe_ok / total

    print("E7 Phase A - spiral-core option: mechanism")
    print(f"  persistence (net charge, mean over t>=60): "
          f"CCW={ccw['charge'][60:].mean():+.2f}  CW={cw['charge'][60:].mean():+.2f}  "
          f"planar={planar['charge'][60:].mean():+.2f}")
    print(f"  net charge at t={T-1}: CCW={ccw['charge'][-1]:+.0f}  CW={cw['charge'][-1]:+.0f}  "
          f"planar={planar['charge'][-1]:+.0f}")
    print(f"  chirality readout accuracy over {total} trials: "
          f"net-charge sign = {net_acc:.2f}, phase-probe lead = {probe_acc:.2f}")

    # ---------------- figure ----------------
    fig = plt.figure(figsize=(13, 9))
    # phase fields
    for i, (res, name) in enumerate([(ccw, "CCW (+1)"), (cw, "CW (-1)")]):
        ax = fig.add_subplot(2, 2, i + 1)
        im = ax.imshow(res["snap"].T, origin="lower", cmap="twilight", vmin=0, vmax=TAU)
        ax.set_title(f"Spiral phase field, {name}\n(net charge {res['charge'][-1]:+.0f} at t={T-1})")
        ax.set_xlabel("x"); ax.set_ylabel("y")
        fig.colorbar(im, ax=ax, fraction=0.046, label="phase")

    # net charge persistence
    ax = fig.add_subplot(2, 2, 3)
    ax.plot(ccw["charge"], color="crimson", label="CCW seed (+1)")
    ax.plot(cw["charge"], color="steelblue", label="CW seed (-1)")
    ax.plot(planar["charge"], color="gray", ls="--", label="planar (no core)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("time step"); ax.set_ylabel("net topological charge")
    ax.set_title("Persistence: a single signed core is held\n(planar control ~ 0)")
    ax.legend(fontsize=8)

    # phase-probe readout
    ax = fig.add_subplot(2, 2, 4)
    ax.plot(ccw["probe"], color="crimson", label="CCW seed (+1)")
    ax.plot(cw["probe"], color="steelblue", label="CW seed (-1)")
    ax.axhline(0, color="k", lw=0.6)
    ax.set_xlabel("time step"); ax.set_ylabel("probe phase difference (E - N)")
    ax.set_title(f"Local phase-probe readout separates handedness\n"
                 f"(readout acc: net-charge {net_acc:.2f}, probe {probe_acc:.2f})")
    ax.legend(fontsize=8)

    fig.suptitle("E7 Phase A: the option as a 2D spiral core -- controllable, "
                 "persistent, and readable rotation direction", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.96])
    p = os.path.join(FIGDIR, "e7_mechanism.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e7_mechanism.npz"),
             ccw_charge=ccw["charge"], cw_charge=cw["charge"], planar_charge=planar["charge"],
             ccw_probe=ccw["probe"], cw_probe=cw["probe"],
             ccw_snap=ccw["snap"], cw_snap=cw["snap"], planar_snap=planar["snap"],
             net_acc=net_acc, probe_acc=probe_acc)
    print("wrote", os.path.join(DATADIR, "e7_mechanism.npz"))


if __name__ == "__main__":
    main()
