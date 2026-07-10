"""
C5 - Is do(chirality) fat-handed?  (topological invariant vs realization)

The C-series showed that intervening on a constituted scalar aggregate W = f(S)
is FAT-HANDED (C2): many micro-states S realize the same W, so do(W=w) does not
pin behaviour -- the realization policy does. C5 asks the same question for the
E7 spiral's ROTATION DIRECTION (chirality), which is different in kind: chirality
is a TOPOLOGICAL INVARIANT (a winding number), realization-invariant by
definition. So does do(chirality) escape C2's fat-handedness, or does the way the
behaviour READS the chirality re-introduce it?

We force the phase field to a target winding chi = +1 (or -1) under a family of
realization policies -- the same topological charge, different micro-states:
  * centered   : core at the lattice centre (the canonical E7 seed)
  * displaced-d: core shifted by d nodes (the key realization freedom)
  * pitch      : an added radial phase gradient (same winding, looser/tighter arm)
  * noisy      : phase jitter on a fraction of nodes (same core, different micro-state)
and read the chirality back three ways:
  * center winding : local winding on a fixed loop at the lattice centre (what the
                     E7 router uses)
  * tracked winding: local winding around the *tracked* core (found first)
  * global charge  : sign of the net topological charge over the whole field

Behaviour proxy B = whether the readout recovers the true chirality (the E7 router
routes deterministically on that bit, so its band tracks the decode band). We
report the achievable BAND of B across realizations per readout, in sigma units
(as in C2/C3), plus a routing confirmation on a frozen E7 router.

Prediction / discriminator:
  Hypothesis A (fat-handed like C2): the fixed-centre readout gives a WIDE band --
    displaced-core realizations carry the same chi but read wrong.
  Hypothesis B (topology escapes C2): the tracked / global readouts give a NARROW
    band -- chirality is a better-posed wave variable than a scalar aggregate.
  The band's dependence on readout LOCALITY is the result: chirality is well-posed
  iff read topologically (tracked/global), fat-handed iff read at a fixed locus.

Outputs
-------
docs/figures/c5_do_chirality.png
result/c5/c5_data.npz
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
from ghca_net import Network, lattice2d
import e7_spiral_option as sp
import e7_learning as e7

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c5")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

L = sp.L
TAU = sp.TAU
C = L // 2
STEPS = 4                      # readout window (as in the E7 router decode)


def make_spiral():
    net = Network(lattice2d(L, r=sp.RANGE, periodic=False), act=sp.ACT, pas=sp.PAS,
                  theta=sp.THETA, p_s=0.0, seed=0)
    net.tau[:] = TAU
    return net


def realize(net, chir, policy, rng, d=0, pitch=0.0, jitter=0.0):
    """Inject a phase field with winding = chir under a realization policy.
    All policies produce the SAME topological charge; they differ in the
    (realization-free) details: core location, arm pitch, micro-noise."""
    if policy == "centered":
        cx = cy = C
    elif policy == "displaced":
        theta = rng.uniform(0, 2 * np.pi)
        cx = int(round(C + d * np.cos(theta)))
        cy = int(round(C + d * np.sin(theta)))
    else:
        cx = cy = C
    xs = np.arange(L)
    X, Y = np.meshgrid(xs, xs, indexing="ij")
    ang = np.arctan2(Y - cy, X - cx)
    frac = ((chir * ang) / (2 * np.pi)) % 1.0
    if policy == "pitch":
        r = np.hypot(X - cx, Y - cy)
        frac = (frac + pitch * r / L) % 1.0        # radial gradient: 0 net winding
    phi = np.clip(np.floor(frac * (TAU + 1)).astype(np.int64), 0, TAU)
    if policy == "noisy" and jitter > 0:
        flip = rng.random((L, L)) < jitter
        phi[flip] = rng.integers(0, TAU + 1, size=int(flip.sum()))
    net.phi = phi.reshape(-1).copy()
    return cx, cy


def read_center(net):
    return sp.local_winding(net.phi, C, C)


def read_tracked(net):
    _, _, _, core = sp.signed_charge(net.phi)
    cx, cy = core
    if not np.isfinite(cx):
        return 0.0
    return sp.local_winding(net.phi, int(round(cx)), int(round(cy)))


def read_global(net):
    _, _, netq, _ = sp.signed_charge(net.phi)
    return float(netq)


READERS = {"center": read_center, "tracked": read_tracked, "global": read_global}


def decode(net, reader, steps=STEPS):
    """Accumulate a readout over `steps` and return the decoded chirality sign."""
    s = 0.0
    for _ in range(steps):
        net.step(None)
        s += READERS[reader](net)
    return 0 if abs(s) < 0.5 else int(np.sign(s))


def measure_bands(n_trials=40):
    """For each readout, decode accuracy per policy; the band across policies is the
    fat-handedness measure. Displacement is swept as the main realization axis."""
    policies = [("centered", {}),
                ("displaced", {"d": 4}), ("displaced", {"d": 8}), ("displaced", {"d": 12}),
                ("pitch", {"pitch": 1.0}), ("pitch", {"pitch": 2.0}),
                ("noisy", {"jitter": 0.08})]
    labels = ["centered", "disp-4", "disp-8", "disp-12", "pitch-1", "pitch-2", "noisy"]
    net = make_spiral()
    # acc[reader][policy] = mean decode-correct; also keep per-trial for sigma
    acc = {r: np.zeros(len(policies)) for r in READERS}
    per_trial = {r: [np.zeros(n_trials) for _ in policies] for r in READERS}
    for pi, (policy, kw) in enumerate(policies):
        for r in READERS:
            for k in range(n_trials):
                rng = np.random.default_rng(1000 * pi + k)
                chir = +1 if k % 2 == 0 else -1
                realize(net, chir, policy, rng, **kw)
                g = decode(net, r)
                ok = float(g == chir)
                per_trial[r][pi][k] = ok
            acc[r][pi] = per_trial[r][pi].mean()
    # band in sigma units: (max-min across policies) / pooled within-policy std
    bands = {}
    for r in READERS:
        pooled = np.sqrt(np.mean([per_trial[r][pi].var() for pi in range(len(policies))]))
        rng_band = acc[r].max() - acc[r].min()
        bands[r] = rng_band / (pooled + 1e-6)
    return labels, acc, bands


def routing_confirmation(n_seeds=3):
    """Behavioural band: freeze an E7 router (trained on the centred spiral), then
    route on displaced realizations decoded by the fixed-centre vs tracked readout.
    B = routing accuracy; a wide center-vs-tracked gap = fat-handedness in behaviour."""
    res = {"center": [], "tracked": []}
    for s in range(n_seeds):
        rnet, roles = e7.make_router(s)
        rng = np.random.default_rng(s + 5)
        spnet = make_spiral()
        # train the router on the centred spiral (correct context) for both rules
        for t in range(1200):
            rule_g = t % 2
            chir = +1 if rule_g == 0 else -1
            sp.seed_spiral(spnet, chir, jitter=0.02, rng=rng)
            gl, _, _ = e7.spiral_decode(spnet)
            g_ctx = gl if gl in (0, 1) else int(rng.integers(2))
            e7.trial(rnet, roles, g_ctx, int(rng.integers(2)), rule_g, rng, learn=True)
        # evaluate on DISPLACED realizations under each readout
        for reader in ("center", "tracked"):
            ok = 0
            for k in range(120):
                r2 = np.random.default_rng(9000 + k)
                rule_g = k % 2
                chir = +1 if rule_g == 0 else -1
                realize(spnet, chir, "displaced", r2, d=8)
                g = decode(spnet, reader)
                g_ctx = (0 if g > 0 else 1) if g != 0 else int(r2.integers(2))
                x = int(r2.integers(2))
                a = e7.trial(rnet, roles, g_ctx, x, rule_g, r2, learn=False)
                ok += a
            res[reader].append(ok / 120)
    return {k: (np.mean(v), np.std(v) / np.sqrt(n_seeds)) for k, v in res.items()}


def main():
    print("C5: decode bands across realization policies ...")
    labels, acc, bands = measure_bands()
    for r in READERS:
        print(f"  readout={r:8s} decode acc by policy: "
              + " ".join(f"{lab}={acc[r][i]:.2f}" for i, lab in enumerate(labels)))
    print("  fat-handedness band (sigma units; low = well-posed):")
    for r in READERS:
        print(f"    {r:8s}: {bands[r]:.1f} sigma")

    print("C5: routing confirmation on displaced realizations (frozen E7 router) ...")
    rc = routing_confirmation()
    print(f"  routing acc @ displaced core: center-readout={rc['center'][0]:.2f}  "
          f"tracked-readout={rc['tracked'][0]:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    colors = {"center": "crimson", "tracked": "seagreen", "global": "steelblue"}
    x = np.arange(len(labels))
    for r in READERS:
        axes[0].plot(x, acc[r], "o-", color=colors[r], label=f"{r} readout")
    axes[0].axhline(0.5, ls=":", color="k", alpha=0.5)
    axes[0].set_xticks(x); axes[0].set_xticklabels(labels, rotation=30, fontsize=8)
    axes[0].set_ylabel("decode accuracy (recover true χ)"); axes[0].set_ylim(0, 1.05)
    axes[0].set_title("do(χ) realizations: the fixed-centre readout fails\n"
                      "on displaced cores; tracked/global do not")
    axes[0].legend(fontsize=8)

    axes[1].bar(list(READERS), [bands[r] for r in READERS],
                color=[colors[r] for r in READERS])
    axes[1].set_ylabel("fat-handedness band (σ)")
    axes[1].set_title("Band of the causal verdict across realizations\n"
                      "(low = well-posed; high = fat-handed)")

    axes[2].bar(["center", "tracked"],
                [rc["center"][0], rc["tracked"][0]],
                yerr=[rc["center"][1], rc["tracked"][1]],
                color=["crimson", "seagreen"], capsize=4)
    axes[2].axhline(0.5, ls=":", color="k", alpha=0.5)
    axes[2].set_ylim(0, 1); axes[2].set_ylabel("routing accuracy (frozen router)")
    axes[2].set_title("Behavioural band @ displaced core:\nreadout locality decides do(χ)")

    fig.suptitle("C5: do(chirality) is fat-handed only when read at a fixed locus — "
                 "topological chirality escapes C2's fat-handedness if read topologically",
                 fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "c5_do_chirality.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c5_data.npz"),
             labels=np.array(labels),
             **{f"acc_{r}": acc[r] for r in READERS},
             **{f"band_{r}": bands[r] for r in READERS},
             routing_center=rc["center"], routing_tracked=rc["tracked"])
    print("wrote", os.path.join(DATADIR, "c5_data.npz"))


if __name__ == "__main__":
    main()
