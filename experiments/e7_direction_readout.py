"""
E7 / Track 1b - A LEARNED direction-selective readout of rotation direction.

E7 reads the spiral's chirality with `local_winding` -- a god's-eye topological
integral of the phase gradient around a loop centred on the core. It is accurate
but (i) *computed*, not learned, and (ii) *fixed-locus*, so it collapses when the
core drifts away from the readout centre (the C5 readout-locality result). Track
1b (see docs/next_steps.md) replaces it with a small population of biologically
plausible **direction-selective cells** whose pooling weights are **learned**:

  * Elementary motion detectors (EMD, Hassenstein-Reichardt / motion-energy):
    for each site and each of the 4 axis directions, the delayed coincidence
    `active[t-1] at i AND active[t] at i+d`, accumulated over a short window and
    coarse-grained -- a local motion-energy field. This is what the *substrate*
    supplies (local, translation-covariant, no core needed).
  * A linear readout over that field, **trained** (logistic regression on the
    known chirality) to output CCW vs CW -- this is what is *learned*, replacing
    the hand-coded winding integral.

Two results:
  A. RECOVERY + ROBUSTNESS. The learned population recovers chirality as well as
     `local_winding` on a centred core, and -- because it pools local detectors
     across the whole field rather than reading one locus -- stays accurate when
     the core is DISPLACED, exactly where the fixed-locus winding reader fails
     (the C5 problem, escaped by a distributed learned reader).
  B. ROUTING. Dropped into E7's switching task in place of `spiral_decode`, the
     learned readout drives the (stimulus x context) routing to the same
     switching accuracy as the computed readout -- rotation direction is read off
     the wave by a learned population and still carries the rule.

Honest boundary: the EMD *feature* is hand-specified (though local and
biologically motivated); what is *learned* is the pooling readout. Training is
supervised on the chirality label (not reward-driven) -- a reward/self-supervised
variant is deferred. So 1b moves the chirality readout from "computed god's-eye
winding" to "learned population of local motion detectors", and removes the
fixed-locus fragility, but the motion primitive itself is still designed.

Outputs
-------
docs/figures/e7_direction_readout.png
result/e7/e7_direction_readout.npz
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
import ghca_stats as st
import e7_spiral_option as sp
import e7_learning as e7

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e7")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

L, TAU, ACT, PAS, THETA, RANGE = sp.L, sp.TAU, sp.ACT, sp.PAS, sp.THETA, sp.RANGE
FRAMES = 6                                   # motion-integration window
BLOCK = 4                                    # coarse-grain block size
DIRS = [(1, 0), (-1, 0), (0, 1), (0, -1)]    # 4 axis motion directions
DISPS = [0, 4, 8, 12]


def make_spiral(seed=0, theta=THETA):
    net = Network(lattice2d(L, r=RANGE, periodic=False), act=ACT, pas=PAS,
                  theta=theta, p_s=0.0, seed=seed)
    net.tau[:] = TAU
    return net


# ---------------------------------------------------------------------------
# The direction-selective feature: a coarse-grained EMD motion field.
# ---------------------------------------------------------------------------

def emd_field(net, frames=FRAMES):
    """Accumulate the 4-direction elementary-motion-detector field over `frames`
    steps and coarse-grain BLOCKxBLOCK. Returns a flat feature vector."""
    prev = net.active_mask().reshape(L, L).astype(float)
    acc = np.zeros((L, L, 4))
    for _ in range(frames):
        net.step(None)
        cur = net.active_mask().reshape(L, L).astype(float)
        for k, (dx, dy) in enumerate(DIRS):
            acc[:, :, k] += prev * np.roll(np.roll(cur, dx, 0), dy, 1)
        prev = cur
    cg = acc.reshape(L // BLOCK, BLOCK, L // BLOCK, BLOCK, 4).sum((1, 3))
    return cg.reshape(-1)


def _seed_spiral_at(net, chir, disp, rng):
    cx = cy = L // 2
    if disp:
        th = rng.uniform(0, 2 * np.pi)
        cx = int(round(L // 2 + disp * np.cos(th)))
        cy = int(round(L // 2 + disp * np.sin(th)))
    sp.seed_spiral(net, chir, cx=cx, cy=cy, jitter=0.02, rng=rng)


def dataset(n, base_seed, disp=0, max_age=0):
    """n labelled EMD examples (alternating CCW/CW), core displaced by `disp`.

    `max_age` > 0 samples the spiral at a random age in [0, max_age] steps after
    seeding (before the EMD window), so the readout is trained on the same
    spiral-age distribution it meets when deployed per-trial in the switching task
    (where the core is nucleated once per block and then ages across trials).
    Part A keeps max_age=0 (fresh cores) — it isolates *spatial* robustness."""
    X, y = [], []
    for i in range(n):
        chir = +1 if i % 2 == 0 else -1
        rng = np.random.default_rng(base_seed + i)
        net = make_spiral(seed=base_seed + i)
        _seed_spiral_at(net, chir, disp, rng)
        for _ in range(int(rng.integers(0, max_age + 1))):
            net.step(None)
        X.append(emd_field(net))
        y.append(0 if chir > 0 else 1)          # 0 = CCW, 1 = CW
    return np.array(X), np.array(y)


# ---------------------------------------------------------------------------
# The learned readout: seeded numpy logistic regression (no global RNG).
# ---------------------------------------------------------------------------

class DirReadout:
    def __init__(self, l2=1.0, iters=300, lr=0.5):
        self.l2, self.iters, self.lr = l2, iters, lr

    def fit(self, X, y):
        Xb = np.hstack([X, np.ones((len(X), 1))])
        self.mu, self.sd = Xb.mean(0), Xb.std(0) + 1e-9
        Xn = (Xb - self.mu) / self.sd
        w = np.zeros(Xn.shape[1])
        for _ in range(self.iters):
            p = 1.0 / (1.0 + np.exp(-Xn @ w))
            w -= self.lr * (Xn.T @ (p - y) / len(y) + self.l2 * w / len(y))
        self.w = w
        return self

    def _prob(self, X):
        Xn = (np.hstack([X, np.ones((len(X), 1))]) - self.mu) / self.sd
        return 1.0 / (1.0 + np.exp(-Xn @ self.w))

    def predict(self, X):
        return (self._prob(X) > 0.5).astype(int)

    def decode_bit(self, net):
        """Chirality bit {0,1} from one spiral's EMD field (the routing hook)."""
        return int(self.predict(emd_field(net)[None, :])[0])


def winding_decode(net, frames=FRAMES):
    """Computed fixed-locus baseline: sign of local_winding at the lattice centre."""
    s = 0.0
    for _ in range(frames):
        net.step(None)
        s += sp.local_winding(net.phi, L // 2, L // 2)
    return 0 if s > 0.5 else (1 if s < -0.5 else -1)


def winding_accuracy(disp, n, base_seed):
    ok = 0
    for i in range(n):
        chir = +1 if i % 2 == 0 else -1
        rng = np.random.default_rng(base_seed + i)
        net = make_spiral(seed=base_seed + i)
        _seed_spiral_at(net, chir, disp, rng)
        ok += (winding_decode(net) == (0 if chir > 0 else 1))
    return ok / n


# ---------------------------------------------------------------------------
# A. recovery + robustness to core displacement.
# ---------------------------------------------------------------------------

def result_A(n_seeds=15, n_train=60, n_test=30):
    learned = {d: np.zeros(n_seeds) for d in DISPS}
    computed = {d: np.zeros(n_seeds) for d in DISPS}
    for s in range(n_seeds):
        Xtr, ytr = dataset(n_train, base_seed=10_000 * s)      # centred training
        rd = DirReadout().fit(Xtr, ytr)
        for d in DISPS:
            Xte, yte = dataset(n_test, base_seed=10_000 * s + 5_000 + d, disp=d)
            learned[d][s] = float((rd.predict(Xte) == yte).mean())
            computed[d][s] = winding_accuracy(d, n_test, base_seed=10_000 * s + 8_000 + d)
        print(f"  A seed {s + 1}/{n_seeds} done", flush=True)
    return learned, computed


# ---------------------------------------------------------------------------
# B. routing confirmation: E7 switching with the learned readout swapped in.
# ---------------------------------------------------------------------------

def run_switching_learned(seed, rd, n_blocks=20, block_len=25):
    """E7 run_switching, but the per-trial context bit comes from the LEARNED
    readout instead of spiral_decode's local_winding."""
    rnet, roles = e7.make_router(seed)
    spnet = make_spiral(seed + 7)
    rng = np.random.default_rng(seed + 5)
    acc = []
    for b in range(n_blocks):
        s = +1 if b % 2 == 0 else -1
        rule_g = 0 if s > 0 else 1
        sp.seed_spiral(spnet, s, jitter=0.02, rng=rng)
        for i in range(block_len):
            g = rd.decode_bit(spnet)
            g_ctx = g if g in (0, 1) else int(rng.integers(2))
            x = int(rng.integers(2))
            r = e7.trial(rnet, roles, g_ctx, x, rule_g, rng, learn=True)
            acc.append(r)
    return np.array(acc)


def result_B(n_seeds=5, n_blocks=20, block_len=25):
    learned_sw, computed_sw = np.zeros(n_seeds), np.zeros(n_seeds)
    tail = block_len * 4
    # the core is nucleated once per block and ages across trials, so the readout
    # must decode aged cores: pre-train across the within-block age distribution.
    max_age = block_len * FRAMES
    for s in range(n_seeds):
        Xtr, ytr = dataset(120, base_seed=700_000 + 10_000 * s, max_age=max_age)
        rd = DirReadout().fit(Xtr, ytr)
        learned_sw[s] = run_switching_learned(s, rd, n_blocks, block_len)[-tail:].mean()
        computed_sw[s] = e7.run_switching(s, ablate=False,
                                          n_blocks=n_blocks, block_len=block_len)[0][-tail:].mean()
        print(f"  B seed {s + 1}/{n_seeds} done", flush=True)
    return learned_sw, computed_sw


def main():
    print("=== A. chirality recovery vs core displacement ===", flush=True)
    learned, computed = result_A()
    print("  displacement:   " + "  ".join(f"d={d}" for d in DISPS))
    for name, res in [("learned", learned), ("computed(winding)", computed)]:
        row = "  ".join(f"{res[d].mean():.2f}" for d in DISPS)
        print(f"  {name:18s} {row}")

    print("\n=== B. routing (switching) with learned readout ===", flush=True)
    learned_sw, computed_sw = result_B()
    print(f"  switching: learned={learned_sw.mean():.2f}  computed={computed_sw.mean():.2f}")

    # figure
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    lm = [learned[d].mean() for d in DISPS]
    cm = [computed[d].mean() for d in DISPS]
    le = [1.96 * learned[d].std() / np.sqrt(len(learned[d])) for d in DISPS]
    ce = [1.96 * computed[d].std() / np.sqrt(len(computed[d])) for d in DISPS]
    ax.errorbar(DISPS, lm, yerr=le, fmt="o-", color="seagreen", capsize=4,
                label="learned population (EMD)")
    ax.errorbar(DISPS, cm, yerr=ce, fmt="s--", color="slategray", capsize=4,
                label="computed local_winding (E7)")
    ax.axhline(0.5, ls=":", color="k", alpha=0.5, label="chance")
    ax.set_xlabel("core displacement from centre (nodes)")
    ax.set_ylabel("chirality-recovery accuracy")
    ax.set_ylim(0, 1.05)
    ax.set_title("A. Learned reader is robust to core meander;\nfixed-locus winding collapses (the C5 problem)")
    ax.legend(fontsize=8)

    ax = axes[1]
    ax.bar(["learned", "computed"], [learned_sw.mean(), computed_sw.mean()],
           yerr=[1.96 * learned_sw.std() / np.sqrt(len(learned_sw)),
                 1.96 * computed_sw.std() / np.sqrt(len(computed_sw))],
           color=["seagreen", "slategray"], capsize=5)
    ax.axhline(0.5, ls=":", color="k", alpha=0.5)
    ax.set_ylabel("E7 switching accuracy (last 4 blocks)")
    ax.set_ylim(0, 1.0)
    ax.set_title("B. Learned readout drives routing\nas well as the computed one")

    fig.suptitle("E7/1b: rotation direction read by a LEARNED population of local "
                 "motion detectors — retires the computed-winding readout", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.93])
    p = os.path.join(FIGDIR, "e7_direction_readout.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e7_direction_readout.npz"),
             disps=np.array(DISPS),
             learned=np.array([learned[d] for d in DISPS]),
             computed=np.array([computed[d] for d in DISPS]),
             switch_learned=learned_sw, switch_computed=computed_sw)
    print("wrote", os.path.join(DATADIR, "e7_direction_readout.npz"))


if __name__ == "__main__":
    main()
