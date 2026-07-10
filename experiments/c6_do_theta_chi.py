"""
C6 - do(theta_chi): the nucleation handle is well-posed; the persistent core is necessary.

C5 showed do(chirality) is fat-handed when the wave is read at a fixed locus
(center band 6.2 sigma) and well-posed only for a topology-aware reader (tracked
1.0 sigma). C6 is the C3-analog: intervene on the GENERATIVE parameter instead of
the wave. The spiral's rotation direction is set by the NUCLEATION seed sign -- a
do(theta) with no realization freedom: it always produces the canonical centred
core, so it sidesteps C5's fat-handedness entirely, and is well-posed even for the
fat-handed fixed-centre reader.

Two results:
  A. Well-posedness. do(theta_chi) (centred nucleation) decoded across seeds gives a
     ~0 sigma band for EVERY readout -- contrast C5's do(chi) bands (6.2 / 1.0 / 2.6).
     "Drive the parameters (nucleate)" is clean; "set the wave (inject a chi
     realization)" is fat-handed.
  B. Necessity (do-ablation). Removing the core's PERSISTENCE (raise the lattice
     threshold so the seeded core dies ~L steps after nucleation) collapses
     switching while sparing the single-rule control (re-cued each trial) -- the
     persistent spiral is necessary for holding the rule across a block, not for the
     routing itself. (Reuses the E7 mechanism as an explicit do-intervention.)

Outputs
-------
docs/figures/c6_do_theta_chi.png
result/c6/c6_data.npz
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
import e7_learning as e7
import c5_do_chirality as c5

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c6")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)


def theta_chi_bands(n_trials=40):
    """do(theta_chi) = the generative nucleation handle: always a centred core.
    Decode across seeds via each readout; the band is the C5-style fat-handedness
    measure. Expect ~0 for every readout (no realization freedom)."""
    net = c5.make_spiral()
    acc = {r: 0.0 for r in c5.READERS}
    per = {r: np.zeros(n_trials) for r in c5.READERS}
    for r in c5.READERS:
        for k in range(n_trials):
            rng = np.random.default_rng(4000 + k)
            chir = +1 if k % 2 == 0 else -1
            c5.realize(net, chir, "centered", rng)     # generative handle: centred core
            g = c5.decode(net, r)
            per[r][k] = float(g == chir)
        acc[r] = per[r].mean()
    bands = {}
    for r in c5.READERS:
        pooled = np.sqrt(per[r].var())
        bands[r] = 0.0 / (pooled + 1e-6) if per[r].max() == per[r].min() else \
            (per[r].max() - per[r].min()) / (pooled + 1e-6)
    return acc, bands


def necessity(n_seeds=5, n_blocks=24, block_len=25):
    """do-ablation of the core's persistence: switching vs single-rule, intact vs
    ablated. Reuses the E7 learning machinery."""
    out = {"switch_intact": [], "switch_ablate": [],
           "single_intact": [], "single_ablate": []}
    for s in range(n_seeds):
        ai, _, _, _ = e7.run_switching(s, ablate=False, n_blocks=n_blocks, block_len=block_len)
        aa, _, _, _ = e7.run_switching(s, ablate=True, n_blocks=n_blocks, block_len=block_len)
        si = e7.run_single_rule(s, ablate=False, n_trials=600)
        sa = e7.run_single_rule(s, ablate=True, n_trials=600)
        out["switch_intact"].append(ai[-block_len * 4:].mean())
        out["switch_ablate"].append(aa[-block_len * 4:].mean())
        out["single_intact"].append(si[-120:].mean())
        out["single_ablate"].append(sa[-120:].mean())
    return {k: (np.mean(v), np.std(v) / np.sqrt(n_seeds)) for k, v in out.items()}


def main():
    print("C6.A: do(theta_chi) well-posedness (centred nucleation) ...")
    acc, tbands = theta_chi_bands()
    for r in c5.READERS:
        print(f"  readout={r:8s} decode acc={acc[r]:.2f}  do(theta_chi) band={tbands[r]:.1f} σ")

    # load C5's do(chi) bands for the contrast
    c5d = np.load(os.path.join(ROOT, "result", "c5", "c5_data.npz"), allow_pickle=True)
    chi_bands = {r: float(c5d[f"band_{r}"]) for r in c5.READERS}
    print("  contrast -- C5 do(chi) bands:",
          " ".join(f"{r}={chi_bands[r]:.1f}σ" for r in c5.READERS))

    print("C6.B: necessity (do-ablation of core persistence) ...")
    nec = necessity()
    print(f"  switching:   intact={nec['switch_intact'][0]:.2f}  ablated={nec['switch_ablate'][0]:.2f}")
    print(f"  single-rule: intact={nec['single_intact'][0]:.2f}  ablated={nec['single_ablate'][0]:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    readers = list(c5.READERS)
    x = np.arange(len(readers)); w = 0.36
    axes[0].bar(x - w / 2, [chi_bands[r] for r in readers], w, color="crimson",
                label="do(χ)  (set the wave)")
    axes[0].bar(x + w / 2, [tbands[r] for r in readers], w, color="seagreen",
                label="do(θ_χ)  (nucleate)")
    axes[0].set_xticks(x); axes[0].set_xticklabels(readers)
    axes[0].set_ylabel("fat-handedness band (σ)")
    axes[0].set_title("A. The generative handle is well-posed for EVERY reader\n"
                      "(do(θ_χ) ≈ 0σ; do(χ) fat-handed at a fixed locus)")
    axes[0].legend(fontsize=8)

    groups = ["switching\n(flexible)", "single-rule\n(re-cued)"]
    xg = np.arange(2)
    intact = [nec["switch_intact"][0], nec["single_intact"][0]]
    ablate = [nec["switch_ablate"][0], nec["single_ablate"][0]]
    ei = [nec["switch_intact"][1], nec["single_intact"][1]]
    ea = [nec["switch_ablate"][1], nec["single_ablate"][1]]
    axes[1].bar(xg - w / 2, intact, w, yerr=ei, capsize=4, color="seagreen", label="intact core")
    axes[1].bar(xg + w / 2, ablate, w, yerr=ea, capsize=4, color="crimson", label="ablated (no persistence)")
    axes[1].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[1].set_xticks(xg); axes[1].set_xticklabels(groups); axes[1].set_ylim(0, 1)
    axes[1].set_ylabel("accuracy")
    axes[1].set_title("B. Necessity: the persistent core is needed for\nswitching, not single-rule routing")
    axes[1].legend(fontsize=8)

    fig.suptitle("C6: do(θ_χ) (nucleation) is the well-posed handle, and the "
                 "persistent spiral core is necessary for switching", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "c6_do_theta_chi.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c6_data.npz"),
             readers=np.array(readers),
             theta_chi_band=np.array([tbands[r] for r in readers]),
             chi_band=np.array([chi_bands[r] for r in readers]),
             **{k: np.array(v) for k, v in nec.items()})
    print("wrote", os.path.join(DATADIR, "c6_data.npz"))


if __name__ == "__main__":
    main()
