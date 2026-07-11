"""
E8.5 - Nested / two-timescale regularities: a slow context gates fast prediction.

E8 (e8_predictive.py) showed prediction of single-timescale sequences and that the
history window is a do(tau) property. E8.5 adds a SECOND, slow timescale: the tone
sequence is piecewise-deterministic under a slowly-switching CONTEXT (regime), and
predicting the next tone needs BOTH the current tone (fast) AND the current regime
(slow). This is the multi-timescale prediction the spec's E8.5 calls for, and it
reuses the E5 "option" idea: a slow persistent trace supplies the regime while the
fast traces predict within it.

Task: an i.i.d. carrier tone stream (which carries NO information about the regime),
plus a regime that is CUED only at each switch and determines a transform of the
current tone -- target = (tone + regime*SHIFT) mod M. Predicting the target needs the
current tone (fast) AND the regime (slow, cued at the switch and then unobserved). The
carrier is i.i.d. so the regime cannot be re-inferred from the tone stream: only a
slow context that PERSISTS from the switch cue carries it. (A deterministic per-regime
sequence rule instead leaks the regime into the recent tones -- then the fast reservoir
infers it and the slow context is not needed; the i.i.d. carrier removes that leak.)

Substrate: fast per-channel traces (tau=14, as E8) + a slow CONTEXT trace bank whose
persistence is set by tau_ctx. A context cue fires the current regime's context node
at each switch; with tau_ctx large it PERSISTS across the regime block (the E5 option
holding the rule), with tau_ctx small it decays and the regime is forgotten mid-block.
A linear next-tone readout reads [fast traces + current tone + context nodes].

Results:
  * do(tau_ctx): nested accuracy rises as the CONTEXT timescale crosses the regime
    duration (tau_ctx >~ T_ctx * ISI) -- the slow option must out-last the block.
  * Discriminator: ablating the slow context (short tau_ctx) collapses nested
    (context-conditioned) prediction to chance-between-regimes, while a
    within-context control (single regime, no switching) is spared.

Outputs
-------
docs/figures/e8_nested.png
result/e8/e8_nested.npz
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
import e8_predictive as e8

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e8")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

M = e8.M
ISI = e8.ISI
TAU_FAST = 14
N_REG = 2
T_CTX = 10                 # regime length in tones (block duration = T_CTX * ISI steps)
SHIFT = 3                  # regime r transforms the current tone by r*SHIFT


def seq_nested(n, rng, switching=True):
    """i.i.d. carrier tones + a switch-cued regime; target = (tone + regime*SHIFT)%M.
    The carrier carries no regime info, so the regime must be HELD from the cue."""
    tones = rng.integers(0, M, size=n)
    regime = np.array([(t // T_CTX) % N_REG if switching else 0 for t in range(n)])
    switches = np.array([True] + [regime[t] != regime[t - 1] for t in range(1, n)])
    targets = (tones + regime * SHIFT) % M
    return tones, regime, switches, targets


def make_net(tau_ctx):
    net = Network(np.zeros((M + N_REG, M + N_REG)), act=e8.ACT, pas=0,
                  theta=99.0, p_s=0.0, seed=0)
    net.tau[:M] = TAU_FAST
    net.tau[M:] = tau_ctx
    return net


def features_targets(tones, regime, switches, targets, tau_ctx):
    net = make_net(tau_ctx)
    feats, targ = [], []
    for t in range(len(tones)):
        c = int(tones[t])
        for s in range(ISI):
            if s == 0:
                net.phi[c] = 1                              # fast: present carrier tone
                if switches[t]:
                    net.phi[M:M + N_REG] = 0                # reset context (E5 ignite_block style)
                    net.phi[M + int(regime[t])] = 1         # slow: context cue at switch
            net.step(None)
        fast = net.phi[:M].astype(float) / TAU_FAST
        ctx = (net.phi[M:] > 0).astype(float)      # which regime is held (binary), not its phase
        onehot = np.zeros(M); onehot[c] = 1.0
        # (current tone x held context) CONJUNCTION -- the fast x slow interaction a
        # linear readout needs to apply a context-dependent transform (cf. E5's
        # conjunction cells). When the slow context is ablated (ctx ~ 0) this term
        # vanishes, so the regime-conditioned prediction is lost.
        interaction = np.outer(onehot, ctx).ravel()
        feats.append(np.concatenate([fast, onehot, ctx, interaction, [1.0]]))
        targ.append(int(targets[t]))
    return np.array(feats), np.array(targ)


def eval_nested(tau_ctx, switching=True, n=6000, seed=0):
    rng = np.random.default_rng(seed)
    tones, reg, sw, tg = seq_nested(n, rng, switching=switching)
    X, y = features_targets(tones, reg, sw, tg, tau_ctx)
    k = len(X) // 2
    W = e8.ridge_fit(X[:k], y[:k])
    pred = (X[k:] @ W).argmax(1)
    return float(np.mean(pred == y[k:]))


def context_available(tau_ctx, seed=0, n=2000):
    """Fraction of tones at which the slow context node for the current regime is
    non-rested (i.e. the regime is actually represented) -- the persistence check."""
    rng = np.random.default_rng(seed)
    tones, reg, sw, _ = seq_nested(n, rng, switching=True)
    net = make_net(tau_ctx)
    hit = tot = 0
    for t in range(len(tones)):
        c = int(tones[t])
        for s in range(ISI):
            if s == 0:
                net.phi[c] = 1
                if sw[t]:
                    net.phi[M:M + N_REG] = 0
                    net.phi[M + int(reg[t])] = 1
            net.step(None)
        hit += int(net.phi[M + int(reg[t])] > 0)
        tot += 1
    return hit / tot


def main():
    taus = [6, 18, 30, 45, 70]
    block_steps = T_CTX * ISI
    nested_acc = [np.mean([eval_nested(tc, True, seed=s) for s in range(3)]) for tc in taus]
    ctx_avail = [np.mean([context_available(tc, seed=s) for s in range(3)]) for tc in taus]
    print(f"E8.5 do(tau_ctx): regime block = {block_steps} steps; chance = {1/M:.3f}")
    for tc, a, ca in zip(taus, nested_acc, ctx_avail):
        print(f"  tau_ctx={tc:2d}: nested_acc={a:.2f}  context_available={ca:.2f}")

    # discriminator: nested vs within-context, intact (long tau_ctx) vs ablated (short)
    tc_intact, tc_ablate = 70, 6
    nested_intact = np.mean([eval_nested(tc_intact, True, seed=s) for s in range(3)])
    nested_ablate = np.mean([eval_nested(tc_ablate, True, seed=s) for s in range(3)])
    within_intact = np.mean([eval_nested(tc_intact, False, seed=s) for s in range(3)])
    within_ablate = np.mean([eval_nested(tc_ablate, False, seed=s) for s in range(3)])
    print("E8.5 discriminator:")
    print(f"  nested (switching):    intact={nested_intact:.2f}  ablated={nested_ablate:.2f}")
    print(f"  within-context (fixed):intact={within_intact:.2f}  ablated={within_ablate:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 4.8))
    ax = axes[0]
    ax.plot(taus, nested_acc, "o-", color="seagreen", label="nested accuracy")
    ax.plot(taus, ctx_avail, "s--", color="steelblue", alpha=0.8, label="context available")
    ax.axvline(block_steps, ls=":", color="k", alpha=0.6, label=f"regime block = {block_steps} steps")
    ax.axhline(1 / M, ls="--", color="gray", alpha=0.5)
    ax.set_xlabel("slow-context timescale τ_ctx  (do(τ_ctx))"); ax.set_ylabel("fraction / accuracy")
    ax.set_ylim(0, 1.05); ax.legend(fontsize=8)
    ax.set_title("E8.5: the slow option must out-last the regime\n"
                 "(nested prediction needs τ_ctx ≳ block duration)")

    groups = ["nested\n(switching)", "within-context\n(fixed regime)"]
    xg = np.arange(2); w = 0.36
    axes[1].bar(xg - w / 2, [nested_intact, within_intact], w, color="seagreen", label="intact slow context")
    axes[1].bar(xg + w / 2, [nested_ablate, within_ablate], w, color="crimson", label="ablated (short τ_ctx)")
    axes[1].axhline(1 / M, ls="--", color="gray", alpha=0.5, label=f"chance={1/M:.2f}")
    axes[1].set_xticks(xg); axes[1].set_xticklabels(groups); axes[1].set_ylim(0, 1.05)
    axes[1].set_ylabel("next-tone accuracy")
    axes[1].set_title("Discriminator: slow context is needed for\nnested prediction, not within-context")
    axes[1].legend(fontsize=8)

    fig.suptitle("E8.5: nested two-timescale prediction — a slow persistent context "
                 "(the E5 option) gates fast next-tone prediction", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.94])
    p = os.path.join(FIGDIR, "e8_nested.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e8_nested.npz"),
             taus=np.array(taus), nested_acc=np.array(nested_acc),
             ctx_avail=np.array(ctx_avail), block_steps=block_steps,
             nested_intact=nested_intact, nested_ablate=nested_ablate,
             within_intact=within_intact, within_ablate=within_ablate)
    print("wrote", os.path.join(DATADIR, "e8_nested.npz"))


if __name__ == "__main__":
    main()
