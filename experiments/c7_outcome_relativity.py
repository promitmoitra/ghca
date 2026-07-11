"""
C7 - Outcome-relativity, mediation, and synthesis for the spiral option.

The C4-analog for the spiral, closing the spiral causal arc. Two results on a
frozen E7 router (trained on the intact spiral):

  A. Outcome-relativity. do(chirality) is causal for the RULE outcome (does the
     router apply identity vs reversal) and EPIPHENOMENAL for the stimulus-CONTENT
     outcome (is the stimulus recoverable from the action) -- while do(g_route)
     (the routing weights) controls content. The same variable (chirality) is
     causal for one outcome and not another: there is no context-free answer to
     "is rotation direction causal?" without naming the outcome.
  B. Mediation / screening. chirality mediates the seed -> behaviour path
     (theta_seed -> chi -> B): intercepting and re-setting chi (do(chi)) overrides
     the nucleation seed, so B follows the injected chirality, not the seed --
     chirality screens off its own generator.

Outcomes (measured as action statistics over stimuli):
  O_rule    = mean_x P(action == x)              (1 = identity applied, 0 = reversal)
  O_content = |P(a=1|x=1) - P(a=1|x=0)|          (1 = stimulus recoverable, 0 = not)

Outputs
-------
docs/figures/c7_outcome_relativity.png
result/c7/c7_data.npz
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
import e7_spiral_option as sp
import e7_learning as e7

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "c7")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)


def train_router(seed, ntrain=1400):
    rnet, roles = e7.make_router(seed)
    rng = np.random.default_rng(seed + 5)
    spnet = e7.make_spiral(seed + 7, ablate=False)
    for t in range(ntrain):
        rule_g = t % 2
        chir = +1 if rule_g == 0 else -1
        sp.seed_spiral(spnet, chir, jitter=0.02, rng=rng)
        gl, _, _ = e7.spiral_decode(spnet)
        g = gl if gl in (0, 1) else int(rng.integers(2))
        e7.trial(rnet, roles, g, int(rng.integers(2)), rule_g, rng, learn=True)
    return rnet, roles, spnet


def route_action(rnet, roles, g_ctx, x, rng):
    rnet.reset_traces()
    cn = roles["ctx"][g_ctx][:e7.N_CTX_ACTIVE] if g_ctx in (0, 1) else None

    def clamp():
        if cn is not None:
            rnet.phi[cn] = 1
    for _ in range(e7.SETTLE):
        clamp(); rnet.step_learn(None)
    for _ in range(e7.CUE):
        clamp(); rnet.step_learn(rnet.sensory_drive(x))
    sc = np.zeros(e7.A)
    for _ in range(e7.WWIN):
        clamp(); rnet.step_learn(None); sc += rnet.motor_scores()
    return int(np.argmax(sc + 1e-6 * rng.standard_normal(e7.A))) if sc.sum() > 0 else -1


def outcomes_for_context(rnet, roles, g, rng, n=80):
    pa = {}
    for x in (0, 1):
        pa[x] = np.mean([route_action(rnet, roles, g, x, rng) == 1 for _ in range(n)])
    O_rule = 0.5 * ((1 - pa[0]) + pa[1])       # P(a==x) averaged over x
    O_content = abs(pa[1] - pa[0])             # stimulus recoverability
    return O_rule, O_content


def decode_context(spnet, chir, rng):
    sp.seed_spiral(spnet, chir, jitter=0.02, rng=rng)
    gl, _, _ = e7.spiral_decode(spnet)
    return gl if gl in (0, 1) else int(rng.integers(2))


def main():
    n_seeds = 3
    mat = np.zeros((2, 2, n_seeds))          # [handle, outcome, seed]; handle 0=do(chi),1=do(g_route)
    screen = np.zeros((2, 2, n_seeds))       # [seed_sign, injected_chi] -> O_rule
    for s in range(n_seeds):
        rnet, roles, spnet = train_router(s)
        rng = np.random.default_rng(s + 31)

        # --- A. outcome matrix ---
        # do(chi): sweep chirality (routing intact)
        gp = decode_context(spnet, +1, rng); Op_rule, Op_con = outcomes_for_context(rnet, roles, gp, rng)
        gm = decode_context(spnet, -1, rng); Om_rule, Om_con = outcomes_for_context(rnet, roles, gm, rng)
        d_chi_rule = abs(Op_rule - Om_rule)
        d_chi_con = abs(Op_con - Om_con)
        # do(g_route): intact vs erased routing (chi fixed +1)
        g = decode_context(spnet, +1, rng)
        Oi_rule, Oi_con = outcomes_for_context(rnet, roles, g, rng)
        Wsave = rnet.W.copy()
        rnet.W[rnet.plastic] = rnet.W[rnet.plastic].mean()      # erase routing selectivity
        Oe_rule, Oe_con = outcomes_for_context(rnet, roles, g, rng)
        rnet.W[:] = Wsave
        d_g_rule = abs(Oi_rule - Oe_rule)
        d_g_con = abs(Oi_con - Oe_con)
        mat[0, 0, s], mat[0, 1, s] = d_chi_rule, d_chi_con
        mat[1, 0, s], mat[1, 1, s] = d_g_rule, d_g_con

        # --- B. mediation / screening: inject do(chi=c) after nucleating seed=si ---
        for i, si in enumerate((+1, -1)):
            for j, cj in enumerate((+1, -1)):
                sp.seed_spiral(spnet, si, jitter=0.02, rng=rng)   # nucleate seed
                for _ in range(6):
                    spnet.step(None)                              # let seed settle
                g = decode_context(spnet, cj, rng)                # do(chi=c): re-set chirality
                O_rule, _ = outcomes_for_context(rnet, roles, g, rng, n=60)
                screen[i, j, s] = O_rule

    M = mat.mean(2)
    # normalize each outcome column by its max range across handles (C4 convention)
    Mn = M / (M.max(0, keepdims=True) + 1e-9)
    scr = screen.mean(2)

    print("C7.A outcome matrix (normalized range of each outcome per handle):")
    print("               O_rule   O_content")
    print(f"  do(chi)      {Mn[0,0]:.2f}     {Mn[0,1]:.2f}")
    print(f"  do(g_route)  {Mn[1,0]:.2f}     {Mn[1,1]:.2f}")
    print("C7.B screening (O_rule; rows = nucleation seed, cols = injected do(chi)):")
    print(f"  seed=+ : inject+={scr[0,0]:.2f}  inject-={scr[0,1]:.2f}")
    print(f"  seed=- : inject+={scr[1,0]:.2f}  inject-={scr[1,1]:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 2, figsize=(12, 5))
    im = axes[0].imshow(Mn, cmap="Greens", vmin=0, vmax=1)
    axes[0].set_xticks([0, 1]); axes[0].set_xticklabels(["O_rule", "O_content"])
    axes[0].set_yticks([0, 1]); axes[0].set_yticklabels(["do(χ)", "do(g_route)"])
    for i in range(2):
        for j in range(2):
            axes[0].text(j, i, f"{Mn[i,j]:.2f}", ha="center", va="center",
                         color="k" if Mn[i, j] < 0.6 else "w", fontsize=13)
    axes[0].set_title("A. Outcome-relativity:\ndo(χ) → rule, not content")
    fig.colorbar(im, ax=axes[0], fraction=0.046)

    im2 = axes[1].imshow(scr, cmap="RdBu_r", vmin=0, vmax=1)
    axes[1].set_xticks([0, 1]); axes[1].set_xticklabels(["inject χ=+", "inject χ=−"])
    axes[1].set_yticks([0, 1]); axes[1].set_yticklabels(["seed +", "seed −"])
    for i in range(2):
        for j in range(2):
            axes[1].text(j, i, f"{scr[i,j]:.2f}", ha="center", va="center", fontsize=13)
    axes[1].set_title("B. Mediation: O_rule follows the INJECTED χ,\n"
                      "not the seed — χ screens its generator")
    fig.colorbar(im2, ax=axes[1], fraction=0.046, label="O_rule (1=identity, 0=reversal)")

    fig.suptitle("C7: rotation direction is a causal MEDIATOR — outcome-relative "
                 "(causal for the rule, not content), not epiphenomenal", fontsize=11)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(FIGDIR, "c7_outcome_relativity.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "c7_data.npz"),
             matrix=M, matrix_norm=Mn, screening=scr)
    print("wrote", os.path.join(DATADIR, "c7_data.npz"))


if __name__ == "__main__":
    main()
