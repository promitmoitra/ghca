"""
E7 (Phase B) - Rotation direction as the learned RULE.

Phase A (experiments/e7_spiral_option.py) showed the substrate holds a single 2D
spiral core whose rotation direction (chirality) is controllable, persistent, and
readable. Phase B makes that rotation direction the RULE of an executive-control /
task-switching problem, exactly as E5 did with a 1D ring -- but here the option is
a genuine 2D spiral and the context is read off the wave:

  * CCW spiral (+charge) -> rule 0 = identity  (action = x)
  * CW  spiral (-charge) -> rule 1 = reversal  (action = 1 - x)
    => action = x XOR rule.

Each trial a direction-selective readout recovers the spiral's rotation direction
from the wave (the sign of the net topological charge; a purely local phase-probe
is reported alongside as the on-thesis alternative). That decoded direction drives
the context input of an E5-style hidden layer of (stimulus x context) conjunction
cells behind a hard-coincidence gate; reward-driven Line A learns the H->M routing.
Nothing is taught per node -- only the scalar reward r = 1[action == x XOR rule].

Predictions (mirroring E5, now carried by a real spiral):
  * switching accuracy is high while the spiral persists; switch cost consolidates.
  * DISCRIMINATOR: ablating the spiral's persistence (raise the lattice threshold so
    the core dies after nucleation) makes the mid-block direction readout chance ->
    switching collapses, while a single-rule control that RE-NUCLEATES the spiral
    each trial (needs no persistence) is spared.
  * the rule is decodable from rotation direction (the fMRI "rotation direction
    classifies the task" result), globally at ~1.0 and locally (phase probe) ~0.9.

Outputs
-------
docs/figures/e7_learning.png
result/e7/e7_learning.npz
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
from ghca_learn import GHLearner
import e7_spiral_option as sp

FIGDIR = os.path.join(ROOT, "docs", "figures")
DATADIR = os.path.join(ROOT, "result", "e7")
os.makedirs(FIGDIR, exist_ok=True)
os.makedirs(DATADIR, exist_ok=True)

# --- router (task network) constants: E5's calibrated hard-coincidence gate ---
K = A = 2
N_S = 8
N_CTX = 16
N_CTX_ACTIVE = 3                 # clamped context nodes (calibrated boost margin)
N_H = 120
N_M = 8
ACT_R = 2
SETTLE, CUE, WWIN = 6, 3, 6
# --- spiral option ---
SPIRAL_STEPS_PER_TRIAL = 4
THETA_LIVE = 4.0                 # spiral persists (E0 band)
THETA_DEAD = 5.0                 # ablation: core dies after nucleation


def build_router(seed=0, channel_bias=0.85, w_sh=0.7, w_ch=0.7, w_hm=0.5,
                 theta_h=3.0, theta_m=1.0, fanin_sh=4, fanin_hm=14, jitter=0.15):
    """S + 2 clamped context groups + hidden (stimulus x context conjunction) + M.
    Same hard-coincidence gate as E5 (stimulus-alone & context-alone subthreshold,
    only their sum fires a hidden unit), but the context is delivered by clamping a
    context group active rather than by an internal ring pulse."""
    rng = np.random.default_rng(seed)
    s0 = 0
    c0 = K * N_S
    h0 = c0 + 2 * N_CTX
    m0 = h0 + N_H
    N = m0 + A * N_M
    sensory = [np.arange(s0 + c * N_S, s0 + (c + 1) * N_S) for c in range(K)]
    ctx = [np.arange(c0 + g * N_CTX, c0 + (g + 1) * N_CTX) for g in range(2)]
    hidden = np.arange(h0, h0 + N_H)
    motor = [np.arange(m0 + c * N_M, m0 + (c + 1) * N_M) for c in range(A)]
    allsen = np.concatenate(sensory)
    W = np.zeros((N, N))
    plastic = np.zeros((N, N), dtype=bool)
    theta = np.zeros(N)
    theta[allsen] = 1.0
    theta[np.concatenate(ctx)] = 99.0            # context nodes only ever clamped
    pref_s = rng.integers(K, size=N_H)
    pref_g = rng.integers(2, size=N_H)
    for i, h in enumerate(hidden):
        spool = sensory[pref_s[i]]
        ns = int(round(channel_bias * fanin_sh))
        src = np.unique(np.concatenate([
            rng.choice(spool, min(ns, len(spool)), replace=False),
            rng.choice(allsen, max(fanin_sh - ns, 0), replace=False)]))
        W[h, src] = w_sh * (1 + jitter * rng.standard_normal(len(src)))
        cg = ctx[pref_g[i]]
        W[h, cg] = w_ch * (1 + jitter * rng.standard_normal(N_CTX))
        theta[h] = theta_h
    for c in range(A):
        for m in motor[c]:
            src = rng.choice(hidden, min(fanin_hm, N_H), replace=False)
            W[m, src] = w_hm * (1 + jitter * rng.standard_normal(len(src)))
            plastic[m, src] = True
            theta[m] = theta_m
    roles = {"sensory": sensory, "ctx": ctx, "hidden": hidden, "motor": motor,
             "N": N, "K": K, "A": A, "n_m": N_M}
    return np.clip(W, 0, None), plastic, roles, theta


def make_router(seed, eta_w=0.06, p_s=3e-3):
    W, plastic, roles, theta = build_router(seed=seed)
    ps_mask = np.zeros(roles["N"], dtype=bool)
    ps_mask[roles["hidden"]] = True
    for c in range(roles["A"]):
        ps_mask[roles["motor"][c]] = True
    net = GHLearner(W, plastic, roles, line="A", act=ACT_R, pas=0, theta=theta,
                    p_s=p_s, p_s_mask=ps_mask, eta_w=eta_w, lam=0.8, w_max=4.0,
                    seed=seed + 41)
    net.tau[:] = ACT_R                            # gapless relay everywhere
    return net, roles


def make_spiral(seed, ablate=False):
    net = Network(lattice2d(sp.L, r=sp.RANGE, periodic=False), act=sp.ACT, pas=sp.PAS,
                  theta=(THETA_DEAD if ablate else THETA_LIVE), p_s=0.0, seed=seed)
    net.tau[:] = sp.TAU
    return net


def spiral_decode(spnet, nsteps=SPIRAL_STEPS_PER_TRIAL):
    """Step the spiral and read its rotation direction from the LOCAL WINDING around
    the (central) core -- a direction-selective readout that reflects the seeded
    handedness immediately (no settling transient) and reads ~0 once the core dies.
    Returns (g_local, g_charge, winding_sum): g in {0,1} (rule 0 = CCW/+, 1 = CW/-),
    or -1 if undetermined. g_charge (global net topological charge) is kept as a
    cross-check."""
    c = sp.L // 2
    s = 0.0
    q = 0.0
    for _ in range(nsteps):
        spnet.step(None)
        s += sp.local_winding(spnet.phi, c, c)
        _, _, netq, _ = sp.signed_charge(spnet.phi)
        q += netq
    g_local = 0 if s > 0.5 else (1 if s < -0.5 else -1)
    g_charge = 0 if q > 0.5 else (1 if q < -0.5 else -1)
    return g_local, g_charge, s


def trial(rnet, roles, g_ctx, x, rule_g, rng, learn=True):
    rnet.reset_traces()
    cnodes = roles["ctx"][g_ctx][:N_CTX_ACTIVE] if g_ctx in (0, 1) else None

    def clamp():
        if cnodes is not None:
            rnet.phi[cnodes] = 1

    for _ in range(SETTLE):
        clamp(); rnet.step_learn(None)
    feats = rnet.features(); V = rnet.value()
    for _ in range(CUE):
        clamp(); rnet.step_learn(rnet.sensory_drive(x))
    sc = np.zeros(A)
    for _ in range(WWIN):
        clamp(); rnet.step_learn(None); sc += rnet.motor_scores()
    action = int(np.argmax(sc + 1e-6 * rng.standard_normal(A))) if sc.sum() > 0 else -1
    target = x ^ rule_g
    r = 1.0 if action == target else 0.0
    if learn:
        rnet.learn(r - V); rnet.update_critic(r - V, feats)
    return r


def run_switching(seed, ablate=False, n_blocks=20, block_len=25):
    rnet, roles = make_router(seed)
    spnet = make_spiral(seed + 7, ablate=ablate)
    rng = np.random.default_rng(seed + 5)
    acc, switch, ctx_ok, charge_ok = [], [], [], []
    for b in range(n_blocks):
        s = +1 if b % 2 == 0 else -1                 # block chirality
        rule_g = 0 if s > 0 else 1
        sp.seed_spiral(spnet, s, jitter=0.02, rng=rng)  # nucleate once per block
        for i in range(block_len):
            gl, gc, _ = spiral_decode(spnet)
            g_ctx = gl if gl in (0, 1) else int(rng.integers(2))
            x = int(rng.integers(2))
            r = trial(rnet, roles, g_ctx, x, rule_g, rng, learn=True)
            acc.append(r); switch.append(1 if i == 0 else 0)
            ctx_ok.append(int(g_ctx == rule_g))
            charge_ok.append(int(gc == rule_g) if gc in (0, 1) else 0)
    return (np.array(acc), np.array(switch), np.array(ctx_ok), np.array(charge_ok))


def run_single_rule(seed, ablate=False, rule_g=0, n_trials=400):
    """Control: RE-NUCLEATE the spiral each trial, so persistence is not needed and
    the direction readout succeeds in both the intact and ablated substrates."""
    rnet, roles = make_router(seed)
    spnet = make_spiral(seed + 7, ablate=ablate)
    rng = np.random.default_rng(seed + 9)
    s = +1 if rule_g == 0 else -1
    acc = []
    for t in range(n_trials):
        sp.seed_spiral(spnet, s, jitter=0.02, rng=rng)      # re-cue each trial
        gl, _, _ = spiral_decode(spnet)
        g_ctx = gl if gl in (0, 1) else int(rng.integers(2))
        x = int(rng.integers(2))
        acc.append(trial(rnet, roles, g_ctx, x, rule_g, rng, learn=True))
    return np.array(acc)


def switch_cost(acc, n_blocks, block_len, post=4):
    return acc.reshape(n_blocks, block_len)[:, :post].mean(1)


def moving_avg(x, w=25):
    return np.convolve(x, np.ones(w) / w, mode="valid")


def main():
    n_seeds = 5
    n_blocks, block_len = 30, 25
    nt = n_blocks * block_len

    print("E7 Phase B: switching (intact) ...")
    ai = np.zeros((n_seeds, nt)); cxi = np.zeros((n_seeds, nt)); qi = np.zeros((n_seeds, nt))
    for s in range(n_seeds):
        a, sw, cx, q = run_switching(s, ablate=False, n_blocks=n_blocks, block_len=block_len)
        ai[s], cxi[s], qi[s] = a, cx, q
    print("E7 Phase B: switching (ablated spiral) ...")
    aa = np.zeros((n_seeds, nt)); cxa = np.zeros((n_seeds, nt))
    for s in range(n_seeds):
        a, sw, cx, q = run_switching(s, ablate=True, n_blocks=n_blocks, block_len=block_len)
        aa[s], cxa[s] = a, cx
    print("E7 Phase B: single-rule (intact / ablated) ...")
    si = np.array([run_single_rule(s, ablate=False, n_trials=600) for s in range(n_seeds)])
    sa = np.array([run_single_rule(s, ablate=True, n_trials=600) for s in range(n_seeds)])

    sci = np.array([switch_cost(ai[s], n_blocks, block_len) for s in range(n_seeds)])

    final_i = ai[:, -block_len * 4:].mean()
    final_a = aa[:, -block_len * 4:].mean()
    ps_early = sci[:, :5].mean(); ps_late = sci[:, -5:].mean()
    single_i = si[:, -120:].mean(); single_a = sa[:, -120:].mean()
    ctx_acc_i = cxi.mean(); ctx_acc_a = cxa.mean()
    charge_acc_i = qi.mean()

    print(f"  switching final:   intact={final_i:.2f}  ablated={final_a:.2f}")
    print(f"  post-switch acc:   early={ps_early:.2f}  late={ps_late:.2f} (intact)")
    print(f"  single-rule final: intact={single_i:.2f}  ablated={single_a:.2f}")
    print(f"  rotation->rule decode acc (local winding): intact={ctx_acc_i:.2f}  ablated={ctx_acc_a:.2f}")
    print(f"  rotation->rule decode acc (global charge, x-check, intact): {charge_acc_i:.2f}")

    # ---------------- figure ----------------
    fig, axes = plt.subplots(1, 3, figsize=(16, 4.8))
    axes[0].plot(np.arange(len(moving_avg(ai.mean(0)))), moving_avg(ai.mean(0)),
                 color="seagreen", label="intact spiral")
    axes[0].plot(np.arange(len(moving_avg(aa.mean(0)))), moving_avg(aa.mean(0)),
                 color="crimson", alpha=0.8, label="ablated spiral")
    for b in range(1, n_blocks):
        axes[0].axvline(b * block_len, color="k", alpha=0.05)
    axes[0].axhline(0.5, ls="--", color="k", alpha=0.4)
    axes[0].set_xlabel("trial (25-trial moving avg)"); axes[0].set_ylabel("accuracy")
    axes[0].set_ylim(0, 1); axes[0].legend(fontsize=8)
    axes[0].set_title("E7: switching driven by spiral rotation direction")

    axes[1].plot(np.arange(n_blocks), sci.mean(0), "o-", color="seagreen")
    axes[1].axhline(0.5, ls="--", color="k", alpha=0.4, label="chance")
    axes[1].set_xlabel("block"); axes[1].set_ylabel("post-switch accuracy (first 4)")
    axes[1].set_ylim(0, 1); axes[1].legend(fontsize=8)
    axes[1].set_title(f"Switch cost consolidates\n({ps_early:.2f} -> {ps_late:.2f})")

    groups = ["switching\n(flexible)", "single-rule\n(re-cued)"]
    xpos = np.arange(2); w = 0.36
    axes[2].bar(xpos - w / 2, [final_i, single_i], w, color="seagreen", label="intact",
                yerr=[ai[:, -block_len*4:].mean(1).std()/np.sqrt(n_seeds),
                      si[:, -120:].mean(1).std()/np.sqrt(n_seeds)], capsize=4)
    axes[2].bar(xpos + w / 2, [final_a, single_a], w, color="crimson", label="ablated",
                yerr=[aa[:, -block_len*4:].mean(1).std()/np.sqrt(n_seeds),
                      sa[:, -120:].mean(1).std()/np.sqrt(n_seeds)], capsize=4)
    axes[2].axhline(0.5, ls="--", color="k", alpha=0.5)
    axes[2].set_xticks(xpos); axes[2].set_xticklabels(groups)
    axes[2].set_ylim(0, 1); axes[2].set_ylabel("final accuracy"); axes[2].legend(fontsize=8)
    axes[2].set_title("Discriminator: spiral needed for\nswitching, not single-rule")

    fig.suptitle("E7 Phase B: rotation direction of a 2D spiral core as the "
                 "executive-control rule", fontsize=12)
    fig.tight_layout(rect=[0, 0, 1, 0.95])
    p = os.path.join(FIGDIR, "e7_learning.png")
    fig.savefig(p, dpi=110)
    print("wrote", p)

    np.savez(os.path.join(DATADIR, "e7_learning.npz"),
             acc_intact=ai, acc_ablated=aa, single_intact=si, single_ablated=sa,
             switchcost_intact=sci, ctx_acc_intact=cxi, ctx_acc_ablated=cxa,
             charge_acc_intact=qi, n_blocks=n_blocks, block_len=block_len)
    print("wrote", os.path.join(DATADIR, "e7_learning.npz"))


if __name__ == "__main__":
    main()
